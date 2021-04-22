/*
 * Worker.c
 *
 *  Created on: Mar 15, 2021
 *      Author: galvez
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <immintrin.h>

#include "Worker.h"
#include "Utilities.h"
#include "NWSW-LSA-512.h"
#include "GlobalData.h"
#include "JobTable.h"

#ifndef KNC
	#include "Vectorization128.h"
	#include "Vectorization256.h"
	#include "Vectorization512.h"
#else
	#include "VectorizationKNC.h"
#endif

int processJob_CISC(struct GlobalData *gd, struct Job *job, struct Node *retFragX, struct Node *retFragY);

/* MAIN FUNCTION OF A WORKER */
/* Enters into a loop to take new jobs to fulfill until there is no more jobs to fulfill.
 * It does use the synchronization functions encapsulated into GlobalData.
 */
void * threadWorker(void * arg) {
	// Takes the globalData as argument
	struct GlobalData *gd = (struct GlobalData *) arg;
	// Sets the function to do the process depending on vectorization support and user preferences
	int (* processJob) (struct GlobalData *, struct Job *, struct Node *, struct Node *);
#ifndef KNC
	switch (gd->vectorization) {
		case SSE3 : processJob = processJob_128; break;
		case AVX2 : processJob = processJob_256; break;
		case AVX512 : processJob = processJob_512; break;
		default: processJob = processJob_CISC;
	}
#else
	processJob = processJob_KNC;
#endif
	// Prepare the best score of all the jobs calculated by this thread in case of S/W
	int bestScore = -1;
	struct Job *bestJob = NULL;
	do {
		struct Job * volatile job = waitForAvailableJob(gd);
		if (job == NULL) break;
		struct Node fragX[job->realSize_X + 1];
		struct Node fragY[job->realSize_Y + 1];
//		printf("Beginning with job (%d,%d)... \n", job->x, job->y);
//		for(int i=0; i<job->realSize_X + 1; i++)
//			displayNode(job->ptrRow+i);
//		printf("\n");
//		for(int j=0; j<job->realSize_Y + 1; j++)
//			displayNode(job->ptrColumn+j);
//		printf("\n");

		/***************************/
		/* MAIN CALL TO PROCESSING */
		/***************************/
		int score = processJob(gd, job, fragX, fragY);
		if (gd->algorithm == SMITH_WATERMAN && score > bestScore) {
			bestScore = score;
			if (gd->pass == FULL_ALIGNMENT) bestJob = job;
		}
//		printf("Returning...\n");
//		for(int i=0; i<job->realSize_X + 1; i++)
//			displayNode(fragX+i);
//		printf("\n");
//		for(int j=0; j<job->realSize_Y + 1; j++)
//			displayNode(fragY+j);
//		printf("\n");
		informFinishedJob(gd, job->x, job->y, fragX, fragY);
	} while(1);
	if (gd->algorithm == SMITH_WATERMAN) {
		saveMyBestScore(gd, bestScore, bestJob);
	}
	pthread_exit(NULL);
}

/* Calculates the  scores in a job */
/* It does the job using NON-VECTOR operations */
int processJob_CISC(struct GlobalData *gd, struct Job *job, struct Node *retFragX, struct Node *retFragY) {
	if (gd->verbose) { block(gd); fprintf(stdout, "Executing non-vectorial job (%d,%d) of size (%dx%d).\n", job->x, job->y, job->realSize_X, job->realSize_Y); unblock(gd); }
	// Prepare the best score of all the jobs calculated by this thread in case of S/W
	int bestScore = -1;
	// Calculus of the initial positions of the sequences to compare
	unsigned queryIdx = job->x * gd->jobTable.fragmentSize_X;
	unsigned subjectIdx = job->y * gd->jobTable.fragmentSize_Y;
	// Create the columns where the nodes will be calculated
	struct Node unus[job->realSize_Y + 1];
	struct Node duo[job->realSize_Y + 1];
	struct Node *prev = unus;
	struct Node *current = duo;
	struct Node *aux;
	// Initializes the first column
	memcpy(prev, job->ptrColumn, sizeof(struct Node)*(job->realSize_Y + 1));
	// Calculates and interchanges columns
	retFragX[0] = prev[job->realSize_Y]; // The first node in the resulting row is already at the end of the starting column
	for(int i=1; i<job->realSize_X+1; i++) {
		current[0] = job->ptrRow[i]; // Initializes the top element (node) of the current row
		for(int j=1; j<job->realSize_Y+1; j++) {
//			printf("top [t=%d, u=%d, s=%d]\n", current[j-1].t, current[j-1].u, current[j-1].s);
//			printf("corner [t=%d, u=%d, s=%d]\n", prev[j-1].t, prev[j-1].u, prev[j-1].s);
//			printf("left [t=%d, u=%d, s=%d]\n", prev[j].t, prev[j].u, prev[j].s);
			current[j].t = max(prev[j].t, prev[j].s - gd->delete) - gd->gapExtend;
			current[j].u = max(current[j-1].u, current[j-1].s - gd->insert) - gd->gapExtend;
			int cellScore = gd->scoreMatrix.matrix[gd->query.dataCoded[queryIdx+i-1]][gd->subject.dataCoded[subjectIdx+j-1]];
			current[j].s = max( max(current[j].t, current[j].u),  prev[j-1].s + cellScore);
			if (gd->algorithm == SMITH_WATERMAN) {
				current[j].s = max(0, current[j].s);
				if (current[j].s > bestScore) bestScore = current[j].s;
			}
//			printf("j=%d [t=%d, u=%d, s=%d] cellScore=%d\n\n", j, current[j].t, current[j].u, current[j].s, cellScore);
		}
		retFragX[i] = current[job->realSize_Y];
		aux = prev;
		prev = current;
		current = aux;
	}
	// Copy the last column into the output parameter
	memcpy(retFragY, prev, sizeof(struct Node)*(job->realSize_Y + 1));
	return (gd->algorithm == SMITH_WATERMAN)? bestScore : 0;
}



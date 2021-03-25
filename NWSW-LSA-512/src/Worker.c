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


int processJob_CISC(struct GlobalData *gd, struct Job *job, struct Node *retFragX, struct Node *retFragY);
int processJob_256(struct GlobalData *gd, struct Job *job, struct Node *retFragX, struct Node *retFragY);

/* MAIN FUNCTION OF A WORKER */
/* Enters into a loop to take new jobs to fulfill until there is no more jobs to fulfill.
 * It does use the synchronization functions encapsulated into GlobalData.
 */
void * threadWorker(void * arg) {
	// Takes the globalData as argument
	struct GlobalData *gd = (struct GlobalData *) arg;
	// Prepare the best score of all the jobs calculated by this thread in case of S/W
	int bestScore = -1;
	struct Job *bestJob = NULL;
	do {
		struct Job *job = waitForAvailableJob(gd);
		if (job == NULL) break;
		struct Node fragX[job->realSize_X + 1];
		struct Node fragY[job->realSize_Y + 1];
		int score = processJob_CISC(gd, job, fragX, fragY);
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
	pthread_exit(NULL);
	if (gd->algorithm == SMITH_WATERMAN) {
		saveMyBestScore(gd, bestScore, bestJob);
	}
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

struct Node_256 {
	__m256i t; // End in horizontal gap
	__m256i u; // End in vertical gap
	__m256i s; // Best value
};
#define lengthVector 8


void initializeVectors_256(int j, struct Node *retFragX, struct GlobalData *gd, struct Job *job, struct Node_256 *arriba, struct Node_256 *izquierda, struct Node_256 *esquina, __m256i *deltaScore);
void calculateAndAdvanceTopLeftDiagonal_256(int j, struct Node *retFragX, int i, struct GlobalData *gd, struct Job *job, struct Node_256 *arriba, struct Node_256 *izquierda, struct Node_256 *esquina, struct Node_256 * resultado, __m256i *deltaScore, int *bestScore);
void calculateAndAdvanceBody_256(int j, struct Node *retFragX, int i, struct GlobalData *gd, struct Job *job, struct Node_256 *arriba, struct Node_256 *izquierda, struct Node_256 *esquina, struct Node_256 * resultado, __m256i *deltaScore, int *bestScore);
void calculateAndAdvanceBottomRightDiagonal_256(int j, struct Node *retFragX, int i, struct GlobalData *gd, struct Job *job, struct Node_256 *arriba, struct Node_256 *izquierda, struct Node_256 *esquina, struct Node_256 * resultado, struct Node *retFragY, __m256i *deltaScore, int *bestScore);
void calculate_256(struct GlobalData *gd, struct Node_256 *arriba, struct Node_256 *izquierda, struct Node_256 *esquina, struct Node_256 * resultado, __m256i *deltaScore, int *bestScore);
void setPos(struct Node_256 *ptrNodeDest_256, int i, struct Node *ptrNodeSrc);
void savePos(struct Node *ptrNodeDest, struct Node_256 *ptrNodeSrc_256, int i);
int reduceMax_256(__m256i *ptrVector);

/* DEBUG FUNCTIONS */
void displayNode_256(struct Node_256 *node256);
void displayVector_256(__m256i *vector);

/* Calculates the  scores in a job */
/* It does the job using TRUE-VECTOR operations */
int processJob_256(struct GlobalData *gd, struct Job *job, struct Node *retFragX, struct Node *retFragY) {
	if (job->realSize_Y % lengthVector != 0 || job->realSize_X <= 2*lengthVector ) return processJob_CISC(gd, job, retFragX, retFragY);
	// Prepare the best score of all the jobs calculated by this thread in case of S/W
	int bestScore = -1;
	// Variable to store values of the score matrix
	__m256i deltaScore;
	// Nodes of vector to work with
	struct Node_256 arriba, izquierda, esquina, resultado;
	// Initializes the top right of the resulting job
	retFragY[0] = job->ptrRow[job->realSize_X];
	// Initializes retFragX. retFragX changes in each pass of vectorization
	memcpy(retFragX, job->ptrRow, sizeof(struct Node)*(job->realSize_X + 1));
	for(int j=0; j<job->realSize_Y - 1; j+=lengthVector) {
		initializeVectors_256(j, retFragX, gd, job, &arriba, &izquierda, &esquina, &deltaScore);
		int i=1;
		for(; i < lengthVector; i++) // Step 1
			calculateAndAdvanceTopLeftDiagonal_256(j, retFragX, i, gd, job, &arriba, &izquierda, &esquina, &resultado, &deltaScore, &bestScore);
		// Update progressively retFragX
		retFragX[0] = job->ptrColumn[j+lengthVector];
		for(; i < job->realSize_X; i++) // Main step
			calculateAndAdvanceBody_256(j, retFragX, i, gd, job, &arriba, &izquierda, &esquina, &resultado, &deltaScore, &bestScore);
		for(; i < job->realSize_X + lengthVector; i++) // Step 3
			calculateAndAdvanceBottomRightDiagonal_256(j, retFragX, i, gd, job, &arriba, &izquierda, &esquina, &resultado, retFragY, &deltaScore, &bestScore);
	}
	return (gd->algorithm == SMITH_WATERMAN)? bestScore : 0;
}

/* Initialization of vectors for an horizontal round of vectorization */
//
void initializeVectors_256(int j, struct Node *retFragX, struct GlobalData *gd, struct Job *job, struct Node_256 *arriba, struct Node_256 *izquierda, struct Node_256 *esquina, __m256i *deltaScore) {
	// Fill first node of esquina with retFragX[0]
	// The rest of nodes may be garbage
	esquina->t = _mm256_maskz_set1_epi32(0x01, retFragX[0].t); // Superfluous
	esquina->u = _mm256_maskz_set1_epi32(0x01, retFragX[0].u); // Superfluous
	esquina->s = _mm256_maskz_set1_epi32(0x01, retFragX[0].s);
	// Fill first node of izquierda with job->ptrColumn[j+1]
	// The rest of nodes may be garbage
	izquierda->t = _mm256_maskz_set1_epi32(0x01, job->ptrColumn[j+1].t);
	izquierda->u = _mm256_maskz_set1_epi32(0x01, job->ptrColumn[j+1].u); // Superfluous
	izquierda->s = _mm256_maskz_set1_epi32(0x01, job->ptrColumn[j+1].s);
	// Fill first node of arriba with retFragX[1]; and the second one with job->ptrColumn[j+1]
	// The rest of nodes may be garbage
	// First
	arriba->t = _mm256_maskz_set1_epi32(0x01, retFragX[1].t); // Superfluous
	arriba->u = _mm256_maskz_set1_epi32(0x01, retFragX[1].u);
	arriba->s = _mm256_maskz_set1_epi32(0x01, retFragX[1].s);
	// Second
	arriba->t = _mm256_mask_set1_epi32 (arriba->t, 0x02, job->ptrColumn[j+1].t);
	arriba->u = _mm256_mask_set1_epi32 (arriba->u, 0x02, job->ptrColumn[j+1].u);
	arriba->s = _mm256_mask_set1_epi32 (arriba->s, 0x02, job->ptrColumn[j+1].s);
	// Fill first node of deltaScore with data from the score matrix
	// The rest of nodes may be garbage
	// Calculus of the initial positions of the sequences to compare
	unsigned queryIdx = job->x * gd->jobTable.fragmentSize_X;
	unsigned subjectIdx = job->y * gd->jobTable.fragmentSize_Y + j;
	// We have to remember that the sequences begin at position 0
	int cellScore = gd->scoreMatrix.matrix[gd->query.dataCoded[queryIdx]][gd->subject.dataCoded[subjectIdx]];
	*deltaScore = _mm256_maskz_set1_epi32(0x01, cellScore);
}

/* Calculus of a single vector in diagonal */
//
void calculate_256(struct GlobalData *gd, struct Node_256 *arriba, struct Node_256 *izquierda, struct Node_256 *esquina, struct Node_256 * resultado, __m256i *deltaScore, int *bestScore) {
	// Creation of vectorized constants
	__m256i gapExtend_256 = _mm256_maskz_set1_epi32(0xFF, gd->gapExtend);
	__m256i openGap_256 = _mm256_maskz_set1_epi32(0xFF, gd->delete + gd->gapExtend);
	// Deletion (horizontal gap). Calculus of resultado->t
	__m256i aux1 = _mm256_maskz_sub_epi32 (0xFF, izquierda->t, gapExtend_256);
	__m256i aux2 = _mm256_maskz_sub_epi32 (0xFF, izquierda->s, openGap_256);
	resultado->t = _mm256_maskz_max_epi32 (0xFF, aux1, aux2);
	// Insertion (vertical gap). Calculus of resultado->u
	aux1 = _mm256_maskz_sub_epi32 (0xFF, arriba->u, gapExtend_256);
	aux2 = _mm256_maskz_sub_epi32 (0xFF, arriba->s, openGap_256);
	resultado->u = _mm256_maskz_max_epi32 (0xFF, aux1, aux2);
	// Continue in diagonal. Calculus of resultado->s
	aux1 = _mm256_maskz_add_epi32 (0xFF, esquina->s, *deltaScore);
	aux1 = _mm256_maskz_max_epi32 (0xFF, aux1, resultado->t);
	resultado->s = _mm256_maskz_max_epi32 (0xFF, aux1, resultado->u);
	if (gd->algorithm == SMITH_WATERMAN) {
		__m256i zeroes = _mm256_maskz_set1_epi32(0xFF, 0);
		resultado->s = _mm256_maskz_max_epi32 (0xFF, zeroes, resultado->s);
		int aux = reduceMax_256(&resultado->s);
		if (aux > *bestScore) *bestScore = aux;
	}
}

void calculateAndAdvanceTopLeftDiagonal_256(int j, struct Node *retFragX, int i, struct GlobalData *gd, struct Job *job, struct Node_256 *arriba, struct Node_256 *izquierda, struct Node_256 *esquina, struct Node_256 * resultado, __m256i *deltaScore, int *bestScore) {
	calculate_256(gd, arriba, izquierda, esquina, resultado, deltaScore, bestScore);
	// ADVANCES esquina-corner
		esquina->t = arriba->t;
		esquina->u = arriba->u;
		esquina->s = arriba->s;
	// ADVANCES izquierda-left: shifts and set pos i+1
		izquierda->t = resultado->t;
		izquierda->u = resultado->u;
		izquierda->s = resultado->s;
		setPos(izquierda, i, &job->ptrColumn[j+i+1]);
	// ADVANCES arriba-top: shifts, sets pos 0 and pos i+1
		// Prepare shuffle
		int arribaShuffle[] __attribute__((aligned (32))) = {0,0,1,2,3,4,5,6} ;
		__m256i arribaMask = _mm256_maskz_load_epi32  (0xFF, arribaShuffle);
		// Shift t, u and s
		arriba->t = _mm256_permutexvar_epi32 (arribaMask, resultado->t);
		arriba->u = _mm256_permutexvar_epi32 (arribaMask, resultado->u);
		arriba->s = _mm256_permutexvar_epi32 (arribaMask, resultado->s);
		setPos(arriba, 0, &(retFragX[i+1]));
		setPos(arriba, i+1, &job->ptrColumn[j+i+1]);
	// Recalculates deltaScore
		// Calculus of the initial positions of the sequences to compare
		unsigned queryIdx = job->x * gd->jobTable.fragmentSize_X + i;
		unsigned subjectIdx = job->y * gd->jobTable.fragmentSize_Y + j;
		int dS[lengthVector] __attribute__((aligned (32)));
		for(int offset=0; offset <= i; offset++) {
			int cellScore = gd->scoreMatrix.matrix[gd->query.dataCoded[queryIdx - offset]][gd->subject.dataCoded[subjectIdx + offset]];
			dS[offset] = cellScore;
		}
		*deltaScore = _mm256_maskz_load_epi32  (0xFF, dS);

}

void calculateAndAdvanceBody_256(int j, struct Node *retFragX, int i, struct GlobalData *gd, struct Job *job, struct Node_256 *arriba, struct Node_256 *izquierda, struct Node_256 *esquina, struct Node_256 * resultado, __m256i *deltaScore, int *bestScore) {
	calculate_256(gd, arriba, izquierda, esquina, resultado, deltaScore, bestScore);
	// Saves last item of resultado into retFragX
		savePos(&retFragX[i - lengthVector + 1], resultado, lengthVector - 1);
	// ADVANCES esquina-corner
		esquina->t = arriba->t;
		esquina->u = arriba->u;
		esquina->s = arriba->s;
	// ADVANCES izquierda-left: shifts
		izquierda->t = resultado->t;
		izquierda->u = resultado->u;
		izquierda->s = resultado->s;
	// ADVANCES arriba-top: shifts and sets pos 0.
		// Prepare shuffle
		int arribaShuffle[] __attribute__((aligned (32))) = {0,0,1,2,3,4,5,6} ;
		__m256i arribaMask = _mm256_maskz_load_epi32  (0xFF, arribaShuffle);
		// Shift t, u and s
		arriba->t = _mm256_permutexvar_epi32 (arribaMask, resultado->t);
		arriba->u = _mm256_permutexvar_epi32 (arribaMask, resultado->u);
		arriba->s = _mm256_permutexvar_epi32 (arribaMask, resultado->s);
		setPos(arriba, 0, &(retFragX[i+1]));
	// Recalculates deltaScore
		// Calculus of the initial positions of the sequences to compare
		unsigned queryIdx = job->x * gd->jobTable.fragmentSize_X + i;
		unsigned subjectIdx = job->y * gd->jobTable.fragmentSize_Y + j;
		int dS[lengthVector] __attribute__((aligned (32)));
		for(int offset=0; offset < lengthVector; offset++) {
			int cellScore = gd->scoreMatrix.matrix[gd->query.dataCoded[queryIdx - offset]][gd->subject.dataCoded[subjectIdx + offset]];
			dS[offset] = cellScore;
		}
		*deltaScore = _mm256_maskz_load_epi32  (0xFF, dS);

}

void calculateAndAdvanceBottomRightDiagonal_256(int j, struct Node *retFragX, int i, struct GlobalData *gd, struct Job *job, struct Node_256 *arriba, struct Node_256 *izquierda, struct Node_256 *esquina, struct Node_256 * resultado, struct Node *retFragY, __m256i *deltaScore, int *bestScore) {
	calculate_256(gd, arriba, izquierda, esquina, resultado, deltaScore, bestScore);
	int progress = i - job->realSize_X;
	// Saves last valid item of resultado into retFragX
		savePos(&retFragX[i - lengthVector + 1], resultado, lengthVector - 1);
	// Saves first valid item of resultado into retFragY
		savePos(&retFragY[j + progress + 1], resultado, progress);
//	if (i==18) {
//		printf("---------------i=%d:\n", i);
//		printf("---------------Arriba:\n");
//		displayNode_256(arriba);
//		printf("---------------Esquina:\n");
//		displayNode_256(esquina);
//		printf("---------------Izquierda:\n");
//		displayNode_256(izquierda);
//		printf("---------------Resultado:\n");
//		displayNode_256(resultado);
//		printf("---------------Fragment temporal X:\n");
//		for(int x=0; x<=job->realSize_X; x++)
//			displayNode(retFragX+x);
//		printf("\n");
//		printf("---------------Fragment temporal Y:\n");
//		for(int y=0; y<=job->realSize_Y; y++)
//			displayNode(retFragY+y);
//		printf("\n");
//		displayVector_256(deltaScore);
//		exit(0);
//	}
	// ADVANCES esquina-corner
		esquina->t = arriba->t;
		esquina->u = arriba->u;
		esquina->s = arriba->s;
	// ADVANCES izquierda-left: shifts
		izquierda->t = resultado->t;
		izquierda->u = resultado->u;
		izquierda->s = resultado->s;
	// ADVANCES arriba-top: only shifts
		// Prepare shuffle
		int arribaShuffle[] __attribute__((aligned (32))) = {0,0,1,2,3,4,5,6} ;
		__m256i arribaMask = _mm256_maskz_load_epi32  (0xFF, arribaShuffle);
		// Shift t, u and s
		arriba->t = _mm256_permutexvar_epi32 (arribaMask, resultado->t);
		arriba->u = _mm256_permutexvar_epi32 (arribaMask, resultado->u);
		arriba->s = _mm256_permutexvar_epi32 (arribaMask, resultado->s);
	// Recalculates deltaScore
		// Calculus of the initial positions of the sequences to compare
		unsigned queryIdx = job->x * gd->jobTable.fragmentSize_X + i;
		unsigned subjectIdx = job->y * gd->jobTable.fragmentSize_Y + j;
		int dS[lengthVector] __attribute__((aligned (32)));
		for(int offset=progress; offset < lengthVector; offset++) {
			int cellScore = gd->scoreMatrix.matrix[gd->query.dataCoded[queryIdx - offset]][gd->subject.dataCoded[subjectIdx + offset]];
			dS[offset] = cellScore;
		}
		*deltaScore = _mm256_maskz_load_epi32  (0xFF, dS);

}

/* Puts a node into a vector of nodes. Puts it at position 'i' */
//
void setPos(struct Node_256 *ptrNodeDest_256, int i, struct Node *ptrNodeSrc) {
	__mmask8 mask = 0x01 << i;
	ptrNodeDest_256->t = _mm256_mask_set1_epi32(ptrNodeDest_256->t, mask, ptrNodeSrc->t);
	ptrNodeDest_256->u = _mm256_mask_set1_epi32(ptrNodeDest_256->u, mask, ptrNodeSrc->u);
	ptrNodeDest_256->s = _mm256_mask_set1_epi32(ptrNodeDest_256->s, mask, ptrNodeSrc->s);
}

/* Gets a node from a vector of nodes. Gets the node at position 'i' */
//
void savePos(struct Node *ptrNodeDest, struct Node_256 *ptrNodeSrc_256, int i) {
	int dump[lengthVector] __attribute__((aligned (32)));
	_mm256_store_epi64(dump, ptrNodeSrc_256->t);
	ptrNodeDest->t = dump[i];
	_mm256_store_epi64(dump, ptrNodeSrc_256->u);
	ptrNodeDest->u = dump[i];
	_mm256_store_epi64(dump, ptrNodeSrc_256->s);
	ptrNodeDest->s = dump[i];

}

/* Returns the max value stored in a vector */
//
int reduceMax_256(__m256i *ptrVector) {
	int dump[lengthVector] __attribute__((aligned (32)));
	_mm256_store_epi64(dump, *ptrVector);
	int ret = dump[0];
	for(int i=1; i<lengthVector - 1; i++)
		if (dump[i] > ret) ret = dump[i];
	return ret;
}

// Displays a struct Node_256 in a formated way
void displayNode_256(struct Node_256 *node256) {
	int t[lengthVector] __attribute__((aligned (32)));
	_mm256_store_epi64(t, node256->t);
	int u[lengthVector] __attribute__((aligned (32)));
	_mm256_store_epi64(u, node256->u);
	int s[lengthVector] __attribute__((aligned (32)));
	_mm256_store_epi64(s, node256->s);

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wmisleading-indentation"
	for(int i=0; i<lengthVector; i++)
		printf("| t=%8d \t", t[i]); printf("\n");
	for(int i=0; i<lengthVector; i++)
		printf("| u=%8d \t", u[i]); printf("\n");
	for(int i=0; i<lengthVector; i++)
		printf("| s=%8d \t", s[i]); printf("\n");
#pragma GCC diagnostic pop
}

void displayVector_256(__m256i *vector) {
	int v[lengthVector] __attribute__((aligned (32)));
	_mm256_store_epi64(v, *vector);
	for(int i=0; i<lengthVector; i++)
		printf("| v=%8d \t", v[i]);
	printf("\n");
}



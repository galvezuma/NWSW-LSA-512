/*
 * GlobalData.c
 *
 *  Created on: Mar 11, 2021
 *      Author: galvez
 */

#include <string.h>
#include <unistd.h>
#include <sys/syscall.h>
#include "GlobalData.h"
#include "JobTable.h"

void copyUserParameters(struct GlobalData *gd, struct UserParameters *up) {
	gd->pass = up->pass;
	gd->tree = up->tree;
	gd->algorithm = up->algorithm;
//	gd->info = up->info;
	gd->insert = up->insert;
	gd->delete = up->delete;
	gd->matchReplace = up->matchReplace;
	gd->gapExtend = up->gapExtend;
	gd->threads = up->threads;
	gd->verbose = up->verbose;
}

/* Function of synchronization */
// Waits until a new job is available to be processed
struct Job * volatile waitForAvailableJob(struct GlobalData *gd) {
	struct Job * volatile ret = NULL;
    pthread_mutex_lock(&gd->globalDataAccess_mutex);
        while(isEmptyStack(&gd->stackOfJobs) && ! gd->jobTableFulfilled)
            pthread_cond_wait(&gd->jobAvailable_condition, &gd->globalDataAccess_mutex);
        if (! isEmptyStack(&gd->stackOfJobs)) {
        	ret = pop(&gd->stackOfJobs);
//    		printf("Working with job (%d,%d)\n", ret->x, ret->y);
//    		displayJob(ret);
        }
    pthread_mutex_unlock(&gd->globalDataAccess_mutex);
    return ret;
}

// Updates the job table, inserts new jobs in the stack when needed and notifies other threads if the stack has changed.
// If the inserted job is the last one, sets the variable 'jobTableFulfilled'
void informFinishedJob(struct GlobalData *gd, unsigned x, unsigned y, struct Node *resultingFragmentX, struct Node *resultingFragmentY) {
//    pthread_mutex_lock(&gd->globalDataAccess_mutex);
//    fprintf(stdout, "Finishing job %d,%d\n", x, y);
//    pthread_mutex_unlock(&gd->globalDataAccess_mutex);
	struct Job *jobDone = getJob(&gd->jobTable, x, y);
	if (x == gd->jobTable.numFragments_X - 1 && y == gd->jobTable.numFragments_Y - 1) { // The just calculated  job is the last one (right-bottom corner)
	    pthread_mutex_lock(&gd->globalDataAccess_mutex);
	    if (gd->algorithm == NEEDLEMAN_WUNSCH) gd->bestScore = resultingFragmentX[jobDone->realSize_X].s;
		gd->jobTableFulfilled = 1;
		pthread_cond_broadcast(&gd->jobAvailable_condition);
	    pthread_mutex_unlock(&gd->globalDataAccess_mutex);
	}
//	printf("(%d,%d) Corner score: %d\n", x, y, resultingFragmentX[jobDone->realSize_X].s);
//	if (x==1 && y==1) exit(0);

	// Let's put the calculated row in the job below
	if (y < gd->jobTable.numFragments_Y - 1) {
		struct Job *jobBelow = getJob(&gd->jobTable, x, y+1);
		// The fragment X is reused in case of ONLY_SCORE; otherwise, we point to a new fragment already allocated
		struct Node * jobBelow_ptrRow_aux;
		jobBelow_ptrRow_aux = (gd->pass == ONLY_SCORE)? jobDone->ptrRow : gd->jobTable.rows[y+1] + x * gd->jobTable.fragmentSize_X;
		memcpy(jobBelow_ptrRow_aux, resultingFragmentX, sizeof(struct Node)*(1+jobBelow->realSize_X));
		//
		// Modification and check condition must be done atomically
	    pthread_mutex_lock(&gd->globalDataAccess_mutex);
		jobBelow->ptrRow = jobBelow_ptrRow_aux;
		if (jobBelow->ptrColumn != NULL) { // If the job below is ready to be fulfilled
			push(&gd->stackOfJobs, jobBelow); // Insert into the stack in mutual exclusion
			pthread_cond_signal(&gd->jobAvailable_condition); // Notifies one thread that a new job is ready to work with
		}
	    pthread_mutex_unlock(&gd->globalDataAccess_mutex);
	}
	//
	// Let's put the calculated column in the job at the right
	if (x < gd->jobTable.numFragments_X - 1) {
		struct Job *jobRight = getJob(&gd->jobTable, x+1, y);
		// The fragment Y is reused in case of ONLY_SCORE; otherwise, We point to a new fragment already allocated
		struct Node * jobRight_ptrColumn_aux;
		jobRight_ptrColumn_aux = (gd->pass == ONLY_SCORE)? jobDone->ptrColumn : gd->jobTable.columns[x+1] + y * gd->jobTable.fragmentSize_Y;
		memcpy(jobRight_ptrColumn_aux, resultingFragmentY, sizeof(struct Node)*(1+jobRight->realSize_Y));
		//
		// Modification and check condition must be done atomically
	    pthread_mutex_lock(&gd->globalDataAccess_mutex);
		jobRight->ptrColumn = jobRight_ptrColumn_aux;
		if (jobRight->ptrRow != NULL) { // If the job at the right is ready to be fulfilled
			push(&gd->stackOfJobs, jobRight); // Insert into the stack in mutual exclusion
			pthread_cond_signal(&gd->jobAvailable_condition); // Notifies one thread that a new job is ready to work with
		}
	    pthread_mutex_unlock(&gd->globalDataAccess_mutex);
	}
}

void saveMyBestScore(struct GlobalData *gd, int score, struct Job *bestJob) {
    pthread_mutex_lock(&gd->globalDataAccess_mutex);
    if (gd->bestScore < score) {
    	gd->bestScore = score;
		if (gd->pass == FULL_ALIGNMENT) gd->bestJob =  bestJob;
    }
    pthread_mutex_unlock(&gd->globalDataAccess_mutex);
}

void block(struct GlobalData *gd) {
    pthread_mutex_lock(&gd->verboseStdOut_mutex);
}

void unblock(struct GlobalData *gd) {
    pthread_mutex_unlock(&gd->verboseStdOut_mutex);
}

void checkSupport(struct GlobalData *gd) {
	if (__builtin_cpu_supports ("avx512f")) {
		gd->lengthVector = 16;
		if (gd->verbose) fprintf(stdout, "AVX512 vectorization supported.\n");
	} else if (__builtin_cpu_supports ("avx2")) {
		gd->lengthVector = 8;
		if (gd->verbose) fprintf(stdout, "AVX2 vectorization supported.\n");
	} else if (__builtin_cpu_supports ("avx")) {
		gd->lengthVector = 4;
		if (gd->verbose) fprintf(stdout, "AVX vectorization supported.\n");
	}
}

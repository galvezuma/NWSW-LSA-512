/*
 * Alignment.c
 *
 *  Created on: Apr 7, 2021
 *      Author: galvez
 */

#include <stdlib.h>
#include <string.h>
#include <stdbool.h>

#include "SingleAlignment.h"
#include "Utilities.h"
#include "NWSW-LSA-512.h"
#include "GlobalData.h"
#include "ScoreMatrix.h"
#include "Fasta.h"
#include "Worker.h"
#include "BackwardsPairwise.h"


// Global data accesible by all threadAligner
// This is used only for synchronization purposes
// There is a single instance of this structure per multiple alignment task (at all)
//
struct AlignerSynchronization {
	// Mutex to access to this structure in mutual exclusion
	pthread_mutex_t alignerSynchronizationAccess_mutex;
	pthread_cond_t alignerAvailable_condition;
	// Data to be accessed in mutual exclusion
	int availableAligners;
	int scoreSize;
	int * scores;
};

int executePairwise(struct GlobalData *gd);


/********************************************************/
/* Main Function to calculate Single Pairwise alignment */
/********************************************************/
//
int singlePairwise(struct UserParameters *ptrUserParams) {
	/* LOAD FILES */
	// The next structure is used throughout the functions of the program
	struct GlobalData globalData;
	// Loads globally the User parameters
	copyUserParameters(&globalData, ptrUserParams);

	// Reads the Score Matrix file or loads the default one
	if (! strcmp(ptrUserParams->matrixFilename, "")) {
		loadDefaultMatrix(&globalData);
	} else {
		loadUserMatrix(ptrUserParams->matrixFilename, &globalData);
	}

	// Reads the Query and Subject sequence Fasta files
	readSingleFastaFile(ptrUserParams->queryFilename, &(globalData.query));
	readSingleFastaFile(ptrUserParams->subjectFilename, &(globalData.subject));

	// VERBOSE
	if (ptrUserParams->verbose) {
		fprintf(stdout, "nwsw_lsa v%.1f. Universidad de Málaga. Spain.\n", VERSION);
		displayUserParameters(ptrUserParams, stdout);
		fprintf(stdout, "****** INPUT SEQUENCES ******\n");
		fprintf(stdout, "Query (%ld nt): %s\n", strlen(globalData.query.data), globalData.query.name);
		fprintf(stdout, "Subject (%ld nt) %s\n", strlen(globalData.subject.data), globalData.subject.name);
		fprintf(stdout, "****** SCORE MATRIX ******\n");
		displayScoreMatrix(stdout, &(globalData.scoreMatrix));
	}

	/* INITIALIZATION */
	// Checks the letters in the sequences and put them in Uppercase
	int invalidLettersQuery, invalidLettersSubject;
	// Checks the Query
	invalidLettersQuery = toUpperCodeAndCheck(	&(globalData.query),
												globalData.scoreMatrix.horizontalAlphabet,
												globalData.scoreMatrix.horizontalCodification
											 );
	if (invalidLettersQuery != 0) fprintf(stderr, "Query has %d letters not found in the alphabet.\n", invalidLettersQuery);

	// Checks the subject
	invalidLettersSubject = toUpperCodeAndCheck(&(globalData.subject),
												globalData.scoreMatrix.verticalAlphabet,
												globalData.scoreMatrix.verticalCodification
											   );
	if (invalidLettersSubject != 0) fprintf(stderr, "Subject has %d letters not found in the alphabet.\n", invalidLettersSubject);

	/***********************************/
	/* MAIN CALL TO PAIRWISE ALIGNMENT */
	/***********************************/
	int score = executePairwise(&globalData);
	fprintf(stdout, "The score is: %d\n", score);
	/***********************************/
	/***********************************/

	/* RELEASE MEMORY */
	freeScoreMatrixStruct(&(globalData.scoreMatrix));
	freeSequenceStruct(&(globalData.query));
	freeSequenceStruct(&(globalData.subject));
	return(EXIT_SUCCESS);
}

/* AUXILIARY FUNCTION TO PERFORM ONE PAIRWISE ALIGNMENT */
/*    This is called from Single and Multi pairwise     */
//
int executePairwise(struct GlobalData *gd) {
	/* CHECKS SUPPORT FOR VECTORIZATION */
	enum Vectorization userParamVectorization = gd->vectorization;
	checkSupport(gd);
	if (gd->vectorization < userParamVectorization)
		fprintf(stdout, "%s not supported. Changing to %s.\n", enumVectorizationToString(userParamVectorization), enumVectorizationToString(gd->vectorization));
	else
		gd->vectorization = userParamVectorization;
	if (gd->verbose) { fprintf(stdout, "Using vectorization: %s.\n", enumVectorizationToString(gd->vectorization)); }
	/* ESTIMATES THE FRAGMENTS SIZE AND CREATES THE TABLE OF JOBS */
	createJobTable(gd);
	/* INITIALIZES THE STACK OF JOBS */
	gd->stackOfJobs = createStack(gd, max(gd->jobTable.numFragments_X, gd->jobTable.numFragments_Y));
	push(&gd->stackOfJobs, getJob(&gd->jobTable, 0, 0)); // The job at position 0,0 is the only available by now

	/* MAIN PROCESSING */
//	displayJob(getJob(&gd->jobTable, 0, 0));
//	displayJob(getJob(&globalData.jobTable, 1, 0));
//	printf("%d, %d\n", globalData.jobTable.numFragments_X, globalData.jobTable.numFragments_Y);
		long initTime = myClock();
		gd->bestScore = -INFINITE;
		// PREPARE SYNCHRONIZATION AMONG WORKERS
		int error;
		gd->jobTableFulfilled = false;
		if (gd->verbose) {
			error = pthread_mutex_init(&gd->verboseStdOut_mutex, NULL);
			if (error) fatalError0("Unable to initialize the verbose mutex.\n");
		}
		error = pthread_mutex_init(&gd->globalDataAccess_mutex, NULL);
		if (error) fatalError0("Unable to initialize the worker's mutex.\n");
		error = pthread_cond_init(&gd->jobAvailable_condition, NULL);
		if (error) fatalError0("Unable to initialize the worker's synchronization condition.\n");
		// CREATE AND LAUNCH THREADS
		if (gd->verbose) printf("Creating %d threads.\n", gd->threads);
		pthread_t * threads= (pthread_t *) internalMalloc(sizeof(pthread_t) * gd->threads);
		for(int i=0; i<gd->threads; i++){
			error = pthread_create(&(threads[i]), NULL, threadWorker, (void *)gd);
			if (error) fatalError1("Unable to create the thread num. %d\n", i);
		}
		// WAIT FOR THREADS
		for(int i=0; i<gd->threads; i++){
			error = pthread_join(threads[i], NULL);
			if (error) fatalError1("Unable to join to the thread num. %d\n", i);
		}
		// DESTROY THREADS AND SYNCHRONIZATION ELEMENTS
		error = pthread_cond_destroy(&gd->jobAvailable_condition);
		if (error) fatalError0("Unable to destroy the worker's synchronization condition.\n");
		error = pthread_mutex_destroy(&gd->globalDataAccess_mutex);
		if (error) fatalError0("Unable to destroy the worker's mutex.\n");
		if (gd->verbose) {
			error = pthread_mutex_destroy(&gd->verboseStdOut_mutex);
			if (error) fatalError0("Unable to destroy the verbose mutex.\n");
		}
		internalFree((void **) &threads);
		// VERBOSE DISPLAY TIME
		long endTime = myClock();
		if (gd->verbose)
			fprintf(stdout, "Time execution: %.3f\n", (float)(endTime-initTime)/1000);

		if (gd->pass == FULL_ALIGNMENT) {
			initTime = myClock();
			struct FastaPairwiseAlignment pairAlign = getFastaAlignment(gd);
			endTime = myClock();
			if (gd->verbose)
				fprintf(stdout, "Time execution of backwards stage: %.3f\n", (float)(endTime-initTime)/1000);
//			fprintf(stdout, "%s\n", pairAlign.sequence[0]->name);
//			fprintf(stdout, "%s\n", pairAlign.alignment[0]);
//			fprintf(stdout, "%s\n", pairAlign.sequence[1]->name);
//			fprintf(stdout, "%s\n", pairAlign.alignment[1]);
		}

	/* RELEASE MEMORY */
	freeStack(&gd->stackOfJobs);
	freeJobTableStruct(&gd->jobTable, gd->pass);
	return gd->bestScore;
}


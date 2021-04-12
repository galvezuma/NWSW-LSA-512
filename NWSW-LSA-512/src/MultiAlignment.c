/*
 * MultiAlignment.c
 *
 *  Created on: Apr 7, 2021
 *      Author: galvez
 */

#include <stdlib.h>
#include <string.h>

#include "SingleAlignment.h"
#include "MultiAlignment.h"
#include "Utilities.h"
#include "NWSW-LSA-512.h"
#include "GlobalData.h"
#include "ScoreMatrix.h"
#include "Fasta.h"
#include "Worker.h"
#include "nj/mainNJ.h"

struct AlignerParam {
	int x, y; // Positions to set in the score array.
	struct GlobalData globalData;
};

// Global data accesible by all threadAligner
// This is used only for synchronization purposes
// There is a single instance of this structure per multiple alignment task (at all)
//
struct AlignerSynchronization {
	// Mutex to access to this structure in mutual exclusion
	pthread_mutex_t alignerSynchronizationAccess_mutex;
	pthread_cond_t alignmentAvailable_condition;
	pthread_cond_t stackNotFull_condition;
	// Data to be accessed in mutual exclusion
	volatile unsigned char noMoreAlignments;
	int numAligners; // This is a copy of the parallel user parameter
	int numSequences; // This is to access exact positions in the square array 'scores'
	int scoreSize;
	int * scores;
	struct {
		struct AlignerParam volatile * ptrArray; // Maximum size is the number of parallel aligners
		int volatile top;
	} stack;
};

void initializeAlignerSynchronization(struct AlignerSynchronization *synchro, struct UserParameters *ptrUserParams, int numSequences);
void freeAlignerSynchronization(struct AlignerSynchronization *as);
void pushParallelPairwise(struct AlignerSynchronization *as, int i, int j, struct GlobalData *gd);
void setNoMoreAlignments(struct AlignerSynchronization *as);
void * threadAligner(void * arg);

/********************************************************/
/* Main Function to calculate Multi Pairwise alignment  */
/* and the Newick tree                                  */
/********************************************************/
//
int multiplePairwise(struct UserParameters *ptrUserParams) {
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
	// Checks that in multi-pairwise the horizontal and vertical alphabets are used interchangeably and, therefore,
	// they should be the same
	if (strcmp(globalData.scoreMatrix.horizontalAlphabet, globalData.scoreMatrix.verticalAlphabet))
		fatalError2("In multi-pairwise alignments, both horizontal and vertical alphabets in matrix should be the same: %s versus %s",
					globalData.scoreMatrix.horizontalAlphabet,
					globalData.scoreMatrix.verticalAlphabet);

	// Reads a file with many fasta sequences
	struct NodeListSequence *ptrNodesFasta = readMultipleFastaFile(ptrUserParams->multifastaFilename);
	if (ptrNodesFasta == NULL || ptrNodesFasta->ptrNext == NULL)
		fatalError1("Multifasta file %s must have two sequences at least.\n", ptrUserParams->multifastaFilename);

	// VERBOSE
	if (ptrUserParams->verbose) {
		fprintf(stdout, "nwsw_lsa v%.1f. Universidad de Málaga. Spain.\n", VERSION);
		displayUserParameters(ptrUserParams, stdout);
		// Display sequences in multifasta file
		fprintf(stdout, "****** INPUT SEQUENCES ******\n");
		int i = 0;
		for (struct NodeListSequence *ptr = ptrNodesFasta; ptr != NULL; ptr = ptr->ptrNext)
			fprintf(stdout, "Sequence %d (%ld nt): %s\n", i++, strlen(ptr->sequence.data), ptr->sequence.name);
		fprintf(stdout, "****** SCORE MATRIX ******\n");
		// Display matrix content
		displayScoreMatrix(stdout, &(globalData.scoreMatrix));
	}

	/* INITIALIZATION */
	// Checks the letters in the sequences and put them in Uppercase
	// Check is made only with horizontal codification
	int i = 0;
	for (struct NodeListSequence *ptr = ptrNodesFasta; ptr != NULL; ptr = ptr->ptrNext, i++) {
		int numInvalidLettersInSequence;
		// Checks the Query
		numInvalidLettersInSequence = toUpperCodeAndCheck(	&ptr->sequence,
													globalData.scoreMatrix.horizontalAlphabet,
													globalData.scoreMatrix.horizontalCodification
												 );
		if (numInvalidLettersInSequence != 0) fprintf(stderr, "Sequence %d has %d letters not found in the alphabet.\n", i, numInvalidLettersInSequence);
	}
	const int numSeq = i;
	// Object to perform synchronization and initializes it
	struct AlignerSynchronization as;
	initializeAlignerSynchronization(&as, ptrUserParams, numSeq);

	/* CREATION OF THREADS */ /* Initially, they will start waiting */
	pthread_t threads[as.numAligners];
	for(int i=0; i<as.numAligners; i++){
		int error = pthread_create(&(threads[i]), NULL, threadAligner, (void *) &as);
		if (error) fatalError1("Unable to create the aligner num. %d\n", i);
	}


	/*****************************************/
	/* MAIN CALL TO MULTI-PAIRWISE ALIGNMENT */
	/*****************************************/
	i = 0;
	for (struct NodeListSequence *ptrI = ptrNodesFasta; ptrI->ptrNext != NULL; ptrI = ptrI->ptrNext, i++) {
		int j = i;
		for (struct NodeListSequence *ptrJ = ptrI->ptrNext; ptrJ != NULL; ptrJ = ptrJ->ptrNext, j++) {
			globalData.query = ptrI->sequence;
			globalData.subject = ptrJ->sequence;
			/* MAIN CALL TO EXECUTE PAIRWISE ALIGNMENT IN PARALLEL: IT PUTS THE TASK IN A STACK AND RETURNS */
			pushParallelPairwise(&as, i, j, &globalData);
		}
	}
	setNoMoreAlignments(&as);

	// WAIT FOR THREADS TO FINISH TO FULFILL THE SCORE ARRAY
	for(int i=0; i<as.numAligners; i++){
		int error = pthread_join(threads[i], NULL);
		if (error) fatalError1("Unable to join to the aligner num. %d\n", i);
	}

	// VERBOSE: Display diagonal matrix
	if (ptrUserParams->verbose) {
		for (int j = 0; j < as.numSequences - 1; j++)
			fprintf(stdout, "\t(%d)", j);
		fprintf(stdout, "\n");
		for (int j = 0; j < as.numSequences - 1; j++) {
			fprintf(stdout, "(%d)", j+1);
			for (int i = 0; i <= j ; i++)
				fprintf(stdout, "\t%d", as.scores[i * (as.numSequences - 1) + j]);
			fprintf(stdout, "\n");
		}
	}

	if (ptrUserParams->tree == NEIGHBOUR_JOIN) {
		char * arraySpeciesName[as.numSequences];
		int i = 0;
		for (struct NodeListSequence *ptr = ptrNodesFasta; ptr != NULL; ptr = ptr->ptrNext, i++)
			arraySpeciesName[i] = ptr->sequence.name;
		// We convert scores into positive distances
		int maxScore = -INFINITE;
		// Calculate maximum value
		for (int j = 0; j < as.numSequences - 1; j++)
			for (int i = 0; i <= j ; i++)
				maxScore = max(maxScore, as.scores[i * (as.numSequences - 1) + j]);
		// Replace scores by positive distances
		for (int j = 0; j < as.numSequences - 1; j++)
			for (int i = 0; i <= j ; i++)
				as.scores[j * (as.numSequences - 1) + i] = as.scores[i * (as.numSequences - 1) + j] = maxScore - as.scores[i * (as.numSequences - 1) + j];
		//
		mainNJ(as.numSequences, as.scores, arraySpeciesName, ptrUserParams->verbose);
	}

	/* RELEASE MEMORY */
	freeAlignerSynchronization(&as); // Also destroys the synchronizationElements
	//
	freeNodeListSequenceStruct(&ptrNodesFasta);
	freeScoreMatrixStruct(&(globalData.scoreMatrix));
	return(EXIT_SUCCESS);

}

/* Initialize mutex and condition to synchronize */
/* Initialize stack and score */
void initializeAlignerSynchronization(struct AlignerSynchronization *as, struct UserParameters *ptrUserParams, int numSequences) {
	// PREPARE SYNCHRONIZATION AMONG ALIGNERS
	int error;
	error = pthread_mutex_init(&as->alignerSynchronizationAccess_mutex, NULL);
	if (error) fatalError0("Unable to initialize the aligner's mutex.\n");
	error = pthread_cond_init(&as->alignmentAvailable_condition, NULL);
	if (error) fatalError0("Unable to initialize the alignment's synchronization condition.\n");
	error = pthread_cond_init(&as->stackNotFull_condition, NULL);
	if (error) fatalError0("Unable to initialize the aligner's synchronization condition.\n");
	// Initialize the number of parallel aligners
	as->numAligners = ptrUserParams->parallel;
	as->numSequences = numSequences;
	// PREPARE STRUCTURE TO STORE SCORES
	as->scoreSize = numSequences-1;
	as->scores = (int *) internalMalloc(sizeof(int)*(numSequences-1)*(numSequences-1));
	as->stack.ptrArray = (struct AlignerParam *) internalMalloc(sizeof(struct AlignerParam) * as->numAligners);
	as->stack.top = 0;
	as->noMoreAlignments = 0;
}

// Frees the memory allocated and destroy the synchronization elements
void freeAlignerSynchronization(struct AlignerSynchronization *as) {
	// Destroy synchronization elements
	int error;
	error = pthread_cond_destroy(&as->stackNotFull_condition);
	if (error) fatalError0("Unable to destroy the alignment's synchronization condition.\n");
	error = pthread_cond_destroy(&as->alignmentAvailable_condition);
	if (error) fatalError0("Unable to destroy the aligner's synchronization condition.\n");
	error = pthread_mutex_destroy(&as->alignerSynchronizationAccess_mutex);
	if (error) fatalError0("Unable to destroy the aligner's mutex.\n");
	// Frees allocated memory
	internalFree((void **) &as->scores);
	internalFree((void **) &as->stack.ptrArray);
}

// Puts a new alinment in a stack in order to be executed when an aligner becomes idle
// GlobalData is passed by value so each pairwise use a different GlobalData object
//
void pushParallelPairwise(struct AlignerSynchronization *as, int i, int j, struct GlobalData *gd) {
    pthread_mutex_lock(&as->alignerSynchronizationAccess_mutex);
        while(as->stack.top >= as->numAligners) // While stack is full
            pthread_cond_wait(&as->stackNotFull_condition, &as->alignerSynchronizationAccess_mutex);
        // Puts a new set of params on top of the stack in order to be taken by a threadAligner
        // GlobalData is passed by value so each pairwise use a different GlobalData object
        as->stack.ptrArray[as->stack.top].x = i;
        as->stack.ptrArray[as->stack.top].y = j;
        as->stack.ptrArray[as->stack.top].globalData = *gd;
        as->stack.top++;
		pthread_cond_broadcast(&as->alignmentAvailable_condition);
	pthread_mutex_unlock(&as->alignerSynchronizationAccess_mutex);
}

// Notifies that no more alignments will be added to the stack
//
void setNoMoreAlignments(struct AlignerSynchronization *as) {
	pthread_mutex_lock(&as->alignerSynchronizationAccess_mutex);
		as->noMoreAlignments = 1;
		pthread_cond_broadcast(&as->alignmentAvailable_condition);
	pthread_mutex_unlock(&as->alignerSynchronizationAccess_mutex);
}

/* MAIN Thread of ALIGNERS */
// It pops alignments until the stack becomes empty and no more alignments will be pushed
//
void * threadAligner(void * arg) {
	struct AlignerSynchronization *as = (struct AlignerSynchronization *) arg;
	int finished = 0;
	while (! finished) {
		pthread_mutex_lock(&as->alignerSynchronizationAccess_mutex);
			while(as->stack.top == 0 && ! as->noMoreAlignments)  { // While stack is empty and more alignments will be included
				pthread_cond_wait(&as->alignmentAvailable_condition, &as->alignerSynchronizationAccess_mutex);
			}
			if (as->stack.top > 0) {
					as->stack.top--;
					int i = as->stack.ptrArray[as->stack.top].x;
					int j = as->stack.ptrArray[as->stack.top].y;
					struct GlobalData gd = as->stack.ptrArray[as->stack.top].globalData; // Copy by value
					pthread_cond_broadcast(&as->stackNotFull_condition);
				pthread_mutex_unlock(&as->alignerSynchronizationAccess_mutex);
				int score = executePairwise(&gd); // Executed out of locks
				pthread_mutex_lock(&as->alignerSynchronizationAccess_mutex);
				as->scores[i * (as->numSequences - 1) + j] = as->scores[j * (as->numSequences - 1) + i] = score;
				fprintf(stdout, "The score %dx%d is: %d\n", i, j, score);
			} else {
				finished = 1;
			}
		pthread_mutex_unlock(&as->alignerSynchronizationAccess_mutex);
	}
	pthread_exit(NULL);
}

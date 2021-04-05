/*
 * GlobalData.h
 *
 *  Created on: Mar 10, 2021
 *      Author: galvez
 */

#ifndef GLOBALDATA_H_
#define GLOBALDATA_H_

#include <stdint.h>
#include <pthread.h>

struct GlobalData; // This pre-declaration is required to avoid cyclic declarations

#include "Sequence.h"
#include "ScoreMatrix.h"
#include "JobTable.h"
#include "Stack.h"
#include "NWSW-LSA-512.h"

/* Global data accesible to any thread.
 * This tries to represent an Object.
 */
// Note that unsigned or uint32_t should not be used. This is because mixing
// int and unsigned int ¡converts negative int into positive!
struct GlobalData {

	// CLI parameters
	enum Pass pass;
	enum Tree tree;
	enum Algorithm algorithm;
//	enum Info info;
	int32_t insert;
	int32_t delete;
	int32_t matchReplace __attribute__((aligned (4)));
	int32_t gapExtend __attribute__((aligned (4)));
	int16_t threads;
	unsigned char verbose;

	// Loaded files
	struct Sequence query;
	struct Sequence subject;
	struct ScoreMatrix scoreMatrix;


	// Job table
	struct JobTable jobTable;
	// Stack of jobs. It must be accessed in Mutual Exclusion
	struct Stack stackOfJobs;

	// Structures for mutual exclusion
	pthread_cond_t jobAvailable_condition;
	pthread_mutex_t globalDataAccess_mutex;
	volatile unsigned char jobTableFulfilled;
	pthread_mutex_t verboseStdOut_mutex;

	// Vectorization availability
	int lengthVector;

	// Results
	volatile int bestScore;
	struct Job * volatile bestJob;
};

void copyUserParameters(struct GlobalData *gd, struct UserParameters *up);
struct Job * volatile waitForAvailableJob(struct GlobalData *gd);
void informFinishedJob(struct GlobalData *gd, unsigned x, unsigned y, struct Node *resultingFragmentX, struct Node *resultingFragmentY);
void saveMyBestScore(struct GlobalData *gd, int score, struct Job *bestJob);
void checkSupport(struct GlobalData *gd);

// Functions to block and unblock while displaying a message in verbose mode.
void block(struct GlobalData *gd);
void unblock(struct GlobalData *gd);

#endif /* GLOBALDATA_H_ */

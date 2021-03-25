/*
 * JobTable.c
 *
 *  Created on: Mar 12, 2021
 *      Author: galvez
 */

#include <string.h>
#include <limits.h>

#include "Utilities.h"
#include "JobTable.h"
#include "GlobalData.h"
#include "NWSW-LSA-512.h"

void allocateFragmentsAndJobs(struct GlobalData *gd);
void initializeFragments(struct GlobalData *gd);

/* Initializes the Job table stores inside a GlobalData struct */
/* Initialization depends on the algorithm (NW or SW) and if the user wants only the score or the alignment too */
void createJobTable(struct GlobalData *gd) {
	gd->jobTable.fragmentSize_X = DEFAULT_FRAGMENTSIZE_X;
	gd->jobTable.fragmentSize_Y = DEFAULT_FRAGMENTSIZE_Y;
	gd->jobTable.numFragments_X = ceilDivision(strlen(gd->query.data), gd->jobTable.fragmentSize_X);
	gd->jobTable.numFragments_Y = ceilDivision(strlen(gd->subject.data), gd->jobTable.fragmentSize_Y);
	allocateFragmentsAndJobs(gd);
	initializeFragments(gd);
}

/* Initializes the rows and columns to store the LSA matrix */
/* To save space, if only the score is needed then not all the fragments are needed because there will not be a backward stage */
/* In case of ONLY_SCORE, only the fragments needed to continue the calculus are needed */
/* Note: all the memory is allocated at the beginning to avoid an Overflow after a long time of program execution */
void allocateFragmentsAndJobs(struct GlobalData *gd) {
	size_t querySize = strlen(gd->query.data);
	size_t subjectSize = strlen(gd->subject.data);
	// THIS ALLOCATES THE FRAGMENTS saving a lot of space in case of ONLY_SCORE.
	if (gd->pass == ONLY_SCORE) { // A single row and column are created
		gd->jobTable.rows = (struct Node **) internalMalloc(sizeof(struct Node *)); // a single pointer because we have a single row
		gd->jobTable.rows[0] = (struct Node*) internalMalloc(sizeof(struct Node) * (1 + querySize));
		//
		gd->jobTable.columns = (struct Node **) internalMalloc(sizeof(struct Node *)); // a single pointer because we have a single column
		gd->jobTable.columns[0] = (struct Node*) internalMalloc(sizeof(struct Node) * (1 + subjectSize));
	} else { // A grid of rows and columns is created
		gd->jobTable.rows = (struct Node **) internalMalloc(sizeof(struct Node *) * gd->jobTable.numFragments_Y); // a single pointer because we have a single row
		for (int i=0; i<gd->jobTable.numFragments_Y; i++)
			gd->jobTable.rows[i] = (struct Node*) internalMalloc(sizeof(struct Node) * (1 + querySize));
		//
		gd->jobTable.columns = (struct Node **) internalMalloc(sizeof(struct Node *) * gd->jobTable.numFragments_X); // a single pointer because we have a single column
		for (int i=0; i<gd->jobTable.numFragments_X; i++)
			gd->jobTable.columns[i] = (struct Node*) internalMalloc(sizeof(struct Node) * (1 + subjectSize));
	}
	//
	// THIS ALLOCATES ALL THE JOBS AND INITIALIZES THEM
	// Jobs are stored linearly although they are managed as a grid
	// Firstly is allocated
	gd->jobTable.jobs = (struct Job *) internalMalloc(sizeof(struct Job) * gd->jobTable.numFragments_X * gd->jobTable.numFragments_Y);
	// Secondly it is initialized, for each Job stored
	for (int y=0; y<gd->jobTable.numFragments_Y; y++)
		for (int x=0; x<gd->jobTable.numFragments_X; x++) {
			// For each Job, we take it ...
			struct Job *job = getJob(&gd->jobTable, x, y);
			// ... and initializes it.
			job->x = x;
			job->y = y;
			// We set its fragments size. The last fragment has a different size
			job->realSize_X = (x!=gd->jobTable.numFragments_X-1)? gd->jobTable.fragmentSize_X : restOfLengthSequence(querySize, gd->jobTable.fragmentSize_X);
			job->realSize_Y = (y!=gd->jobTable.numFragments_Y-1)? gd->jobTable.fragmentSize_Y : restOfLengthSequence(subjectSize, gd->jobTable.fragmentSize_Y);
			// We make the fragments to point to the previously allocated memory, x.e., rows and columns
			// job->ptrRow = (gd->pass == ONLY_SCORE && y > 0)? NULL : gd->jobTable.rows[y] + x;
			// job->ptrColumn = (gd->pass == ONLY_SCORE && x > 0)? NULL : gd->jobTable.columns[x] + y;
			// A job contains two pointers. Each of them has a different from NULL value when
			// it points to an already calculated row or column.
			job->ptrRow = (y > 0)? NULL : gd->jobTable.rows[y] + x * gd->jobTable.fragmentSize_X;
			job->ptrColumn = (x > 0)? NULL : gd->jobTable.columns[x] + y * gd->jobTable.fragmentSize_Y;
		}
}

/* Frees the memory previously allocated */
/* Note: the pass parameter must be sent to know if there is a complete grid to be freed or a single row and columns */
void freeJobTableStruct(struct JobTable * jobTable, enum Pass pass) {
	internalFree((void **) &jobTable->jobs);
	if (pass == ONLY_SCORE) { // A single row and column should be freed
		internalFree((void **) &jobTable->rows[0]);
		internalFree((void **)&jobTable->rows);
		internalFree((void **)&jobTable->columns[0]);
		internalFree((void **)&jobTable->columns);
	} else { // A grid of rows and columns is freed
		for (int i=0; i<jobTable->numFragments_Y; i++)
			internalFree((void **)&jobTable->rows[i]);
		internalFree((void **)&jobTable->rows);
		for (int i=0; i<jobTable->numFragments_X; i++)
			internalFree((void **)&jobTable->columns[i]);
		internalFree((void **)&jobTable->columns);
	}
}

void initializeStructNode(struct Node *node, int t, int u, int s) {
	node->t = t;
	node->u = u;
	node->s = s;
}

/* Initialize the fragments of the first row and column depending on NW or SW */
void initializeFragments(struct GlobalData *gd) {
		// Initializes first row
		size_t querySize = strlen(gd->query.data);
		initializeStructNode(&gd->jobTable.rows[0][0], 0, 0, 0); // Initialize first cell of the row
		for(int i=1; i<querySize+1; i++) // Initialize the rest of the cells of the row
			initializeStructNode(&gd->jobTable.rows[0][i],
								-(gd->delete + i*gd->gapExtend),
								-INFINITE,
								(gd->algorithm == NEEDLEMAN_WUNSCH)? -(gd->delete + i*gd->gapExtend) : 0
								);
		// Initializes first column
		size_t subjectSize = strlen(gd->subject.data);
		initializeStructNode(&gd->jobTable.columns[0][0], 0, 0, 0); // Initialize first cell of the column
		for(int j=1; j<subjectSize+1; j++) // Initialize the rest of the cells of the column
			initializeStructNode(&gd->jobTable.columns[0][j],
								-INFINITE,
								-(gd->insert + j*gd->gapExtend),
								(gd->algorithm == NEEDLEMAN_WUNSCH)? -(gd->insert + j*gd->gapExtend) : 0
								);
}


struct Job *getJob(struct JobTable *jobTable, int posX, int posY) {
	return jobTable->jobs + posX * jobTable->numFragments_X + posY;
}

/* DEBUG FUNCTIONS */
// Display the indicated node
void displayNode(struct Node *node) {
	printf("[t=%d,u=%d,s=%d]-", node->t, node->u, node->s);
}

// Displays the indicated Job
void displayJobFromJobTable(struct JobTable *jobTable, int x, int y) {
	if (x>jobTable->numFragments_X || y>jobTable->numFragments_Y) {
		printf("Inexistent node\n");
		return;
	}
	displayJob(getJob(jobTable,x,y));
}
void displayJob(struct Job *job) {
	printf("(%d,%d)\n", job->x, job->y);
	printf("Fragment X size: %d\n", job->realSize_X);
	printf("Fragment Y size: %d\n", job->realSize_Y);
	if (job->ptrRow != NULL) for(int i=0; i<=job->realSize_X; i++)
								displayNode(job->ptrRow+i);
	else printf("NULL");
	printf("\n");
	if (job->ptrColumn != NULL) for(int j=0; j<=job->realSize_Y; j++)
									displayNode(job->ptrColumn+j);
	else printf("NULL");
	printf("\n");
}

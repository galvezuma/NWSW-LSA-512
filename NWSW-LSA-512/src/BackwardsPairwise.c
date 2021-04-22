/*
 * BackwardsPairwise.c
 *
 *  Created on: Apr 12, 2021
 *      Author: galvez
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <assert.h>

#include "BackwardsPairwise.h"
#include "Utilities.h"

// This internal structure is used to return the sub-alignment from a single Job
// The complete alignment is made by composing several of these structures
struct FastaJobAlignment {
	char * fastaSubAlign[2];
	int upperLeft_x;
	int upperLeft_y;
	bool swFinished; // To register if the SW backwards has finished (reached a 0 score)
};

struct FastaJobAlignment fulfillJob(struct GlobalData * gd, struct Job *job, int starting_x, int starting_y, bool initialJob);
void inline freeFastaJobAlignmentStruct(struct FastaJobAlignment *subAlign);
void reconcatStr(char **ptrDst, char *src);
void reverseStr(char *str);

/*  Returns a Fasta alignment in a 'struct FastaPairwiseAlignment'.
 *  The starting point if a JobTable already fulfilled.
 */
struct FastaPairwiseAlignment getFastaAlignment(struct GlobalData * gd) {
	struct FastaPairwiseAlignment ret;
	if (gd->pass != FULL_ALIGNMENT) fatalError1("Final alignment can not be obtained because %s has not been specified as option.\n", enumPassToString(FULL_ALIGNMENT));
	if (! gd->jobTableFulfilled) fatalError0("Filling of Dynamic Programming Table has not been carried out.\n");

	// Initialize the result
	ret.sequence[0] = &gd->query;
	ret.sequence[1] = &gd->subject;
	// Initialize with the first slice of the alignment, taken from the lower-right job of the Job table.
	int currentFragment_X, currentFragment_Y;
	int starting_x, starting_y;
	struct Job * lowerRightJob;
	if (gd->algorithm == NEEDLEMAN_WUNSCH) {
		currentFragment_X = gd->jobTable.numFragments_X - 1;
		currentFragment_Y = gd->jobTable.numFragments_Y - 1;
		lowerRightJob = getJob(&gd->jobTable, currentFragment_X, currentFragment_Y);
	} else { // gd->algorithm == SMITH_WATERMAN
		lowerRightJob = gd->bestJob;
		currentFragment_X = lowerRightJob->x;
		currentFragment_Y = lowerRightJob->y;
	}
	starting_x = lowerRightJob->realSize_X;
	starting_y = lowerRightJob->realSize_Y;
	/********************************/
	/* Let's do the initial fulfill */
	/********************************/
	//
	// Initially, the alignments are obtained in REVERSE ORDER
	//
	struct FastaJobAlignment subAlign = fulfillJob(gd, lowerRightJob, starting_x, starting_y, true);
	ret.alignment[0] = internalMalloc(strlen(ret.sequence[0]->data) + strlen(ret.sequence[1]->data) + 1); // Maximum possible length
	ret.alignment[1] = internalMalloc(strlen(ret.sequence[0]->data) + strlen(ret.sequence[1]->data) + 1); // Maximum possible length
	strcpy(ret.alignment[0], subAlign.fastaSubAlign[0]);
	strcpy(ret.alignment[1], subAlign.fastaSubAlign[1]);
	int currentLength = strlen(ret.alignment[0]); // Lengths of ret.alignment[0] and ret.alignment[1] are the same.

	freeFastaJobAlignmentStruct(&subAlign);


	// Once initialized with the lower-right Job, let's continue with the rest of jobs
	// depending on where the previous Job (backwards stage) has finished.
	do {
		// Let's locate the next Job to calculate
		if (subAlign.upperLeft_x == 0) currentFragment_X--;
		if (subAlign.upperLeft_y == 0) currentFragment_Y--;
		// We check if the loop must end
		if (currentFragment_X < 0 || currentFragment_Y < 0) break; // We have reached one wall of the Job table.
		if (subAlign.swFinished) break; // We have reached the beginning of a local alignment.

		struct Job * currentJob = getJob(&gd->jobTable, currentFragment_X, currentFragment_Y);
		// Let's locate the position from where to start
		starting_x = (subAlign.upperLeft_x == 0)? currentJob->realSize_X : subAlign.upperLeft_x;
		starting_y = (subAlign.upperLeft_y == 0)? currentJob->realSize_Y: subAlign.upperLeft_y;
		/************************/
		/* Let's do the fulfill */
		/************************/
		subAlign = fulfillJob(gd, currentJob, starting_x, starting_y, false);
		strcpy(ret.alignment[0] + currentLength, subAlign.fastaSubAlign[0]);
		strcpy(ret.alignment[1] + currentLength, subAlign.fastaSubAlign[1]);
		currentLength += strlen(subAlign.fastaSubAlign[0]); // Lengths of subAlign.fastaSubAlign[0] and subAlign.fastaSubAlign[1] are the same.
		freeFastaJobAlignmentStruct(&subAlign);
	} while(true);

	// Finally, a likely initial gap in the alignment must be appended.
	// We must put it in reverse order, like the other slices/subalignments.
	// If the algorithm is Smith-Waterman, this process is not done.
	if (gd->algorithm == NEEDLEMAN_WUNSCH) {
		if (currentFragment_X < 0) { // Check for going up up to the top of the Job Table
			int numChars = subAlign.upperLeft_y + (currentFragment_Y) * gd->jobTable.fragmentSize_Y;
			if (numChars > 0) { // There is a gap at the beginning
				char sliceQuery[numChars + 1];
				char sliceSubject[numChars + 1];
				for(int i=0; i<numChars; i++) {
					sliceQuery[i] = '-';
					sliceSubject[i] = gd->subject.data[numChars - i - 1];
				}
				sliceQuery[numChars] = sliceSubject[numChars] = '\0';
				strcpy(ret.alignment[0] + currentLength, sliceQuery);
				strcpy(ret.alignment[1] + currentLength, sliceSubject);
				currentLength += numChars;
			}
		} else if (currentFragment_Y < 0) { // Check for going up up to the top of the Job Table
			int numChars = subAlign.upperLeft_x + (currentFragment_X) * gd->jobTable.fragmentSize_X;
			if (numChars > 0) { // There is a gap at the beginning
				char sliceQuery[numChars + 1];
				char sliceSubject[numChars + 1];
				for(int i=0; i<numChars; i++) {
					sliceQuery[i] = gd->query.data[numChars - i - 1];
					sliceSubject[i] = '-';
				}
				sliceQuery[numChars] = sliceSubject[numChars] = '\0';
				strcpy(ret.alignment[0] + currentLength, sliceQuery);
				strcpy(ret.alignment[1] + currentLength, sliceSubject);
				currentLength += numChars;
			}
		}
	}

	// Alignments are in reverse order and must be reoriented
	reverseStr(ret.alignment[0]);
	reverseStr(ret.alignment[1]);
	ret.alignment[0] = realloc(ret.alignment[0], currentLength + 1);
	ret.alignment[1] = realloc(ret.alignment[1], currentLength + 1);
	return ret;
}

/* Populates a Job completely (every cell) and calculates the alignment of such a Job.
 * Actually, it fulfills only the amount of rows and columns indicated by 'starting_x' and 'starting_y'.
 * These two parameters represent the lower-right cell from where the backward alignment has to start.
 * It does the calculus using NON-VECTOR operations
 */
struct FastaJobAlignment fulfillJob(struct GlobalData * gd, struct Job *job, int const starting_x, int const starting_y, bool initialJob) {
	// FIRST PASS: Fulfill the Job
		if (gd->verbose) { fprintf(stdout, "Executing non-vectorial full job (%d,%d) of size (%dx%d).\n", job->x, job->y, job->realSize_X, job->realSize_Y); }
		// Calculus of the initial positions of the sequences to compare
		unsigned const queryIdx = job->x * gd->jobTable.fragmentSize_X;
		unsigned const subjectIdx = job->y * gd->jobTable.fragmentSize_Y;
		// Create the array where the nodes will be calculated
		// Equivalent to struct Node array[1 + starting_x][1 + starting_y];
		struct Node * array = (struct Node *) internalMalloc((1 + starting_x) * (1 + starting_y) * sizeof(struct Node));
#define array_1(i) (array + (i) * (1 + starting_y))
#define array_2(i, j) (array + (i) * (1 + starting_y) + (j))
		// Copy the row ...
		for(int i=0; i <= starting_x; i++)
			*array_2(i, 0) = job->ptrRow[i]; // array[i][0] = job->ptrRow[i];
		// ...and copy the column
		for(int j=0; j <= starting_y; j++)
			*array_2(0, j) = job->ptrColumn[j]; // array[0][j] = job->ptrColumn[j];
		//
		/* This fulfills entirely the matrix of the Job with Nodes */
		//
		int bestX = starting_x;
		int bestY = starting_y;
		int bestScore = -1;
		for(int i=1; i <= starting_x; i++) {
			struct Node *prev = array_1(i-1); // ... = array[i-1];
			struct Node *current = array_1(i); // ... = array[i];
			for(int j=1; j <= starting_y; j++) {
				current[j].t = max(prev[j].t, prev[j].s - gd->delete) - gd->gapExtend;
				current[j].u = max(current[j-1].u, current[j-1].s - gd->insert) - gd->gapExtend;
				int cellScore = gd->scoreMatrix.matrix[gd->query.dataCoded[queryIdx+i-1]][gd->subject.dataCoded[subjectIdx+j-1]];
				current[j].s = max( max(current[j].t, current[j].u),  prev[j-1].s + cellScore);
				if (gd->algorithm == SMITH_WATERMAN) {
					current[j].s = max(0, current[j].s);
					if (initialJob && current[j].s > bestScore) {
						bestScore = current[j].s;
						bestX = i;
						bestY = j;
					}
				}
			}
		}
		if (gd->algorithm == SMITH_WATERMAN && initialJob) assert(bestScore == gd->bestScore);

	// SECOND PASS: Let's obtain the alignment in REVERSE ORDER
		struct FastaJobAlignment ret;
		ret.fastaSubAlign[0] = (char *) internalMalloc(starting_x + starting_y + 1);
		ret.fastaSubAlign[1] = (char *) internalMalloc(starting_x + starting_y + 1);
		int posSubAlign = 0;
		int i = (gd->algorithm == SMITH_WATERMAN)? bestX : starting_x;
		int j = (gd->algorithm == SMITH_WATERMAN)? bestY : starting_y;
		ret.swFinished = false;
		while(i > 0 && j > 0 && ! ret.swFinished) {
			if (gd->algorithm == SMITH_WATERMAN && array_2(i, j)->s == 0) ret.swFinished = true;
			if (array_2(i, j)->s == array_2(i, j)->t) { // Comes from the left // if (array[i][j].s == array[i][j].t)
				ret.fastaSubAlign[0][posSubAlign] = gd->query.data[queryIdx+i-1];
				ret.fastaSubAlign[1][posSubAlign] = '-';
				i--;
			} else if (array_2(i, j)->s == array_2(i, j)->u) { // Comes from the top // if (array[i][j].s == array[i][j].u)
				ret.fastaSubAlign[0][posSubAlign] = '-';
				ret.fastaSubAlign[1][posSubAlign] = gd->subject.data[subjectIdx+j-1];
				j--;
			} else { // Comes from the corner
				ret.fastaSubAlign[0][posSubAlign] = gd->query.data[queryIdx+i-1];
				ret.fastaSubAlign[1][posSubAlign] = gd->subject.data[subjectIdx+j-1];
				i--;
				j--;
			}
			posSubAlign++;
		}
#undef array_1
#undef array_2
	internalFree((void **) &array);
	ret.fastaSubAlign[0][posSubAlign] = ret.fastaSubAlign[1][posSubAlign] = '\0';
	ret.upperLeft_x = i;
	ret.upperLeft_y = j;
	return ret;
}

void inline freeFastaJobAlignmentStruct(struct FastaJobAlignment *subAlign) {
		internalFree((void **) &subAlign->fastaSubAlign[0]);
		internalFree((void **) &subAlign->fastaSubAlign[1]);
}

void freeFastaPairwiseAlignmentStruct(struct FastaPairwiseAlignment *pairAlign) {
	if (pairAlign->alignment[0] != NULL) internalFree((void **) &pairAlign->alignment[0]);
	if (pairAlign->alignment[1] != NULL) internalFree((void **) &pairAlign->alignment[1]);
}

// UNUSED
// Concat a string to another string by allocating memory.
void reconcatStr(char **ptrDst, char *src) {
	char *ret;
	ret = (char *) realloc(*ptrDst, strlen(*ptrDst)+strlen(src)+1);
	if (ret == *ptrDst) { // The realloc is successful and returns to the same memory address than dst
		strcat(*ptrDst, src);
	} else if (ret != NULL) { // The realloc is successful but returns a different address
		*ptrDst = ret;
		strcat(*ptrDst, src);
	} else { // The realloc is not successful
		// Let's try with an usual malloc
		ret = (char *) internalMalloc(strlen(*ptrDst)+strlen(src)+1);
		strcpy(ret, *ptrDst);
		internalFree((void **) ptrDst);
		*ptrDst = ret;
		strcat(*ptrDst, src);
	}
}

// Reverse a string by swapping its chars.
void reverseStr(char *str) {
  int last = strlen(str);
  for (int i=0; i < last/2; i++)  {
    char ch = str[i];
    str[i] = str[last - i - 1];
    str[last - i - 1] = ch;
  }
}

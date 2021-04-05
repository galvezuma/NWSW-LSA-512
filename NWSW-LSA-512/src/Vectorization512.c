/*
 * Vectorization512.c
 *
 *  Created on: Mar 25, 2021
 *      Author: galvez
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <immintrin.h>
#include "Utilities.h"
#include "Vectorization512.h"

#define lengthVector 16

struct Node_512 {
	__m512i t; // End in horizontal gap
	__m512i u; // End in vertical gap
	__m512i s; // Best value
};
// Vector constants
__m512i arribaMask_512;
__m512i gapExtend_512;
__m512i openGapDeletion_512;
__m512i openGapInsertion_512;
__m512i zeroes_512;

extern int processJob_CISC(struct GlobalData *gd, struct Job *job, struct Node *retFragX, struct Node *retFragY);
void initializeVectors_512(int j, struct Node *retFragX, struct GlobalData *gd, struct Job *job, struct Node_512 *arriba, struct Node_512 *izquierda, struct Node_512 *esquina, __m512i *deltaScore);
void calculateAndAdvanceTopLeftDiagonal_512(int j, struct Node *retFragX, int i, struct GlobalData *gd, struct Job *job, struct Node_512 *arriba, struct Node_512 *izquierda, struct Node_512 *esquina, struct Node_512 * resultado, __m512i *deltaScore, int *bestScore);
void calculateAndAdvanceBody_512(int j, struct Node *retFragX, int i, struct GlobalData *gd, struct Job *job, struct Node_512 *arriba, struct Node_512 *izquierda, struct Node_512 *esquina, struct Node_512 * resultado, __m512i *deltaScore, int *bestScore);
void calculateAndAdvanceBottomRightDiagonal_512(int j, struct Node *retFragX, int i, struct GlobalData *gd, struct Job *job, struct Node_512 *arriba, struct Node_512 *izquierda, struct Node_512 *esquina, struct Node_512 * resultado, struct Node *retFragY, __m512i *deltaScore, int *bestScore);
inline void calculate_512(struct GlobalData *gd, struct Node_512 *arriba, struct Node_512 *izquierda, struct Node_512 *esquina, struct Node_512 * resultado, __m512i *deltaScore, int *bestScore);
inline void setPos(struct Node_512 *ptrNodeDest_512, int i, struct Node *ptrNodeSrc);
inline void savePos(struct Node *ptrNodeDest, struct Node_512 *ptrNodeSrc_512, int i);

/* DEBUG FUNCTIONS */
void displayNode_512(struct Node_512 *node512);
void displayVector_512(__m512i *vector);

/* Calculates the  scores in a job */
/* It does the job using TRUE-VECTOR operations */
int processJob_512(struct GlobalData *gd, struct Job *job, struct Node *retFragX, struct Node *retFragY) {
	if (job->realSize_Y % lengthVector != 0 || job->realSize_X <= 2*lengthVector ) return processJob_CISC(gd, job, retFragX, retFragY);
	// Load vector constants
	// Prepare shuffle
	const int arribaShuffle_512[] __attribute__((aligned (64))) = {0,0,1,2,3,4,5,6,7,8,9,10,11,12,13,14};
	arribaMask_512 = _mm512_load_epi32 (arribaShuffle_512);
	// Prepare vectorized constants to extend gap, open Deletion and open insertion
	gapExtend_512 = _mm512_set1_epi32(gd->gapExtend);
	openGapDeletion_512 = _mm512_set1_epi32(gd->delete + gd->gapExtend);
	openGapInsertion_512 = _mm512_set1_epi32(gd->delete + gd->gapExtend);
	// Prepare vector of zeroes
	zeroes_512 = _mm512_setzero_epi32();
	//
	// Prepare the best score of all the jobs calculated by this thread in case of S/W
	int bestScore = -1;
	// Variable to store values of the score matrix
	__m512i deltaScore;
	// Nodes of vector to work with
	struct Node_512 arriba, izquierda, esquina, resultado;
	// Initializes the top right of the resulting job
	retFragY[0] = job->ptrRow[job->realSize_X];
	// Initializes retFragX. retFragX changes in each pass of vectorization
	memcpy(retFragX, job->ptrRow, sizeof(struct Node)*(job->realSize_X + 1));
	for(int j=0; j<job->realSize_Y - 1; j+=lengthVector) {
		initializeVectors_512(j, retFragX, gd, job, &arriba, &izquierda, &esquina, &deltaScore);
		int i=1;
		for(; i < lengthVector; i++) // Step 1
			calculateAndAdvanceTopLeftDiagonal_512(j, retFragX, i, gd, job, &arriba, &izquierda, &esquina, &resultado, &deltaScore, &bestScore);
		// Update progressively retFragX
		retFragX[0] = job->ptrColumn[j+lengthVector];
		for(; i < job->realSize_X; i++) // Main step
			calculateAndAdvanceBody_512(j, retFragX, i, gd, job, &arriba, &izquierda, &esquina, &resultado, &deltaScore, &bestScore);
		for(; i < job->realSize_X + lengthVector; i++) // Step 3
			calculateAndAdvanceBottomRightDiagonal_512(j, retFragX, i, gd, job, &arriba, &izquierda, &esquina, &resultado, retFragY, &deltaScore, &bestScore);
	}
	return (gd->algorithm == SMITH_WATERMAN)? bestScore : 0;
}

/* Initialization of vectors for an horizontal round of vectorization */
//
void initializeVectors_512(int j, struct Node *retFragX, struct GlobalData *gd, struct Job *job, struct Node_512 *arriba, struct Node_512 *izquierda, struct Node_512 *esquina, __m512i *deltaScore) {
	// Fill first node of esquina with retFragX[0]
	// The rest of nodes may be garbage
	esquina->t = _mm512_maskz_set1_epi32(0x0001, retFragX[0].t); // Superfluous
	esquina->u = _mm512_maskz_set1_epi32(0x0001, retFragX[0].u); // Superfluous
	esquina->s = _mm512_maskz_set1_epi32(0x0001, retFragX[0].s);
	// Fill first node of izquierda with job->ptrColumn[j+1]
	// The rest of nodes may be garbage
	izquierda->t = _mm512_maskz_set1_epi32(0x0001, job->ptrColumn[j+1].t);
	izquierda->u = _mm512_maskz_set1_epi32(0x0001, job->ptrColumn[j+1].u); // Superfluous
	izquierda->s = _mm512_maskz_set1_epi32(0x0001, job->ptrColumn[j+1].s);
	// Fill first node of arriba with retFragX[1]; and the second one with job->ptrColumn[j+1]
	// The rest of nodes may be garbage
	// First
	arriba->t = _mm512_maskz_set1_epi32(0x0001, retFragX[1].t); // Superfluous
	arriba->u = _mm512_maskz_set1_epi32(0x0001, retFragX[1].u);
	arriba->s = _mm512_maskz_set1_epi32(0x0001, retFragX[1].s);
	// Second
	arriba->t = _mm512_mask_set1_epi32 (arriba->t, 0x0002, job->ptrColumn[j+1].t);
	arriba->u = _mm512_mask_set1_epi32 (arriba->u, 0x0002, job->ptrColumn[j+1].u);
	arriba->s = _mm512_mask_set1_epi32 (arriba->s, 0x0002, job->ptrColumn[j+1].s);
	// Fill first node of deltaScore with data from the score matrix
	// The rest of nodes may be garbage
	// Calculus of the initial positions of the sequences to compare
	unsigned queryIdx = job->x * gd->jobTable.fragmentSize_X;
	unsigned subjectIdx = job->y * gd->jobTable.fragmentSize_Y + j;
	// We have to remember that the sequences begin at position 0
	int cellScore = gd->scoreMatrix.matrix[gd->query.dataCoded[queryIdx]][gd->subject.dataCoded[subjectIdx]];
	*deltaScore = _mm512_maskz_set1_epi32(0x0001, cellScore);
}

/* Calculus of a single vector in diagonal */
//
inline void calculate_512(struct GlobalData *gd, struct Node_512 *arriba, struct Node_512 *izquierda, struct Node_512 *esquina, struct Node_512 * resultado, __m512i *deltaScore, int *bestScore) {
	// Deletion (horizontal gap). Calculus of resultado->t
	__m512i aux1 = _mm512_sub_epi32 (izquierda->t, gapExtend_512);
	__m512i aux2 = _mm512_sub_epi32 (izquierda->s, openGapDeletion_512);
	resultado->t = _mm512_max_epi32 (aux1, aux2);
	// Insertion (vertical gap). Calculus of resultado->u
	aux1 = _mm512_sub_epi32 (arriba->u, gapExtend_512);
	aux2 = _mm512_sub_epi32 (arriba->s, openGapInsertion_512);
	resultado->u = _mm512_max_epi32 (aux1, aux2);
	// Continue in diagonal. Calculus of resultado->s
	aux1 = _mm512_add_epi32 (esquina->s, *deltaScore);
	aux1 = _mm512_max_epi32 (aux1, resultado->t);
	resultado->s = _mm512_max_epi32 (aux1, resultado->u);
	if (gd->algorithm == SMITH_WATERMAN) {
		resultado->s = _mm512_max_epi32 (zeroes_512, resultado->s);
		int aux = _mm512_reduce_max_epi32(resultado->s);
		if (aux > *bestScore) *bestScore = aux;
	}
}

void calculateAndAdvanceTopLeftDiagonal_512(int j, struct Node *retFragX, int i, struct GlobalData *gd, struct Job *job, struct Node_512 *arriba, struct Node_512 *izquierda, struct Node_512 *esquina, struct Node_512 * resultado, __m512i *deltaScore, int *bestScore) {
	calculate_512(gd, arriba, izquierda, esquina, resultado, deltaScore, bestScore);
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
		// Shift t, u and s
		arriba->t = _mm512_permutexvar_epi32 (arribaMask_512, resultado->t);
		arriba->u = _mm512_permutexvar_epi32 (arribaMask_512, resultado->u);
		arriba->s = _mm512_permutexvar_epi32 (arribaMask_512, resultado->s);
		setPos(arriba, 0, &(retFragX[i+1]));
		setPos(arriba, i+1, &job->ptrColumn[j+i+1]);
	// Recalculates deltaScore
		// Calculus of the initial positions of the sequences to compare
		unsigned queryIdx = job->x * gd->jobTable.fragmentSize_X + i;
		unsigned subjectIdx = job->y * gd->jobTable.fragmentSize_Y + j;
		int dS[lengthVector] __attribute__((aligned (64)));
		for(int offset=0; offset <= i; offset++) {
			int cellScore = gd->scoreMatrix.matrix[gd->query.dataCoded[queryIdx - offset]][gd->subject.dataCoded[subjectIdx + offset]];
			dS[offset] = cellScore;
		}
		// Sets to zero last cells to reset unused values of the vector
		for(int offset=i+1; offset<lengthVector; offset++)
			dS[offset] = -INFINITE;
		*deltaScore = _mm512_load_epi32 (dS);
}

void calculateAndAdvanceBody_512(int j, struct Node *retFragX, int i, struct GlobalData *gd, struct Job *job, struct Node_512 *arriba, struct Node_512 *izquierda, struct Node_512 *esquina, struct Node_512 * resultado, __m512i *deltaScore, int *bestScore) {
	calculate_512(gd, arriba, izquierda, esquina, resultado, deltaScore, bestScore);
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
		// Shift t, u and s
		arriba->t = _mm512_permutexvar_epi32 (arribaMask_512, resultado->t);
		arriba->u = _mm512_permutexvar_epi32 (arribaMask_512, resultado->u);
		arriba->s = _mm512_permutexvar_epi32 (arribaMask_512, resultado->s);
		setPos(arriba, 0, &(retFragX[i+1]));
	// Recalculates deltaScore
		// Calculus of the initial positions of the sequences to compare
		unsigned queryIdx = job->x * gd->jobTable.fragmentSize_X + i;
		unsigned subjectIdx = job->y * gd->jobTable.fragmentSize_Y + j;
		int dS[lengthVector] __attribute__((aligned (64)));
		for(int offset=0; offset < lengthVector; offset++) {
			int cellScore = gd->scoreMatrix.matrix[gd->query.dataCoded[queryIdx - offset]][gd->subject.dataCoded[subjectIdx + offset]];
			dS[offset] = cellScore;
		}
		*deltaScore = _mm512_load_epi32 (dS);
}

void calculateAndAdvanceBottomRightDiagonal_512(int j, struct Node *retFragX, int i, struct GlobalData *gd, struct Job *job, struct Node_512 *arriba, struct Node_512 *izquierda, struct Node_512 *esquina, struct Node_512 * resultado, struct Node *retFragY, __m512i *deltaScore, int *bestScore) {
	calculate_512(gd, arriba, izquierda, esquina, resultado, deltaScore, bestScore);
	int progress = i - job->realSize_X;
	// Saves last valid item of resultado into retFragX
		savePos(&retFragX[i - lengthVector + 1], resultado, lengthVector - 1);
	// Saves first valid item of resultado into retFragY
		savePos(&retFragY[j + progress + 1], resultado, progress);
//	if (i==18) {
//		printf("---------------i=%d:\n", i);
//		printf("---------------Arriba:\n");
//		displayNode_512(arriba);
//		printf("---------------Esquina:\n");
//		displayNode_512(esquina);
//		printf("---------------Izquierda:\n");
//		displayNode_512(izquierda);
//		printf("---------------Resultado:\n");
//		displayNode_512(resultado);
//		printf("---------------Fragment temporal X:\n");
//		for(int x=0; x<=job->realSize_X; x++)
//			displayNode(retFragX+x);
//		printf("\n");
//		printf("---------------Fragment temporal Y:\n");
//		for(int y=0; y<=job->realSize_Y; y++)
//			displayNode(retFragY+y);
//		printf("\n");
//		displayVector_512(deltaScore);
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
		// Shift t, u and s
		arriba->t = _mm512_permutexvar_epi32 (arribaMask_512, resultado->t);
		arriba->u = _mm512_permutexvar_epi32 (arribaMask_512, resultado->u);
		arriba->s = _mm512_permutexvar_epi32 (arribaMask_512, resultado->s);
	// Recalculates deltaScore
		// Calculus of the initial positions of the sequences to compare
		unsigned queryIdx = job->x * gd->jobTable.fragmentSize_X + i;
		unsigned subjectIdx = job->y * gd->jobTable.fragmentSize_Y + j;
		int dS[lengthVector] __attribute__((aligned (64)));
		// Sets to zero former cells to reset unused values of the vector
		for(int offset=0; offset < progress; offset++)
			dS[offset] = -INFINITE;
		for(int offset=progress; offset < lengthVector; offset++) {
			int cellScore = gd->scoreMatrix.matrix[gd->query.dataCoded[queryIdx - offset]][gd->subject.dataCoded[subjectIdx + offset]];
			dS[offset] = cellScore;
		}
		*deltaScore = _mm512_load_epi32 (dS);
}

/* Puts a node into a vector of nodes. Puts it at position 'i' */
//
inline void setPos(struct Node_512 *ptrNodeDest_512, int i, struct Node *ptrNodeSrc) {
	__mmask16 mask = 0x0001 << i;
	ptrNodeDest_512->t = _mm512_mask_set1_epi32(ptrNodeDest_512->t, mask, ptrNodeSrc->t);
	ptrNodeDest_512->u = _mm512_mask_set1_epi32(ptrNodeDest_512->u, mask, ptrNodeSrc->u);
	ptrNodeDest_512->s = _mm512_mask_set1_epi32(ptrNodeDest_512->s, mask, ptrNodeSrc->s);
}

/* Gets a node from a vector of nodes. Gets the node at position 'i' */
//
inline void savePos(struct Node *ptrNodeDest, struct Node_512 *ptrNodeSrc_512, int i) {
	int dump[lengthVector] __attribute__((aligned (64)));
	_mm512_store_epi32(dump, ptrNodeSrc_512->t);
	ptrNodeDest->t = dump[i];
	_mm512_store_epi32(dump, ptrNodeSrc_512->u);
	ptrNodeDest->u = dump[i];
	_mm512_store_epi32(dump, ptrNodeSrc_512->s);
	ptrNodeDest->s = dump[i];

}

// Displays a struct Node_512 in a formated way
void displayNode_512(struct Node_512 *node512) {
	int t[lengthVector] __attribute__((aligned (64)));
	_mm512_store_epi32(t, node512->t);
	int u[lengthVector] __attribute__((aligned (64)));
	_mm512_store_epi32(u, node512->u);
	int s[lengthVector] __attribute__((aligned (64)));
	_mm512_store_epi32(s, node512->s);

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

void displayVector_512(__m512i *vector) {
	int v[lengthVector] __attribute__((aligned (64)));
	_mm512_store_epi32(v, *vector);
	for(int i=0; i<lengthVector; i++)
		printf("| v=%8d \t", v[i]);
	printf("\n");
}

#undef lengthVector


/*
 * Vectorization256.c
 *
 *  Created on: Mar 25, 2021
 *      Author: galvez
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <immintrin.h>
#include "Utilities.h"
#include "Vectorization256.h"

#define lengthVector 8

struct Node_256 {
	__m256i t; // End in horizontal gap
	__m256i u; // End in vertical gap
	__m256i s; // Best value
};

// Declaration of vectorized constants
__m256i gapExtend_256;
__m256i openGap_256;
__m256i zeroes_256;
// Prepare shuffle
int arribaShuffle[] __attribute__((aligned (32))) = {0,0,1,2,3,4,5,6} ;
__m256i arribaMask_256;

extern int processJob_CISC(struct GlobalData *gd, struct Job *job, struct Node *retFragX, struct Node *retFragY);
void initializeVectors_256(int j, struct Node *retFragX, struct GlobalData *gd, struct Job *job, struct Node_256 *arriba, struct Node_256 *izquierda, struct Node_256 *esquina, __m256i *deltaScore);
void calculateAndAdvanceTopLeftDiagonal_256(int j, struct Node *retFragX, int i, struct GlobalData *gd, struct Job *job, struct Node_256 *arriba, struct Node_256 *izquierda, struct Node_256 *esquina, struct Node_256 * resultado, __m256i *deltaScore, int *bestScore);
void calculateAndAdvanceBody_256(int j, struct Node *retFragX, int i, struct GlobalData *gd, struct Job *job, struct Node_256 *arriba, struct Node_256 *izquierda, struct Node_256 *esquina, struct Node_256 * resultado, __m256i *deltaScore, int *bestScore);
void calculateAndAdvanceBottomRightDiagonal_256(int j, struct Node *retFragX, int i, struct GlobalData *gd, struct Job *job, struct Node_256 *arriba, struct Node_256 *izquierda, struct Node_256 *esquina, struct Node_256 * resultado, struct Node *retFragY, __m256i *deltaScore, int *bestScore);
inline void calculate_256(struct GlobalData *gd, struct Node_256 *arriba, struct Node_256 *izquierda, struct Node_256 *esquina, struct Node_256 * resultado, __m256i *deltaScore, int *bestScore);
inline void setPos(struct Node_256 *ptrNodeDest_256, int i, struct Node *ptrNodeSrc);
inline void savePos(struct Node *ptrNodeDest, struct Node_256 *ptrNodeSrc_256, int i);
inline int reduceMax_256(__m256i *ptrVector);
inline __m256i my_mm256_max_epi32(__m256i *a, __m256i *b);

/* DEBUG FUNCTIONS */
void displayNode_256(struct Node_256 *node256);
void displayVector_256(__m256i *vector);

/* Calculates the  scores in a job */
/* It does the job using TRUE-VECTOR operations */
int processJob_256(struct GlobalData *gd, struct Job *job, struct Node *retFragX, struct Node *retFragY) {
	if (job->realSize_Y % lengthVector != 0 || job->realSize_X <= 2*lengthVector ) return processJob_CISC(gd, job, retFragX, retFragY);
	// Creation of vectorized constants
	gapExtend_256 = _mm256_maskz_set1_epi32(0xFF, gd->gapExtend);
	openGap_256 = _mm256_maskz_set1_epi32(0xFF, gd->delete + gd->gapExtend);
	zeroes_256 = _mm256_maskz_set1_epi32(0xFF, 0);
	arribaMask_256 = _mm256_maskz_load_epi32  (0xFF, arribaShuffle);
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
inline void calculate_256(struct GlobalData *gd, struct Node_256 *arriba, struct Node_256 *izquierda, struct Node_256 *esquina, struct Node_256 * resultado, __m256i *deltaScore, int *bestScore) {
	// Deletion (horizontal gap). Calculus of resultado->t
	__m256i aux1 = _mm256_maskz_sub_epi32 (0xFF, izquierda->t, gapExtend_256);
	__m256i aux2 = _mm256_maskz_sub_epi32 (0xFF, izquierda->s, openGap_256);
	resultado->t = my_mm256_max_epi32 (&aux1, &aux2);
	// Insertion (vertical gap). Calculus of resultado->u
	aux1 = _mm256_maskz_sub_epi32 (0xFF, arriba->u, gapExtend_256);
	aux2 = _mm256_maskz_sub_epi32 (0xFF, arriba->s, openGap_256);
	resultado->u = my_mm256_max_epi32 (&aux1, &aux2);
	// Continue in diagonal. Calculus of resultado->s
	aux1 = _mm256_maskz_add_epi32 (0xFF, esquina->s, *deltaScore);
	aux1 = my_mm256_max_epi32 (&aux1, &resultado->t);
	resultado->s = my_mm256_max_epi32 (&aux1, &resultado->u);
	if (gd->algorithm == SMITH_WATERMAN) {
		resultado->s = my_mm256_max_epi32 (&zeroes_256, &resultado->s);
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
		// Shift t, u and s
		arriba->t = _mm256_permutexvar_epi32 (arribaMask_256, resultado->t);
		arriba->u = _mm256_permutexvar_epi32 (arribaMask_256, resultado->u);
		arriba->s = _mm256_permutexvar_epi32 (arribaMask_256, resultado->s);
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
		// Sets to zero last cells to reset unused values of the vector
		for(int offset=i+1; offset<lengthVector; offset++)
			dS[offset] = -INFINITE;
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
		// Shift t, u and s
		arriba->t = _mm256_permutexvar_epi32 (arribaMask_256, resultado->t);
		arriba->u = _mm256_permutexvar_epi32 (arribaMask_256, resultado->u);
		arriba->s = _mm256_permutexvar_epi32 (arribaMask_256, resultado->s);
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
		// Shift t, u and s
		arriba->t = _mm256_permutexvar_epi32 (arribaMask_256, resultado->t);
		arriba->u = _mm256_permutexvar_epi32 (arribaMask_256, resultado->u);
		arriba->s = _mm256_permutexvar_epi32 (arribaMask_256, resultado->s);
	// Recalculates deltaScore
		// Calculus of the initial positions of the sequences to compare
		unsigned queryIdx = job->x * gd->jobTable.fragmentSize_X + i;
		unsigned subjectIdx = job->y * gd->jobTable.fragmentSize_Y + j;
		int dS[lengthVector] __attribute__((aligned (32)));
		// Sets to zero former cells to reset unused values of the vector
		for(int offset=0; offset < progress; offset++)
			dS[offset] = -INFINITE;
		for(int offset=progress; offset < lengthVector; offset++) {
			int cellScore = gd->scoreMatrix.matrix[gd->query.dataCoded[queryIdx - offset]][gd->subject.dataCoded[subjectIdx + offset]];
			dS[offset] = cellScore;
		}
		*deltaScore = _mm256_maskz_load_epi32  (0xFF, dS);

}

/* Puts a node into a vector of nodes. Puts it at position 'i' */
//
inline void setPos(struct Node_256 *ptrNodeDest_256, int i, struct Node *ptrNodeSrc) {
	__mmask8 mask = 0x01 << i;
	ptrNodeDest_256->t = _mm256_mask_set1_epi32(ptrNodeDest_256->t, mask, ptrNodeSrc->t);
	ptrNodeDest_256->u = _mm256_mask_set1_epi32(ptrNodeDest_256->u, mask, ptrNodeSrc->u);
	ptrNodeDest_256->s = _mm256_mask_set1_epi32(ptrNodeDest_256->s, mask, ptrNodeSrc->s);
}

/* Gets a node from a vector of nodes. Gets the node at position 'i' */
//
inline void savePos(struct Node *ptrNodeDest, struct Node_256 *ptrNodeSrc_256, int i) {
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
inline int reduceMax_256(__m256i *ptrVector) {
	int dump[lengthVector] __attribute__((aligned (32)));
	_mm256_store_epi64(dump, *ptrVector);
	int ret = dump[0];
	for(int i=1; i<lengthVector; i++)
		if (dump[i] > ret) ret = dump[i];
	return ret;
}


/* AVX2 has no _mm_max_epi32 so we have to implement it by our own: https://fgiesen.wordpress.com/2016/04/03/sse-mind-the-gap/
*/
inline __m256i my_mm256_max_epi32(__m256i *a, __m256i *b) {
	__m256i cond = _mm256_cmpgt_epi32 (*a, *b);
	return _mm256_or_si256(_mm256_and_si256(*a, cond), _mm256_andnot_si256(cond, *b));

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

#undef lengthVector

/*
 * Vectorization128.c
 *
 *  Created on: Mar 25, 2021
 *      Author: galvez
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <immintrin.h>
#include "Utilities.h"
#include "Vectorization128.h"

#define lengthVector 4

struct Node_128 {
	__m128i t; // End in horizontal gap
	__m128i u; // End in vertical gap
	__m128i s; // Best value
};
// Vectorized constants
__m128i gapExtend_128;
__m128i openGap_128;
__m128i zeroes_128;

extern int processJob_CISC(struct GlobalData *gd, struct Job *job, struct Node *retFragX, struct Node *retFragY);
void initializeVectors_128(int j, struct Node *retFragX, struct GlobalData *gd, struct Job *job, struct Node_128 *arriba, struct Node_128 *izquierda, struct Node_128 *esquina, __m128i *deltaScore);
void calculateAndAdvanceTopLeftDiagonal_128(int j, struct Node *retFragX, int i, struct GlobalData *gd, struct Job *job, struct Node_128 *arriba, struct Node_128 *izquierda, struct Node_128 *esquina, struct Node_128 * resultado, __m128i *deltaScore, int *bestScore);
void calculateAndAdvanceBody_128(int j, struct Node *retFragX, int i, struct GlobalData *gd, struct Job *job, struct Node_128 *arriba, struct Node_128 *izquierda, struct Node_128 *esquina, struct Node_128 * resultado, __m128i *deltaScore, int *bestScore);
void calculateAndAdvanceBottomRightDiagonal_128(int j, struct Node *retFragX, int i, struct GlobalData *gd, struct Job *job, struct Node_128 *arriba, struct Node_128 *izquierda, struct Node_128 *esquina, struct Node_128 * resultado, struct Node *retFragY, __m128i *deltaScore, int *bestScore);
inline void calculate_128(struct GlobalData *gd, struct Node_128 *arriba, struct Node_128 *izquierda, struct Node_128 *esquina, struct Node_128 * resultado, __m128i *deltaScore, int *bestScore);
inline void setPos(struct Node_128 *ptrNodeDest_128, int i, struct Node *ptrNodeSrc);
inline void savePos(struct Node *ptrNodeDest, struct Node_128 *ptrNodeSrc_128, int i);
inline int reduceMax_128(__m128i *ptrVector);
inline __m128i my_mm128_max_epi32(__m128i *a, __m128i *b);

/* DEBUG FUNCTIONS */
void displayNode_128(struct Node_128 *node128);
void displayVector_128(__m128i *vector);

/* Calculates the  scores in a job */
/* It does the job using TRUE-VECTOR operations */
int processJob_128(struct GlobalData *gd, struct Job *job, struct Node *retFragX, struct Node *retFragY) {
	if (job->realSize_Y % lengthVector != 0 || job->realSize_X <= 2*lengthVector ) return processJob_CISC(gd, job, retFragX, retFragY);
	// Creation of vectorized constants
	gapExtend_128 = _mm_set1_epi32(gd->gapExtend);
	openGap_128 = _mm_set1_epi32(gd->delete + gd->gapExtend);
	zeroes_128 = _mm_set1_epi32(0);
	// Prepare the best score of all the jobs calculated by this thread in case of S/W
	int bestScore = -1;
	// Variable to store values of the score matrix
	__m128i deltaScore;
	// Nodes of vector to work with
	struct Node_128 arriba, izquierda, esquina, resultado;
	// Initializes the top right of the resulting job
	retFragY[0] = job->ptrRow[job->realSize_X];
	// Initializes retFragX. retFragX changes in each pass of vectorization
	memcpy(retFragX, job->ptrRow, sizeof(struct Node)*(job->realSize_X + 1));
	for(int j=0; j<job->realSize_Y - 1; j+=lengthVector) {
		initializeVectors_128(j, retFragX, gd, job, &arriba, &izquierda, &esquina, &deltaScore);
		int i=1;
		for(; i < lengthVector; i++) // Step 1
			calculateAndAdvanceTopLeftDiagonal_128(j, retFragX, i, gd, job, &arriba, &izquierda, &esquina, &resultado, &deltaScore, &bestScore);
		// Update progressively retFragX
		retFragX[0] = job->ptrColumn[j+lengthVector];
		for(; i < job->realSize_X; i++) // Main step
			calculateAndAdvanceBody_128(j, retFragX, i, gd, job, &arriba, &izquierda, &esquina, &resultado, &deltaScore, &bestScore);
		for(; i < job->realSize_X + lengthVector; i++) // Step 3
			calculateAndAdvanceBottomRightDiagonal_128(j, retFragX, i, gd, job, &arriba, &izquierda, &esquina, &resultado, retFragY, &deltaScore, &bestScore);
	}
	return (gd->algorithm == SMITH_WATERMAN)? bestScore : 0;
}

/* Initialization of vectors for an horizontal round of vectorization */
//
void initializeVectors_128(int j, struct Node *retFragX, struct GlobalData *gd, struct Job *job, struct Node_128 *arriba, struct Node_128 *izquierda, struct Node_128 *esquina, __m128i *deltaScore) {
	// Fill first node of esquina with retFragX[0]
	// The rest of nodes may be garbage
	esquina->t = _mm_set1_epi32(retFragX[0].t); // Superfluous
	esquina->u = _mm_set1_epi32(retFragX[0].u); // Superfluous
	esquina->s = _mm_set1_epi32(retFragX[0].s);
	// Fill first node of izquierda with job->ptrColumn[j+1]
	// The rest of nodes may be garbage
	izquierda->t = _mm_set1_epi32(job->ptrColumn[j+1].t);
	izquierda->u = _mm_set1_epi32(job->ptrColumn[j+1].u); // Superfluous
	izquierda->s = _mm_set1_epi32(job->ptrColumn[j+1].s);
	// Fill first node of arriba with retFragX[1]; and the second one with job->ptrColumn[j+1]
	// The rest of nodes may be garbage
	arriba->t = _mm_setr_epi32(0, 0, job->ptrColumn[j+1].t, retFragX[1].t);
	arriba->u = _mm_setr_epi32(0, 0, job->ptrColumn[j+1].u, retFragX[1].u);
	arriba->s = _mm_setr_epi32(0, 0, job->ptrColumn[j+1].s, retFragX[1].s);
	// Fill first node of deltaScore with data from the score matrix
	// The rest of nodes may be garbage
	// Calculus of the initial positions of the sequences to compare
	unsigned queryIdx = job->x * gd->jobTable.fragmentSize_X;
	unsigned subjectIdx = job->y * gd->jobTable.fragmentSize_Y + j;
	// We have to remember that the sequences begin at position 0
	int cellScore = gd->scoreMatrix.matrix[gd->query.dataCoded[queryIdx]][gd->subject.dataCoded[subjectIdx]];
	*deltaScore = _mm_set1_epi32(cellScore);
}

/* Calculus of a single vector in diagonal */
//
inline void calculate_128(struct GlobalData *gd, struct Node_128 *arriba, struct Node_128 *izquierda, struct Node_128 *esquina, struct Node_128 * resultado, __m128i *deltaScore, int *bestScore) {
	// Deletion (horizontal gap). Calculus of resultado->t
	__m128i aux1 = _mm_sub_epi32 (izquierda->t, gapExtend_128);
	__m128i aux2 = _mm_sub_epi32 (izquierda->s, openGap_128);
	resultado->t = my_mm128_max_epi32 (&aux1, &aux2);
	// Insertion (vertical gap). Calculus of resultado->u
	aux1 = _mm_sub_epi32 (arriba->u, gapExtend_128);
	aux2 = _mm_sub_epi32 (arriba->s, openGap_128);
	resultado->u = my_mm128_max_epi32 (&aux1, &aux2);
	// Continue in diagonal. Calculus of resultado->s
	aux1 = _mm_add_epi32 (esquina->s, *deltaScore);
	aux1 = my_mm128_max_epi32 (&aux1, &(resultado->t));
	resultado->s = my_mm128_max_epi32 (&aux1, &(resultado->u));
	if (gd->algorithm == SMITH_WATERMAN) {
		resultado->s = my_mm128_max_epi32 (&zeroes_128, &(resultado->s));
		int aux = reduceMax_128(&resultado->s);
		if (aux > *bestScore) *bestScore = aux;
	}
}

void calculateAndAdvanceTopLeftDiagonal_128(int j, struct Node *retFragX, int i, struct GlobalData *gd, struct Job *job, struct Node_128 *arriba, struct Node_128 *izquierda, struct Node_128 *esquina, struct Node_128 * resultado, __m128i *deltaScore, int *bestScore) {
	calculate_128(gd, arriba, izquierda, esquina, resultado, deltaScore, bestScore);
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
		// int arribaShuffle[] __attribute__((aligned (16))) = {0,0,1,2} ;
		const int imm8 = 0b10010000;
		// Shift t, u and s
		arriba->t = _mm_shuffle_epi32 (resultado->t, imm8);
		arriba->u = _mm_shuffle_epi32 (resultado->u, imm8);
		arriba->s = _mm_shuffle_epi32 (resultado->s, imm8);
		setPos(arriba, 0, &(retFragX[i+1]));
		setPos(arriba, i+1, &job->ptrColumn[j+i+1]);
	// Recalculates deltaScore
		// Calculus of the initial positions of the sequences to compare
		unsigned queryIdx = job->x * gd->jobTable.fragmentSize_X + i;
		unsigned subjectIdx = job->y * gd->jobTable.fragmentSize_Y + j;
		int dS[lengthVector] __attribute__((aligned (16)));
		for(int offset=0; offset <= i; offset++) {
			int cellScore = gd->scoreMatrix.matrix[gd->query.dataCoded[queryIdx - offset]][gd->subject.dataCoded[subjectIdx + offset]];
			dS[offset] = cellScore;
		}
		// Sets to zero last cells to reset unused values of the vector
		for(int offset=i+1; offset<lengthVector; offset++)
			dS[offset] = -INFINITE;
		*deltaScore = _mm_load_si128  ((__m128i *) dS);
}

void calculateAndAdvanceBody_128(int j, struct Node *retFragX, int i, struct GlobalData *gd, struct Job *job, struct Node_128 *arriba, struct Node_128 *izquierda, struct Node_128 *esquina, struct Node_128 * resultado, __m128i *deltaScore, int *bestScore) {
	calculate_128(gd, arriba, izquierda, esquina, resultado, deltaScore, bestScore);
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
		// int arribaShuffle[] __attribute__((aligned (16))) = {0,0,1,2} ;
		const int imm8 = 0b10010000;
		// Shift t, u and s
		arriba->t = _mm_shuffle_epi32 (resultado->t, imm8);
		arriba->u = _mm_shuffle_epi32 (resultado->u, imm8);
		arriba->s = _mm_shuffle_epi32 (resultado->s, imm8);
		setPos(arriba, 0, &(retFragX[i+1]));
	// Recalculates deltaScore
		// Calculus of the initial positions of the sequences to compare
		unsigned queryIdx = job->x * gd->jobTable.fragmentSize_X + i;
		unsigned subjectIdx = job->y * gd->jobTable.fragmentSize_Y + j;
		int dS[lengthVector] __attribute__((aligned (16)));
		for(int offset=0; offset < lengthVector; offset++) {
			int cellScore = gd->scoreMatrix.matrix[gd->query.dataCoded[queryIdx - offset]][gd->subject.dataCoded[subjectIdx + offset]];
			dS[offset] = cellScore;
		}
		*deltaScore = _mm_load_si128 ((__m128i *) dS);
}

void calculateAndAdvanceBottomRightDiagonal_128(int j, struct Node *retFragX, int i, struct GlobalData *gd, struct Job *job, struct Node_128 *arriba, struct Node_128 *izquierda, struct Node_128 *esquina, struct Node_128 * resultado, struct Node *retFragY, __m128i *deltaScore, int *bestScore) {
	calculate_128(gd, arriba, izquierda, esquina, resultado, deltaScore, bestScore);
	int progress = i - job->realSize_X;
	// Saves last valid item of resultado into retFragX
		savePos(&retFragX[i - lengthVector + 1], resultado, lengthVector - 1);
	// Saves first valid item of resultado into retFragY
		savePos(&retFragY[j + progress + 1], resultado, progress);
//	if (i==18) {
//		printf("---------------i=%d:\n", i);
//		printf("---------------Arriba:\n");
//		displayNode_128(arriba);
//		printf("---------------Esquina:\n");
//		displayNode_128(esquina);
//		printf("---------------Izquierda:\n");
//		displayNode_128(izquierda);
//		printf("---------------Resultado:\n");
//		displayNode_128(resultado);
//		printf("---------------Fragment temporal X:\n");
//		for(int x=0; x<=job->realSize_X; x++)
//			displayNode(retFragX+x);
//		printf("\n");
//		printf("---------------Fragment temporal Y:\n");
//		for(int y=0; y<=job->realSize_Y; y++)
//			displayNode(retFragY+y);
//		printf("\n");
//		displayVector_128(deltaScore);
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
		// int arribaShuffle[] __attribute__((aligned (16))) = {0,0,1,2} ;
		const int imm8 = 0b10010000;
		// Shift t, u and s
		arriba->t = _mm_shuffle_epi32 (resultado->t, imm8);
		arriba->u = _mm_shuffle_epi32 (resultado->u, imm8);
		arriba->s = _mm_shuffle_epi32 (resultado->s, imm8);
	// Recalculates deltaScore
		// Calculus of the initial positions of the sequences to compare
		unsigned queryIdx = job->x * gd->jobTable.fragmentSize_X + i;
		unsigned subjectIdx = job->y * gd->jobTable.fragmentSize_Y + j;
		int dS[lengthVector] __attribute__((aligned (16)));
		// Sets to zero former cells to reset unused values of the vector
		for(int offset=0; offset < progress; offset++)
			dS[offset] = -INFINITE;
		for(int offset=progress; offset < lengthVector; offset++) {
			int cellScore = gd->scoreMatrix.matrix[gd->query.dataCoded[queryIdx - offset]][gd->subject.dataCoded[subjectIdx + offset]];
			dS[offset] = cellScore;
		}
		*deltaScore = _mm_load_si128 ((__m128i *) dS);
}

/* Puts a node into a vector of nodes. Puts it at position 'i' */
//
inline void setPos(struct Node_128 *ptrNodeDest_128, int i, struct Node *ptrNodeSrc) {
	int dump[lengthVector] __attribute__((aligned (16)));
	_mm_store_si128((__m128i *) dump, ptrNodeDest_128->t);
	dump[i] = ptrNodeSrc->t;
	ptrNodeDest_128->t = _mm_load_si128((__m128i *) dump);
	_mm_store_si128((__m128i *) dump, ptrNodeDest_128->u);
	dump[i] = ptrNodeSrc->u;
	ptrNodeDest_128->u = _mm_load_si128((__m128i *) dump);
	_mm_store_si128((__m128i *) dump, ptrNodeDest_128->s);
	dump[i] = ptrNodeSrc->s;
	ptrNodeDest_128->s = _mm_load_si128((__m128i *) dump);
}

/* Gets a node from a vector of nodes. Gets the node at position 'i' */
//
inline void savePos(struct Node *ptrNodeDest, struct Node_128 *ptrNodeSrc_128, int i) {
	int dump[lengthVector] __attribute__((aligned (16)));
	_mm_store_si128((__m128i *) dump, ptrNodeSrc_128->t);
	ptrNodeDest->t = dump[i];
	_mm_store_si128((__m128i *) dump, ptrNodeSrc_128->u);
	ptrNodeDest->u = dump[i];
	_mm_store_si128((__m128i *) dump, ptrNodeSrc_128->s);
	ptrNodeDest->s = dump[i];

}

/* Returns the max value stored in a vector */
//
inline int reduceMax_128(__m128i *ptrVector) {
	int dump[lengthVector] __attribute__((aligned (16)));
	_mm_store_si128((__m128i *) dump, *ptrVector);
	int ret = dump[0];
	for(int i=1; i<lengthVector; i++)
		if (dump[i] > ret) ret = dump[i];
	return ret;
}

/* SSE3 has no _mm_max_epi32 so we have to implement it by our own: https://fgiesen.wordpress.com/2016/04/03/sse-mind-the-gap/
*/
inline __m128i my_mm128_max_epi32(__m128i *a, __m128i *b) {
	__m128i cond = _mm_cmpgt_epi32 (*a, *b);
	return _mm_or_si128(_mm_and_si128(*a, cond), _mm_andnot_si128(cond, *b));

}

// Displays a struct Node_128 in a formated way
void displayNode_128(struct Node_128 *node128) {
	int t[lengthVector] __attribute__((aligned (16)));
	_mm_store_si128((__m128i *) t, node128->t);
	int u[lengthVector] __attribute__((aligned (16)));
	_mm_store_si128((__m128i *) u, node128->u);
	int s[lengthVector] __attribute__((aligned (16)));
	_mm_store_si128((__m128i *) s, node128->s);

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

void displayVector_128(__m128i *vector) {
	int v[lengthVector] __attribute__((aligned (16)));
	_mm_store_si128((__m128i *) v, *vector);
	for(int i=0; i<lengthVector; i++)
		printf("| v=%8d \t", v[i]);
	printf("\n");
}

#undef lengthVector


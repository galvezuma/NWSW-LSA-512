/*
 * ScoreMatrix.h
 *
 *  Created on: Mar 11, 2021
 *      Author: galvez
 */

#ifndef SCOREMATRIX_H_
#define SCOREMATRIX_H_

#include <stdint.h>
#include <stdio.h>


#define MAX_LENGTH_ALPHABET 128
#define OUT_OF_ALPHABET 126

/* Definition of the score matrix. This stores i) the data itself, ii) the alphabets horizontal and vertical (just in case they are different)
 * related to columns and rows of the matrix (their lengths delimit the matrix size)
 * and the bytes [0-128] to code the sequences into numbers to address indirectly the matrix. */
struct GlobalData;
struct ScoreMatrix {
   int **matrix; // Store the numbers of the matrix
   char horizontalAlphabet[MAX_LENGTH_ALPHABET];
   char verticalAlphabet[MAX_LENGTH_ALPHABET];
   uint8_t horizontalCodification[MAX_LENGTH_ALPHABET];
   uint8_t verticalCodification[MAX_LENGTH_ALPHABET];
   // int id; // Just in case we should use an adaptive matrix
};

void loadDefaultMatrix(struct GlobalData *gd);
void loadUserMatrix(const char *filename, struct GlobalData *gd);
void displayScoreMatrix(FILE *out, struct ScoreMatrix *matrix);
void freeScoreMatrixStruct(struct ScoreMatrix *matrix);

#endif /* SCOREMATRIX_H_ */

/*
 * GlobalData.h
 *
 *  Created on: Mar 10, 2021
 *      Author: galvez
 */

#ifndef GLOBALDATA_H_
#define GLOBALDATA_H_

#include <stdint.h>

#define MAX_SIZEHEADER 512 // Maximum length of the header in a Fasta file
#define MAX_SIZELINE 256 // Maximum length of any data line in a Fasta file

/* Structure to store a Fasta file */
struct Sequence {
	char name[MAX_SIZEHEADER]; // Header of the file.
	char * data; // Content of the file, i.e., sequence of consecutive nucleotides
};

/* Global data accesible to any thread.
 * This tries to represent an Object.
 */
struct GlobalData {
	struct Sequence query;
	struct Sequence subject;
};

/* Definition of the score matrix. This stores i) the data itself, ii) the alphabets horizontal and vertical (just in case they are different)
 * related to columns and rows of the matrix (their lengths delimit the matrix size)
 * and the bytes [0-128] to code the sequences into numbers to address indirectly the matrix. */
 typedef struct ScoreMatrix ScoreMatrix;

 struct ScoreMatrix {
   int** matrix; // Store the numbers of the matrix
   char* horizontalAlphabet;
   char* verticalAlphabet;
   uint8_t * horizontalCodification;
   uint8_t * verticalCodification;
   // int id; // Just in case we should use an adaptive matrix
 };

#endif /* GLOBALDATA_H_ */

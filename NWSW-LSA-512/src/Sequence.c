/*
 * Sequence.c
 *
 *  Created on: Mar 11, 2021
 *      Author: galvez
 */

#include <ctype.h>
#include <stdint.h>
#include <string.h>
#include "Utilities.h"
#include "ScoreMatrix.h"

/* Puts the data of a sequence in Uppercase.
 * In addition, checks how many of the letters in the sequence (once in uppercase)
 * are not in the alphabet. Returns such a number.
 */
int toUpperCodeAndCheck(struct Sequence *seq, char *alphabet, uint8_t *codification, int **matrix) {
	// Auxiliar array to check quickly if a letter is valid or not.
	uint8_t validLetters[MAX_LENGTH_ALPHABET];
	memset(validLetters, 0, sizeof(validLetters)); // Sets to FALSE
	// Sets to TRUE only valid letters of the alphabet
	// We disable temporarily the warning message: array subscript has type ‘char’
	#pragma GCC diagnostic push
		#pragma GCC diagnostic ignored "-Wchar-subscripts"
		for(int i=0;i<strlen(alphabet); i++)
			validLetters[alphabet[i]] = 1;

		// Let's fulfill the dataCoded field of the sequence at the same time that we check the letters
		seq->dataCoded = (uint8_t *) internalMalloc(strlen(seq->data) * sizeof(uint8_t));

		// This loop iterates through the data of the sequence to:
		// 1.- Put each letter into Uppercase
		// 2.- Populate the dataCoded with the index
		// 3.- Check if the letter is valid (in the alphabet) or not
		// Calculate best score against itself
		int countInvalidLetters = 0;
		seq->againstItselfScore = 0;
		int const strlenSeqData = strlen(seq->data);
		for (int i=0; i<strlenSeqData; i++) {
			seq->data[i] = toupper(seq->data[i]);
			seq->dataCoded[i] = codification[seq->data[i]];
			if (! validLetters[seq->data[i]])
				countInvalidLetters++;
			seq->againstItselfScore += matrix[seq->dataCoded[i]][seq->dataCoded[i]];
		}
	#pragma GCC diagnostic pop
	return countInvalidLetters;
}

void freeSequenceStruct(struct Sequence *seq) {
	internalFree((void **) &(seq->dataCoded));
	internalFree((void **) &(seq->data));
}

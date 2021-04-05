/*
 * Sequence.h
 *
 *  Created on: Mar 11, 2021
 *      Author: galvez
 */

#ifndef SEQUENCE_H_
#define SEQUENCE_H_

#include <stdint.h>

#define MAX_SIZEHEADER 512 // Maximum length of the header in a Fasta file

/* Structure to store a Fasta file */
struct Sequence {
	char name[MAX_SIZEHEADER]; // Header of the file.
	char * data; // Content of the file, i.e., sequence of consecutive nucleotides
	uint8_t *dataCoded; // Translation of the letters in the file into positions in an alphabet
};

struct NodeListSequence {
	struct Sequence sequence;
	struct NodeListSequence * ptrNext;
};

int toUpperCodeAndCheck(struct Sequence *seq, char *alphabet, uint8_t *codification);
void freeSequenceStruct(struct Sequence *seq);

#endif /* SEQUENCE_H_ */

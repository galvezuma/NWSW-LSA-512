/*
 * Fasta.h
 *
 *  Created on: Apr 2, 2021
 *      Author: galvez
 */

#ifndef FASTA_H_
#define FASTA_H_

#include "Sequence.h"

void readSingleFastaFile(char* filename, struct Sequence *seq);
struct NodeListSequence * readMultipleFastaFile(char* filename);
void freeNodeListSequenceStruct(struct NodeListSequence **ptrList);
void displayListSequence(struct NodeListSequence *ptrList);

#endif /* FASTA_H_ */

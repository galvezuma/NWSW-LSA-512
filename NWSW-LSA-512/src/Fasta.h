/*
 * Fasta.h
 *
 *  Created on: Apr 2, 2021
 *      Author: galvez
 */

#ifndef FASTA_H_
#define FASTA_H_

#include "Sequence.h"
#include "BackwardsPairwise.h"

void readSingleFastaFile(char* filename, struct Sequence *seq);
struct NodeListSequence * readMultipleFastaFile(char* filename);
void saveFastaPairwiseAlignment(char *filename, struct FastaPairwiseAlignment *fpa);
void freeNodeListSequenceStruct(struct NodeListSequence **ptrList);
void displayListSequence(struct NodeListSequence *ptrList);

#endif /* FASTA_H_ */

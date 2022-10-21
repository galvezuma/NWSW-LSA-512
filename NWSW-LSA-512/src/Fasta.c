/*
 * Fasta.c
 *
 *  Created on: Apr 2, 2021
 *      Author: galvez
 */

#include <ctype.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>
#include "Fasta.h"
#include "Utilities.h"
#include "Sequence.h"

/** Reads a Sequence from a Fasta file. Takes the header and the content. */
void readSingleFastaFile(char* filename, struct Sequence *seq) {
	char* fastaHeader = seq->name;
	char** fastaContent = &(seq->data);
	char * ok;
	FILE *file;
	if((file = fopen(filename, "r")) == NULL) {
    	fatalError1("Unable to open the input Fasta file: %s\n", filename); // This HALTS
  	}

  	// The first line is the header (fastaHeader FASTA)
  	ok = fgets(fastaHeader, MAX_SIZEHEADER, file);
  	if ((! ok) || (fastaHeader[0]!='>')) { // If it is empty o does not starts with '>', then it is not Fasta
    	fatalError1("Fasta file %s should begin with >.\n", filename); // This HALTS
  	}
  	int sizeFirstLine = strlen(fastaHeader);
  	// Remove the CR and LF at the end of the header
  	while(fastaHeader[strlen(fastaHeader)-1] == '\n' || fastaHeader[strlen(fastaHeader)-1] == '\r')
  		fastaHeader[strlen(fastaHeader)-1] = '\0';

  	// Let's allocate space for the Fasta content (the sequence in nucleotides)
  	// Actually, we allocate more space than required because we store also de CR and LF. See later realloc
  	fseek(file, 0, SEEK_END); // Let's go to the end of the file to estimate the maximum size of the Fasta sequence
  	long size = ftell(file) - sizeFirstLine;
  	*fastaContent = (char *) internalMalloc(sizeof(char) * (size + 1));

  	// We will iterate through the lines of the file
  	int offset = 0; //Initially the allocated memory is empty
  	fseek(file, sizeFirstLine, SEEK_SET); // We start reading after the header of the Fasta file
  	while (fgets((*fastaContent)+offset, MAX_SIZELINE, file) != NULL) { // Data is directly read over the allocated memory
  		int tempOffset = offset;
  		if (((*fastaContent)+offset)[0] == '>') fatalError2("Fasta file %s should contain a single sequence.\nAnother sequence has been found %s\n", filename, (*fastaContent)+offset); // This HALTS
  		offset += strlen((*fastaContent)+offset);
  		while ((tempOffset < offset) && ( (*fastaContent)[offset-1]=='\n' || (*fastaContent)[offset-1]=='\r' ) ){ // CR and LF are removed
  			offset--;
		}
  	}
	if (offset == 0) fatalError1("Sequence %s has length 0.\n", fastaHeader);
  	(*fastaContent)[offset] = '\0'; // The Fasta string must be finished (delimited by \0)
  	// A realloc is needed to take exclusively the required amount of data
  	*fastaContent = (char *) realloc((*fastaContent), sizeof(char) * (strlen(*fastaContent) + 1));

	fclose(file); // The file should be closed
}

// Inserts a node at the end of a list of nodes
// Implementation is recursive
void insertNodeSequenceList(struct NodeListSequence **ptrList, char *sequenceHeader, char *sequenceContent) {
	if ((*ptrList) != NULL) { // Recursive case
		insertNodeSequenceList(&((*ptrList)->ptrNext), sequenceHeader, sequenceContent);
	} else { // Base case
		if (strlen(sequenceContent) == 0) fatalError1("Sequence %s has length 0.\n", sequenceHeader);
		(*ptrList) = (struct NodeListSequence *) malloc(sizeof(struct NodeListSequence));
		strcpy((*ptrList)->sequence.name, sequenceHeader);
		(*ptrList)->sequence.data = strdup(sequenceContent);
		(*ptrList)->ptrNext = NULL;
	}
}

/** Reads a Sequence from a Fasta file. Takes the header and the content. */
struct NodeListSequence * readMultipleFastaFile(char* filename) {
	struct NodeListSequence * ret = NULL;
	//
	char * ok;
	FILE *file;
	if((file = fopen(filename, "r")) == NULL) {
    	fatalError1("Unable to open the input Fasta file: %s\n", filename); // This HALTS
  	}
  	// Initially, the memory for each sequence is the length of the file.
  	// Later, we will malloc it to the exact amount of memory required
  	fseek(file, 0, SEEK_END); // Let's go to the end of the file to estimate the maximum size of the Fasta sequence
  	long maxSizeOfSequence = ftell(file);
  	fseek(file, 0, SEEK_SET); // We start from the beginning of the Fasta file
  	//
  	// LET'S READ SEQUENCE BY SEQUENCE
	char sequenceHeader[MAX_SIZEHEADER];
	char sequenceContent[maxSizeOfSequence];
  	int offset; // To put new lines read
	do {
		//
	  	// The first line is the header (sequenceHeader FASTA)
	  	ok = fgets(sequenceHeader, MAX_SIZEHEADER, file);
	  	if ((! ok) || (sequenceHeader[0]!='>')) // If it is empty o does not starts with '>', then it is not Fasta
	    	fatalError1("Fasta file %s should begin with >.\n", filename); // This HALTS
	  	// Remove the CR and LF at the end of the header
	  	while(sequenceHeader[strlen(sequenceHeader)-1] == '\n' || sequenceHeader[strlen(sequenceHeader)-1] == '\r')
	  		sequenceHeader[strlen(sequenceHeader)-1] = '\0';
	  	//
	  	// The next lines are the sequence content
	  	offset = 0; //Initially the allocated memory is empty
	  	long toRewindToPreviousLine = ftell(file);
	  	while (fgets(sequenceContent+offset, MAX_SIZELINE, file) != NULL) { // Data is directly read over the allocated memory
	  		int tempOffset = offset;
	  		if (sequenceContent[offset] == '>') { // A whole sequence has been read and it must be stored into a node
	  			sequenceContent[offset] = '\0'; // The Fasta string must be finished (delimited by \0)
	  			insertNodeSequenceList(&ret, sequenceHeader, sequenceContent);
	  			fseek(file, toRewindToPreviousLine, SEEK_SET);
	  			break;
	  		}
	  		offset += strlen(sequenceContent+offset);
	  		while ((tempOffset < offset) && ( sequenceContent[offset-1]=='\n' || sequenceContent[offset-1]=='\r' ) ){ // CR and LF are removed
	  			offset--;
			}
		  	toRewindToPreviousLine = ftell(file);
	  	}
	} while(! feof(file));
  	sequenceContent[offset] = '\0'; // The Fasta string must be finished (delimited by \0)
	insertNodeSequenceList(&ret, sequenceHeader, sequenceContent);
	return ret;
}

// Saves a pairwise alignment result in text Fasta format
//
void saveFastaPairwiseAlignment(char *filename, struct FastaPairwiseAlignment *fpa) {
	FILE *file;
	if((file = fopen(filename, "wt")) == NULL) {
    	fatalError1("Unable to open the output Fasta alignment file: %s\n", filename); // This HALTS
  	}
	// Saves one aligned sequence after the other
	fprintf(file, "%s\n", fpa->sequence[0]->name);
	fprintf(file, "%s\n", fpa->alignment[0]);
	fprintf(file, "%s\n", fpa->sequence[1]->name);
	fprintf(file, "%s\n", fpa->alignment[1]);
	//
	fclose(file);
}

// Frees the memory allocated for a list of nodes of sequences
// Implementation is recursive
void freeNodeListSequenceStruct(struct NodeListSequence **ptrList) {
	if ((*ptrList) != NULL) { // Recursive case
		freeNodeListSequenceStruct(&((*ptrList)->ptrNext)); // Recursive call
		freeSequenceStruct(&((*ptrList)->sequence));
		internalFree((void **)ptrList);
	} // Base case is empty
}

void displayListSequence(struct NodeListSequence *list) {
	while(list != NULL) {
		fprintf(stdout, "Header: %s\n", list->sequence.name);
		fprintf(stdout, "Content length: %zu: %s\n", strlen(list->sequence.data), list->sequence.data);
		list = list->ptrNext;
	}
}


/*
 * Utilities.h
 *
 *  Created on: Mar 10, 2021
 *      Author: galvez
 */

#ifndef UTILITIES_H_
#define UTILITIES_H_

#include <unistd.h>
#include <stdlib.h>
#include "Sequence.h"

// Definition of our infinite. It could be also INT_MAX
#define INFINITE 999999

/* Displays a message and HALTS */
#define fatalError0(message)  { fprintf(stderr, message); exit(EXIT_FAILURE); }
#define fatalError1(message, a) { fprintf(stderr, message, a); exit(EXIT_FAILURE); }
#define fatalError2(message, a, b) { fprintf(stderr, message, a, b); exit(EXIT_FAILURE); }
//max is taken from https://stackoverflow.com/questions/3437404/min-and-max-in-c
#define max(a,b) \
   ({ __typeof__ (a) _a = (a); \
       __typeof__ (b) _b = (b); \
     _a > _b ? _a : _b; })

void * internalMalloc(size_t numBytes);
void internalFree(void **ptr);
char * toString(char c);
int ceilDivision(int a, int b);
int restOfLengthSequence(int a, int b);
int getNumberOfCores();

void readFastaFile(char* filename, struct Sequence *seq);

#endif /* UTILITIES_H_ */

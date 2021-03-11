/*
 * Utilities.h
 *
 *  Created on: Mar 10, 2021
 *      Author: galvez
 */

#ifndef UTILITIES_H_
#define UTILITIES_H_

#include "GlobalData.h"

/* Displays a message and HALTS */
#define fatalError0(message)  { fprintf(stderr, message); exit(EXIT_FAILURE); }
#define fatalError1(message, a) { fprintf(stderr, message, a); exit(EXIT_FAILURE); }
#define fatalError2(message, a, b) { fprintf(stderr, message, a, b); exit(EXIT_FAILURE); }

void * internalMalloc(size_t numBytes);

void readFastaFile(char* filename, struct Sequence *seq);

#endif /* UTILITIES_H_ */

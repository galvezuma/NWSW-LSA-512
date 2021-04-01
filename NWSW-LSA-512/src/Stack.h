/*
 * Stack.h
 *
 *  Created on: Mar 15, 2021
 *      Author: galvez
 */

#ifndef STACK_H_
#define STACK_H_

#include <stdio.h>
#include <stdlib.h>
#include "JobTable.h"

struct Stack {
	struct Job * volatile *ptrArray; // Array of pointers to Job
	unsigned maxSize;
	volatile unsigned size;
};

struct Stack createStack(struct GlobalData *gd, unsigned size);
unsigned volatile isEmptyStack(struct Stack *ptrStack);
void push(struct Stack *ptrStack, struct Job *job);
struct Job * volatile pop(struct Stack *ptrStack);
void freeStack(struct Stack *ptrStack);

#endif /* STACK_H_ */

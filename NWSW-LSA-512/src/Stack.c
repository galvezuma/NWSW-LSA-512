/*
 * Stack.c
 *
 *  Created on: Mar 15, 2021
 *      Author: galvez
 */

#include <stdio.h>
#include <stdlib.h>
#include "Stack.h"
#include "JobTable.h"
#include "Utilities.h"

struct Stack createStack(unsigned size) {
	struct Stack ret;
	ret.ptrArray = (struct Job **) internalMalloc(size * sizeof(struct Job *));
	ret.maxSize = size;
	ret.size = 0;
	return ret;
}

unsigned isEmptyStack(struct Stack *ptrStack) {
	return ptrStack->size == 0;
}

void push(struct Stack *ptrStack, struct Job *job) {
	if (ptrStack->size == ptrStack->maxSize) fatalError0("Stack of jobs overflow\n");
	ptrStack->ptrArray[ptrStack->size] = job;
	ptrStack->size++;
}

struct Job * pop(struct Stack *ptrStack) {
	if (isEmptyStack(ptrStack )) fatalError0("Stack of jobs underflow\n");
	ptrStack->size--;
	return ptrStack->ptrArray[ptrStack->size];
}

void freeStack(struct Stack *ptrStack) {
	internalFree((void **)  &ptrStack->ptrArray);
}

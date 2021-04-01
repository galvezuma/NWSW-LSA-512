/*
 * Stack.c
 *
 *  Created on: Mar 15, 2021
 *      Author: galvez
 */

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <sys/syscall.h>
#include "Stack.h"
#include "JobTable.h"
#include "Utilities.h"
#include "GlobalData.h"

struct Stack createStack(struct GlobalData *gd, unsigned size) {
	struct Stack ret;
	ret.ptrArray = (struct Job **) internalMalloc(size * sizeof(struct Job *));
	ret.maxSize = size;
	ret.size = 0;
	if (gd->verbose) {
		fprintf(stdout, "Created a stack of job with maximum size: %d\n", ret.maxSize);
	}
	return ret;
}

unsigned volatile isEmptyStack(struct Stack *ptrStack) {
	return ptrStack->size == 0;
}

void push(struct Stack *ptrStack, struct Job *job) {
	if (ptrStack->size == ptrStack->maxSize) {
		printf("Stack size: %d\n", ptrStack->size);
		for (int i=0; i<ptrStack->maxSize; i++)
			printf("[i=%d (%d,%d)]", i, ptrStack->ptrArray[i]->x, ptrStack->ptrArray[i]->y);
		printf("\nInserting: %d,%d\n", job->x, job->y);
		fatalError0("Stack of jobs overflow\n");
	}
	ptrStack->ptrArray[ptrStack->size] = job;
	ptrStack->size++;
	pid_t id = syscall(__NR_gettid);
	printf("\tInserting job=(%d,%d) ID=%d; Size=%d\n", ptrStack->ptrArray[ptrStack->size-1]->x, ptrStack->ptrArray[ptrStack->size-1]->y, id, ptrStack->size);
}

struct Job * volatile pop(struct Stack *ptrStack) {
	if (isEmptyStack(ptrStack )) fatalError0("Stack of jobs underflow\n");
	ptrStack->size--;
	pid_t id = syscall(__NR_gettid);
	printf("Taking job=(%d,%d) ID=%d; Size=%d\n", ptrStack->ptrArray[ptrStack->size]->x, ptrStack->ptrArray[ptrStack->size]->y, id, ptrStack->size);
	return ptrStack->ptrArray[ptrStack->size];
}

void freeStack(struct Stack *ptrStack) {
	internalFree((void **)  &ptrStack->ptrArray);
}

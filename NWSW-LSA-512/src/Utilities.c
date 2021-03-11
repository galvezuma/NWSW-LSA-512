/*
 * Utilities.c
 *
 *  Created on: Mar 10, 2021
 *      Author: galvez
 */

#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include "Utilities.h"

/* Calls to malloc. I any error happens then HALTS */
void * internalMalloc(size_t numBytes) {
	void *temp = malloc(numBytes);
	if (temp == NULL) fatalError0("Memory overflow. malloc returns NULL.\n");
	return temp;
}

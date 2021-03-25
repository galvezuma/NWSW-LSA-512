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

#ifdef _WIN32
#include <windows.h>
#elif MACOS
#include <sys/param.h>
#include <sys/sysctl.h>
#else
#include <unistd.h>
#endif

/* Calls to malloc. I any error happens then HALTS */
void * internalMalloc(size_t numBytes) {
	void *temp = malloc(numBytes);
	if (temp == NULL) fatalError0("Memory overflow. malloc returns NULL.\n");
	return temp;
}

/* Calls to free and resets the pointer to NULL */
void internalFree(void **ptr) {
	free(*ptr);
	*ptr = NULL;
}

/* Converts a char into a string (char *) to be printf-ed with %s */
static char internalText[]=" ";
char * toString(char c) {
	internalText[0] = c;
	return internalText;
}

/* Calculates a/b and return the nearest greater integer */
int ceilDivision(int a, int b) {
	if (a%b == 0) return a/b;
	return 1+(a/b);
}

int restOfLengthSequence(int a, int b) {
	return (a%b == 0)? b : a%b;
}

/* Returns the number of logical cores in the system */
int getNumberOfCores() {
#ifdef WIN32
    SYSTEM_INFO sysinfo;
    GetSystemInfo(&sysinfo);
    return sysinfo.dwNumberOfProcessors;
#elif MACOS
    int nm[2];
    size_t len = 4;
    uint32_t count;

    nm[0] = CTL_HW; nm[1] = HW_AVAILCPU;
    sysctl(nm, 2, &count, &len, NULL, 0);

    if(count < 1) {
    nm[1] = HW_NCPU;
    sysctl(nm, 2, &count, &len, NULL, 0);
    if(count < 1) { count = 1; }
    }
    return count;
#else
    return sysconf(_SC_NPROCESSORS_ONLN);
#endif
}

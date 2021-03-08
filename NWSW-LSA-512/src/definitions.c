/*
 * definitions.c
 *
 *  Created on: Mar 8, 2021
 *      Author: galvez
 */

/* Implementación de la función malloc_safe_agr248() necesitaria para evadir un bug de la Tilera. */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "definitions.h"

long my_clock(){
	struct timeval t;
	gettimeofday(&t, NULL);
	// printf("Ejecutando clock %ld %ld\n", t.tv_sec, t.tv_usec);
	return (t.tv_sec&0x0000FFFFFFFFFFFF) * 1000l + t.tv_usec / 1000l;
}


 void * malloc_safe_agr248(int size) {
	void * memory = malloc(size);
	if (memory == NULL) // Tile64 bug in ranges 130961-131016 (1024*128 - [111,56]) and 65425-65480 (1024*64 - [111,56]): need to malloc again
		memory = malloc(size);
	return memory;
 }

  void * malloc_shared(int size) {
	void * memory = malloc(size);
	return memory;
 }

 void * memdup(void * ptr, int size){
	void * ptr_aux = malloc(size);
	memcpy(ptr_aux, ptr, size);
	return ptr_aux;
}

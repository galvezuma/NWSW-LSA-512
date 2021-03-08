/*
 * MyLib.c
 *
 *  Created on: Mar 8, 2021
 *      Author: galvez
 */

#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <pthread.h>
#include <unistd.h>
#include <stdint.h>
#include <immintrin.h>
#include "MyLib.h"
#include "definitions.h"

void mylib_init() {
	;
}

int mylib_proc_remaining(){
	return INT_MAX;
}

#pragma GCC diagnostic ignored "-Wformat-security"
void mylib_die(char * msg, ...){
	fprintf(stderr, msg);
};
#pragma GCC diagnostic warning "-Wformat-security"

pthread_t * mylib_proc_spawn(int num_hilos, mylibProcParam param) {
	int rc;
	mylibThreadParam * ptr_t_param;
	pthread_t * threads= (pthread_t *)malloc(sizeof(pthread_t)*num_hilos);
	for(int i=0; i<num_hilos; i++){
		ptr_t_param = memdup(&(param.t_params), sizeof(param.t_params));
		ptr_t_param->num=i;
		rc = pthread_create(&(threads[i]), NULL, param.nameOfFunction, (void *)ptr_t_param);
		if (rc) return NULL;
	}
	printf("Hilos creados: %d.\n", num_hilos);
	return threads;
}

void mylib_wait_proc_vanishes(pthread_t * set, int num_hilos){
	for(int i=0; i<num_hilos; i++){
		pthread_join(set[i], NULL);
		printf("Proceso %d finalizado\n", i);
	}
	free(set);
}

// Functions for AVX512
inline int my_mm512_extract_epi32(__m512i v, int pos) {
	int32_t address[16] __attribute__((aligned (64)));
	_mm512_store_epi32 (&address, v);
	return v[pos];
}


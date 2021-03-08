/*
 * MyLib.h
 *
 *  Created on: Mar 8, 2021
 *      Author: galvez
 */

#ifndef MYLIB_H_
#define MYLIB_H_

#include <pthread.h>

#define MYLIB_SUCCESS 0

typedef struct {
	int num;
	char estado;
} mylibThreadParam;

typedef struct {
	void * (*nameOfFunction) (void *);
	mylibThreadParam t_params;
} mylibProcParam;

extern void mylib_init();
#define mylib_finish()	printf("\tCerrando\n"); pthread_exit(NULL);
extern int mylib_proc_remaining();
extern void mylib_die(char * msg, ...);
extern pthread_t * mylib_proc_spawn(int num_hilos, mylibProcParam param);
extern void mylib_wait_proc_vanishes(pthread_t * set, int num);

#endif /* MYLIB_H_ */

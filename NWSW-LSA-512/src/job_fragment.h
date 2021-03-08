/*
 * job_fragment.h
 *
 *  Created on: Mar 8, 2021
 *      Author: galvez
 */

#ifndef JOB_FRAGMENT_H_
#define JOB_FRAGMENT_H_

#include "definitions.h"

/* Definición de un fragmento de trabajo, que es el cuanto de filas o columnas con datos con los que el
 * algoritmo Needleman-Wunsch paralelo y sus subtrabajos trabajarán. */

 extern struct Nodo_fragmento* crea_fragmento_shared(int size);

 extern struct Nodo_fragmento* crea_fragmento(int size);

 extern void copia_fragmento(struct Nodo_fragmento** copia, struct Nodo_fragmento* original);

 extern void copia_fragmento_shared(struct Nodo_fragmento** copia, struct Nodo_fragmento* original);

 extern void libera_fragmento(struct Nodo_fragmento** ptr_fragmento);

 extern void inicializa_fragmento(struct Nodo_fragmento* fragmento, int horizontal, long valor_inicial);

 extern int set_elemento_fragmento(struct Nodo_fragmento* fragmento, int i, struct Nodo_valor *e);

 extern int get_elemento_fragmento(struct Nodo_fragmento* fragmento, int i, struct Nodo_valor *e);

 extern long get_max_fragmento(struct Nodo_fragmento* fragmento, int* posicion);

 extern int get_size(struct Nodo_fragmento* fragmento);

 extern void print_fragmento();

 /** Macro para copiar  los valores del nodo de B en A. */
 #define copia_valor_nodo(a,b) { \
 		Nodo_valor *nodo_temp1; \
 		Nodo_valor *nodo_temp2; \
		nodo_temp1 = &(b); \
      	nodo_temp2 = &(a); \
    	nodo_temp2->valor = nodo_temp1->valor; \
    	nodo_temp2->vertical = nodo_temp1->vertical; \
    	nodo_temp2->horizontal = nodo_temp1->horizontal; \
 }

#endif /* JOB_FRAGMENT_H_ */

/*
 * job_fragment.c
 *
 *  Created on: Mar 8, 2021
 *      Author: galvez
 */

#include <malloc.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "job_fragment.h"

/* Implementación del módulo de tratamiento de fragmentos de trabajo, que son el cuanto de filas o columnas con datos con los que el
 * algoritmo Needleman-Wunsch paralelo y sus subtrabajos trabajarán. */

 /* Crea y devuelve un nuevo fragmento de tamaño de array indicado en MEMORIA COMPARTIDA. */
 struct Nodo_fragmento* crea_fragmento_shared(int size){
 	//Liberamos memoria para el fragmento (estructura)
 	struct Nodo_fragmento* frag = (struct Nodo_fragmento *)malloc_shared(sizeof (struct Nodo_fragmento));
 	//Liberamos memoria para el array del tamaño determinado
 	struct Nodo_valor* array = (struct Nodo_valor *)malloc_shared(sizeof(struct Nodo_valor) * size);

 	if (frag == NULL || array == NULL) {
 		printf("ERROR: Al crear un fragmento en memoria compartida, la tarjeta se ha quedado sin memoria.\n");
 		imprime_estado_memoria();
 		exit(-1);
 	}

 	//Asignamos y devolvemos
 	frag->contenido = array;
 	frag->longitud = size;
 	return frag;
 }

 /* Crea y devuelve un nuevo fragmento de tamaño de array indicado en MEMORIA NO COMPARTIDA. */
 struct Nodo_fragmento* crea_fragmento(int size){
 	//Liberamos memoria para el fragmento (estructura)
 	struct Nodo_fragmento* frag = (struct Nodo_fragmento *)malloc_safe_agr248(sizeof (struct Nodo_fragmento));
 	//Liberamos memoria para el array del tamaño determinado
 	struct Nodo_valor* array = (struct Nodo_valor *)malloc_safe_agr248(sizeof(struct Nodo_valor) * size);

 	if (frag == NULL || array == NULL) {
 		printf("ERROR: Al crear un fragmento en memoria local, la tarjeta se ha quedado sin memoria.\n");
 		imprime_estado_memoria();
 		exit(-1);
 	}

 	//Asignamos y devolvemos
 	frag->contenido = array;
 	frag->longitud = size;
 	return frag;
 }

 /* Crea una copia del fragmento en MEMORIA NO COMPARTIDA. */
 void copia_fragmento(struct Nodo_fragmento** copia, struct Nodo_fragmento* original){
 	//Liberamos memoria para el fragmento (estructura)
 	struct Nodo_fragmento* frag = (struct Nodo_fragmento *)malloc_safe_agr248(sizeof (struct Nodo_fragmento));
 	//Liberamos memoria para el array del tamaño determinado y copiamos datos
 	int size = original->longitud;
 	struct Nodo_valor* array = (struct Nodo_valor *)malloc_safe_agr248(sizeof(struct Nodo_valor) * size);

 	if (frag == NULL || array == NULL) {
 		printf("ERROR: Al realizar la copia de un fragmento a memoria local, la tarjeta se ha quedado sin memoria.\n");
 		imprime_estado_memoria();
 		exit(-1);
 	}

 	memcpy(array, original->contenido, sizeof(struct Nodo_valor) * size);

 	//Asignamos y devolvemos
 	frag->contenido = array;
 	frag->longitud = size;

 	//Apuntamos a la copia
 	*copia = frag;
 }

 /* Crea una copia del fragmento en MEMORIA COMPARTIDA. */
 void copia_fragmento_shared(struct Nodo_fragmento** copia, struct Nodo_fragmento* original){
 	//Liberamos memoria para el fragmento (estructura)
 	struct Nodo_fragmento* frag = (struct Nodo_fragmento *)malloc_shared(sizeof (struct Nodo_fragmento));
 	//Liberamos memoria para el array del tamaño determinado y copiamos datos
 	int size = original->longitud;
 	struct Nodo_valor* array = (struct Nodo_valor *)malloc_shared(sizeof(struct Nodo_valor) * size);

 	if (frag == NULL || array == NULL) {
 		printf("ERROR: Al realizar la copia de un fragmento a memoria compartida, la tarjeta se ha quedado sin memoria.\n");
 		imprime_estado_memoria();
 		exit(-1);
 	}

 	memcpy(array, original->contenido, sizeof(struct Nodo_valor) * size);

 	//Asignamos y devolvemos
 	frag->contenido = array;
 	frag->longitud = size;

 	//Apuntamos a la copia
 	*copia = frag;
 }

 /* Libera la memoria de un fragmento. */
 void libera_fragmento(struct Nodo_fragmento** ptr_fragmento) {
	// Liberamos el array
	free((*ptr_fragmento)->contenido);
	// Liberamos el propio puntero al fragmento
	free(*ptr_fragmento);
	*ptr_fragmento = NULL;
 }

 /* Inicializa un fragmento, dando los valores a las celdas de las tres matrices.
  * Optimización para alineamientos que comienzan con desplazamientos (huecos en una) sin
  * penalización: los fragmentos horizontales toman [0,0,valor] y los verticales [0,valor,0] */
 void inicializa_fragmento(struct Nodo_fragmento* fragmento, int horizontal, long valor_inicial){
 	int size = fragmento->longitud;
 	Nodo_valor *temp;
	for(int i=0; i<size; i++){
		temp = &(fragmento->contenido[i]);

		// Valor 0 para permitir desplazamientos con penalización 0
		temp->valor = 0;

        // Valor horizontal y vertical a 0 en sus casos, en otro, toman valor como inic.
        // Versión que almacena diferencias entre valor y matrices auxiliares: valor (0) - aux (valor_inicial) = - valor_inicial
        if (horizontal) {
        	temp->horizontal = 0;
#ifdef CLUSTALW
        	temp->vertical = valor_inicial; // Valor tal cual (es negativo)
#else
        	temp->vertical = - valor_inicial; // Diferencia (es positivo)
#endif
        } else {
#ifdef CLUSTALW
        	temp->horizontal = valor_inicial; // Valor tal cual (es negativo)
#else
        	temp->horizontal = - valor_inicial; // Diferencia (es positivo)
#endif
        	temp->vertical = 0;
        }
    }
 }

 /* Sustituye el elemento de indice i dentro del fragmento por el indicado en e.
  * Devuelve 0 si se ha podido susituir sin problemas, -1 en otro caso (posible fuera de rango).
  * NORMALMENTE NO UTILIZAREMOS ESTA FUNCIÓN porque lo hace todo más engorroso, con el control de
  * límites del array, pero ahí está. */
 int set_elemento_fragmento(struct Nodo_fragmento* fragmento, int i, struct Nodo_valor *e){
 	if (fragmento->longitud < i)
 		return -1;
 	fragmento->contenido[i] = *e;
 	return 0;
 }

 /* Obtiene el elemento de indice i dentro del fragmento, asignándolo a e.
  * Devuelve 0 si se ha encontrado sin problemas, -1 en otro caso (posible fuera de rango).
  * NORMALMENTE NO UTILIZAREMOS ESTA FUNCIÓN porque lo hace todo más engorroso, con el control de
  * límites del array, pero ahí está. */
 int get_elemento_fragmento(struct Nodo_fragmento* fragmento, int i, struct Nodo_valor *e){
 	if (fragmento->longitud < i)
 		return -1;
 	*e = fragmento->contenido[i];
 	return 0;
 }

 /* Obtiene el máximo valor de todo el fragmento (en campo valor), devoldiendo el índice en que lo encontró también */
 long get_max_fragmento(struct Nodo_fragmento* fragmento, int * posicion){
 	long max=LONG_MIN;
 	for (unsigned int i = 0; i<fragmento->longitud; i++){
 		long valor = fragmento->contenido[i].valor;
 		if (valor >= max) {
 			max = valor;
 			*posicion = i;
 		}
 	}
 	return max;
 }

 /* Devuelve el tamaño del fragmento. */
 int get_size(struct Nodo_fragmento* fragmento){
 	return fragmento->longitud;
 }

 /* Escribe por pantalla el contenido de un fragmento de trabajo. */
 void print_fragmento(struct Nodo_fragmento* fragmento){
  	printf("[");
  	for (unsigned int i = 0; i<fragmento->longitud; i++){
  		struct Nodo_valor valor = fragmento->contenido[i];
  		printf("(%d,%d,%d),", valor.valor, valor.horizontal, valor.vertical);
  	}
  	printf("]\n");
 }


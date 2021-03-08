/*
 * priority_queue.c
 *
 *  Created on: Mar 8, 2021
 *      Author: galvez
 */

#include <stdlib.h>
#include "priority_queue.h"

#define maximo(a, b) ((a>b)?a:b)
#define minimo(a, b) ((a>b)?b:a)

// La posición 0 no se usa en un heap implementado con arrays.
struct QueuePriorityRecord {
  int capacity;
  int size;
  IndiceCola *array;
};

int isEmptyQueuePrty(QueuePrty q) {
  return q->size == 0;
}

int isFullQueuePrty(QueuePrty q) {
  return q->size == q->capacity;
}

int getQueuePrtySize(QueuePrty q) {
  return q->size;
}

QueuePrty createQueuePrty(int maxElements) {
  QueuePrty q;

  if (maxElements < MinQueueSize) {
    Error("CreateQueue Error: Queue size is too small.");
  }

  q = malloc (sizeof(struct QueuePriorityRecord));
  if (q == NULL) {
    FatalError("CreateQueue Error: Unable to allocate more memory.");
  }

  // La posición 0 no se usa en un heap implementado con arrays.
  // Es preferible desperdiciarla a complicar el algoritmo.
  q->array = malloc( sizeof(IndiceCola) * (maxElements+1) );
  if (q->array == NULL) {
    FatalError("CreateQueue Error: Unable to allocate more memory.");
  }

  q->capacity = maxElements;
  makeEmptyQueuePrty(q);

  return q;
}

void makeEmptyQueuePrty(QueuePrty q) {
  q->size = 0;
}

void disposeQueuePrty(QueuePrty q) {
  if (q != NULL) {
    free(q->array);
    free(q);
  }
}

// Para calcular la altura de un heap basta con calcular la longitud de la rama más a la izquierda
int alturaQueuePrty(QueuePrty q, int posicion){
	int resultado=0;
	while(posicion<=q->size){ resultado++; posicion*=2; }
	return resultado;
}

// Un heap es lleno si la altura de su rama más a la izquierda coincide con la de su rama más a la derecha
// ya que siempre es completo.
int esLlenoQueuePrty(QueuePrty q, int posicion){
	int resultado=0;
	int alturaIzq = alturaQueuePrty(q, posicion);
	while(posicion<=q->size){ resultado++; posicion=posicion*2+1; }
	return alturaIzq==resultado;
}

void pushQueuePrty(QueuePrty q, IndiceCola x) {
  	//printf("Insertado %d - ", x);
  if (isFullQueuePrty(q)) {
    Error("Enqueue Error: The queue is full.");
  } else {
  	int actual=1;
  	while(actual<=q->size){
  		IndiceCola raiz=q->array[actual];
  		if (x < raiz){ // Se intercambian raíz y x
  			q->array[actual]=x;
  			x=raiz;
  		}
  		// Averiguar por donde debe seguir la inserción
  		if ( (alturaQueuePrty(q, actual*2) > alturaQueuePrty(q, actual*2+1)) && ! esLlenoQueuePrty(q, actual*2)){ // Insertar a la izquierda
  			actual*=2;
  		}else{ // Insertar por la derecha
  			actual=actual*2+1;
  		}
  	}
  	// En este momento actual es igual a size+1
  	q->array[++(q->size)]=x;
  	//printQueuePrty(q);
  }
}

IndiceCola headerQueuePrty(QueuePrty q) {
  if (!isEmptyQueuePrty(q)) {
    return q->array[1]; // Primer elemento del Heap
  }
  Error("Front Error: The queue is empty.");

  /* Return value to avoid warnings from the compiler */
  return 0;

}

void deleteHeaderQueuePrty(QueuePrty q) {
  if (isEmptyQueuePrty(q)) {
    Error("Dequeue Error: The queue is empty.");
  } else {
  	int actual=1;
  	IndiceCola raizActual;
  	// Ponemos como raiz el último elemento del heap
  	q->array[actual] = q->array[q->size--];
  	// Y ahora lo hundimos.
  	while(actual < q->size){
  		IndiceCola raizIzq, raizDch;
  		raizActual = q->array[actual];
  		if (actual*2 > q->size){
  			break;
  		}else{ // Hay hijo izquierdo ...
  			raizIzq=q->array[actual*2];
  			if (actual*2+1 > q->size){ // ... pero no hay hijo derecho.
  				// Dejamos el mas grande arriba y listo.
  				q->array[actual]=minimo(raizActual, raizIzq);
  				q->array[actual*2]=maximo(raizActual, raizIzq);
  				actual=actual*2;
  				break;
  			} else { //... y también hay hijo derecho
  				raizDch=q->array[actual*2+1];
  				// Si el más chico es la raíz se acabó
  				if (raizActual <= minimo(raizIzq, raizDch)){
  					break;
  				} else if (raizIzq <= raizDch) { // Tiramos por la izquierda
  					q->array[actual]=raizIzq;
  					actual=actual*2;
  					q->array[actual]=raizActual;
  				} else { // Tiramos por la derecha
  					q->array[actual]=raizDch;
  					actual=actual*2+1;
  					q->array[actual]=raizActual;
  				}
  			}
  		}
  	}
  }

}

IndiceCola popQueuePrty(QueuePrty q) {
  IndiceCola X = 0;

  if (isEmptyQueuePrty(q)) {
    Error("FrontAndDequeue Error: The queue is empty.");
  } else {
  	X = headerQueuePrty(q);
  	//printf("Cabeza: %d ", X);
  	deleteHeaderQueuePrty(q);
  	//printQueuePrty(q);
  }
  return X;
}


void printQueuePrty(QueuePrty q) {
  if (isEmptyQueuePrty(q)) {
    printf("The queue is empty.\n");
  } else {
  	printf("Size: %d ~ ", q->size);
  	printf("Queue: [");
  	int size = q->size;
  	for (int i = 1; i <= size; i++){
		printf("%d,",q->array[i]);
  	}
    printf("]\n");
  }
}


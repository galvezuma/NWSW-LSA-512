/*
 * priority_queue.h
 *
 *  Created on: Mar 8, 2021
 *      Author: galvez
 */

#ifndef PRIORITY_QUEUE_H_
#define PRIORITY_QUEUE_H_


/* COLA DE PRIORIDAD EN ARRAY QUE ALMACENARÁ LOS ÍNDICES DE LOS TRABAJOS PENDIENTES
 * (los trabajos del tipo [x,y] serán el trabajo y*MaxX + x). */
#include <stdio.h>
#include <stdlib.h>

#define Error(Str)        FatalError(Str)
#define FatalError(Str)   fprintf(stderr, "%s\n", Str), exit(1)

#define MinQueueSize (2) //Tamaño mínimo sin el cual no funcionará la cola (por los punteros a front y rear)

typedef int IndiceCola;

struct QueuePriorityRecord;
typedef struct QueuePriorityRecord *QueuePrty;
typedef QueuePrty Queue; //Para nosotros, con esta implementación, una cola será una cola de prioridad

// Todas las funciones de colas serán las de priority_queue.h
#define isEmptyQueuePrty is_empty
#define isFullQueuePrty is_full
#define create_queue createQueuePrty
#define dispose_queue disposeQueuePrty
#define make_empty makeEmptyQueuePrty
#define push pushQueuePrty
#define header headerQueuePrty
#define delete_header deleteHeaderQueuePrty
#define pop popQueuePrty
#define print_queue printQueuePrty

  int         isEmptyQueuePrty(QueuePrty Q);
  int         isFullQueuePrty(QueuePrty Q);
  int         getQueuePrtySize(QueuePrty Q);
  QueuePrty   createQueuePrty(int MaxElements);
  void        disposeQueuePrty(QueuePrty Q);
  void        makeEmptyQueuePrty(QueuePrty Q);
  void        pushQueuePrty(QueuePrty Q, IndiceCola X);
  IndiceCola  headerQueuePrty(QueuePrty Q);
  void        deleteHeaderQueuePrty(QueuePrty Q);
  IndiceCola  popQueuePrty(QueuePrty Q);
  void        printQueuePrty(QueuePrty Q);

#endif /* PRIORITY_QUEUE_H_ */

/*
 * job_table.h
 *
 *  Created on: Mar 8, 2021
 *      Author: galvez
 */

#ifndef JOB_TABLE_H_
#define JOB_TABLE_H_


#include "job_fragment.h"
#include "NWSW_LSA_512_controller.h"

/* Módulo de manipulación y tratamiento de la tabla de trabajos del sistema, de dimensiones:
 * [floor(long.sec A / SIZE_SUBTRABAJO_H), floor(long.secB % SIZE_SUBTRABAJO_V)]
 * Será única en el sistema, pero ya que las funciones serán usadas por todos los núcleos,
 * no podemos crear una instancia única estática en el módulo. */

/* Cada posición de la tabla tendrá un registro llamado "trabajo" compuesto por dos punteros a "fragmento_trabajo",
 * el primero será al fragmento horizontal del subtrabajo y el segundo al fragmento vertical. */

#define nombreFichero "TempTile"

struct args_salvarFila {
	 struct Nodo_trabajo ** tablaTrabajos;
	 int y;
	 int tabla_trabajos_MAX_X;
	 int cierra_conexion;
	 pthread_t * hilo_actual;
} ;

 extern struct Nodo_trabajo** crea_tabla_trabajos(int anchura, int altura);

 extern struct Nodo_trabajo** crea_tabla_trabajos_compartida(int anchura, int altura);

 extern void libera_tabla_trabajos(Nodo_trabajo *** ptr_tabla_trabajos, int anchura, int altura);

 extern struct Nodo_fragmento* is_fragmento_fila(struct Nodo_trabajo ** tabla_trabajos, int x, int y);

 extern struct Nodo_fragmento* get_fragmento_fila(struct Nodo_trabajo ** tabla_trabajos, int x, int y, int tabla_trabajos_MAX_X, int tabla_trabajos_MAX_Y);

 extern struct Nodo_fragmento* is_fragmento_columna(struct Nodo_trabajo ** tabla_trabajos, int x, int y);

 extern struct Nodo_fragmento* get_fragmento_columna(struct Nodo_trabajo ** tabla_trabajos, int x, int y, int tabla_trabajos_MAX_X, int tabla_trabajos_MAX_Y);

 extern void set_fragmento_fila(struct Nodo_trabajo ** tabla_trabajos, int x, int y, struct Nodo_fragmento* fragmento_fila, int tabla_trabajos_MAX_X);

 extern void set_fragmento_columna(struct Nodo_trabajo ** tabla_trabajos, int x, int y, struct Nodo_fragmento* fragmento_columna, int tabla_trabajos_MAX_X);

 extern  void * pthread_salvarFila(void * a);

 extern void salvarFila(struct Nodo_trabajo ** tablaTrabajos, int y, int tabla_trabajos_MAX_X, int cierra_conexion);

 extern void recuperarFila(struct Nodo_trabajo ** tablaTrabajos, int y, int tabla_trabajos_MAX_X);

 extern void liberarFila(struct Nodo_trabajo ** tablaTrabajos, int y, int tabla_trabajos_MAX_X);

 extern void imprime_estado_tabla_trabajos(struct Nodo_trabajo ** tabla_trabajos, int anchura, int altura);

 extern void escribe_tabla_fichero(char* file, struct Nodo_trabajo ** tabla_trabajos, int anchura, int altura);

 extern void imprime_fragmento_fichero(FILE** fichero, Nodo_fragmento* fragmento);

#endif /* JOB_TABLE_H_ */

/*
 * NWSW_LSA_512_controller.h
 *
 *  Created on: Mar 8, 2021
 *      Author: galvez
 */

#ifndef NWSW_LSA_512_CONTROLLER_H_
#define NWSW_LSA_512_CONTROLLER_H_


#include "job_table.h"
#include "priority_queue.h"
#include "alphabet.h"


#pragma GCC diagnostic ignored "-Wunused-variable"

 /** Única instancia estática de la tabla de trabajos del sistema (sólo hay un controlador). */
 static struct Nodo_trabajo ** tabla_trabajos;

 /** Única instancia estática de la cola del sistema para la gestión de tiles disponibles. */
 static Queue cola_procesos;
 /** Única instancia estática de la cola del sistema para la gestión de trabajos disponibles.
  * Un trabajo [x,y] será un entero del tipo y*MaxX + x. */
 static QueuePrty cola_trabajos;

 /** Variables globales estáticas para tamaños horizontal y vertical de la tabla de trabajos. */
 static int tabla_trabajos_MAX_X;
 static int tabla_trabajos_MAX_Y;

 /** Variable de estado del algoritmo: primera fase (FORWARD) o segunda fase retroceso (BACKWARD) */
 estado_algoritmo estado;
 /** Variable del tipo de ejecución del algoritmo: si global (NW) o local (SW). */
 static int global;
 /** Nombre del fichero de salida (será el prefijo de los archivos de la tabla temporal) */
 char * prefijo_fichero_temporal;

 /** Salvado del máximo del algoritmo (en cualquier punto en SW, en cualquiera de las últimas filas o columnas en NW) */
 long maximo_valor; //Valor máximo para alineamiento
 int maximo_indice_x; //Índice X del subtrabajo donde se encontró el máximo
 int maximo_indice_y; //Índice Y del subtrabajo donde se encontró el máximo
 int maximo_posicion_x; //Posición X donde se encontró el máximo en el subtrabajo
 int maximo_posicion_y; //Posición Y donde se encontró el máximo en el subtrabajo

 /** Variable global estática del grupo de trabajadores dependientes del controlador. */
 // static ilibGroup grupo_trabajadores;

 /** Variables globales estáticas que referencian a las secuencias a alinear (están en memoria compartida). */
 static char* secuencia_horizontal;
 static char* secuencia_vertical;

 static byte direccion_proxima;

 /** Variables globales estáticas propias del Needleman-Wunsch (en realidad son constantes sólo inicializadas una vez). */
 static short INSERT_COST;
 static short DELETE_COST;
 static long MATCHREPLACE_COST;
 static short GAPEXTEND_COST;

 /** Variables globales estáticas propias del algoritmo paralelo. */
 static int SIZE_SUBTRABAJO_H;
 static int SIZE_SUBTRABAJO_V;
 static int NUMERO_TILES;

 /* Variable global estática para control de la tabla de puntuaciones a usar (y sus alfabetos). */
 static struct Matriz_puntuaciones* tabla_puntuacion;

 /** Lanza un algoritmo Needleman-Wunsch paralelo, haciendo uso del proceso controlador del Needleman-Wunsch, que inicializará
 *  la matriz de resultados, colocará los trabajos y lanzará y coordinará a los núcleos, devolviendo las cadenas de alineamiento. */
 void MC64NWSW_Controlador(char * prefijo, int global_bool, estado_algoritmo stage, char* horizontal, char* vertical, Nodo_resultado * resultado,
	short valor_insert, short valor_delete, int valor_matchreplace, short valor_gapextend, char* matriz_puntuaciones,
	int valor_k_horizontal, int valor_k_vertical, int num_tiles);

 /** Encola un trabajo de coordenadas [x,y] a la cola de trabajos pendientes. */
 extern void encolar_trabajo(int x, int y);
 /** Desencola un trabajo de la cola de trabajos pendientes. */
 extern void desencolar_trabajo(int* x, int* y);
 /** Inicializa la tabla de trabajos y encola el primer trabajo, el [0,0]. */
 extern void inicializar_trabajos();
 /** Lanza al resto de procesadores y les asigna el código binario de "NW_worker". */
 extern void inicializar_hilos();
 /** Define el canal sumidero (sink) y establece como emisores del mismo a todos los procesos del grupo grupo_trabajadores. */
 extern void define_canal_sumidero();
 /** Bucle principal de la fase 1 que entrega y recibe trabajos entre los procesadores hasta que se reciba el último trabajo. */
 extern void coordinar_tiles();
 /** Envía el próximo trabajo disponible al primer núcleo de la cola de trabajos. */
 extern void enviar_trabajo(int limite_x, int limite_y);
 /** Envía un trabajo de índice x=-1, indicando así al núcleo que el algoritmo ha terminado y puede finalizar. */
 extern void enviar_finalizacion(int rango);
 /** Envía la señal de finalización a todos los núcleos del sistema. */
 extern void finalizar_hilos();
 /** Recibe un trabajo ya finalizado por parte de alguno de los núcleos, o espera a que éstos envíen dicho trabajo. */
 extern int recibir_trabajo();
 /** Bucle principal de la fase 2 que calcula las secuencias alineadas de atrás hacia delante coordinando a los núcleos. */
 extern void controller_calcula_alineamiento(char ** al_horizontal, char ** al_vertical, long * align_s1_start, long * align_s1_end, long * align_s2_start, long * align_s2_end);
 /** Recibe un fragmento de alineamiento ya calculado por parte de alguno de los núcleos, o espera a que éstos envíen dicho trabajo. */
 extern int recibir_alineamiento(char ** frag_aln_horizontal, char ** frag_aln_vertical, int * pos_x_final, int * pos_y_final, long * fin_alineamiento_x, long * fin_alineamiento_y);

 extern void inicializar_mutexes(int num_tiles);
 extern void finalizar_mutexes(int num_tiles);

 /** Escribe los datos de la tabla de trabajos en el fichero seleccionado. */
 extern void escribe_tabla_trabajos(char * fichero);

 /** Función auxiliar para el cálculo del mínimo entre dos números. */
 extern int min(int a, int b);
 /** Función auxiliar para el cálculo del máximo entre dos números. */
 extern int max(int a, int b);
 /** Función auxiliar que modifica una cadena y la convierte en su inversa. */
 extern void reverse(char *t);
 /** Función auxiliar que calcula ceil(k/total) en enteros para sacar el número de trabajosa partir
 *  del tamaño total y el subtamaño de cada uno. */
 extern int total_subtrabajos(int k, int total);

 /** Escribe los valores temporales necesarios para la segunda fase. */
 extern void salva_maximos_fichero(char* file_name, long maximo_valor, int maximo_indice_x, int maximo_indice_y, int maximo_posicion_x, int maximo_posicion_y, clock_t tiempo_inicio, clock_t tiempo_intermedio);
 /** Lee los valores temporales necesarios para la segunda fase. */
 void recupera_maximos_fichero(char * file_name, long *maximo_valor, int *maximo_indice_x, int *maximo_indice_y, int *maximo_posicion_x, int *maximo_posicion_y, clock_t *tiempo_inicio, clock_t *tiempo_intermedio);


#pragma GCC diagnostic warning "-Wunused-variable"

#endif /* NWSW_LSA_512_CONTROLLER_H_ */

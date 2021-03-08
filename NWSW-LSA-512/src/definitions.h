/*
 * definitions.h
 *
 *  Created on: Mar 8, 2021
 *      Author: galvez
 */

#ifndef DEFINITIONS_H_
#define DEFINITIONS_H_

#include <limits.h>
#include "MyLib.h"
#include <sys/time.h>
#include <sys/types.h>
#include <pthread.h>
//#include <tmc/cmem.h> //Sólo para la versión 2.0 del MDE

// Gestion del tiempo
// Para permitir contar hasta milisegundos
// nos quedamos con la primera parte del long y luego le ponemos los milisegundos
long my_clock();

//Esto lo haré en la línea del compilador gcc con -D SMITHWATERMAN
#define SMITHWATERMAN 1 //Para alineamiento local: cambios en NW_worker por SW_worker: recompilamos con esta directivaç

//#define CLUSTALW 1 // Si activamos la flag, se aplicarán los cambios de la versión ClustalW: uso de 12 bytes por registro en lugar de 6
// Se usarán 3 long en lugar de 1 long y 2 bytes, ya que la matriz de puntuaciones de ClustalW tiene valores > 128. ¿Resultado? Algoritmo más pesado y lento

/* Tamaño fragmentos de la matriz para dividir en subtrabajos. */
#define SIZE_SUBTRABAJO_H_DEFAULT 1000 //En un principio, cuadrados de 1000x1000
#define SIZE_SUBTRABAJO_V_DEFAULT 1000

/* Parametrización de constantes del algoritmo. */
#define NUMERO_TILES_DEFAULT 59 //Número de procesadores paralelos a utilizar como trabajadores (sin contar el controlador). Max=59 con G64wo: 3 shared, 1 dedicated)
//#define SHARED_TILES 3 // Constante asociada al bootrom del G64wo: 3 shared y 1 dedicated (no lo devuelve proc_remaining() QUE NO SE PUEDEN USAR (sólo 59 como máximo)
#define MENSAJE_TAG 210 //Tag de mensaje para el paso de mensajes entre trabajadores y controladors

#define SIZELINE_ALIGNMENT 50 //Longitud de la línea de un alineamiento a la hora de escribir
#define MAX_SIZEHEADER 400 //Longitud máxima de la cabecera de un fichero FASTA
#define MAX_SIZELINE 200 //Longitud máxima de una línea de secuencia de un fichero FASTA

/* Parámetros y elementos auxiliares constantes para el algoritmo Needleman-Wunsch. */
#define INSERT_COST_DEFAULT 4 // Coste de la operación de apertura de una inserción (hueco en la secuencia original [horizontal])
#define DELETE_COST_DEFAULT 4 // Coste de la operación de apertura de una delección (hueco en la secuencia con la que comparar [vertical])
#define MATCHREPLACE_COST_DEFAULT 1 // Puntuación por defecto en caso de que ambos nucleótidos coincidan, pero no estén en la matriz
#define GAPEXTEND_COST_DEFAULT 1 // Coste de la extensión de un hueco una vez abierta
//Eliminado para mayor optimización en acceso a matriz de puntuación: suponemos que ambos valdrán lo mismo
//#define REPLACE_COST 1 // Puntuación por defecto en caso de que ambos nucleótidos NO coincidan, pero no estén en la matriz

#define DEBUG 0 //Para escritura de tabla temporal y debug
//DEPRECATED: Debug guardaba la tabla, al estar ahora salvada en sistema de ficheros, sólo guarda la primera fila
#define TRAZA 0 // Para escritura de la traza de resultados para el Applet Java
#define TIEMPOS 0 // Para obtener los tiempos al comienzo, fin de primera fase y final
#define DIMENSION_TARJETA 8
#define CLOCKS_PER_CENTSEC 10 //(CLOCKS_PER_SEC/100), asumiendo que CLOCKS_PER_SEC = 1.000
#define MILLIS_PER_SEC 1000
#define byte unsigned char //Hacemos este reemplazo para usar un byte (char en C)

/* Cambios para arreglar los warnings de la versión del MDE 2.0 */
// Definición de la nueva función de malloc de la librería <tmc/cmem.h>
//#undef malloc_shared
//#define malloc_shared tmc_cmem_malloc

// Definición del FINCAD como caracter '\0' (el compilador dice que NULL era un puntero a void, no caracter)
#define FINCAD '\0'

#define DIRECCION_DIAGONAL 0
#define DIRECCION_HORIZONTAL 1
#define DIRECCION_VERTICAL 2

/* Definición del enumerado que determina la fase en la que se encuentra el algoritmo. */
 typedef enum {FORWARD, BACKWARD} estado_algoritmo;


//Definición del canal de la conversión
enum Channels {
	CANAL_SINK
};
//Macros del raw channel sink para puertos de entrada y salida
//ILIB_RAW_SEND_PORT(puerto_envio, 0); //Definimos el puerto de envío de cada proceso
//ILIB_RAW_RECEIVE_PORT(puerto_recepcion, 0); //Definimos el puerto de recepción del proceso 0

 /* Resultado del algoritmo NWunsch de ejecución completa: secuencias originales alineadas,
  * cada una en un array de char diferente. */
 typedef struct Alineamiento_resultado Alineamiento_resultado;

 struct Alineamiento_resultado {
   char * alineamiento_query;
   char * alineamiento_subject;
   long align_s1_start;
   long align_s1_end;
   long align_s2_start;
   long align_s2_end;
 };

 /* Estructura con los resultados del Needleman-Wunsch:
  * · 1pass y 2pass: Valor long de la puntuación total del alineamiento entre las secuencias.
  * · 2pass: Estructura con los alineamientos resultantes de las dos secuencias. */
 typedef struct Nodo_datos_resultado Nodo_datos_resultado;

 struct Nodo_datos_resultado {
   Alineamiento_resultado alineamientos;
   long puntuacion_total;
 };

/* Resultado final: compuesto de la estructura resultante (alineamiento y puntuación) y datos estadísticos. */
 typedef struct Nodo_resultado Nodo_resultado;

 struct Nodo_resultado {
   Nodo_datos_resultado datos;
   time_t inicio, intermedio, final;
 };

/* Definición de una matriz de puntuaciones, que almacenará dicha matriz, así como los alfabetos
 * horizontal y vertical que dan sentido a la misma (sus longitudes delimitarán el tamaño de la matriz)
 * y los bytes [0-128] que codificarán las secuencias a números para direccionamiento indirecto en la matriz. */
 typedef struct Matriz_puntuaciones Matriz_puntuaciones;

 struct Matriz_puntuaciones {
   int** matriz;
   char* alfabeto_horizontal;
   char* alfabeto_vertical;
   unsigned short * codificacion_horizontal;
   unsigned short * codificacion_vertical;
   int id; // Para uso en caso de que sea adaptativa
 };

 /* Definición de un valor, una estructura que almacenará el valor de las 3 matrices para NW:
 * · valor: valor long de la matriz natural de Needleman-Wunsch
 * · vertical: valor short de la versión AFFINE GAP de NW (valores verticales, matriz E en BioJava)
 * · horizontal: valor short de la versión AFFINE GAP de NW (valores horizontales, matriz D en BioJava)
 * */
 typedef struct Nodo_valor Nodo_valor;

/* OPTIMIZACIÓN DE MEMORIA: struct en C se alinea a 8 bytes (múltiplo del mayor, long, 4), desperdicio de un 33% de memoria
 * en todo momento: tabla de trabajos en memoria compartida, sistema de ficheros a la hora de salvarse.
 * Utilizando la macro pragma, obligamos a alinearlo a 6 bytes, ganamos un 33% de memoria pero algunas operaciones serán más lentas. */
//#ifdef CLUSTALW
 // SIN PRAGMA, ya que son longs
 struct Nodo_valor {
   int32_t valor;
   int32_t vertical; //Para la primera fase de ClustalW, cambiar a long
   int32_t horizontal; //Para la primera fase de ClustalW, cambiar a long
 };
//#endif
/*
 // antigua definici�n para Smith-Waterman. Ahora da igual que sea local o global
 #pragma pack(push,1)
 struct Nodo_valor {
   long valor;
   byte vertical;
   byte horizontal;
 };
 // FIN PRAGMA
 #pragma pack(pop)
*/

/* Definición de un fragmento, que será una estructura con significado semántico propio que almacenará:
 * · contenido: un array de enteros, el fragmento en sí mismo
 * · longitud: longitud del fragmento, ya que aunque será de un tamaño estándar la mayoría de las ocasiones, hay excepciones
 * */
 typedef struct Nodo_fragmento Nodo_fragmento;

 struct Nodo_fragmento {
   struct Nodo_valor* contenido;
   int longitud;
 };

 /* Definición de un trabajo, compuesto por dos punteros a "fragmento_trabajo",
  * el primero será al fragmento horizontal del subtrabajo y el segundo al fragmento vertical. */
 typedef struct Nodo_trabajo Nodo_trabajo;

 struct Nodo_trabajo {
   struct Nodo_fragmento* fragmento_horizontal;
   struct Nodo_fragmento* fragmento_vertical;
 };

 /* Definición de un nodo de trabajo que contiene información del trabajo al que pertenece en la tabla y
  * el worker que se encarga de dicho trabajo, además del trabajo en sí y las cadenas a alinear. */
 typedef struct Nodo_trabajo_completo Nodo_trabajo_completo;

 struct Nodo_trabajo_completo {
   struct Nodo_fragmento* fragmento_horizontal;
   struct Nodo_fragmento* fragmento_vertical;
   int x;
   int y;
   int indice_worker;
   char* secuencia_horizontal;
   char* secuencia_vertical;
   struct Matriz_puntuaciones* tabla_puntuacion;
   int limite_x;
   int limite_y;
   short INSERT_COST;
   short DELETE_COST;
   short GAPEXTEND_COST;
   int SIZE_SUBTRABAJO_H;
   int SIZE_SUBTRABAJO_V;
   long valor_maximo; //Valor máximo para alineamiento local
   int valor_maximo_x; //Posición X donde se encontró: uso en SW
   int valor_maximo_y; //Posición Y donde se encontró: uso en SW
   byte direccion; // De uso en fase 2
 };

 /* Definición de un nodo de alineamiento, que contiene fragmentos del
  * alineamiento para las secuencias horizontal y vertical. */
 typedef struct Nodo_alineamiento Nodo_alineamiento;

 struct Nodo_alineamiento {
   char * fragmento_horizontal;
   char * fragmento_vertical;
   int x;
   int y;
   int indice_worker;
   int limite_x;
   int limite_y;
   int valor_minimo_x; //Posición X donde termina el alineamiento: uso en SW
   int valor_minimo_y; //Posición Y donde termina el alineamiento: uso en SW
   byte direccion; // De uso sólo en fase 2
 };

 /* Definicion de mecanismos de comunicacion
  * En lugar de canales se usa memoria compartida
  */
pthread_mutex_t sumidero_mutex;
struct {
	union {
		Nodo_alineamiento* alineamiento;
		Nodo_trabajo_completo * trabajo;
	} datos;
	int sender;
	unsigned char datos_disponibles; // booleano
} sumidero;
pthread_cond_t sumidero_datos_disponibles_cond;
pthread_cond_t sumidero_libre_cond;

pthread_mutex_t *peer_to_peer_mutex; // Array de "canales" punto a punto
typedef struct {
	union {
		Nodo_trabajo_completo trabajo;
	} datos;
	int sender;
	unsigned char datos_disponibles; // booleano
} peer_to_peer_struct;
peer_to_peer_struct * peer_to_peer; // Array
pthread_cond_t * peer_to_peer_datos_disponibles_cond; // Array

 /** Para depuraciones y gestión memoria en caso de error: escribe la utilización de la memoria. */
 #define imprime_estado_memoria() { \
    malloc_stats(); \
    struct mallinfo info = mallinfo(); \
    unsigned long memoria_total = info.arena + info.usmblks; \
    printf("arena: %lu, ordblks: %u, uordblks: %u, fordblks: %u, keepcost: %u\n", \
    	memoria_total, info.ordblks, info.uordblks, info.fordblks, info.keepcost); \
    printf("Memory in use: %u bytes\n", \
        info.usmblks + info.uordblks); \
    printf("Total heap size: %u bytes\n", info.arena); \
 }

/* Arreglo de bug del malloc en Tilera no fixeado al hacer mallocs de memoria de tamaño cercano a 64k y 128k. No ocurre con 32, 256 y 512 ni con malloc_shared() */
 extern void * malloc_safe_agr248(int size);
 extern void * malloc_shared(int size);
 extern void * memdup(void * ptr, int size);

#endif /* DEFINITIONS_H_ */

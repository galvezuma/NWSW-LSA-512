/*
 * NWSW_LSA_512_worker.c
 *
 *  Created on: Mar 8, 2021
 *      Author: galvez
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <immintrin.h>
#include "MyLib.h"
#include "job_fragment.h"
#include "alphabet.h"

#define RANGO_PADRE 0 //Ya que es el único proceso de ILIB_GROUP_PARENT

#define VECTOR_16
#define VECTORIZE
#ifdef VECTOR_16
	#define int_per_simdreg 16 // Número de enteros que cabe en un registro vector
#endif
#ifdef VECTOR_8
	#define int_per_simdreg 8 // Número de enteros que cabe en un registro vector
#endif
#ifdef VECTOR_4
	#define int_per_simdreg 4 // Número de enteros que cabe en un registro vector
#endif

/* Parámetros del algoritmo NW (en realidad serán constantes a todas las ejecuciones) */
static short INSERT_COST;
static short DELETE_COST;
static short GAPEXTEND_COST;
/* Parámetros del algoritmo paralelo (en realidad serán constantes a todas las ejecuciones) */
static int SIZE_SUBTRABAJO_H;
static int SIZE_SUBTRABAJO_V;

/* Parametros del estado del algoritmo */
static estado_algoritmo estado;

/* NOTA IMPORTANTE:
 * No estamos trabajando con orientacion a objetos.
 * A esto se une que trabajamos con memoria compartida y no con procesos independientes
 * Ello nos lleva a que todos los threads comparten las mismas variables, a diferencia de lo que sucedia
 * con Tilera, donde las variables globales eran, en realidad, locales a cada proceso.
 * Para facilitar el trabajo con threads sin cambiar demasiado la estructura del programa
 * se ha decidido crear una estructura que engloba a todas las variables de cada proceso, de manera que cada
 * thread hara un malloc de dicha estructura y trabajara con la suya propia.
 */
typedef struct {

	/* Fragmentos de trabajos iniciales y finales: todos en memoria local para mayor velocidad acceso */
	struct Nodo_fragmento* trozoSubColumnaPartida;
	struct Nodo_fragmento* trozoSubFilaPartida;
	struct Nodo_fragmento* trozoSubColumnaTrabajo;
	struct Nodo_fragmento* trozoSubFilaTrabajo;
	int indice_trabajo_x;
	int indice_trabajo_y;
	int limite_x;
	int limite_y;
	long valor_maximo_inicial; // Guardamos el máximo
	long valor_maximo;
	char* secuencia_horizontal; //Guardamos referencias a las secuencias para poder hacer malloc/free
	char* secuencia_vertical;
	struct Matriz_puntuaciones* tabla_puntuacion;
	int tabla_puntuacion_id; // ID de la tabla actual. SI cambia al recibir un trabajo, es que es adaptativa
	char* alineamiento_horizontal;
	char* alineamiento_vertical;
	byte direccion_inicial;

	#ifdef SMITHWATERMAN
	/* Posiciones del subtrabajo donde se encuentra el último máximo (sólo para SW) y el punto de finalización del alineamiento. */
	int pos_max_x;
	int pos_max_y;
	int pos_min_x;
	int pos_min_y;
	#endif
} estado_del_Objeto;

void init_estado_del_Objeto(estado_del_Objeto * Obj){
	Obj->trozoSubColumnaPartida=NULL;
	Obj->trozoSubFilaPartida=NULL;
	Obj->trozoSubColumnaTrabajo=NULL;
	Obj->trozoSubFilaTrabajo=NULL;
	Obj->indice_trabajo_x=0;
	Obj->indice_trabajo_y=0;
	Obj->limite_x=0;
	Obj->limite_y=0;
	Obj->valor_maximo_inicial=0; // Guardamos el máximo
	Obj->valor_maximo=0;
	Obj->secuencia_horizontal=NULL; //Guardamos referencias a las secuencias para poder hacer malloc/free
	Obj->secuencia_vertical=NULL;
	Obj->tabla_puntuacion=NULL;
	Obj->tabla_puntuacion_id=0; // ID de la tabla actual. SI cambia al recibir un trabajo, es que es adaptativa
	Obj->alineamiento_horizontal=NULL;
	Obj->alineamiento_vertical=NULL;
	Obj->direccion_inicial=0;
	#ifdef SMITHWATERMAN
		Obj->pos_max_x=-1;
		Obj->pos_max_y=-1;
		Obj->pos_min_x=-1;
		Obj->pos_min_y=-1;
	#endif
}

/** Recepción del trabajo a realizar por paso de mensajes desde el controlador. */
int recepcion_trabajo(int rango_tile, estado_del_Objeto * Obj);

/** Entrega del trabajo ya finalizado al controlador por rawchannel en sumidero. */
void entregar_trabajo(int rango_tile, estado_del_Objeto * Obj);

/** Realiza el procesamiento de un trabajo: algoritmo NW en una submatriz utilizando sólo una fila y una columna. */
void realiza_procesamiento(estado_del_Objeto * Obj);
/** Procesamiento de un trabajo por columnas. */
void realiza_procesamiento_columnas(estado_del_Objeto * Obj);
void realiza_procesamiento_columnasSIMD(estado_del_Objeto * Obj);
/** Hace el cálculo de una columna de ancho int_per_simdreg, de forma vectorizada */
void vectorizacionColumnas(estado_del_Objeto * Obj, struct Nodo_valor* contenidoSubColumnaPartida, struct Nodo_valor* contenidoSubColumnaTrabajo, struct Nodo_valor * r_diag_a, struct Nodo_valor * r_diag_b, char * letraSecuenciaHorizontal);
/** Procesamiento de un trabajo por filas. */
void realiza_procesamiento_filas(estado_del_Objeto * Obj);

/** Cálculo de la submatriz para la fase 2: en función de los límites será calculada parcialmente. */
struct Nodo_valor** calcula_submatriz(estado_del_Objeto * Obj, int anchura, int altura);

/** Liberación de memoria de la submatriz para la fase 2. */
void libera_submatriz(struct Nodo_valor*** ptr_matriz, int anchura);

/** Cálculo de las secuencias de alineamiento a partir de una submatriz y la esquina inferior derecha de la misma. */
void calcula_alineamiento(estado_del_Objeto * Obj, struct Nodo_valor** matriz, int num_filas, int num_columnas);

/** Entrega del alineamiento ya calculado al controlador por rawchannel en sumidero. */
void entregar_alineamiento(int rango_tile, estado_del_Objeto * Obj);

/** Calcula el valor de puntuación para los dos caracteres en función de la matriz de sustituciones. */
inline long calcula_matchReplace(estado_del_Objeto * Obj, unsigned char letraSecuenciaHorizontal, unsigned char letraSecuenciaVertical);

/** AUXILIARES */
/** Devuelve el máximo entre 2 números long. */
inline long max2(long a, long b);
/** Devuelve el máximo entre 3 números long. */
inline long max3(long a, long b, long c);

/** Devuelve el índice que es el máximo entre 3 números. */
inline int max_indice(long a, long b, long c);
/*
//#define calcula_valores(diagonal, arriba, izquierda, maximo) { \
	long result_vert, result_hor; \
	// El que está arriba, hueco en la original -> inserción o sigo haciendo hueco  \
    result_vert = arriba->valor + max2(-arriba->vertical, -INSERT_COST) - GAPEXTEND_COST; \
    // El que está a la izquierda, hueco en la otra -> delección o sigo haciendo hueco  \
    result_hor = izquierda->valor + max2(-izquierda->horizontal, -DELETE_COST) - GAPEXTEND_COST; \
    // Calculamos el valor real entre INSERT/DELETE o REPLACE  \
    maximo->valor = max( \
    	result_vert, \
    	result_hor, \
    	diagonal->valor + calcula_matchReplace(letraSecuenciaHorizontal, letraSecuenciaVertical) \
    ); \
    maximo->vertical = maximo->valor - result_vert; \
    maximo->horizontal = maximo->valor - result_hor; \
}
*/
#define calcula_valores(diagonal, arriba, izquierda, maximo) { \
	int32_t result_vert, result_hor; \
	/* El que está arriba, hueco en la original -> inserción o sigo haciendo hueco */ \
    result_vert = arriba->valor - arriba->vertical; \
    /* El que está a la izquierda, hueco en la otra -> delección o sigo haciendo hueco */ \
    result_hor = izquierda->valor - izquierda->horizontal; \
    /* Calculamos el valor real entre INSERT/DELETE o REPLACE */ \
    maximo->valor = max3( \
    	result_vert, \
    	result_hor, \
    	diagonal->valor + calcula_matchReplace(Obj, letraSecuenciaHorizontal, letraSecuenciaVertical) \
    ); \
    maximo->valor = max2(maximo->valor, 0); \
    maximo->vertical = maximo->valor - (max2(result_vert, maximo->valor - INSERT_COST) - GAPEXTEND_COST); \
    maximo->horizontal = maximo->valor - (max2(result_hor, maximo->valor - DELETE_COST) - GAPEXTEND_COST); \
}
/* MUY IMPORTANTE: Como ClustalW necesita usar longs en lugar de bytes porque las diferencias pueden ser mayor de 128,
 * podemos ahorrarnos el guardar las diferencias y guardaremos el valor en sí, algo mucho más óptimo. Por ello, hay dos macros. */
#ifdef CLUSTALW

#define calcula_valores(diagonal, arriba, izquierda, maximo) { \
    /* Calculamos el valor real entre INSERT/DELETE o REPLACE */ \
    maximo->valor = max3( \
    	arriba->vertical, \
    	izquierda->horizontal, \
    	diagonal->valor + calcula_matchReplace(Obj, letraSecuenciaHorizontal, letraSecuenciaVertical) \
    ); \
    maximo->vertical = max2(arriba->vertical, maximo->valor - INSERT_COST) - GAPEXTEND_COST;  \
    maximo->horizontal = max2(izquierda->horizontal, maximo->valor - DELETE_COST) - GAPEXTEND_COST; \
}

#endif
/*
#ifdef SMITHWATERMAN
#undef calcula_valores
#define calcula_valores(diagonal, arriba, izquierda, maximo) { \
	long result_vert, result_hor; \
	 El que está arriba, hueco en la original -> inserción o sigo haciendo hueco  \
    result_vert = arriba->valor - arriba->vertical; \
     El que está a la izquierda, hueco en la otra -> delección o sigo haciendo hueco  \
    result_hor = izquierda->valor - izquierda->horizontal; \
     Calculamos el valor real entre INSERT/DELETE o REPLACE  \
    maximo->valor = max3( \
    	result_vert, \
    	result_hor, \
    	diagonal->valor + calcula_matchReplace(Obj, letraSecuenciaHorizontal, letraSecuenciaVertical) \
    ); \
     En Smith-Waterman no puede haber negativos  \
    maximo->valor = max2(maximo->valor, 0); \
    maximo->vertical = maximo->valor - (max2(result_vert, maximo->valor - INSERT_COST) - GAPEXTEND_COST); \
    maximo->horizontal = maximo->valor - (max2(result_hor, maximo->valor - DELETE_COST) - GAPEXTEND_COST); \
}

#endif
*/

//NUEVA VECTORIZACION  **********************************************************************


//*************************************************************************************************





/*
void vectorizacionColumnas(	estado_del_Objeto * Obj,
							struct Nodo_valor* contenidoSubColumnaPartida,
							struct Nodo_valor* contenidoSubColumnaTrabajo,
							struct Nodo_valor * r_diag_a,
							struct Nodo_valor * r_diag_b,
							char * letraSecuenciaHorizontal){
	__m128i r_up_valor __attribute__((aligned (16)));
	__m128i r_up_vertical __attribute__((aligned (16)));
	__m128i r_left_valor __attribute__((aligned (16)));
	__m128i r_left_horizontal __attribute__((aligned (16)));
	__m128i r_res_vertical __attribute__((aligned (16)));
	__m128i r_res_horizontal __attribute__((aligned (16)));
	__m128i r_max_valor __attribute__((aligned (16)));
	__m128i r_diag_valor __attribute__((aligned (16)));
	__m128i r_matchReplace __attribute__((aligned (16)));
	__m128i r_res_diag __attribute__((aligned (16)));
	__m128i r_max_vertical __attribute__((aligned (16)));
	__m128i r_max_horizontal __attribute__((aligned (16)));
	__m128i r_const_insertCost __attribute__((aligned (16)));
	__m128i r_const_deleteCost __attribute__((aligned (16)));
	__m128i r_const_gapExtend __attribute__((aligned (16)));
	#ifdef SMITHWATERMAN
		__m128i r_const_zero __attribute__((aligned (16))) = _mm_set_epi32(0,0,0,0);
	#endif
	//Trabajaremos con los arrays para mayor sencillez de acceso
	int size_subcolumna_partida = get_size(Obj->trozoSubColumnaPartida);
	// Inicialización de constantes
	//
	r_const_deleteCost = _mm_set_epi32(DELETE_COST,DELETE_COST,DELETE_COST,DELETE_COST);
	r_const_insertCost = _mm_set_epi32(INSERT_COST,INSERT_COST,INSERT_COST,INSERT_COST);
	r_const_gapExtend = _mm_set_epi32(GAPEXTEND_COST,GAPEXTEND_COST,GAPEXTEND_COST,GAPEXTEND_COST);


	// Inicialización de valores de entrada
	//
	r_up_valor = _mm_set_epi32(r_diag_b[4].valor,r_diag_b[3].valor,r_diag_b[2].valor,r_diag_b[1].valor); // reverse order
	r_up_vertical = _mm_set_epi32(r_diag_b[4].vertical,r_diag_b[3].vertical,r_diag_b[2].vertical,r_diag_b[1].vertical); // reverse order
	//
	r_left_valor = _mm_set_epi32(r_diag_b[3].valor,r_diag_b[2].valor,r_diag_b[1].valor,r_diag_b[0].valor); // reverse order
	r_left_horizontal = _mm_set_epi32(r_diag_b[3].horizontal,r_diag_b[2].horizontal,r_diag_b[1].horizontal,r_diag_b[0].horizontal); // reverse order
	//
	r_diag_valor =  _mm_set_epi32(r_diag_a[3].valor,r_diag_a[2].valor,r_diag_a[1].valor,r_diag_a[0].valor); // reverse order

	// Comienzo del bucle
	for(int y=int_per_simdreg - 1; y<size_subcolumna_partida-1; y++) {

		//A PARTIR DE AQUI HAY MODIFICACIONES************************************************************************************





		//r_res_vertical = _mm_sub_epi32(r_up_valor, r_up_vertical);
		//r_res_horizontal = _mm_sub_epi32(r_left_valor, r_left_horizontal);
		//


		r_max_valor = _mm_max_epi32(r_up_vertical, r_left_horizontal);


		r_matchReplace = _mm_set_epi32(	calcula_matchReplace(Obj, letraSecuenciaHorizontal[3], Obj->secuencia_vertical[y-3]),
										calcula_matchReplace(Obj, letraSecuenciaHorizontal[2], Obj->secuencia_vertical[y-2]),
										calcula_matchReplace(Obj, letraSecuenciaHorizontal[1], Obj->secuencia_vertical[y-1]),
										calcula_matchReplace(Obj, letraSecuenciaHorizontal[0], Obj->secuencia_vertical[y])
				); // reverse order

		r_res_diag = _mm_add_epi32(r_diag_valor, r_matchReplace);
		//

		r_max_valor = _mm_max_epi32(r_max_valor, r_res_diag);




#ifdef DEBUG_VECTOR
		printf("%c %c, %c %c, %c %c, %c %c\n", letraSecuenciaHorizontal[0],Obj->secuencia_vertical[y],letraSecuenciaHorizontal[1],Obj->secuencia_vertical[y-1],letraSecuenciaHorizontal[2],Obj->secuencia_vertical[y-2],letraSecuenciaHorizontal[3],Obj->secuencia_vertical[y-3]);
		print_reg("MatchReplace ", r_matchReplace);
		print_reg("Valor ", r_max_valor);
#endif
		#ifdef SMITHWATERMAN
			r_max_valor = _mm_max_epi32(r_max_valor, r_const_zero);
		#endif
		//
		//
		r_max_vertical = _mm_sub_epi32(r_max_valor, r_const_insertCost);
		r_max_vertical = _mm_max_epi32(r_max_vertical, r_up_vertical);//UP
		r_max_vertical = _mm_sub_epi32(r_max_vertical, r_const_gapExtend);
		//r_max_vertical = _mm_sub_epi32(r_max_valor, r_max_vertical);

		//
		//
		r_max_horizontal = _mm_sub_epi32(r_max_valor, r_const_deleteCost);//maximoValor-delete
		r_max_horizontal = _mm_max_epi32(r_max_horizontal, r_left_horizontal);//MAX(izqda_horizontal,(maximoValor-delete))
		r_max_horizontal = _mm_sub_epi32(r_max_horizontal, r_const_gapExtend);
		//r_max_horizontal = _mm_sub_epi32(r_max_valor, r_max_horizontal);



//		printf("%d\n", y);
//		print_reg("Max valor ", r_max_valor);
//		print_reg("Max horiz ", r_max_horizontal);
//		print_reg("Max vertical ", r_max_vertical);
		//
		//
		//Guardamos los valores extremos
		contenidoSubColumnaTrabajo[y+2 - int_per_simdreg].valor = _mm_extract_epi32(r_max_valor, 3);
		contenidoSubColumnaTrabajo[y+2 - int_per_simdreg].vertical = _mm_extract_epi32(r_max_vertical, 3);
		contenidoSubColumnaTrabajo[y+2 - int_per_simdreg].horizontal = _mm_extract_epi32(r_max_horizontal, 3);

		// *****************************************************************************************
		//HASTA AQUI HECHO, (es hasta donde habia q llegar)

#ifdef DEBUG_VECTOR
		print_nodo("Guardo ", contenidoSubColumnaTrabajo[y+2 - int_per_simdreg], y+2 - int_per_simdreg);
		printf("SubColumna %d: valor:%ld vertical:%d horizontal:%d \n",y+2 - int_per_simdreg,contenidoSubColumnaTrabajo[y+2 - int_per_simdreg].valor,contenidoSubColumnaTrabajo[y+2 - int_per_simdreg].vertical, contenidoSubColumnaTrabajo[y+2 - int_per_simdreg].horizontal);
		printf("%ld\n",sizeof(contenidoSubColumnaTrabajo[y+2 - int_per_simdreg].horizontal));
#endif
		#ifdef SMITHWATERMAN
			// TODO Ver si hemos obtenido un nuevo máximo
			int local_max = horizontal_max_Vec4i(r_max_valor);
			if (local_max > Obj->valor_maximo){
				int valor[int_per_simdreg] __attribute__((aligned (16)));
				_mm_stream_si128 ((__m128i *)valor, r_max_valor);
				for(int diag=1;diag<=int_per_simdreg;diag++){ // Reverse order
					r_diag_b[diag].valor=valor[int_per_simdreg - diag];
				}

			}
		#endif
		// Preparamos al siguiente bucle
		// r_max se convierte en r_up
		r_up_valor = r_max_valor;
		r_up_vertical = r_max_vertical;
		// r_left se convierte en r_diag
		r_diag_valor = r_left_valor;
		// r_max desplazado se convierte en r_left
//		printf("Antes\n");
		// TODO Al finalizar el bucle se hace un acceso ilegar a memoria (indexOutOfbounds)
		r_left_valor = _mm_slli_si128(r_max_valor, sizeof(int));
		r_left_horizontal = _mm_slli_si128(r_max_horizontal, sizeof(int));
		if (y+2<size_subcolumna_partida){
			r_left_valor = _mm_insert_epi32(r_left_valor, contenidoSubColumnaPartida[y+2].valor, 0);
			r_left_horizontal = _mm_insert_epi32(r_left_horizontal, contenidoSubColumnaPartida[y+2].horizontal, 0);
		}
//		printf("Despues\n");
#ifdef DEBUG_VECTOR
		print_nodo("Parto de valor ",contenidoSubColumnaPartida[y+2], y+2);
#endif
	}
	// Al acabar el proceso vectorial debemos reobtener r_diag_a y r_diag_b para que finalice el cálculo en la función llamante
	// Guardar r_max en r_diag_b
	int valor[int_per_simdreg] __attribute__((aligned (16)));
	int vertical[int_per_simdreg] __attribute__((aligned (16)));
	int horizontal[int_per_simdreg] __attribute__((aligned (16)));
	_mm_stream_si128 ((__m128i *)valor, r_max_valor);
	_mm_stream_si128 ((__m128i *)vertical, r_max_vertical);
	_mm_stream_si128 ((__m128i *)horizontal, r_max_horizontal);
	for(int diag=1;diag<=int_per_simdreg;diag++){ // Reverse order
		r_diag_b[diag].valor=valor[diag-1];//int_per_simdreg - diag];
		r_diag_b[diag].horizontal=horizontal[diag-1];//int_per_simdreg - diag];
		r_diag_b[diag].vertical=vertical[diag-1];//int_per_simdreg - diag];
	}
	// En r_diag_a sólo necesitamos el valor y los elementos 1ÃÂº y último tampoco nos hacen falta
	_mm_stream_si128 ((__m128i *)valor, r_diag_valor);
	for(int diag=1;diag<int_per_simdreg;diag++){ // Reverse order
		r_diag_a[diag].valor=valor[diag];// diag-1 //int_per_simdreg - diag];
	}
}


*/







//***************************************************************************************************


//END NUEVA VECTORIZACION  *******************************************************************




/** Cada procesador estará en espera de recibir un trabajo y, una vez recibido, ejecutará el procesamiento.
 *  Posteriormente, se entregará el trabajo resultante al núcleo controlador y se repetirá el ciclo. */
void * worker_main(void * param){
	/* Creacion del estado interno del worker que hay que pasar de funcion a funcion */
	//fprintf(stderr, "Entrando en un hilo.\n");
	estado_del_Objeto Obj;
	init_estado_del_Objeto(&Obj);
	/* Constante de posición del core en la tajeta */
	int rango_tile;
	char param_estado;
	mylibThreadParam * t_params = (mylibThreadParam *)param;
	//printf("Me han creado con un %c y soy %d\n", t_params->estado, t_params->num);
	// Inicializamos la variable global de almacenamiento de posición/rango
	rango_tile = t_params->num;
	param_estado = t_params->estado;
	// Liberamos la memoria de parametrizacion creada ad hoc.
	free(param);
	t_params=param=NULL;
	/* long secondsMIN = 99, usecondsMIN = 999999, secondsMAX = 0, usecondsMAX = 0; */
	mylib_init();


	//Abrimos el canal sumidero por el lado del emisor
  	//ilib_rawchan_open_sender(CANAL_SINK, puerto_envio);

  	// Sacamos los datos de los argumentos
  	if (param_estado == 0x00) {
    	mylib_die("Error al recibir los parámetros en el Worker. El parámetro no puede ser nulo.\n");
	}
  	if (param_estado == '1')
  		estado = FORWARD;
  	else
  		estado = BACKWARD;

	// En función del estado, el tile se encargará de una u otra tarea, y cuando acabe la fase, FINALIZARÃï¿½ SIEMPRE, PARA PODER LIBERAR TILES
	if (estado == FORWARD) {
	/* FASE 1: Mientras recibamos trabajos no finales, calculamos matrices y
	 * guardamos sólo fila y columna final, que entregamos al controlador. */
		while (recepcion_trabajo(rango_tile, &Obj)){
			/* struct timeval tiempo_inicial, tiempo;
			long seconds, useconds;
			gettimeofday(&tiempo_inicial, NULL); */
	//		clock_t inicio, final;

	//		inicio = clock();
	//		printf("1. Worker %d comenzó trabajo [%d, %d].\n", rango_tile, indice_trabajo_x, indice_trabajo_y);

			//Recibo trabajo por paso de mensajes
			realiza_procesamiento(&Obj);

	/*		gettimeofday(&tiempo, NULL);
			seconds  = tiempo.tv_sec  - tiempo_inicial.tv_sec;
	    	useconds = tiempo.tv_usec - tiempo_inicial.tv_usec;
	    	if (seconds*1000000+useconds > secondsMAX*1000000+usecondsMAX) { secondsMAX=seconds; usecondsMAX=useconds; }
	    	if (seconds*1000000+useconds < secondsMIN*1000000+usecondsMIN) { secondsMIN=seconds; usecondsMIN=useconds; } */

	//		printf("2. Worker %d terminó trabajo [%d, %d]. Esperando a enviar\n", rango_tile, indice_trabajo_x, indice_trabajo_y);

			//Envío trabajo por sink
			entregar_trabajo(rango_tile, &Obj);
	//    	final = clock();
	//   		printf("3. Worker %d terminó trabajo [%d, %d]. DUR: %.2f\n", rango_tile, indice_trabajo_x, indice_trabajo_y, (float)(final - inicio)/CLOCKS_PER_SEC);

	/*	 	printf("fin trabajo MAX: %ld.%06ld \n", secondsMAX, usecondsMAX);
		printf("fin trabajo MIN: %ld.%06ld \n", secondsMIN, usecondsMIN); */
		}
	} else { // estado = BACKWARD
	/* FASE 2: Mientras recibamos trabajos no finales, calculamos submatrices y
	 * guardamos todo. Con estas submatrices, calculamos el alineamiento en sentido
	 * inverso y lo enviamos al controlador. */
		struct Nodo_valor** submatriz;
		int num_filas, num_columnas;
		printf("En fase de retroceso y soy %d\n", rango_tile);

		while (recepcion_trabajo(rango_tile, &Obj)){
			//Recibo trabajo por paso de mensajes
			//printf("Hola soy el núcleo %d con trabajo_fase2 [%d, %d]\n", rango_tile, indice_trabajo_x, indice_trabajo_y);

			//Control de límites: si es 0, es completo
			num_columnas = ((Obj.limite_x>0) ? Obj.limite_x : get_size(Obj.trozoSubFilaPartida));
			num_filas = ((Obj.limite_y>0) ? Obj.limite_y : get_size(Obj.trozoSubColumnaPartida));

			submatriz = calcula_submatriz(&Obj, num_columnas, num_filas);
			calcula_alineamiento(&Obj, submatriz, num_filas, num_columnas);
			libera_submatriz(&submatriz, num_columnas);

			//Envío trabajo por sink
			entregar_alineamiento(rango_tile, &Obj);
		}
	}

	if (Obj.tabla_puntuacion != NULL) // En caso de que el worker no haya trabajado nunca (raro, pero puede ocurrir)
		matrix_free(&(Obj.tabla_puntuacion));

	//esto no libera realmente el tile, es inútil: ilib_group_free(ILIB_GROUP_SIBLINGS);
	mylib_finish();

	return 0;
}

/** Recepción del trabajo a realizar por paso de mensajes desde el controlador.
 *  Devuelve 1 si se ha recibido trabajo y 0 en caso de recibir trabajo.x==-1, para finalizar al núcleo. */
int recepcion_trabajo(int rango_tile, estado_del_Objeto * Obj){
	int retorno;
	//Creamos el nodo: no es necesario en memoria compartida puesto que se envía el buffer completo
	Nodo_trabajo_completo * trabajo;
	//ilibStatus status;
	//Recibimos el mensaje con los dos punteros a columna y fila iniciales del proceso padre (grupo PARENTS, rango 0, tag MENSAJE_TAG)
	//Es una recepción bloqueante, por lo que mientras no haya tabajo o espere en la cola, estará detenido
	//if (ilib_msg_receive(ILIB_GROUP_PARENT, RANGO_PADRE, MENSAJE_TAG, &trabajo, sizeof(trabajo), &status) != ILIB_SUCCESS){
	//	ilib_die("Error crítico, el trabajador %d no ha podido recoger el trabajo.", rango_tile);
	//}

	pthread_mutex_lock(peer_to_peer_mutex+rango_tile);
//	 printf("\tBloqueo worker para recibir %d\n", rango_tile);
		while(!peer_to_peer[rango_tile].datos_disponibles){
			pthread_cond_wait(peer_to_peer_datos_disponibles_cond+rango_tile, peer_to_peer_mutex+rango_tile);
		}
	trabajo = &(peer_to_peer[rango_tile].datos.trabajo);
//	 printf("\tHe recibido algo %d\n", rango_tile);
	if (trabajo->x==-1) {
		retorno= 0;
	} else {
		Obj->indice_trabajo_x = trabajo->x;
		Obj->indice_trabajo_y = trabajo->y;
		Obj->limite_x = trabajo->limite_x;
		Obj->limite_y = trabajo->limite_y;
		Obj->valor_maximo_inicial = trabajo->valor_maximo;
		Obj->valor_maximo = Obj->valor_maximo_inicial;

		INSERT_COST = trabajo->INSERT_COST;
		DELETE_COST = trabajo->DELETE_COST;
		GAPEXTEND_COST = trabajo->GAPEXTEND_COST;
		SIZE_SUBTRABAJO_H = trabajo->SIZE_SUBTRABAJO_H;
		SIZE_SUBTRABAJO_V = trabajo->SIZE_SUBTRABAJO_V;
		Obj->direccion_inicial = trabajo->direccion;
		/* Hacemos copias de los fragmentos a memoria local para un acceso óptimo. */
		copia_fragmento(&(Obj->trozoSubColumnaPartida), trabajo->fragmento_vertical);
		copia_fragmento(&(Obj->trozoSubFilaPartida), trabajo->fragmento_horizontal);
		/* Ya que sólo necesitaremos un fragmento de las secuencias originales:
		 * - Horizontal: [indice_trabajo_x*SIZE_SUBTRABAJO_H, indice_trabajo_x*SIZE_SUBTRABAJO_H + tam_fila]
		 * - Vertical: [indice_trabajo_y*SIZE_SUBTRABAJO_V, indice_trabajo_y*SIZE_SUBTRABAJO_V + tam_col]
		 * Y el acceso a memoria compartida es mucho más lento, copiaremos a local los fragmentos necesarios,
		 * acelerando enormemente todos los accesos. */
		int num_columnas = ((Obj->limite_x>0) ? Obj->limite_x : get_size(Obj->trozoSubFilaPartida));
		Obj->secuencia_horizontal = (char *) malloc (sizeof(char) * (num_columnas + 1));
		unsigned int i;
		for (i = 0; i < num_columnas; i++)
			Obj->secuencia_horizontal[i] = trabajo->secuencia_horizontal[Obj->indice_trabajo_x*SIZE_SUBTRABAJO_H + i];
		Obj->secuencia_horizontal[i] = FINCAD; // \000 para poder escribir y calcular longitud

		int num_filas = ((Obj->limite_y>0) ? Obj->limite_y : get_size(Obj->trozoSubColumnaPartida));
		Obj->secuencia_vertical = (char *) malloc (sizeof(char) * (num_filas + 1));
		for (i = 0; i < num_filas; i++)
			Obj->secuencia_vertical[i] = trabajo->secuencia_vertical[Obj->indice_trabajo_y*SIZE_SUBTRABAJO_V + i];
		Obj->secuencia_vertical[i] = FINCAD; // \000 para poder escribir y calcular longitud

		/* Ya que los accesos a memoria compartida son MUCHO más lentos y hacemos un continuo uso de la matriz,
		 * realizaremos un copia de la estructura, pues le sacaremos rentabilidad.
		 * Además, sólo la copiaremos la primera vez o cuando haya cambiado el ID de la matriz (adaptativa). */
		 if ((Obj->tabla_puntuacion == NULL) || (Obj->tabla_puntuacion_id != trabajo->tabla_puntuacion->id)) {
		 	 if (Obj->tabla_puntuacion != NULL)
		 	 	matrix_free(&(Obj->tabla_puntuacion)); // Liberamos la anterior matriz
		 	 matrix_copy(&(Obj->tabla_puntuacion), trabajo->tabla_puntuacion);
		 	 Obj->tabla_puntuacion_id = trabajo->tabla_puntuacion->id;
		 }

		/* Acceso en memoria compartida: mucho más ineficiente
		secuencia_horizontal = trabajo.secuencia_horizontal;
		secuencia_vertical = trabajo.secuencia_vertical; */

		// Uso desde memoria compartida: mucho más ineficiente
		//trozoSubColumnaPartida = trabajo.fragmento_vertical;
		//trozoSubFilaPartida = trabajo.fragmento_horizontal;

		 // Acceso en memoria compartida: mucho más ineficiente
//		tabla_puntuacion = trabajo.tabla_puntuacion;
		retorno= 1;
	}
	peer_to_peer[rango_tile].datos_disponibles=0;
	pthread_mutex_unlock(peer_to_peer_mutex+rango_tile);
//	 printf("\tDesbloqueo worker para recibir %d\n", rango_tile);
	return retorno;
}

/** Entrega un puntero al trabajo ya finalizado (en memoria compartida) por rawchannel sumidero al controlador. */
void entregar_trabajo(int rango_tile, estado_del_Objeto * Obj){
	//Creamos el nodo en memoria compartida para que sea accesible por el controlador
	Nodo_trabajo_completo* trabajo = (struct Nodo_trabajo_completo *)malloc_shared(sizeof (struct Nodo_trabajo_completo));


	//Asignamos la fila y columna final, los índices del trabajo que era y el índice del trabajador
	trabajo->x = Obj->indice_trabajo_x;
	trabajo->y = Obj->indice_trabajo_y;
	trabajo->indice_worker = rango_tile;
	//Si no había un máximo nuevo, devolvemos LONG_MIN (el mínimo long) para que no se guarde que estaba aquí
	trabajo->valor_maximo = ((Obj->valor_maximo == Obj->valor_maximo_inicial) ? LONG_MIN : Obj->valor_maximo);
	#ifdef SMITHWATERMAN
		trabajo->valor_maximo_x = Obj->pos_max_x;
		trabajo->valor_maximo_y = Obj->pos_max_y;
	#endif
	// Realizamos las copias a memoria compartida para entregar como resultado
	copia_fragmento_shared(&(trabajo->fragmento_vertical), Obj->trozoSubColumnaTrabajo);
	copia_fragmento_shared(&(trabajo->fragmento_horizontal), Obj->trozoSubFilaTrabajo);

	/* Limpiamos la memoria de todo aquello que hubiÃÂ©semos copiado a memoria local */
	// Fragmentos de secuencias en memoria local
	free(Obj->secuencia_horizontal);
	free(Obj->secuencia_vertical);

	// Los 4 fragmentos en memoria local
	libera_fragmento(&(Obj->trozoSubColumnaPartida));
 	libera_fragmento(&(Obj->trozoSubFilaPartida));
	libera_fragmento(&(Obj->trozoSubColumnaTrabajo));
 	libera_fragmento(&(Obj->trozoSubFilaTrabajo));

	//Enviamos los resultados al núcleo (es un envío NO bloqueante, si el controlador está ocupado, se queda en su buffer)
	/* Hay que esperar a que el controlador este preparado para recibir por el sumidero */
	pthread_mutex_lock(&sumidero_mutex);
	// printf("\tBloqueo worker para enviar %d\n", rango_tile);
		while(sumidero.datos_disponibles){
			pthread_cond_wait(&sumidero_libre_cond, &sumidero_mutex);
		}
	sumidero.datos.trabajo = trabajo;
	sumidero.datos_disponibles=1;
	pthread_cond_signal(&sumidero_datos_disponibles_cond);
	pthread_mutex_unlock(&sumidero_mutex);
//	 printf("\tDesbloqueo worker para enviar %d\n", rango_tile);
	//ilib_rawchan_send_pointer(puerto_envio, trabajo);
}


/** Realiza el procesamiento de uno de los trabajadores del NW paralelo: a partir de una subfila y una subcolumna de partida,
 *  va realizando operaciones de comparación para ir generando la tabla (sin guardarse) y llegar a la última fila y columna,
 *  que se almacenarán en la tabla de trabajos y abrirán las puertas a otros. */
void realiza_procesamiento(estado_del_Objeto * Obj) {
	//En realidad en esta primera fase siempre se calculan los fragmentos completos, pero
	//damos la posibilidad de que en futuro no fuese así teniendo una implementación para
	//filas y columnas y un control del límites
	if (Obj->limite_x >= 0) {
		//realiza_procesamiento_columnas(Obj);
		realiza_procesamiento_columnasSIMD(Obj);
	} else { //hasta cierta fila
		realiza_procesamiento_filas(Obj);
	}
}



void print_nodo(char *s, struct Nodo_valor r, int pos){
	printf("%s pos(%d) %d %d %d\n", s, pos, r.valor, r.horizontal, r.vertical);

}

/**  EL CALCULO SE REALIZARÁ POR COLUMNAS, calculando en cada iteración una posición de la última subfila y al final la última subcolumna. */
void realiza_procesamiento_columnas(estado_del_Objeto * Obj) {
	//Guardamos los valores enteros para evitar tener que recalcular más veces
	int size_subcolumna_partida = get_size(Obj->trozoSubColumnaPartida);
	int size_subfila_partida = get_size(Obj->trozoSubFilaPartida);
	#ifdef SMITHWATERMAN
		// Inicializamos las posiciones del máximo
		Obj->pos_max_x = Obj->pos_max_y = -1;
	#endif

	// Traducimos las letras para ir más rápidos
		// Inhabilitamos temporalmente el mensaje de advertencia: array subscript has type ‘char’
#pragma GCC diagnostic push
	#pragma GCC diagnostic ignored "-Wchar-subscripts"
	for(int i=0; i < size_subfila_partida; i++)
		Obj->secuencia_horizontal[i] = Obj->tabla_puntuacion->codificacion_horizontal[Obj->secuencia_horizontal[i]];
	for(int j=0; j < size_subcolumna_partida; j++)
		Obj->secuencia_vertical[j] = Obj->tabla_puntuacion->codificacion_vertical[Obj->secuencia_vertical[j]];
#pragma GCC diagnostic pop

	// Creamos los fragmentos de trabajo en memoria local
	Obj->trozoSubColumnaTrabajo = crea_fragmento(size_subcolumna_partida);
	Obj->trozoSubFilaTrabajo = crea_fragmento(size_subfila_partida);

	//Trabajaremos con los arrays para mayor sencillez de acceso
	struct Nodo_valor* contenidoSubColumnaTrabajo = Obj->trozoSubColumnaTrabajo->contenido;
	struct Nodo_valor* contenidoSubFilaTrabajo = Obj->trozoSubFilaTrabajo->contenido;
	struct Nodo_valor* contenidoSubColumnaPartida = Obj->trozoSubColumnaPartida->contenido;
	struct Nodo_valor* contenidoSubFilaPartida = Obj->trozoSubFilaPartida->contenido;

	/* Trabajaremos por columnas: optimizando el número de columnas a calcular si no necesitamos todo */
	int columnas_a_calcular = ((Obj->limite_x>0) ? Obj->limite_x : size_subfila_partida);

 	//Punteros a nodo para mayor optimización: ahorro de acceso a array con offset para cada valor de la estructura
	register struct Nodo_valor *temp, *temp_izq, *temp_arr, *temp_diag, *aux;
    /* Cálculo del contenido */
    //contenidoSubFilaTrabajo[0] = contenidoSubColumnaPartida[size_subcolumna_partida-1];
    copia_valor_nodo(contenidoSubFilaTrabajo[0], contenidoSubColumnaPartida[size_subcolumna_partida-1]);
    for (int x=1; x < columnas_a_calcular; x++){
    	//Leemos las letras de la secuencia H a partir del índice que nos corresponde
    	char letraSecuenciaHorizontal = Obj->secuencia_horizontal[x-1];
//    	printf("Columna %d %c\n",x, letraSecuenciaHorizontal);
    	/* Calcular trozoSubColumnaTrabajo a partir de trozoSubColumnaPartida */
    	// contenidoSubColumnaTrabajo[0] = contenidoSubFilaPartida[x]
    	copia_valor_nodo(contenidoSubColumnaTrabajo[0], contenidoSubFilaPartida[x]);

    	for (int y=1; y < size_subcolumna_partida; y++){
    		//Leemos las letras de la secuencia V a partir del índice que nos corresponde
    		char letraSecuenciaVertical = Obj->secuencia_vertical[y-1];

    		temp = &(contenidoSubColumnaTrabajo[y]);
    		temp_arr = &(contenidoSubColumnaTrabajo[y-1]);
    		temp_izq = &(contenidoSubColumnaPartida[y]);
    		temp_diag = &(contenidoSubColumnaPartida[y-1]);

    		// Calculamos los valores haciendo uso de una macro
    		calcula_valores(temp_diag, temp_arr, temp_izq, temp);

    		//Sólo si estamos en Smith-Waterman: vamos recalculando el máximo
    		#ifdef SMITHWATERMAN
    			if (contenidoSubColumnaTrabajo[y].valor > Obj->valor_maximo) {
    				Obj->valor_maximo = contenidoSubColumnaTrabajo[y].valor;
    				Obj->pos_max_x = x;
    				Obj->pos_max_y = y;
    			}
    		#endif

//    	    	print_nodo("Nodo ",*temp, y);
    	}
    	// contenidoSubFilaTrabajo[x] = contenidoSubColumnaTrabajo[size_subcolumna_partida - 1]
    	copia_valor_nodo(contenidoSubFilaTrabajo[x], contenidoSubColumnaTrabajo[size_subcolumna_partida - 1]);
        /* Nos preparamos para la siguien	te iteración evitando tomar más memoria.
         * Para ello se hace un intercambio de arrays que hay que deshacer tras el último ciclo. */
        aux = contenidoSubColumnaPartida;
        contenidoSubColumnaPartida = contenidoSubColumnaTrabajo;
        contenidoSubColumnaTrabajo = aux;
    }
    // El último cambio hay que deshacerlo (se sale de los límites de la matriz en la última iteración)
    contenidoSubColumnaTrabajo = contenidoSubColumnaPartida;



    //Sólo si estamos en Needleman-Wunsch: el máximo será la esquina inferior derecha
    #ifndef SMITHWATERMAN
    if (contenidoSubColumnaTrabajo[size_subcolumna_partida - 1].valor > Obj->valor_maximo) {
    	Obj->valor_maximo = contenidoSubColumnaTrabajo[size_subcolumna_partida - 1].valor;
    }
    #endif

    // Ya que no garantizamos que apunte realmente a su array (dependiendo de si es par o impar, HAY QUE COPIAR LOS VALORES)
    for (int y=0; y < size_subcolumna_partida; y++){
    	copia_valor_nodo(Obj->trozoSubColumnaTrabajo->contenido[y], contenidoSubColumnaTrabajo[y]);
    }

//	printf("Es normal: Fila\n");

}

#ifdef VECTORIZE

#ifdef VECTOR_4
// Vectorización 4

void print_reg(char *s, __m128i r){
	int valor[int_per_simdreg] __attribute__((aligned (16)));
	_mm_stream_si128 ((__m128i *)valor, r);
	printf("%s", s);
	for(int i=0; i<int_per_simdreg; i++){
		printf("%d ", valor[i]);
	}
	printf("\n");
}


int horizontal_max_Vec4i(__m128i x) {
    __m128i max1 = _mm_shuffle_epi32(x, _MM_SHUFFLE(0,0,3,2));
    __m128i max2 = _mm_max_epi32(x,max1);
    __m128i max3 = _mm_shuffle_epi32(max2, _MM_SHUFFLE(0,0,0,1));
    __m128i max4 = _mm_max_epi32(max2,max3);
    return _mm_cvtsi128_si32(max4);
}




void vectorizacionColumnas(	estado_del_Objeto * Obj,
							struct Nodo_valor* contenidoSubColumnaPartida,
							struct Nodo_valor* contenidoSubColumnaTrabajo,
							struct Nodo_valor * r_diag_a,
							struct Nodo_valor * r_diag_b,
							char * letraSecuenciaHorizontal){
	__m128i r_up_valor __attribute__((aligned (16)));
	__m128i r_up_vertical __attribute__((aligned (16)));
	__m128i r_left_valor __attribute__((aligned (16)));
	__m128i r_left_horizontal __attribute__((aligned (16)));
	__m128i r_res_vertical __attribute__((aligned (16)));
	__m128i r_res_horizontal __attribute__((aligned (16)));
	__m128i r_max_valor __attribute__((aligned (16)));
	__m128i r_diag_valor __attribute__((aligned (16)));
	__m128i r_matchReplace __attribute__((aligned (16)));
	__m128i r_res_diag __attribute__((aligned (16)));
	__m128i r_max_vertical __attribute__((aligned (16)));
	__m128i r_max_horizontal __attribute__((aligned (16)));
	__m128i r_const_insertCost __attribute__((aligned (16)));
	__m128i r_const_deleteCost __attribute__((aligned (16)));
	__m128i r_const_gapExtend __attribute__((aligned (16)));
	__m128i r_const_zero __attribute__((aligned (16))) = _mm_set_epi32(0,0,0,0);
#ifdef SMITHWATERMAN
	__m128i r_maximum __attribute__((aligned (32))) = _mm_set_epi32(0,0,0,0);
#endif

	//Trabajaremos con los arrays para mayor sencillez de acceso
	int size_subcolumna_partida = get_size(Obj->trozoSubColumnaPartida);
	// Inicialización de constantes
	//
	r_const_deleteCost = _mm_set_epi32(DELETE_COST,DELETE_COST,DELETE_COST,DELETE_COST);
	r_const_insertCost = _mm_set_epi32(INSERT_COST,INSERT_COST,INSERT_COST,INSERT_COST);
	r_const_gapExtend = _mm_set_epi32(GAPEXTEND_COST,GAPEXTEND_COST,GAPEXTEND_COST,GAPEXTEND_COST);


	// Inicialización de valores de entrada
	//
	r_up_valor = _mm_set_epi32(r_diag_b[4].valor,r_diag_b[3].valor,r_diag_b[2].valor,r_diag_b[1].valor); // reverse order
	r_up_vertical = _mm_set_epi32(r_diag_b[4].vertical,r_diag_b[3].vertical,r_diag_b[2].vertical,r_diag_b[1].vertical); // reverse order
	//
	r_left_valor = _mm_set_epi32(r_diag_b[3].valor,r_diag_b[2].valor,r_diag_b[1].valor,r_diag_b[0].valor); // reverse order
	r_left_horizontal = _mm_set_epi32(r_diag_b[3].horizontal,r_diag_b[2].horizontal,r_diag_b[1].horizontal,r_diag_b[0].horizontal); // reverse order
	//
	r_diag_valor =  _mm_set_epi32(r_diag_a[3].valor,r_diag_a[2].valor,r_diag_a[1].valor,r_diag_a[0].valor); // reverse order

	// Comienzo del bucle
//	struct timespec tstart={0,0}, tend={0,0};
//	    clock_gettime(CLOCK_MONOTONIC, &tstart);
	for(int y=int_per_simdreg - 1; y<size_subcolumna_partida-1; y++) {

		r_res_vertical = _mm_sub_epi32(r_up_valor, r_up_vertical);
		r_res_horizontal = _mm_sub_epi32(r_left_valor, r_left_horizontal);
		//
		r_max_valor = _mm_max_epi32(r_res_vertical, r_res_horizontal);
		r_matchReplace = _mm_set_epi32(	calcula_matchReplace(Obj, letraSecuenciaHorizontal[3], Obj->secuencia_vertical[y-3]),
										calcula_matchReplace(Obj, letraSecuenciaHorizontal[2], Obj->secuencia_vertical[y-2]),
										calcula_matchReplace(Obj, letraSecuenciaHorizontal[1], Obj->secuencia_vertical[y-1]),
										calcula_matchReplace(Obj, letraSecuenciaHorizontal[0], Obj->secuencia_vertical[y])
				); // reverse order

		r_res_diag = _mm_add_epi32(r_diag_valor, r_matchReplace);
		//
		r_max_valor = _mm_max_epi32(r_max_valor, r_res_diag);
#ifdef DEBUG_VECTOR
		printf("%c %c, %c %c, %c %c, %c %c\n", letraSecuenciaHorizontal[0],Obj->secuencia_vertical[y],letraSecuenciaHorizontal[1],Obj->secuencia_vertical[y-1],letraSecuenciaHorizontal[2],Obj->secuencia_vertical[y-2],letraSecuenciaHorizontal[3],Obj->secuencia_vertical[y-3]);
		print_reg("MatchReplace ", r_matchReplace);
		print_reg("Valor ", r_max_valor);
#endif
		//#ifdef SMITHWATERMAN
			r_max_valor = _mm_max_epi32(r_max_valor, r_const_zero);
		//#endif
		//
		//
		r_max_vertical = _mm_sub_epi32(r_max_valor, r_const_insertCost);
		r_max_vertical = _mm_max_epi32(r_max_vertical, r_res_vertical);
		r_max_vertical = _mm_sub_epi32(r_max_vertical, r_const_gapExtend);
		r_max_vertical = _mm_sub_epi32(r_max_valor, r_max_vertical);
		//
		//
		r_max_horizontal = _mm_sub_epi32(r_max_valor, r_const_deleteCost);
		r_max_horizontal = _mm_max_epi32(r_max_horizontal, r_res_horizontal);
		r_max_horizontal = _mm_sub_epi32(r_max_horizontal, r_const_gapExtend);
		r_max_horizontal = _mm_sub_epi32(r_max_valor, r_max_horizontal);

//		printf("%d\n", y);
//		print_reg("Max valor ", r_max_valor);
//		print_reg("Max horiz ", r_max_horizontal);
//		print_reg("Max vertical ", r_max_vertical);
		//
		//
		//Guardamos los valores extremos
		contenidoSubColumnaTrabajo[y+2 - int_per_simdreg].valor = _mm_extract_epi32(r_max_valor, 3);
		contenidoSubColumnaTrabajo[y+2 - int_per_simdreg].vertical = _mm_extract_epi32(r_max_vertical, 3);
		contenidoSubColumnaTrabajo[y+2 - int_per_simdreg].horizontal = _mm_extract_epi32(r_max_horizontal, 3);
#ifdef DEBUG_VECTOR
		print_nodo("Guardo ", contenidoSubColumnaTrabajo[y+2 - int_per_simdreg], y+2 - int_per_simdreg);
		printf("SubColumna %d: valor:%ld vertical:%d horizontal:%d \n",y+2 - int_per_simdreg,contenidoSubColumnaTrabajo[y+2 - int_per_simdreg].valor,contenidoSubColumnaTrabajo[y+2 - int_per_simdreg].vertical, contenidoSubColumnaTrabajo[y+2 - int_per_simdreg].horizontal);
		printf("%ld\n",sizeof(contenidoSubColumnaTrabajo[y+2 - int_per_simdreg].horizontal));
#endif
#ifdef SMITHWATERMAN
		r_maximum = _mm_max_epi32(r_maximum, r_max_valor);
#endif
		// Preparamos al siguiente bucle
		// r_max se convierte en r_up
		r_up_valor = r_max_valor;
		r_up_vertical = r_max_vertical;
		// r_left se convierte en r_diag
		r_diag_valor = r_left_valor;
		// r_max desplazado se convierte en r_left
//		printf("Antes\n");
		// TODO Al finalizar el bucle se hace un acceso ilegar a memoria (indexOutOfbounds)
		r_left_valor = _mm_slli_si128(r_max_valor, sizeof(int));
		r_left_horizontal = _mm_slli_si128(r_max_horizontal, sizeof(int));
		if (y+2<size_subcolumna_partida){
			r_left_valor = _mm_insert_epi32(r_left_valor, contenidoSubColumnaPartida[y+2].valor, 0);
			r_left_horizontal = _mm_insert_epi32(r_left_horizontal, contenidoSubColumnaPartida[y+2].horizontal, 0);
		}
//		printf("Despues\n");
#ifdef DEBUG_VECTOR
		print_nodo("Parto de valor ",contenidoSubColumnaPartida[y+2], y+2);
#endif
	}

//    clock_gettime(CLOCK_MONOTONIC, &tend);
//    printf("some_long_computation took about %.9f seconds\n",
//           ((double)tend.tv_sec + 1.0e-9*tend.tv_nsec) -
//           ((double)tstart.tv_sec + 1.0e-9*tstart.tv_nsec));

	// Al acabar el proceso vectorial debemos reobtener r_diag_a y r_diag_b para que finalice el cálculo en la función llamante
	// Guardar r_max en r_diag_b
	int valor[int_per_simdreg] __attribute__((aligned (16)));
	int vertical[int_per_simdreg] __attribute__((aligned (16)));
	int horizontal[int_per_simdreg] __attribute__((aligned (16)));
	_mm_stream_si128 ((__m128i *)valor, r_max_valor);
	_mm_stream_si128 ((__m128i *)vertical, r_max_vertical);
	_mm_stream_si128 ((__m128i *)horizontal, r_max_horizontal);
	for(int diag=1;diag<=int_per_simdreg;diag++){ // Reverse order
		r_diag_b[diag].valor=valor[diag-1];//int_per_simdreg - diag];
		r_diag_b[diag].horizontal=horizontal[diag-1];//int_per_simdreg - diag];
		r_diag_b[diag].vertical=vertical[diag-1];//int_per_simdreg - diag];
	}
	// En r_diag_a sólo necesitamos el valor y los elementos 1ÃÂº y último tampoco nos hacen falta
	_mm_stream_si128 ((__m128i *)valor, r_diag_valor);
	for(int diag=1;diag<int_per_simdreg;diag++){ // Reverse order
		r_diag_a[diag].valor=valor[diag];// diag-1 //int_per_simdreg - diag];
	}
#ifdef SMITHWATERMAN
	// Miramos si hay un nuevo máximo
	int local_max = horizontal_max_Vec4i(r_maximum);
	if (local_max > Obj->valor_maximo){
		Obj->valor_maximo = local_max;
	}
#endif
}

// END Vectorización 4
#endif

// Vectorización 8
#ifdef VECTOR_8

void print_reg(char *s, __m256i r){
	int valor[int_per_simdreg] __attribute__((aligned (32)));
	__m256i mascara __attribute__((aligned (32)));
	mascara = _mm256_set_epi32((unsigned int)0xFFFFFFFFFFFFFFFF, (unsigned int)0xFFFFFFFFFFFFFFFF, (unsigned int)0xFFFFFFFFFFFFFFFF, (unsigned int)0xFFFFFFFFFFFFFFFF, (unsigned int)0xFFFFFFFFFFFFFFFF, (unsigned int)0xFFFFFFFFFFFFFFFF, (unsigned int)0xFFFFFFFFFFFFFFFF, (unsigned int)0xFFFFFFFFFFFFFFFF);
	_mm256_maskstore_epi32 (valor, mascara, r);
	printf("%s", s);
	for(int i=0; i<int_per_simdreg; i++){
		printf("%d ", valor[i]);
	}
	printf("\n");
}


int horizontal_max_Vec8i(__m256i x) {
	__m256i max1 = _mm256_shuffle_epi32(x, _MM_SHUFFLE(0,0,3,2));
	__m256i max2 = _mm256_max_epi32(x,max1);
	__m256i max3 = _mm256_shuffle_epi32(max2, _MM_SHUFFLE(0,0,0,1));
	__m256i max4 = _mm256_max_epi32(max2,max3);
	//
	__m256i	max5 = _mm256_permute4x64_epi64(max4, _MM_SHUFFLE(0,0,0,2));
	__m256i max6 = _mm256_max_epi32(max4,max5);
	return _mm256_extract_epi32(max6, 0);
}




void vectorizacionColumnas(	estado_del_Objeto * Obj,
							struct Nodo_valor* contenidoSubColumnaPartida,
							struct Nodo_valor* contenidoSubColumnaTrabajo,
							struct Nodo_valor * r_diag_a,
							struct Nodo_valor * r_diag_b,
							char * letraSecuenciaHorizontal){
	__m256i r_up_valor __attribute__((aligned (32)));
	__m256i r_up_vertical __attribute__((aligned (32)));
	__m256i r_left_valor __attribute__((aligned (32)));
	__m256i r_left_horizontal __attribute__((aligned (32)));
	__m256i r_res_vertical __attribute__((aligned (32)));
	__m256i r_res_horizontal __attribute__((aligned (32)));
	__m256i r_max_valor __attribute__((aligned (32)));
	__m256i r_diag_valor __attribute__((aligned (32)));
	__m256i r_matchReplace __attribute__((aligned (32)));
	__m256i r_res_diag __attribute__((aligned (32)));
	__m256i r_max_vertical __attribute__((aligned (32)));
	__m256i r_max_horizontal __attribute__((aligned (32)));
	__m256i r_const_insertCost __attribute__((aligned (32)));
	__m256i r_const_deleteCost __attribute__((aligned (32)));
	__m256i r_const_gapExtend __attribute__((aligned (32)));
	__m256i r_const_permutacion __attribute__((aligned (32)));
	__m256i r_const_mascara __attribute__((aligned (32)));
	__m256i r_const_zero __attribute__((aligned (32))) = _mm256_set_epi32(0,0,0,0,0,0,0,0);
#ifdef SMITHWATERMAN
	__m256i r_maximum __attribute__((aligned (32))) = _mm256_set_epi32(0,0,0,0,0,0,0,0);
#endif

		//printf("Entrando...%ld\n", letraSecuenciaHorizontal);

	//Trabajaremos con los arrays para mayor sencillez de acceso
	int size_subcolumna_partida = get_size(Obj->trozoSubColumnaPartida);
	// Inicialización de constantes
	//
	r_const_deleteCost = _mm256_set_epi32(DELETE_COST,DELETE_COST,DELETE_COST,DELETE_COST,DELETE_COST,DELETE_COST,DELETE_COST,DELETE_COST);
	r_const_insertCost = _mm256_set_epi32(INSERT_COST,INSERT_COST,INSERT_COST,INSERT_COST,INSERT_COST,INSERT_COST,INSERT_COST,INSERT_COST);
	r_const_gapExtend = _mm256_set_epi32(GAPEXTEND_COST,GAPEXTEND_COST,GAPEXTEND_COST,GAPEXTEND_COST,GAPEXTEND_COST,GAPEXTEND_COST,GAPEXTEND_COST,GAPEXTEND_COST);
	r_const_permutacion = _mm256_set_epi32(6,5,4,3,2,1,0,0);
	r_const_mascara = _mm256_set_epi32((unsigned int)0xFFFFFFFFFFFFFFFF, (unsigned int)0xFFFFFFFFFFFFFFFF, (unsigned int)0xFFFFFFFFFFFFFFFF, (unsigned int)0xFFFFFFFFFFFFFFFF, (unsigned int)0xFFFFFFFFFFFFFFFF, (unsigned int)0xFFFFFFFFFFFFFFFF, (unsigned int)0xFFFFFFFFFFFFFFFF, (unsigned int)0xFFFFFFFFFFFFFFFF);

	// Inicialización de valores de entrada
	//
	r_up_valor = _mm256_set_epi32(r_diag_b[8].valor,r_diag_b[7].valor,r_diag_b[6].valor,r_diag_b[5].valor,r_diag_b[4].valor,r_diag_b[3].valor,r_diag_b[2].valor,r_diag_b[1].valor); // reverse order
	r_up_vertical = _mm256_set_epi32(r_diag_b[8].vertical,r_diag_b[7].vertical,r_diag_b[6].vertical,r_diag_b[5].vertical,r_diag_b[4].vertical,r_diag_b[3].vertical,r_diag_b[2].vertical,r_diag_b[1].vertical); // reverse order
	//
	r_left_valor = _mm256_set_epi32(r_diag_b[7].valor,r_diag_b[6].valor,r_diag_b[5].valor,r_diag_b[4].valor,r_diag_b[3].valor,r_diag_b[2].valor,r_diag_b[1].valor,r_diag_b[0].valor); // reverse order
	r_left_horizontal = _mm256_set_epi32(r_diag_b[7].horizontal,r_diag_b[6].horizontal,r_diag_b[5].horizontal,r_diag_b[4].horizontal,r_diag_b[3].horizontal,r_diag_b[2].horizontal,r_diag_b[1].horizontal,r_diag_b[0].horizontal); // reverse order
	//
	r_diag_valor =  _mm256_set_epi32(r_diag_a[7].valor,r_diag_a[6].valor,r_diag_a[5].valor,r_diag_a[4].valor,r_diag_a[3].valor,r_diag_a[2].valor,r_diag_a[1].valor,r_diag_a[0].valor); // reverse order


	// Comienzo del bucle
//	struct timespec tstart={0,0}, tend={0,0};
//	    clock_gettime(CLOCK_MONOTONIC, &tstart);
	for(int y=int_per_simdreg - 1; y<size_subcolumna_partida-1; y++) {

		r_res_vertical = _mm256_sub_epi32(r_up_valor, r_up_vertical);
		r_res_horizontal = _mm256_sub_epi32(r_left_valor, r_left_horizontal);
		//
		r_max_valor = _mm256_max_epi32(r_res_vertical, r_res_horizontal);
		r_matchReplace = _mm256_set_epi32(
											calcula_matchReplace(Obj, letraSecuenciaHorizontal[7], Obj->secuencia_vertical[y-7]),
											calcula_matchReplace(Obj, letraSecuenciaHorizontal[6], Obj->secuencia_vertical[y-6]),
											calcula_matchReplace(Obj, letraSecuenciaHorizontal[5], Obj->secuencia_vertical[y-5]),
											calcula_matchReplace(Obj, letraSecuenciaHorizontal[4], Obj->secuencia_vertical[y-4]),
											calcula_matchReplace(Obj, letraSecuenciaHorizontal[3], Obj->secuencia_vertical[y-3]),
											calcula_matchReplace(Obj, letraSecuenciaHorizontal[2], Obj->secuencia_vertical[y-2]),
											calcula_matchReplace(Obj, letraSecuenciaHorizontal[1], Obj->secuencia_vertical[y-1]),
											calcula_matchReplace(Obj, letraSecuenciaHorizontal[0], Obj->secuencia_vertical[y])
				); // reverse order

		r_res_diag = _mm256_add_epi32(r_diag_valor, r_matchReplace);
		//
		r_max_valor = _mm256_max_epi32(r_max_valor, r_res_diag);
#ifdef DEBUG_VECTOR
		printf("%c %c, %c %c, %c %c, %c %c, %c %c, %c %c, %c %c, %c %c\n",
				letraSecuenciaHorizontal[0],Obj->secuencia_vertical[y],
				letraSecuenciaHorizontal[1],Obj->secuencia_vertical[y-1],
				letraSecuenciaHorizontal[2],Obj->secuencia_vertical[y-2],
				letraSecuenciaHorizontal[3],Obj->secuencia_vertical[y-3],
				letraSecuenciaHorizontal[4],Obj->secuencia_vertical[y-4],
				letraSecuenciaHorizontal[5],Obj->secuencia_vertical[y-5],
				letraSecuenciaHorizontal[6],Obj->secuencia_vertical[y-6],
				letraSecuenciaHorizontal[7],Obj->secuencia_vertical[y-7]);
		print_reg("MatchReplace ", r_matchReplace);
		print_reg("Valor ", r_max_valor);
#endif
		//#ifdef SMITHWATERMAN
			r_max_valor = _mm256_max_epi32(r_max_valor, r_const_zero);
		//#endif
		//
		//
		r_max_vertical = _mm256_sub_epi32(r_max_valor, r_const_insertCost);
		r_max_vertical = _mm256_max_epi32(r_max_vertical, r_res_vertical);
		r_max_vertical = _mm256_sub_epi32(r_max_vertical, r_const_gapExtend);
		r_max_vertical = _mm256_sub_epi32(r_max_valor, r_max_vertical);
		//
		//
		r_max_horizontal = _mm256_sub_epi32(r_max_valor, r_const_deleteCost);
		r_max_horizontal = _mm256_max_epi32(r_max_horizontal, r_res_horizontal);
		r_max_horizontal = _mm256_sub_epi32(r_max_horizontal, r_const_gapExtend);
		r_max_horizontal = _mm256_sub_epi32(r_max_valor, r_max_horizontal);

//		printf("%d\n", y);
//		print_reg("Max valor ", r_max_valor);
//		print_reg("Max horiz ", r_max_horizontal);
//		print_reg("Max vertical ", r_max_vertical);
		//
		//
		//Guardamos los valores extremos
		// TODO: Esto tarda una barbaridad.
		contenidoSubColumnaTrabajo[y+2 - int_per_simdreg].valor = _mm256_extract_epi32(r_max_valor, 7);
		contenidoSubColumnaTrabajo[y+2 - int_per_simdreg].vertical = _mm256_extract_epi32(r_max_vertical, 7);
		contenidoSubColumnaTrabajo[y+2 - int_per_simdreg].horizontal = _mm256_extract_epi32(r_max_horizontal, 7);
#ifdef DEBUG_VECTOR
		print_nodo("Guardo ", contenidoSubColumnaTrabajo[y+2 - int_per_simdreg], y+2 - int_per_simdreg);
		printf("SubColumna %d: valor:%ld vertical:%d horizontal:%d \n",y+2 - int_per_simdreg,contenidoSubColumnaTrabajo[y+2 - int_per_simdreg].valor,contenidoSubColumnaTrabajo[y+2 - int_per_simdreg].vertical, contenidoSubColumnaTrabajo[y+2 - int_per_simdreg].horizontal);
		printf("%ld\n",sizeof(contenidoSubColumnaTrabajo[y+2 - int_per_simdreg].horizontal));
#endif
#ifdef SMITHWATERMAN
		r_maximum = _mm256_max_epi32(r_maximum, r_max_valor);
#endif
		// Preparamos al siguiente bucle
		// r_max se convierte en r_up
		r_up_valor = r_max_valor;
		r_up_vertical = r_max_vertical;
		// r_left se convierte en r_diag
		r_diag_valor = r_left_valor;
		// r_max desplazado se convierte en r_left
//		printf("Antes\n");
		// TODO Al finalizar el bucle se hace un acceso ilegal a memoria (indexOutOfbounds)
		r_left_valor = _mm256_permutevar8x32_epi32(r_max_valor, r_const_permutacion);
		r_left_horizontal = _mm256_permutevar8x32_epi32(r_max_horizontal, r_const_permutacion);
		if (y+2<size_subcolumna_partida){
			r_left_valor = _mm256_insert_epi32(r_left_valor, contenidoSubColumnaPartida[y+2].valor, 0);
			r_left_horizontal = _mm256_insert_epi32(r_left_horizontal, contenidoSubColumnaPartida[y+2].horizontal, 0);
		}
//		printf("Despues\n");
#ifdef DEBUG_VECTOR
		print_nodo("Parto de valor ",contenidoSubColumnaPartida[y+2], y+2);
#endif
	}
//    clock_gettime(CLOCK_MONOTONIC, &tend);
//    printf("some_long_computation took about %.9f seconds\n",
//           ((double)tend.tv_sec + 1.0e-9*tend.tv_nsec) -
//           ((double)tstart.tv_sec + 1.0e-9*tstart.tv_nsec));

	// Al acabar el proceso vectorial debemos reobtener r_diag_a y r_diag_b para que finalice el cálculo en la función llamante
	// Guardar r_max en r_diag_b
	int valor[int_per_simdreg] __attribute__((aligned (32)));
	int vertical[int_per_simdreg] __attribute__((aligned (32)));
	int horizontal[int_per_simdreg] __attribute__((aligned (32)));
	_mm256_maskstore_epi32 (valor, r_const_mascara, r_max_valor);
	_mm256_maskstore_epi32 (vertical, r_const_mascara, r_max_vertical);
	_mm256_maskstore_epi32 (horizontal, r_const_mascara, r_max_horizontal);
	for(int diag=1;diag<=int_per_simdreg;diag++){ // Reverse order
		r_diag_b[diag].valor=valor[diag-1];//int_per_simdreg - diag];
		r_diag_b[diag].horizontal=horizontal[diag-1];//int_per_simdreg - diag];
		r_diag_b[diag].vertical=vertical[diag-1];//int_per_simdreg - diag];
	}
	// En r_diag_a sólo necesitamos el valor y los elementos 1Âº y último tampoco nos hacen falta
	_mm256_maskstore_epi32 (valor, r_const_mascara, r_diag_valor);
	for(int diag=1;diag<int_per_simdreg;diag++){ // Reverse order
		r_diag_a[diag].valor=valor[diag];// diag-1 //int_per_simdreg - diag];
	}
#ifdef SMITHWATERMAN
	// Miramos si hay un nuevo máximo
	int local_max = horizontal_max_Vec8i(r_maximum);
	if (local_max > Obj->valor_maximo){
		Obj->valor_maximo = local_max;
	}
#endif
}

// END Vectorización 8
#endif


// Vectorización 16
#ifdef VECTOR_16

void print_reg(char *s, __m512i r){
	int valor[16] __attribute__((aligned (16)));
	_mm512_store_epi32 ((__m512i *)valor, r);
	printf("%s", s);
	for(int i=0; i<16; i++){
		printf("%d ", valor[i]);
	}
	printf("\n");
}

int horizontal_max_Vec16i(__m512i x) {
	/*
	__m512i shuf = _mm512_set_epi32(0,0,0,0,0,0,0,0,15,14,13,12,11,10,9,8);
	__m512i max =  _mm512_max_epi32(_mm512_permutexvar_epi32(shuf, x), x);
    shuf = _mm512_set_epi32(0,0,0,0,0,0,0,0,0,0,0,0,7,6,5,4);
    max =  _mm512_max_epi32(_mm512_permutexvar_epi32(shuf, max), max);
    shuf = _mm512_set_epi32(0,0,0,0,0,0,0,0,0,0,0,0,0,0,3,2);
    max =  _mm512_max_epi32(_mm512_permutexvar_epi32(shuf, max), max);
    shuf = _mm512_set_epi32(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1);
    max =  _mm512_max_epi32(_mm512_permutexvar_epi32(shuf, max), max);
	return _mm512_extract_epi32((union __m512i)max, 0);
	*/
	return _mm512_reduce_max_epi32 (x);
}

void vectorizacionColumnas(	estado_del_Objeto * Obj,
							struct Nodo_valor* contenidoSubColumnaPartida,
							struct Nodo_valor* contenidoSubColumnaTrabajo,
							struct Nodo_valor * r_diag_a,
							struct Nodo_valor * r_diag_b,
							char * letraSecuenciaHorizontal){
	__m512i r_up_valor __attribute__((aligned (64)));
	__m512i r_up_vertical __attribute__((aligned (64)));
	__m512i r_left_valor __attribute__((aligned (64)));
	__m512i r_left_horizontal __attribute__((aligned (64)));
	__m512i r_res_vertical __attribute__((aligned (64)));
	__m512i r_res_horizontal __attribute__((aligned (64)));
	__m512i r_max_valor __attribute__((aligned (64)));
	__m512i r_diag_valor __attribute__((aligned (64)));
	__m512i r_matchReplace __attribute__((aligned (64)));
	__m512i r_res_diag __attribute__((aligned (64)));
	__m512i r_max_vertical __attribute__((aligned (64)));
	__m512i r_max_horizontal __attribute__((aligned (64)));
	__m512i r_const_insertCost __attribute__((aligned (64)));
	__m512i r_const_deleteCost __attribute__((aligned (64)));
	__m512i r_const_gapExtend __attribute__((aligned (64)));
	__m512i r_const_permutacion __attribute__((aligned (64)));
	__m512i r_const_zero __attribute__((aligned (64))) = _mm512_setzero_epi32();
	#ifdef SMITHWATERMAN
	__m512i r_maximum __attribute__((aligned (64))) = _mm512_setzero_epi32();
	#endif

		//printf("Entrando...%ld\n", letraSecuenciaHorizontal);

	//Trabajaremos con los arrays para mayor sencillez de acceso
	int size_subcolumna_partida = get_size(Obj->trozoSubColumnaPartida);
	// Inicialización de constantes
	//
	r_const_deleteCost = _mm512_set1_epi32(DELETE_COST);
	r_const_insertCost = _mm512_set1_epi32(INSERT_COST);
	r_const_gapExtend = _mm512_set1_epi32(GAPEXTEND_COST);
	r_const_permutacion = _mm512_set_epi32(14,13,12,11,10,9,8,7,6,5,4,3,2,1,0,0);

	// Inicialización de valores de entrada
	//
	r_up_valor = _mm512_set_epi32(r_diag_b[16].valor,r_diag_b[15].valor,r_diag_b[14].valor,r_diag_b[13].valor,r_diag_b[12].valor,r_diag_b[11].valor,r_diag_b[10].valor,r_diag_b[9].valor,r_diag_b[8].valor,r_diag_b[7].valor,r_diag_b[6].valor,r_diag_b[5].valor,r_diag_b[4].valor,r_diag_b[3].valor,r_diag_b[2].valor,r_diag_b[1].valor); // reverse order
	r_up_vertical = _mm512_set_epi32(r_diag_b[16].vertical,r_diag_b[15].vertical,r_diag_b[14].vertical,r_diag_b[13].vertical,r_diag_b[12].vertical,r_diag_b[11].vertical,r_diag_b[10].vertical,r_diag_b[9].vertical,r_diag_b[8].vertical,r_diag_b[7].vertical,r_diag_b[6].vertical,r_diag_b[5].vertical,r_diag_b[4].vertical,r_diag_b[3].vertical,r_diag_b[2].vertical,r_diag_b[1].vertical); // reverse order
	//
	r_left_valor = _mm512_set_epi32(r_diag_b[15].valor,r_diag_b[14].valor,r_diag_b[13].valor,r_diag_b[12].valor,r_diag_b[11].valor,r_diag_b[10].valor,r_diag_b[9].valor,r_diag_b[8].valor,r_diag_b[7].valor,r_diag_b[6].valor,r_diag_b[5].valor,r_diag_b[4].valor,r_diag_b[3].valor,r_diag_b[2].valor,r_diag_b[1].valor,r_diag_b[0].valor); // reverse order
	r_left_horizontal = _mm512_set_epi32(r_diag_b[15].horizontal,r_diag_b[14].horizontal,r_diag_b[13].horizontal,r_diag_b[12].horizontal,r_diag_b[11].horizontal,r_diag_b[10].horizontal,r_diag_b[9].horizontal,r_diag_b[8].horizontal,r_diag_b[7].horizontal,r_diag_b[6].horizontal,r_diag_b[5].horizontal,r_diag_b[4].horizontal,r_diag_b[3].horizontal,r_diag_b[2].horizontal,r_diag_b[1].horizontal,r_diag_b[0].horizontal); // reverse order
	//
	r_diag_valor =  _mm512_set_epi32(r_diag_a[15].valor,r_diag_a[14].valor,r_diag_a[13].valor,r_diag_a[12].valor,r_diag_a[11].valor,r_diag_a[10].valor,r_diag_a[9].valor,r_diag_a[8].valor,r_diag_a[7].valor,r_diag_a[6].valor,r_diag_a[5].valor,r_diag_a[4].valor,r_diag_a[3].valor,r_diag_a[2].valor,r_diag_a[1].valor,r_diag_a[0].valor); // reverse order


	// Comienzo del bucle
//	struct timespec tstart={0,0}, tend={0,0};
//	    clock_gettime(CLOCK_MONOTONIC, &tstart);
	for(int y=int_per_simdreg - 1; y<size_subcolumna_partida-1; y++) {

		r_res_vertical = _mm512_sub_epi32(r_up_valor, r_up_vertical);
		r_res_horizontal = _mm512_sub_epi32(r_left_valor, r_left_horizontal);
		//
		r_max_valor = _mm512_max_epi32(r_res_vertical, r_res_horizontal);
		r_matchReplace = _mm512_set_epi32(
											calcula_matchReplace(Obj, letraSecuenciaHorizontal[15], Obj->secuencia_vertical[y-15]),
											calcula_matchReplace(Obj, letraSecuenciaHorizontal[14], Obj->secuencia_vertical[y-14]),
											calcula_matchReplace(Obj, letraSecuenciaHorizontal[13], Obj->secuencia_vertical[y-13]),
											calcula_matchReplace(Obj, letraSecuenciaHorizontal[12], Obj->secuencia_vertical[y-12]),
											calcula_matchReplace(Obj, letraSecuenciaHorizontal[11], Obj->secuencia_vertical[y-11]),
											calcula_matchReplace(Obj, letraSecuenciaHorizontal[10], Obj->secuencia_vertical[y-10]),
											calcula_matchReplace(Obj, letraSecuenciaHorizontal[9], Obj->secuencia_vertical[y-9]),
											calcula_matchReplace(Obj, letraSecuenciaHorizontal[8], Obj->secuencia_vertical[y-8]),
											calcula_matchReplace(Obj, letraSecuenciaHorizontal[7], Obj->secuencia_vertical[y-7]),
											calcula_matchReplace(Obj, letraSecuenciaHorizontal[6], Obj->secuencia_vertical[y-6]),
											calcula_matchReplace(Obj, letraSecuenciaHorizontal[5], Obj->secuencia_vertical[y-5]),
											calcula_matchReplace(Obj, letraSecuenciaHorizontal[4], Obj->secuencia_vertical[y-4]),
											calcula_matchReplace(Obj, letraSecuenciaHorizontal[3], Obj->secuencia_vertical[y-3]),
											calcula_matchReplace(Obj, letraSecuenciaHorizontal[2], Obj->secuencia_vertical[y-2]),
											calcula_matchReplace(Obj, letraSecuenciaHorizontal[1], Obj->secuencia_vertical[y-1]),
											calcula_matchReplace(Obj, letraSecuenciaHorizontal[0], Obj->secuencia_vertical[y])
				); // reverse order

		r_res_diag = _mm512_add_epi32(r_diag_valor, r_matchReplace);
		//
		r_max_valor = _mm512_max_epi32(r_max_valor, r_res_diag);
#ifdef DEBUG_VECTOR
		printf("%c %c, %c %c, %c %c, %c %c, %c %c, %c %c, %c %c, %c %c, %c %c, %c %c, %c %c, %c %c, %c %c, %c %c, %c %c, %c %c\n",
				letraSecuenciaHorizontal[0],Obj->secuencia_vertical[y],
				letraSecuenciaHorizontal[1],Obj->secuencia_vertical[y-1],
				letraSecuenciaHorizontal[2],Obj->secuencia_vertical[y-2],
				letraSecuenciaHorizontal[3],Obj->secuencia_vertical[y-3],
				letraSecuenciaHorizontal[4],Obj->secuencia_vertical[y-4],
				letraSecuenciaHorizontal[5],Obj->secuencia_vertical[y-5],
				letraSecuenciaHorizontal[6],Obj->secuencia_vertical[y-6],
				letraSecuenciaHorizontal[7],Obj->secuencia_vertical[y-7],
				letraSecuenciaHorizontal[8],Obj->secuencia_vertical[y-8],
				letraSecuenciaHorizontal[9],Obj->secuencia_vertical[y-9],
				letraSecuenciaHorizontal[10],Obj->secuencia_vertical[y-10],
				letraSecuenciaHorizontal[11],Obj->secuencia_vertical[y-11],
				letraSecuenciaHorizontal[12],Obj->secuencia_vertical[y-12],
				letraSecuenciaHorizontal[13],Obj->secuencia_vertical[y-13],
				letraSecuenciaHorizontal[14],Obj->secuencia_vertical[y-14],
				letraSecuenciaHorizontal[15],Obj->secuencia_vertical[y-15]);
		print_reg("MatchReplace ", r_matchReplace);
		print_reg("Valor ", r_max_valor);
#endif
		//#ifdef SMITHWATERMAN
			r_max_valor = _mm512_max_epi32(r_max_valor, r_const_zero);
		//#endif
		//
		//
		r_max_vertical = _mm512_sub_epi32(r_max_valor, r_const_insertCost);
		r_max_vertical = _mm512_max_epi32(r_max_vertical, r_res_vertical);
		r_max_vertical = _mm512_sub_epi32(r_max_vertical, r_const_gapExtend);
		r_max_vertical = _mm512_sub_epi32(r_max_valor, r_max_vertical);
		//
		//
		r_max_horizontal = _mm512_sub_epi32(r_max_valor, r_const_deleteCost);
		r_max_horizontal = _mm512_max_epi32(r_max_horizontal, r_res_horizontal);
		r_max_horizontal = _mm512_sub_epi32(r_max_horizontal, r_const_gapExtend);
		r_max_horizontal = _mm512_sub_epi32(r_max_valor, r_max_horizontal);

//		printf("%d\n", y);
//		print_reg("Max valor ", r_max_valor);
//		print_reg("Max horiz ", r_max_horizontal);
//		print_reg("Max vertical ", r_max_vertical);
		//
		//
		//Guardamos los valores extremos
		// TODO: Esto tarda una barbaridad.
		contenidoSubColumnaTrabajo[y+2 - int_per_simdreg].valor = my_mm512_extract_epi32(r_max_valor, 15);
		contenidoSubColumnaTrabajo[y+2 - int_per_simdreg].vertical = my_mm512_extract_epi32(r_max_vertical, 15);
		contenidoSubColumnaTrabajo[y+2 - int_per_simdreg].horizontal = my_mm512_extract_epi32(r_max_horizontal, 15);
#ifdef DEBUG_VECTOR
		print_nodo("Guardo ", contenidoSubColumnaTrabajo[y+2 - int_per_simdreg], y+2 - int_per_simdreg);
		printf("SubColumna %d: valor:%ld vertical:%d horizontal:%d \n",y+2 - int_per_simdreg,contenidoSubColumnaTrabajo[y+2 - int_per_simdreg].valor,contenidoSubColumnaTrabajo[y+2 - int_per_simdreg].vertical, contenidoSubColumnaTrabajo[y+2 - int_per_simdreg].horizontal);
		printf("%ld\n",sizeof(contenidoSubColumnaTrabajo[y+2 - int_per_simdreg].horizontal));
#endif
		#ifdef SMITHWATERMAN
			r_maximum = _mm512_max_epi32(r_maximum, r_max_valor);
		#endif
		// Preparamos al siguiente bucle
		// r_max se convierte en r_up
		r_up_valor = r_max_valor;
		r_up_vertical = r_max_vertical;
		// r_left se convierte en r_diag
		r_diag_valor = r_left_valor;
		// r_max desplazado se convierte en r_left
//		printf("Antes\n");
		// TODO Al finalizar el bucle se hace un acceso ilegal a memoria (indexOutOfbounds)
		r_left_valor = _mm512_permutexvar_epi32(r_const_permutacion, r_max_valor);
		r_left_horizontal = _mm512_permutexvar_epi32(r_const_permutacion, r_max_horizontal);
		if (y+2<size_subcolumna_partida){
			r_left_valor = _mm512_mask_i32gather_epi32(r_left_valor, 0x0001, r_const_zero, &(contenidoSubColumnaPartida[y+2].valor), 1);
			r_left_horizontal = _mm512_mask_i32gather_epi32(r_left_horizontal, 0x0001, r_const_zero, &(contenidoSubColumnaPartida[y+2].horizontal), 1);
			//r_left_valor = _mm512_insert_epi32(r_left_valor, contenidoSubColumnaPartida[y+2].valor, 0);
			//r_left_horizontal = _mm512_insert_epi32(r_left_horizontal, contenidoSubColumnaPartida[y+2].horizontal, 0);
		}
//		printf("Despues\n");
#ifdef DEBUG_VECTOR
		print_nodo("Parto de valor ",contenidoSubColumnaPartida[y+2], y+2);
#endif
	}
//    clock_gettime(CLOCK_MONOTONIC, &tend);
//    printf("some_long_computation took about %.9f seconds\n",
//           ((double)tend.tv_sec + 1.0e-9*tend.tv_nsec) -
//           ((double)tstart.tv_sec + 1.0e-9*tstart.tv_nsec));

	// Al acabar el proceso vectorial debemos reobtener r_diag_a y r_diag_b para que finalice el cálculo en la función llamante
	// Guardar r_max en r_diag_b
	int valor[int_per_simdreg] __attribute__((aligned (64)));
	int vertical[int_per_simdreg] __attribute__((aligned (64)));
	int horizontal[int_per_simdreg] __attribute__((aligned (64)));
	_mm512_store_epi32 (valor, r_max_valor);
	_mm512_store_epi32 (vertical, r_max_vertical);
	_mm512_store_epi32 (horizontal, r_max_horizontal);
	for(int diag=1;diag<=int_per_simdreg;diag++){ // Reverse order
		r_diag_b[diag].valor=valor[diag-1];//int_per_simdreg - diag];
		r_diag_b[diag].horizontal=horizontal[diag-1];//int_per_simdreg - diag];
		r_diag_b[diag].vertical=vertical[diag-1];//int_per_simdreg - diag];
	}
	// En r_diag_a sólo necesitamos el valor y los elementos 1Âº y último tampoco nos hacen falta
	_mm512_store_epi32 (valor, r_diag_valor);
	for(int diag=1;diag<int_per_simdreg;diag++){ // Reverse order
		r_diag_a[diag].valor=valor[diag];// diag-1 //int_per_simdreg - diag];
	}
#ifdef SMITHWATERMAN
	// Miramos si hay un nuevo máximo
	int local_max = horizontal_max_Vec16i(r_maximum);
	if (local_max > Obj->valor_maximo){
		Obj->valor_maximo = local_max;
	}
#endif
}

// END Vectorización 16
#endif



//HASTA AQUI VECTORIZACION


void realiza_procesamiento_columnasSIMD(estado_del_Objeto * Obj) {
	//Guardamos los valores enteros para evitar tener que recalcular más veces
	int size_subcolumna_partida = get_size(Obj->trozoSubColumnaPartida);
	int size_subfila_partida = get_size(Obj->trozoSubFilaPartida);
	// CONTROLAR SI MERECE LA PENA VECTORIZAR
	if (size_subcolumna_partida < 4*int_per_simdreg){
		// No merece la pena hacer vectorización
		realiza_procesamiento_columnas(Obj);
		return;
	}
	#ifdef SMITHWATERMAN
		// Inicializamos las posiciones del máximo
		Obj->pos_max_x = Obj->pos_max_y = -1;
	#endif

	#pragma GCC diagnostic push
		#pragma GCC diagnostic ignored "-Wchar-subscripts"
		// Traducimos las letras para ir más rápidos
		for(int i=0; i < size_subfila_partida; i++)
			Obj->secuencia_horizontal[i] = Obj->tabla_puntuacion->codificacion_horizontal[Obj->secuencia_horizontal[i]];
		for(int j=0; j < size_subcolumna_partida; j++)
			Obj->secuencia_vertical[j] = Obj->tabla_puntuacion->codificacion_vertical[Obj->secuencia_vertical[j]];
	#pragma GCC diagnostic pop

	// Creamos los fragmentos de trabajo en memoria local
	Obj->trozoSubColumnaTrabajo = crea_fragmento(size_subcolumna_partida);
	Obj->trozoSubFilaTrabajo = crea_fragmento(size_subfila_partida);

	//Trabajaremos con los arrays para mayor sencillez de acceso
	struct Nodo_valor* contenidoSubColumnaTrabajo = Obj->trozoSubColumnaTrabajo->contenido;
	struct Nodo_valor* contenidoSubFilaTrabajo = Obj->trozoSubFilaTrabajo->contenido;
	struct Nodo_valor* contenidoSubColumnaPartida = Obj->trozoSubColumnaPartida->contenido;
	struct Nodo_valor* contenidoSubFilaPartida = Obj->trozoSubFilaPartida->contenido;

	/* Trabajaremos por columnas: optimizando el número de columnas a calcular si no necesitamos todo */
	int columnas_a_calcular = ((Obj->limite_x>0) ? Obj->limite_x : size_subfila_partida);

 	//Punteros a nodo para mayor optimización: ahorro de acceso a array con offset para cada valor de la estructura
	register struct Nodo_valor *temp, *temp_izq, *temp_arr, *temp_diag, *aux;

	// Diagonales menores para vectorización
	struct Nodo_valor * r_diag_a = (struct Nodo_valor *)malloc_safe_agr248(sizeof(struct Nodo_valor) * (int_per_simdreg + 1));
	struct Nodo_valor * r_diag_b = (struct Nodo_valor *)malloc_safe_agr248(sizeof(struct Nodo_valor) * (int_per_simdreg + 1));
	struct Nodo_valor * r_diag_c = (struct Nodo_valor *)malloc_safe_agr248(sizeof(struct Nodo_valor) * (int_per_simdreg + 1));

	/* Cálculo del contenido */
    //contenidoSubFilaTrabajo[0] = contenidoSubColumnaPartida[size_subcolumna_partida-1];
    copia_valor_nodo(contenidoSubFilaTrabajo[0], contenidoSubColumnaPartida[size_subcolumna_partida-1]);

    int x;
    for (x=1; x+int_per_simdreg-1 < columnas_a_calcular; x=x+int_per_simdreg){
#ifdef DEBUG_VECTOR
    	printf("\nX vale %d\n\n", x);
#endif
    	//Leemos las letras de la secuencia H a partir del Ã­ndice que nos corresponde
    	char * letraSecuenciaHorizontal;
    	// for(int i=0; i<int_per_simdreg; i++)
    		//letraSecuenciaHorizontal[i] = Obj->secuencia_horizontal[x+i-1];
    	letraSecuenciaHorizontal = Obj->secuencia_horizontal+x-1;
    	/* Calcular trozoSubColumnaTrabajo a partir de trozoSubColumnaPartida */
    	copia_valor_nodo(contenidoSubColumnaTrabajo[0], contenidoSubFilaPartida[x+int_per_simdreg-1]);
    	//
    	//
    	// Calculamos la diagonal superior para preparar el cálculo vectorizado
    	//
    	copia_valor_nodo(r_diag_a[0], contenidoSubColumnaPartida[0]);
    	copia_valor_nodo(r_diag_b[0], contenidoSubColumnaPartida[1]);
    	copia_valor_nodo(r_diag_b[1], contenidoSubFilaPartida[x]);
    	for(int diag=2; diag <= int_per_simdreg; diag++){
    		copia_valor_nodo(r_diag_c[0], contenidoSubColumnaPartida[diag]);
    		copia_valor_nodo(r_diag_c[diag], contenidoSubFilaPartida[x + diag - 1]);
    		for(int pos=1; pos<diag; pos++){
    			char letraSecuenciaHorizontal = Obj->secuencia_horizontal[x + pos-2];
    			char letraSecuenciaVertical = Obj->secuencia_vertical[diag - 1 - pos];
    			  temp = &(r_diag_c[pos]);
    			  temp_arr = &(r_diag_b[pos]);
    			  temp_izq = &(r_diag_b[pos - 1]);
    			  temp_diag = &(r_diag_a[pos - 1]);
    			  // Calculamos los valores haciendo uso de una macro
    			  calcula_valores(temp_diag, temp_arr, temp_izq, temp);

          		//Sólo si estamos en Smith-Waterman: vamos recalculando el máximo
          		#ifdef SMITHWATERMAN
          			if (temp->valor > Obj->valor_maximo) {
          				Obj->valor_maximo = temp->valor;
          				Obj->pos_max_x = x + pos-1;
          				Obj->pos_max_y = diag - pos;
          			}
          		#endif
#ifdef DEBUG_VECTOR
          			printf("lHor %c lVert %c Valor %ld vertical %d horizontal %d\n", letraSecuenciaHorizontal, letraSecuenciaVertical, temp->valor, temp->vertical, temp->horizontal);
                	printf("-----------\n");
#endif
    		}
    		//printf("-----------\n");
    		// r_diag_a ya no se va a usar más
    		struct Nodo_valor * aux;
    		aux = r_diag_a;
    		r_diag_a = r_diag_b;
    		r_diag_b = r_diag_c;
    		r_diag_c = aux;
    	}

    	//Sólo si estamos en Smith-Waterman: vamos recalculando el máximo
		#ifdef SMITHWATERMAN
			int prev_valor_maximo = Obj->valor_maximo;
		#endif
    	//
    	// r_diag_a y r_diag_b tienen los datos de partida para vectorizar.
    	vectorizacionColumnas(Obj, contenidoSubColumnaPartida, contenidoSubColumnaTrabajo, r_diag_a, r_diag_b, letraSecuenciaHorizontal);
    	//
    	//Sólo si estamos en Smith-Waterman: estimamos dónde está el máximo
		#ifdef SMITHWATERMAN
			if (Obj->valor_maximo > prev_valor_maximo){
  				Obj->pos_max_x = x+int_per_simdreg-1;
  				Obj->pos_max_y = size_subcolumna_partida - 1;
			}
		#endif

    	// Calculamos la diagonal inferior que no se ha podido vectorizar
    	// Justo antes guardamos en la subfila de trabajo el primer nodo de la última diagonal vectorizada
    	copia_valor_nodo(contenidoSubFilaTrabajo[x], r_diag_b[1]);
    	for(int diag=int_per_simdreg - 2; diag >= 0; diag--){

//    		for(int w=0;w<5;w++) print_nodo("diag_a", r_diag_a[w], w);
//    		for(int w=0;w<5;w++) print_nodo("diag_b", r_diag_b[w], w);
//    		for(int w=0;w<5;w++) print_nodo("diag_c", r_diag_c[w], w);

    		for(int pos=int_per_simdreg - diag; pos<=int_per_simdreg; pos++){
    			char letraSecuenciaHorizontal = Obj->secuencia_horizontal[x+pos-2];
    			char letraSecuenciaVertical = Obj->secuencia_vertical[size_subcolumna_partida - 2 + (int_per_simdreg-diag) - pos];
    			  temp = &(r_diag_c[pos]);
    			  temp_arr = &(r_diag_b[pos]);
    			  temp_izq = &(r_diag_b[pos - 1]);
    			  temp_diag = &(r_diag_a[pos - 1]);
    			  // Calculamos los valores haciendo uso de una macro
    			  calcula_valores(temp_diag, temp_arr, temp_izq, temp);

          		//Sólo si estamos en Smith-Waterman: vamos recalculando el máximo
          		#ifdef SMITHWATERMAN
          			if (temp->valor > Obj->valor_maximo) {
          				Obj->valor_maximo = temp->valor;
          				Obj->pos_max_x = x+pos-1;
          				Obj->pos_max_y = size_subcolumna_partida - 1 + (int_per_simdreg-diag) - pos;
          			}
          		#endif

#ifdef DEBUG_VECTOR
          		printf("diag %d pos %d lHoriz %c lVerti %c\n", diag, pos, letraSecuenciaHorizontal, letraSecuenciaVertical);
#endif
    		}
//    		print_nodo("En columna ", *temp, size_subcolumna_partida - 1 - diag);
    		// El último nodo calculado va a parar a la subcolumna de trabajo ...
        	copia_valor_nodo(contenidoSubColumnaTrabajo[size_subcolumna_partida - 1 - diag], *temp);
        	// ... y el primero va a parar a la subfila de trabajo.
//        	print_nodo("En fila ", r_diag_c[int_per_simdreg - diag], x + int_per_simdreg - diag - 1);
        	copia_valor_nodo(contenidoSubFilaTrabajo[x + int_per_simdreg - diag - 1], r_diag_c[int_per_simdreg - diag]);
    		// r_diag_a ya no se va a usar más
    		struct Nodo_valor * aux;
    		aux = r_diag_a;
    		r_diag_a = r_diag_b;
    		r_diag_b = r_diag_c;
    		r_diag_c = aux;
    	}


        /* Nos preparamos para la siguiente iteración evitando tomar más memoria.
         * Para ello se hace un intercambio de arrays que hay que deshacer tras el último ciclo. */
        aux = contenidoSubColumnaPartida;
        contenidoSubColumnaPartida = contenidoSubColumnaTrabajo;
        contenidoSubColumnaTrabajo = aux;
    }
    // Liberamos las diagonales menores
	free(r_diag_c);
	free(r_diag_b);
	free(r_diag_a);


#ifdef DEBUG_VECTOR
	printf("Max %d\n", Obj->valor_maximo);
#endif
	// Por si queda alguna columna que calcular,
	// se procede a actuar de manera tradicional.
	// CUIDADO: ESTO ES UN COPY & PASTE DEL BLOQUE CENTRAL DE LA FUNCIÃ“N realiza_procesamiento_columnas
    for (/* x=x */; x < columnas_a_calcular; x++){


    	// printf("Un fleco %d\n", x);

    	//Leemos las letras de la secuencia H a partir del índice que nos corresponde
    	char letraSecuenciaHorizontal = Obj->secuencia_horizontal[x-1];
    	/* Calcular trozoSubColumnaTrabajo a partir de trozoSubColumnaPartida */
    	// contenidoSubColumnaTrabajo[0] = contenidoSubFilaPartida[x]
    	copia_valor_nodo(contenidoSubColumnaTrabajo[0], contenidoSubFilaPartida[x]);

    	for (int y=1; y < size_subcolumna_partida; y++){
    		//Leemos las letras de la secuencia V a partir del índice que nos corresponde
    		char letraSecuenciaVertical = Obj->secuencia_vertical[y-1];

    		temp = &(contenidoSubColumnaTrabajo[y]);
    		temp_arr = &(contenidoSubColumnaTrabajo[y-1]);
    		temp_izq = &(contenidoSubColumnaPartida[y]);
    		temp_diag = &(contenidoSubColumnaPartida[y-1]);

    		// Calculamos los valores haciendo uso de una macro
    		calcula_valores(temp_diag, temp_arr, temp_izq, temp);

    		//Sólo si estamos en Smith-Waterman: vamos recalculando el máximo
    		#ifdef SMITHWATERMAN
    			if (contenidoSubColumnaTrabajo[y].valor > Obj->valor_maximo) {
    				Obj->valor_maximo = contenidoSubColumnaTrabajo[y].valor;
    				Obj->pos_max_x = x;
    				Obj->pos_max_y = y;
    			}
    		#endif
    	}
    	// contenidoSubFilaTrabajo[x] = contenidoSubColumnaTrabajo[size_subcolumna_partida - 1]
    	copia_valor_nodo(contenidoSubFilaTrabajo[x], contenidoSubColumnaTrabajo[size_subcolumna_partida - 1]);

        /* Nos preparamos para la siguien	te iteración evitando tomar más memoria.
         * Para ello se hace un intercambio de arrays que hay que deshacer tras el último ciclo. */
        aux = contenidoSubColumnaPartida;
        contenidoSubColumnaPartida = contenidoSubColumnaTrabajo;
        contenidoSubColumnaTrabajo = aux;
    }
#ifdef DEBUG_VECTOR
	printf("Max %d\n", Obj->valor_maximo);
#endif

    // El último cambio hay que deshacerlo (se sale de los límites de la matriz en la última iteración)
    contenidoSubColumnaTrabajo = contenidoSubColumnaPartida;


    //Sólo si estamos en Needleman-Wunsch: el máximo será la esquina inferior derecha
    #ifndef SMITHWATERMAN
    if (contenidoSubColumnaTrabajo[size_subcolumna_partida - 1].valor > Obj->valor_maximo) {
    	Obj->valor_maximo = contenidoSubColumnaTrabajo[size_subcolumna_partida - 1].valor;
    }
    #endif




    // Ya que no garantizamos que apunte realmente a su array (dependiendo de si es par o impar, HAY QUE COPIAR LOS VALORES)
    for (int y=0; y < size_subcolumna_partida; y++){
    	copia_valor_nodo(Obj->trozoSubColumnaTrabajo->contenido[y], contenidoSubColumnaTrabajo[y]);
    }
//    for(int x=0;x<=size_subfila_partida; x++){print_nodo("", Obj->trozoSubFilaTrabajo->contenido[x],x);}
//    printf("Columna---------------------\n");
//    	for (int y=0; y < size_subcolumna_partida; y++){print_nodo("", Obj->trozoSubColumnaTrabajo->contenido[y],y);}
}

#endif

/**  EL CALCULO SE REALIZARÃï¿½ POR FILAS, calculando en cada iteración una posición de la última subcolumna y al final la última subfila. */
void realiza_procesamiento_filas(estado_del_Objeto * Obj) {
	//Guardamos los valores enteros para evitar acceder numerosas veces a memoria compartida (más lento)
	int size_subcolumna_partida = get_size(Obj->trozoSubColumnaPartida);
	int size_subfila_partida = get_size(Obj->trozoSubFilaPartida);

	#ifdef SMITHWATERMAN
		// Inicializamos las posiciones del máximo
		Obj->pos_max_x = Obj->pos_max_y = -1;
	#endif

	// Creamos los fragmentos de trabajo en memoria local
	Obj->trozoSubColumnaTrabajo = crea_fragmento(size_subcolumna_partida);
	Obj->trozoSubFilaTrabajo = crea_fragmento(size_subfila_partida);

	//Trabajaremos con los arrays para mayor sencillez de acceso
	struct Nodo_valor* contenidoSubColumnaTrabajo = Obj->trozoSubColumnaTrabajo->contenido;
	struct Nodo_valor* contenidoSubFilaTrabajo = Obj->trozoSubFilaTrabajo->contenido;
	struct Nodo_valor* contenidoSubColumnaPartida = Obj->trozoSubColumnaPartida->contenido;
	struct Nodo_valor* contenidoSubFilaPartida = Obj->trozoSubFilaPartida->contenido;

	// Trabajaremos por filas: optimizando el número de filas a calcular si no necesitamos todo
	int filas_a_calcular = ((Obj->limite_y>0) ? Obj->limite_y : size_subcolumna_partida);

 	//Punteros a nodo para mayor optimización: ahorro de acceso a array con offset para cada valor de la estructura
	register struct Nodo_valor *temp, *temp_izq, *temp_arr, *temp_diag, *aux;
    // Cálculo del contenido
    //contenidoSubColumnaTrabajo[0] = contenidoSubFilaPartida[size_subfila_partida-1];
    copia_valor_nodo(contenidoSubColumnaTrabajo[0], contenidoSubFilaPartida[size_subfila_partida-1]);
    for (int y=1; y < filas_a_calcular; y++){
    	//Leemos las letras de la secuencia V a partir del índice que nos corresponde
    	char letraSecuenciaVertical = Obj->secuencia_vertical[y-1];
    	// Calcular trozoSubFilaTrabajo a partir de trozoSubFilaPartida
    	// contenidoSubFilaTrabajo[0] = contenidoSubColumnaPartida[y]
    	copia_valor_nodo(contenidoSubFilaTrabajo[0], contenidoSubColumnaPartida[y]);

    	for (int x=1; x < size_subfila_partida; x++){
    		//Leemos las letras de la secuencia H a partir del índice que nos corresponde
    		char letraSecuenciaHorizontal = Obj->secuencia_horizontal[x-1];

    		temp = &(contenidoSubFilaTrabajo[x]);
    		temp_arr = &(contenidoSubFilaPartida[x]);
    		temp_izq = &(contenidoSubFilaTrabajo[x-1]);
    		temp_diag = &(contenidoSubFilaPartida[x-1]);

    		// Calculamos los valores haciendo uso de una macro
    		calcula_valores(temp_diag, temp_arr, temp_izq, temp);

    		//Sólo si estamos en Smith-Waterman: vamos recalculando el máximo
    		#ifdef SMITHWATERMAN
    			if (contenidoSubFilaTrabajo[x].valor > Obj->valor_maximo) {
    				Obj->valor_maximo = contenidoSubFilaTrabajo[x].valor;
    				Obj->pos_max_x = x;
    				Obj->pos_max_y = y;
    			}
    		#endif
    	}
    	// contenidoSubColumnaTrabajo[y] = contenidoSubFilaTrabajo[size_subfila_partida - 1]
    	copia_valor_nodo(contenidoSubColumnaTrabajo[y], contenidoSubFilaTrabajo[size_subfila_partida - 1]);

        //Nos preparamos para la siguiente iteración evitando tomar más memoria.
        //Para ello se hace un intercambio de arrays que hay que deshacer tras el último ciclo.
        aux = contenidoSubFilaPartida;
        contenidoSubFilaPartida = contenidoSubFilaTrabajo;
        contenidoSubFilaTrabajo = aux;
    }
    // El último cambio hay que deshacerlo (se sale de los límites de la matriz en la última iteración)
    contenidoSubFilaTrabajo = contenidoSubFilaPartida;

    //Sólo si estamos en Needleman-Wunsch: el máximo será la esquina inferior derecha
    #ifndef SMITHWATERMAN
    if (contenidoSubFilaTrabajo[size_subfila_partida - 1].valor > Obj->valor_maximo) {
    	Obj->valor_maximo = contenidoSubFilaTrabajo[size_subfila_partida - 1].valor;
    }
    #endif

    // Ya que no garantizamos que apunte realmente a su array (dependiendo de si es par o impar, HAY QUE COPIAR LOS VALORES)
    for (int x=0; x < size_subfila_partida; x++){
    	copia_valor_nodo(Obj->trozoSubFilaTrabajo->contenido[x], contenidoSubFilaTrabajo[x]);
    }
}

/** Calcula una submatriz completa para un trabajo.
 *  Los límites limite_x y limite_y determinan el límite vertical y horizontal respectivamente.
 *  ÃÂ· Si limite_x==0 -> calcular todas las columnas
 *  ÃÂ· Si limite_y==0 -> calcular todas las filas
 *   -> [0,0] = todo */
Nodo_valor** calcula_submatriz(estado_del_Objeto * Obj, int anchura, int altura) {
	char letraSecuenciaHorizontal, letraSecuenciaVertical;
	Nodo_valor** matriz;
	//Punteros a nodo para mayor optimización: ahorro de acceso a array con offset para cada valor de la estructura
	struct Nodo_valor *temp, *temp_izq, *temp_arr, *temp_diag;

	/* CREACIÃâ€œN DE LA ESTRUCTURA */
	matriz = (Nodo_valor **) malloc_safe_agr248(sizeof(Nodo_valor *) * anchura);
	if (matriz == NULL) //No hay memoria suficiente
		return NULL;
    for (unsigned int i = 0; i < anchura; i++) {
       matriz[i] = (Nodo_valor *) malloc_safe_agr248(sizeof(Nodo_valor) * altura);

       if (matriz[i] == NULL) { //No hay memoria suficiente
           /* handle error, free() the previously allocated rows */
           for (unsigned int x = 0; x < i; x++)
           		free(matriz[x]);
           free(matriz);
           return NULL;
       }
    }

	/* INICIALIZACIÃâ€œN: fila y columna iniciales */
    for (int i=0; i<anchura; i++){
    	// matriz[i][0]=trozoSubFilaPartida->contenido[i]
    	copia_valor_nodo(matriz[i][0], Obj->trozoSubFilaPartida->contenido[i]);
    }

    for (int j=1; j<altura; j++){
    	// matriz[0][j]=trozoSubColumnaPartida->contenido[j]
    	copia_valor_nodo(matriz[0][j], Obj->trozoSubColumnaPartida->contenido[j]);
    }

    /* CÃï¿½LCULO NW CLÃï¿½SICO: por columnas */
    for (int i=1; i<anchura; i++){
		letraSecuenciaHorizontal = Obj->secuencia_horizontal[i-1];
        for (int j=1; j<altura; j++){
	    	letraSecuenciaVertical = Obj->secuencia_vertical[j-1];

	    	temp = &(matriz[i][j]);
	    	temp_arr = &(matriz[i][j-1]);
	    	temp_izq = &(matriz[i-1][j]);
	    	temp_diag = &(matriz[i-1][j-1]);

    		// Calculamos los valores haciendo uso de una macro
    		calcula_valores(temp_diag, temp_arr, temp_izq, temp);
        }
    }

    //Imprimimos la tabla
/*
#ifdef CLUSTALW
    //Imprimimos la tabla
   	for (int j = 0; j < altura; j++) {
   		printf("%d: ",j);
    	for (int i = 0; i < anchura; i++)
    	  	printf("(%4ld,%4ld,%4ld) | ", matriz[i][j].valor, matriz[i][j].horizontal, matriz[i][j].vertical);
    	printf("\n");
   	}
#else
    //Imprimimos la tabla
   	for (int j = 0; j < altura; j++) {
   		printf("%d: ",j);
    	for (int i = 0; i < anchura; i++)
    	  	printf("(%4ld,%4ld,%4ld) | ", matriz[i][j].valor, (matriz[i][j].valor - matriz[i][j].horizontal), (matriz[i][j].valor - matriz[i][j].vertical));
    	printf("\n");
   	}
#endif */

    return matriz;
}

/** Libera la memoria de submatriz de anchura determinada. */
void libera_submatriz(struct Nodo_valor*** ptr_matriz, int anchura){
	unsigned int i;
    for (i = 0; i < anchura; i++) {
    	free((*ptr_matriz)[i]);
    }
   	free(*ptr_matriz);
}

#ifndef SMITHWATERMAN // Función para Needleman-Wunsch
/** Cálculo de las secuencias de alineamiento a partir de una submatriz y la esquina inferior derecha de la misma. */
	void calcula_alineamiento(estado_del_Objeto * Obj, struct Nodo_valor** matriz, int num_filas, int num_columnas){
	int i=num_columnas-1, j=num_filas-1; //Ya que la posición tamaÃÂ±o está fuera de los índices
	//Punteros a nodo para mayor optimización: ahorro de acceso a array con offset para cada valor de la estructura
	struct Nodo_valor *temp, *temp_arr, *temp_izq;

	/* Alineamientos EN MEMORIA LOCAL y su contador. */
	Obj->alineamiento_vertical = (char *) malloc_safe_agr248 (sizeof(char) * (num_filas + num_columnas + 1));
	Obj->alineamiento_horizontal = (char *) malloc_safe_agr248 (sizeof(char) * (num_filas + num_columnas + 1));
	int pos_actual=0;

	// Implementación vuelta atrás: no busca máximo, sino de dónde viene realmente el valor

	char letraSecuenciaHorizontal, letraSecuenciaVertical;
	int direccion = Obj->direccion_inicial;
	int temp_valor_arriba, temp_valor_izquierda, temp_arriba, temp_izquierda;

	//Mientras no lleguemos a ningún límite, vamos desplazando y aÃÂ±adiendo letras
	while ((i>0) && (j>0)) {
		letraSecuenciaHorizontal = Obj->secuencia_horizontal[i-1];
    	letraSecuenciaVertical = Obj->secuencia_vertical[j-1];

    	temp = &(matriz[i][j]);
    	temp_arr = &(matriz[i][j-1]);
    	temp_izq = &(matriz[i-1][j]);

        // Como en las tuplas guardamos el vertical y horizontal de las siguientes celdas (ver diagramas Gotoh con diferencias),
        // para ver de dónde hemos venido, tenemos que precalcular y mirar en los horizontal/vertical de anteriores celdas
        // Además, hay que atender a la dirección de la que vengo para controlar los casos:
        // Vertical: ÃÂ¿Vengo realmente del anterior vertical o fue el principio del gap y me voy en diagonal/horizontal?
        // Horizontal: ÃÂ¿Vengo realmente del anterior horizontal o fue el principio del gap y me voy en diagonal/vertical?
        // Diagonal: ÃÂ¿De dónde vine? Comparar con los valores vertical y horizontal regenerados (se guardan diferencias) y si son iguales, es por ahí, sino, diagonal again
    	switch (direccion) {
    		case DIRECCION_VERTICAL:
#ifdef CLUSTALW
    			temp_valor_arriba = temp_arr->vertical - GAPEXTEND_COST; // Valor - GAP_EXTEND
    			temp_arriba = temp->vertical; // Ni diferencia ni ná de ná
#else
    			temp_valor_arriba = temp_arr->valor - (temp_arr->vertical + GAPEXTEND_COST); // Regeneramos el valor desde la diferencia - GAP_EXTEND
    			temp_arriba = temp->valor - temp->vertical; // Regeneramos el valor de la diferencia de esta casilla
#endif
        	// INSERCIÃâ€œN: desplazamiento vertical hacia arriba (vertical letra, horizontal hueco)
    			if (temp_arriba == temp_valor_arriba) // Vengo de arriba, coincide arriba->vert + GAP_EXTEND = yo
					direccion = DIRECCION_VERTICAL;
    			else {
    				// No sigo subiendo, vengo del valor máximo de mi tupla, no de la superior
    				// Lo normal es que venga de la diagonal, pero puede darse el caso X-/-X si mismatch > open_gap x 2
    				// Tenemos que verificar que el valor almacenado (máximo) de esta casilla no vino de la izquierda -> si es así, X-/-X, sino, diagonal
#ifdef CLUSTALW
    				temp_valor_izquierda = temp_izq->horizontal; // Ni diferencia ni ná de ná
#else
    				temp_valor_izquierda = temp_izq->valor - temp_izq->horizontal; // Regeneramos el valor de la diferencia
#endif
	    			if (temp->valor == temp_valor_izquierda) // Vengo de la izquierda, coincide izquierda->horiz = yo
						direccion = DIRECCION_HORIZONTAL;
    				else
    					direccion = DIRECCION_DIAGONAL; // Vengo de la diagonal, he abierto (valor = INSERT + GAPEXTEND)
    			}
    			break;
    		case DIRECCION_HORIZONTAL:
#ifdef CLUSTALW
    			temp_valor_izquierda = temp_izq->horizontal - GAPEXTEND_COST; // Valor - GAP_EXTEND
    			temp_izquierda = temp->horizontal; // Ni diferencia ni ná de ná
#else
    			temp_valor_izquierda = temp_izq->valor - (temp_izq->horizontal + GAPEXTEND_COST); // Regeneramos el valor desde la diferencia - GAP_EXTEND
    			temp_izquierda = temp->valor - temp->horizontal; // Regeneramos el valor de la diferencia de esta casilla
#endif
        	// DELETE: desplazamiento horizontal hacia la izquierda (vertical hueco, horizontal letra)
    			if (temp_izquierda == temp_valor_izquierda) // Vengo de la izquierda, coincide izquierda->horiz + GAP_EXTEND = yo
					direccion = DIRECCION_HORIZONTAL;
    			else {
    				// No sigo desplazándome a izquierda, vengo del valor máximo de mi tupla, no de la izquierda
    				// Lo normal es que venga de la diagonal, pero puede darse el caso -X/X- si mismatch > open_gap x 2
    				// Tenemos que verificar que el valor almacenado (máximo) de esta casilla no vino de arriba -> si es así, -X/X-, sino, diagonal
#ifdef CLUSTALW
    				temp_valor_arriba = temp_arr->vertical; // Ni diferencia ni ná de ná
#else
    				temp_valor_arriba = temp_arr->valor - temp_arr->vertical; // Regeneramos el valor de la diferencia
#endif
	    			if (temp->valor == temp_valor_arriba) // Vengo de arriba, coincide arriba->vert = yo
						direccion = DIRECCION_VERTICAL;
    				else
    					direccion = DIRECCION_DIAGONAL; // Vengo de la diagonal, he abierto (valor = DELETE + GAPEXTEND)
    			}
    			break;
    		default:
#ifdef CLUSTALW
    		 	// No hay diferencia que regenerar (si soy igual a ellos, vengo de ahí, en otro caso, diagonal)
    			temp_valor_izquierda = temp_izq->horizontal;
    			temp_valor_arriba = temp_arr->vertical;
#else
    		 	// Regeneramos los valores desde la diferencia (si soy igual a ellos, vengo de ahí, en otro caso, diagonal)
    			temp_valor_izquierda = temp_izq->valor - temp_izq->horizontal;
    			temp_valor_arriba = temp_arr->valor - temp_arr->vertical;
#endif
        	// MATCH / REPLACE: vertical letra, horizontal letra
    			if (temp->valor == temp_valor_izquierda) // Vengo de la izquierda, coincide izquierda->horiz = yo
					direccion = DIRECCION_HORIZONTAL;
    			else if (temp->valor == temp_valor_arriba) // Vengo de arriba, coincide arriba->vert = yo
					direccion = DIRECCION_VERTICAL;
    			else
    				direccion = DIRECCION_DIAGONAL; // Sigo por la diagonal
    			break;
    	}

    	// Desplazamientos: se harán en función de de dónde vine, no el movimiento actual (que es el anterior)
    	switch (direccion) {
    		case DIRECCION_VERTICAL:
				Obj->alineamiento_horizontal[pos_actual] = '-';
				Obj->alineamiento_vertical[pos_actual] = letraSecuenciaVertical;
				//printf("Estoy en %d/%d, voy en VERTICAL\n", i, j);
				j--;
    		break;
    		case DIRECCION_HORIZONTAL:
				Obj->alineamiento_horizontal[pos_actual] = letraSecuenciaHorizontal;
				Obj->alineamiento_vertical[pos_actual] = '-';
				//printf("Estoy en %d/%d, voy en HORIZONTAL\n", i, j);
				i--;
    			break;
    		default:
				Obj->alineamiento_horizontal[pos_actual] = letraSecuenciaHorizontal;
				Obj->alineamiento_vertical[pos_actual] = letraSecuenciaVertical;
				//printf("Estoy en %d/%d, voy en DIAGONAL\n", i, j);
				i--;
				j--;
    			break;
    	}
		//printf("%c/%c\n", alineamiento_horizontal[pos_actual], alineamiento_vertical[pos_actual]);
		pos_actual++;
	}

	//Control de casos especiales límites: alineamiento que se topa con el borde: comienza con ----
	//Además, puede no ser subtrabajo [0,0] por lo que puede expandirse a otros subtrabajos

	//Límite horizontal == 0 -> subimos para arriba hasta y = 0
	while (Obj->indice_trabajo_x==0 && i==0 && j>0) {
		Obj->alineamiento_horizontal[pos_actual] = '-';
		Obj->alineamiento_vertical[pos_actual] = Obj->secuencia_vertical[j-1];
		j--;
		pos_actual++;
		direccion = DIRECCION_VERTICAL;
	}

	//Límite vertical == 0 -> subimos para la izquierda hasta x = 0
	while (Obj->indice_trabajo_y==0 && j==0 && i>0) {
		Obj->alineamiento_horizontal[pos_actual] = Obj->secuencia_horizontal[i-1];
		Obj->alineamiento_vertical[pos_actual] = '-';
		i--;
		pos_actual++;
		direccion = DIRECCION_HORIZONTAL;
	}

	//Colocamos los finales de línea
	Obj->alineamiento_horizontal[pos_actual]=0x00;
	Obj->alineamiento_vertical[pos_actual]=0x00;

	Obj->direccion_inicial = direccion; // Volvemos a guardar

	Obj->limite_x = (i==0) ? 0 : ++i; //Si ha llegado a 0, se calculan sólo unas pocas de columnas en el siguiente
	Obj->limite_y = (j==0) ? 0 : ++j; //Si ha llegado a 0, se calculan sólo unas pocas de filas en el siguiente
	//El incremento es necesario para indicar el número a calcular [1..] no posición [0..]

/*    // DEPURACIÃâ€œN
	printf("Alineamientos:\nH: %s\nV: %s\n", alineamiento_horizontal, alineamiento_vertical);
	printf("Ãï¿½ndices finales:\nH: %d\nV: %d\n", limite_x, limite_y); */
}
#endif

#ifdef SMITHWATERMAN //Función para Smith-Waterman
/** Cálculo de las secuencias de alineamiento a partir de una submatriz y la esquina inferior derecha de la misma. */
void calcula_alineamiento(estado_del_Objeto * Obj, struct Nodo_valor** matriz, int num_filas, int num_columnas){
	int i=num_columnas-1, j=num_filas-1; //Ya que la posición tamaÃÂ±o está fuera de los índices
	//Punteros a nodo para mayor optimización: ahorro de acceso a array con offset para cada valor de la estructura
	struct Nodo_valor *temp, *temp_arr, *temp_izq;

	/* Alineamientos EN MEMORIA LOCAL y su contador. */
	Obj->alineamiento_vertical = (char *) malloc_safe_agr248 (sizeof(char) * (num_filas + num_columnas + 1));
	Obj->alineamiento_horizontal = (char *) malloc_safe_agr248 (sizeof(char) * (num_filas + num_columnas + 1));
	int pos_actual=0;

	// Implementación vuelta atrás: no busca máximo, sino de dónde viene realmente el valor

	int fin_calculo = 0; //Control del fin del cálculo: llegamos a un 0
	char letraSecuenciaHorizontal, letraSecuenciaVertical;
	int direccion = Obj->direccion_inicial;
	int temp_valor_arriba, temp_valor_izquierda, temp_arriba, temp_izquierda;

	//Mientras no lleguemos a ningún límite, vamos desplazando y aÃÂ±adiendo letras
	while ((i>0) && (j>0) && !fin_calculo) {
		letraSecuenciaHorizontal = Obj->secuencia_horizontal[i-1];
    	letraSecuenciaVertical = Obj->secuencia_vertical[j-1];

    	temp = &(matriz[i][j]);
    	temp_arr = &(matriz[i][j-1]);
    	temp_izq = &(matriz[i-1][j]);

        // Como en las tuplas guardamos el vertical y horizontal de las siguientes celdas (ver diagramas Gotoh con diferencias),
        // para ver de dónde hemos venido, tenemos que precalcular y mirar en los horizontal/vertical de anteriores celdas
        // Además, hay que atender a la dirección de la que vengo para controlar los casos:
        // Vertical: ÃÂ¿Vengo realmente del anterior vertical o fue el principio del gap y me voy en diagonal?
        // Horizontal: ÃÂ¿Vengo realmente del anterior horizontal o fue el principio del gap y me voy en diagonal?
        // Diagonal: ÃÂ¿De dónde vine? Comparar con los valores vertical y horizontal regenerados (se guardan diferencias) y si son iguales, es por ahí, sino, diagonal again
    	switch (direccion) {
    		case DIRECCION_VERTICAL:
#ifdef CLUSTALW
    			temp_valor_arriba = temp_arr->vertical - GAPEXTEND_COST; // Valor - GAP_EXTEND
    			temp_arriba = temp->vertical; // Ni diferencia ni ná de ná
#else
    			temp_valor_arriba = temp_arr->valor - (temp_arr->vertical + GAPEXTEND_COST); // Regeneramos el valor desde la diferencia - GAP_EXTEND
    			temp_arriba = temp->valor - temp->vertical; // Regeneramos el valor de la diferencia de esta casilla
#endif
        	// INSERCIÃâ€œN: desplazamiento vertical hacia arriba (vertical letra, horizontal hueco)
    			if (temp_arriba == temp_valor_arriba) // Vengo de arriba, coincide arriba->vert + GAP_EXTEND = yo
					direccion = DIRECCION_VERTICAL;
    			else {
    				// No sigo subiendo, vengo del valor máximo de mi tupla, no de la superior
    				// Lo normal es que venga de la diagonal, pero puede darse el caso X-/-X si mismatch > open_gap x 2
    				// Tenemos que verificar que el valor almacenado (máximo) de esta casilla no vino de la izquierda -> si es así, X-/-X, sino, diagonal
#ifdef CLUSTALW
    				temp_valor_izquierda = temp_izq->horizontal; // Ni diferencia ni ná de ná
#else
    				temp_valor_izquierda = temp_izq->valor - temp_izq->horizontal; // Regeneramos el valor de la diferencia
#endif
	    			if (temp->valor == temp_valor_izquierda) // Vengo de la izquierda, coincide izquierda->horiz = yo
						direccion = DIRECCION_HORIZONTAL;
    				else
    					direccion = DIRECCION_DIAGONAL; // Vengo de la diagonal, he abierto (valor = INSERT + GAPEXTEND)
    			}
    			break;
    		case DIRECCION_HORIZONTAL:
#ifdef CLUSTALW
    			temp_valor_izquierda = temp_izq->horizontal - GAPEXTEND_COST; // Valor - GAP_EXTEND
    			temp_izquierda = temp->horizontal; // Ni diferencia ni ná de ná
#else
    			temp_valor_izquierda = temp_izq->valor - (temp_izq->horizontal + GAPEXTEND_COST); // Regeneramos el valor desde la diferencia - GAP_EXTEND
    			temp_izquierda = temp->valor - temp->horizontal; // Regeneramos el valor de la diferencia de esta casilla
#endif
        	// DELETE: desplazamiento horizontal hacia la izquierda (vertical hueco, horizontal letra)
    			if (temp_izquierda == temp_valor_izquierda) // Vengo de la izquierda, coincide izquierda->horiz + GAP_EXTEND = yo
					direccion = DIRECCION_HORIZONTAL;
    			else {
    				// No sigo desplazándome a izquierda, vengo del valor máximo de mi tupla, no de la izquierda
    				// Lo normal es que venga de la diagonal, pero puede darse el caso -X/X- si mismatch > open_gap x 2
    				// Tenemos que verificar que el valor almacenado (máximo) de esta casilla no vino de arriba -> si es así, -X/X-, sino, diagonal
#ifdef CLUSTALW
    				temp_valor_arriba = temp_arr->vertical; // Ni diferencia ni ná de ná
#else
    				temp_valor_arriba = temp_arr->valor - temp_arr->vertical; // Regeneramos el valor de la diferencia
#endif
	    			if (temp->valor == temp_valor_arriba) // Vengo de arriba, coincide arriba->vert = yo
						direccion = DIRECCION_VERTICAL;
    				else
    					direccion = DIRECCION_DIAGONAL; // Vengo de la diagonal, he abierto (valor = DELETE + GAPEXTEND)
    			}
    			break;
    		default:
#ifdef CLUSTALW
    		 	// No hay diferencia que regenerar (si soy igual a ellos, vengo de ahí, en otro caso, diagonal)
    			temp_valor_izquierda = temp_izq->horizontal;
    			temp_valor_arriba = temp_arr->vertical;
#else
    		 	// Regeneramos los valores desde la diferencia (si soy igual a ellos, vengo de ahí, en otro caso, diagonal)
    			temp_valor_izquierda = temp_izq->valor - temp_izq->horizontal;
    			temp_valor_arriba = temp_arr->valor - temp_arr->vertical;
#endif
        	// MATCH / REPLACE: vertical letra, horizontal letra
    			if (temp->valor == temp_valor_izquierda) // Vengo de la izquierda, coincide izquierda->horiz = yo
					direccion = DIRECCION_HORIZONTAL;
    			else if (temp->valor == temp_valor_arriba) // Vengo de arriba, coincide arriba->vert = yo
					direccion = DIRECCION_VERTICAL;
    			else
    				direccion = DIRECCION_DIAGONAL; // Sigo por la diagonal
    			break;
    	}

    	// Desplazamientos: se harán en función de de dónde vine, no el movimiento actual (que es el anterior)
    	switch (direccion) {
    		case DIRECCION_VERTICAL:
				Obj->alineamiento_horizontal[pos_actual] = '-';
				Obj->alineamiento_vertical[pos_actual] = letraSecuenciaVertical;
				//printf("Estoy en %d/%d, voy en VERTICAL\n", i, j);
				j--;
    		break;
    		case DIRECCION_HORIZONTAL:
				Obj->alineamiento_horizontal[pos_actual] = letraSecuenciaHorizontal;
				Obj->alineamiento_vertical[pos_actual] = '-';
				//printf("Estoy en %d/%d, voy en HORIZONTAL\n", i, j);
				i--;
    			break;
    		default:
				Obj->alineamiento_horizontal[pos_actual] = letraSecuenciaHorizontal;
				Obj->alineamiento_vertical[pos_actual] = letraSecuenciaVertical;
				//printf("Estoy en %d/%d, voy en DIAGONAL\n", i, j);
				i--;
				j--;
    			break;
    	}
		//printf("%c/%c\n", alineamiento_horizontal[pos_actual], alineamiento_vertical[pos_actual]);
		pos_actual++;

		fin_calculo = (matriz[i][j].valor == 0);
	}

	//No es necesario el control de casos especiales límites: alineamiento que se topa con el borde: comienza con ----
	//siempre es inicializado con [0,0,0] en los bordes iniciales

	Obj->direccion_inicial = direccion; // Volvemos a guardar

	//Colocamos los finales de línea
	Obj->alineamiento_horizontal[pos_actual]=FINCAD;
	Obj->alineamiento_vertical[pos_actual]=FINCAD;

	Obj->limite_x = (i==0) ? 0 : ++i; //Si ha llegado a 0, se calculan sólo unas pocas de columnas en el siguiente
	Obj->limite_y = (j==0) ? 0 : ++j; //Si ha llegado a 0, se calculan sólo unas pocas de filas en el siguiente
	//El incremento es necesario para indicar el número a calcular [1..] no posición [0..]

	// En cambio si hemos llegado a fin_calculo, pasamos -1
	if (fin_calculo) {
		Obj->limite_x = -1;
	}
}
#endif

/** Entrega un puntero al alineamiento ya calculado (en memoria compartida) por rawchannel sumidero al controlador. */
void entregar_alineamiento(int rango_tile, estado_del_Objeto * Obj){
	//Creamos el nodo en memoria compartida para que sea accesible por el controlador
	Nodo_alineamiento* alineamiento = (struct Nodo_alineamiento *)malloc_shared(sizeof (struct Nodo_alineamiento));

	//Asignamos los alineamientos y los datos que necesita el controlador
	//Los alineamientos los copiamos a memoria compartida
	int size_horizontal = sizeof(char) * (strlen(Obj->alineamiento_horizontal) + 1);
	alineamiento->fragmento_horizontal = (char*) malloc_shared (size_horizontal);
	memcpy(alineamiento->fragmento_horizontal, Obj->alineamiento_horizontal, size_horizontal);
	int size_vertical = sizeof(char) * (strlen(Obj->alineamiento_vertical) + 1);
	alineamiento->fragmento_vertical = (char*) malloc_shared (size_vertical);
	memcpy(alineamiento->fragmento_vertical, Obj->alineamiento_vertical, size_vertical);

	//alineamiento->fragmento_horizontal=alineamiento_horizontal;
	//alineamiento->fragmento_vertical=alineamiento_vertical;

	alineamiento->x = Obj->indice_trabajo_x;
	alineamiento->y = Obj->indice_trabajo_y;
	alineamiento->limite_x = Obj->limite_x;
	alineamiento->limite_y = Obj->limite_y;
	alineamiento->indice_worker = rango_tile;
	alineamiento->direccion = Obj->direccion_inicial;
	#ifdef SMITHWATERMAN
		alineamiento->valor_minimo_x = Obj->pos_min_x;
		alineamiento->valor_minimo_y = Obj->pos_min_y;
	#endif

	/* Limpiamos la memoria de todo aquello que hubiÃÂ©semos copiado a memoria local */
	// Fragmentos de secuencias en memoria local
	free(Obj->secuencia_horizontal);
	free(Obj->secuencia_vertical);

	// Los 2 fragmentos en memoria local (los de trabajo no se utilizan)
	libera_fragmento(&(Obj->trozoSubColumnaPartida));
 	libera_fragmento(&(Obj->trozoSubFilaPartida));
 	// Los alineamientos
 	free(Obj->alineamiento_horizontal);
 	free(Obj->alineamiento_vertical);

	// Tabla de puntuaciones con su alfabeto y codificaciones
	matrix_free(&(Obj->tabla_puntuacion));

	//Enviamos los resultados al núcleo (es un envío NO bloqueante, si el controlador está ocupado, se queda en su buffer)
	pthread_mutex_lock(&sumidero_mutex);
		while(sumidero.datos_disponibles){
			pthread_cond_wait(&sumidero_libre_cond, &sumidero_mutex);
		}
	sumidero.datos.alineamiento = alineamiento;
	sumidero.datos_disponibles=1;
	pthread_cond_signal(&sumidero_datos_disponibles_cond);
	pthread_mutex_unlock(&sumidero_mutex);
}

/** Calcula el valor de puntuación para los dos caracteres en función de la matriz de sustituciones. */

long calcula_matchReplace(estado_del_Objeto * Obj, unsigned char letraSecuenciaHorizontal, unsigned char letraSecuenciaVertical){
	//Direccionamiento indirecto: codificacion[letra]=índice de dicha letra en la matriz
	//Las 3 líneas son lo mismo que el return pero en 3 líneas más claras pero menos óptimo
//	unsigned short x = tabla_puntuacion->codificacion_horizontal[letraSecuenciaHorizontal];
//	unsigned short y = tabla_puntuacion->codificacion_vertical[letraSecuenciaVertical];
//	return tabla_puntuacion->matriz[x][y];
	// DONE: Esta función era un cuello de botella. La transformación de las letras se hace al recibir cada trabajo a procesar
	// con lo que se reduce el número de traducciones drásticamente.
//	return Obj->tabla_puntuacion->matriz
//		[Obj->tabla_puntuacion->codificacion_horizontal[letraSecuenciaHorizontal]]
//		[Obj->tabla_puntuacion->codificacion_vertical[letraSecuenciaVertical]];
	return Obj->tabla_puntuacion->matriz
		[letraSecuenciaHorizontal]
		[letraSecuenciaVertical];
	// Versión directa en memoria: 12secs para 4k
//	return tabla_puntuacion_default[indice(alfabeto_default, letraSecuenciaHorizontal)][indice(alfabeto_default, letraSecuenciaVertical)];
}


/** Devuelve el máximo entre 2 números. */
long max2(long a, long b){
	return (a>b)?a:b;
}

/** Devuelve el máximo entre 3 números. */
long max3(long a, long b, long c){
	#ifndef SMITHWATERMAN
		return (a>b)?((a>c)?a:c):((b>c)?b:c);
	#endif
	/** Devuelve el máximo entre 3 números y un 0. */
	#ifdef SMITHWATERMAN
		return (a>b)?((a>c)?((a>0)?a:0):((c>0)?c:0)):((b>c)?((b>0)?b:0):((c>0)?c:0));
	#endif
}

/** Devuelve el índice que es el máximo entre 3 números. */
inline int max_indice(long a, long b, long c){
    return (a>b)?((a>c)?1:3):((b>c)?2:3);
}


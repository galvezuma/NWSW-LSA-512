/*
 * job_table.c
 *
 *  Created on: Mar 8, 2021
 *      Author: galvez
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <sys/types.h>
#include <malloc.h>
#include <pthread.h>
#include "job_table.h"
#include "scif.h"

extern char * prefijo_fichero_temporal; // La necesito de controller.h
static pthread_mutex_t host_comm_mutex;
static pthread_t * ultimo_hilo_no_liberado = NULL;

/* Implementación del módulo de la tabla de trabajos. */

/** Crea una nueva tabla de trabajos en memoria, reservando memoria física para ella.
 * Devuelve un puntero a la tabla cuya memoria ya ha sido reservada. */
 struct Nodo_trabajo ** crea_tabla_trabajos(int anchura, int altura){
 	 struct Nodo_trabajo ** tabla_trabajos = (struct Nodo_trabajo **) malloc_safe_agr248(sizeof(struct Nodo_trabajo *) * anchura);
	 if (tabla_trabajos == NULL) //No hay memoria suficiente
		return NULL;
     for (unsigned int i = 0; i < anchura; i++) {
       	tabla_trabajos[i] = (struct Nodo_trabajo *) malloc_safe_agr248(sizeof(struct Nodo_trabajo) * altura);

       	if (tabla_trabajos[i] == NULL) { //No hay memoria suficiente
        	/* handle error, free() the previously allocated rows */
           	for (unsigned int x = 0; x < i; x++)
           		free(tabla_trabajos[x]);
           	free(tabla_trabajos);
           	return NULL;
       	} else {
       		//Nos aseguramos de poner los punteros a nulos por si las moscas
       		for (unsigned int j = 0; j < altura; j++) {
           		tabla_trabajos[i][j].fragmento_horizontal = NULL;
           		tabla_trabajos[i][j].fragmento_vertical = NULL;
       		}
       	}
   	 }
     pthread_mutex_init(&host_comm_mutex, NULL);
     return tabla_trabajos;
 }

/** Crea una nueva tabla de trabajos en memoria compartida, reservando memoria física para ella en memoria compartida.
 * Devuelve un puntero a la tabla cuya memoria ya ha sido reservada. */
 // Esta funci󮠹a no se usa.
 struct Nodo_trabajo ** crea_tabla_trabajos_compartida(int anchura, int altura){
	struct Nodo_trabajo ** tabla_trabajos = (struct Nodo_trabajo **) malloc_shared(sizeof(struct Nodo_trabajo *) * anchura);
	if (tabla_trabajos == NULL) //No hay memoria suficiente
		return NULL;
    for (unsigned int i = 0; i < anchura; i++) {
        tabla_trabajos[i] = (struct Nodo_trabajo *) malloc_shared(sizeof(struct Nodo_trabajo) * altura);

        if (tabla_trabajos[i] == NULL) { //No hay memoria suficiente
            /* handle error, free() the previously allocated rows */
            for (unsigned int x = 0; x < i; x++)
           		free(tabla_trabajos[x]);
            free(tabla_trabajos);
            return NULL;
       	} else {
       		//Nos aseguramos de poner los punteros a nulos por si las moscas
       		for (unsigned int j = 0; j < altura; j++) {
           		tabla_trabajos[i][j].fragmento_horizontal = NULL;
           		tabla_trabajos[i][j].fragmento_vertical = NULL;
       		}
       	}
    }
    return tabla_trabajos;
 }

/** Libera la memoria una tabla de trabajos del tamaño especificado. */
 void libera_tabla_trabajos(Nodo_trabajo *** ptr_tabla_trabajos, int anchura, int altura){
 	struct Nodo_trabajo ** tabla_trabajos = *ptr_tabla_trabajos;
    for (unsigned int i = 0; i < anchura; i++) {
       //Si se tienen datos, se eliminan
       for (unsigned int j = 0; j < altura; j++) {
       		if (tabla_trabajos[i][j].fragmento_horizontal!=NULL) {
       			libera_fragmento(&(tabla_trabajos[i][j].fragmento_horizontal));
       		}
       		if (tabla_trabajos[i][j].fragmento_vertical!=NULL) {
       			libera_fragmento(&(tabla_trabajos[i][j].fragmento_vertical));
       		}
       }
       free(tabla_trabajos[i]);
    }
    free(tabla_trabajos);
    pthread_mutex_destroy(&host_comm_mutex);
 }

/** Consulta si el fragmento de trabajo horizontal del subtrabajo de coordenadas [x,y] es NULO */
 inline struct Nodo_fragmento* is_fragmento_fila(struct Nodo_trabajo ** tabla_trabajos, int x, int y){
 	return tabla_trabajos[x][y].fragmento_horizontal;
 }


/** Obtiene un puntero al fragmento de trabajo horizontal del subtrabajo de coordenadas [x,y] */
 inline struct Nodo_fragmento* get_fragmento_fila(struct Nodo_trabajo ** tabla_trabajos, int x, int y, int tabla_trabajos_MAX_X, int tabla_trabajos_MAX_Y){
 	if (tabla_trabajos[x][y].fragmento_horizontal == NULL){
 		recuperarFila(tabla_trabajos, y, tabla_trabajos_MAX_X);
 		// Sólo en el caso de las filas y libero una única vez
 		if (y < tabla_trabajos_MAX_Y-1 && tabla_trabajos[x][y+1].fragmento_horizontal!=NULL){
 			liberarFila(tabla_trabajos, y+1, tabla_trabajos_MAX_X);
 		}
 	}
 	return tabla_trabajos[x][y].fragmento_horizontal;
 }

/** Consulta si el fragmento de trabajo vertical del subtrabajo de coordenadas [x,y] es NULO */
 inline struct Nodo_fragmento* is_fragmento_columna(struct Nodo_trabajo ** tabla_trabajos, int x, int y){
		// printf("Asking for fragment %d %d\n",x, y); //, tabla_trabajos, tabla_trabajos[x][y].fragmento_vertical);
 	return tabla_trabajos[x][y].fragmento_vertical;
 }

/** Obtiene un puntero al fragmento de trabajo vertical del subtrabajo de coordenadas [x,y] */
 inline struct Nodo_fragmento* get_fragmento_columna(struct Nodo_trabajo ** tabla_trabajos, int x, int y, int tabla_trabajos_MAX_X, int tabla_trabajos_MAX_Y){
 	if (tabla_trabajos[x][y].fragmento_vertical==NULL){
		recuperarFila(tabla_trabajos, y, tabla_trabajos_MAX_X);
 	}
 	return tabla_trabajos[x][y].fragmento_vertical;
 }

 void * pthread_salvarFila(void * a){
	 struct args_salvarFila * args = (struct args_salvarFila *) a;
	 pthread_mutex_lock(&host_comm_mutex);
	 if (ultimo_hilo_no_liberado != NULL){
		 pthread_join(*ultimo_hilo_no_liberado, NULL);
		 free(ultimo_hilo_no_liberado);
	 }
	 ultimo_hilo_no_liberado = args->hilo_actual; // Current thread
	 salvarFila(args->tablaTrabajos, args->y, args->tabla_trabajos_MAX_X, args->cierra_conexion);
	 pthread_mutex_unlock(&host_comm_mutex);
	 free(args);
	 return NULL;
 }

/** Asigna el fragmento de trabajo vertical al subtrabajo de coordenadas [x,y] */
 inline void set_fragmento_fila(struct Nodo_trabajo ** tabla_trabajos, int x, int y, struct Nodo_fragmento* fragmento_fila, int tabla_trabajos_MAX_X){
 	tabla_trabajos[x][y].fragmento_horizontal=fragmento_fila;
 	if (x == tabla_trabajos_MAX_X-1 && y>0){
		//printf("Saving row %d\n", y-1);
 		pthread_t *s = (pthread_t *)malloc(sizeof(pthread_t));
 		struct args_salvarFila * args = (struct args_salvarFila *) malloc(sizeof(struct args_salvarFila));
 		args->tablaTrabajos = tabla_trabajos;
 		args->y = y-1;
 		args->tabla_trabajos_MAX_X = tabla_trabajos_MAX_X;
 		args->cierra_conexion = 0;
 		args->hilo_actual = s;
 		pthread_create(s, NULL, pthread_salvarFila, (void *)args);
 		// No pthread_exit() is required.
 		// El valor de s lo libera el siguiente hilo.
 		// El último hilo lo lanza el controlador, no es un puntero y no es necesario liberarlo.
 		//salvarFila(tabla_trabajos, y-1, tabla_trabajos_MAX_X, 0);
 	}
 	// printf("y values after saving row: %d\n", y);
 }

/** Asigna el fragmento de trabajo horizontal al subtrabajo de coordenadas [x,y] */
 inline void set_fragmento_columna(struct Nodo_trabajo ** tabla_trabajos, int x, int y, struct Nodo_fragmento* fragmento_columna, int tabla_trabajos_MAX_X){
 	tabla_trabajos[x][y].fragmento_vertical=fragmento_columna;
 }

 static const int PUERTO = 9070;
 static const int NODO_HOST = 0;
 static const int LONG_BUFFER = 64 * 1024;
 static const int PAGE_SIZE = 8192;
 static const int TRANSFER_DONE_FLAG = 1;
 static const int ERROR_RESULT = -1;
 static const int SLEEP_MICROSECONDS = 250000;

 // Esta funci󮠰rocede de https://software.intel.com/en-us/forums/topic/530366
 static void* allocRegisteredMem(scif_epd_t endpoint, size_t num_bytes) {
 	// allocate message buffer (registered memory)
 	void *ptr;
 	int result;
 	ptr = memalign(PAGE_SIZE, num_bytes);
 	if (ptr == NULL) {
 		fprintf(stderr, "failed to allocate message buffer (posix_memalign)\n");
 		exit(EXIT_FAILURE);
 	}

 	// clear message buffer
 	memset(ptr, 0, num_bytes);

 	// register message buffer memory
 	result = scif_register(endpoint, ptr, num_bytes,
 			(off_t)ptr, SCIF_PROT_READ|SCIF_PROT_WRITE,
 			SCIF_MAP_FIXED);
 	if (result == SCIF_REGISTER_FAILED) {
 		fprintf(stderr, "scif_register() failed");
 		exit(EXIT_FAILURE);
 	}

 	return ptr;
 }

 // Salva una fila completa de la tabla de trabajos enviᮤola al host por SCIF.
 void salvarFila(struct Nodo_trabajo ** tablaTrabajos, int y, int tabla_trabajos_MAX_X, int cierra_conexion){
	int i;
	int32_t ok;
	static scif_epd_t mi_end_point_final;
	static void* reg_mem = NULL;
	int longitud, offset;
	const int CONTINUAR = 0;
	const int ACABAR = 1;
	if (reg_mem == NULL) {
		// Hay que iniciar los canales de comunicaci󮊉	int tam_buffer;
		struct scif_portID mi_peer;
		int bytes_enviados;
		int bytes_recibidos;

		mi_end_point_final = scif_open();
		if (mi_end_point_final == SCIF_OPEN_FAILED){
			fprintf(stderr, "Error irrecuperable: no puede abrirse el endpoint.\n");
			return;
		}

		mi_peer.node = NODO_HOST;
		mi_peer.port = PUERTO;
		// Quiz᳠esto deba ser un bucle de conexi󮠨ver documentaci󮠤e la funci󮠳cif_connect).
		ok = scif_connect(mi_end_point_final, &mi_peer);
		if (ok == ERROR_RESULT){
			fprintf(stderr, "Error irrecuperable: problema al conectarse al nodo %d puerto %d. Errno: %d\n", mi_peer.node, mi_peer.port, errno);
			return;
		}
		// Ubica y registra memoria para transferencias RMA
		// Lo siguiente es el tama񯠣orrecto.
		tam_buffer = sizeof(struct Nodo_valor) * (tablaTrabajos[0][y].fragmento_horizontal->longitud +
												  tablaTrabajos[0][y].fragmento_vertical->longitud)
				   + sizeof(int) * 2;
		// Lo siguiente es para cortar por lo sano.
		tam_buffer = LONG_BUFFER;
		reg_mem = allocRegisteredMem(mi_end_point_final, tam_buffer);
		// Env�la direcci󮠬ocal al host
		bytes_enviados = scif_send(mi_end_point_final, &reg_mem, sizeof(void *), SCIF_SEND_BLOCK);
		if (bytes_enviados < sizeof(void *)) {
			fprintf(stderr, "Error irrecuperable: problema al hacer scif_send().\n");
			return;
		}
	}
	for(i=0; i<tabla_trabajos_MAX_X; i++){
		offset = 0;

 		longitud=tablaTrabajos[i][y].fragmento_horizontal->longitud;
		memcpy(reg_mem+offset, &longitud, sizeof(longitud));
		offset += sizeof(longitud);

		memcpy(reg_mem+offset, tablaTrabajos[i][y].fragmento_horizontal->contenido, longitud*sizeof(struct Nodo_valor));
		offset += longitud*sizeof(struct Nodo_valor);

 		longitud=tablaTrabajos[i][y].fragmento_vertical->longitud;
		memcpy(reg_mem+offset, &longitud, sizeof(longitud));
		offset += sizeof(longitud);

		memcpy(reg_mem+offset, tablaTrabajos[i][y].fragmento_vertical->contenido, longitud*sizeof(struct Nodo_valor));
		offset += longitud*sizeof(struct Nodo_valor);

		ok = CONTINUAR;
		// Env�se񡬠de que los datos pueden recogerse
		scif_send(mi_end_point_final, (void *)&ok, sizeof(ok), SCIF_SEND_BLOCK);

		// Espera se񡬠bloqueante de que los datos se han recogido.
		scif_recv(mi_end_point_final, &ok, sizeof(ok), SCIF_RECV_BLOCK);

		// Los datos ya se han enviado y la memoria puede liberarse.
 		libera_fragmento(&(tablaTrabajos[i][y].fragmento_horizontal));
 		libera_fragmento(&(tablaTrabajos[i][y].fragmento_vertical));
 		tablaTrabajos[i][y].fragmento_horizontal = NULL;
 		tablaTrabajos[i][y].fragmento_vertical = NULL;
	}
	if (cierra_conexion){
		// Esto tambi鮠desregistra.
		// Env�se񡬠de que la comunicaci󮠨a acabado.
		ok = ACABAR;
		scif_send(mi_end_point_final, (void *)&ok, sizeof(ok), SCIF_SEND_BLOCK);
		scif_close(mi_end_point_final);
		if (reg_mem != NULL)
			free(reg_mem);
		reg_mem = NULL;
	}

	// printf("Row saved %d\n", y);

	if (TRAZA) {
		printf("%08ld,F,%d,%d\n", my_clock()/CLOCKS_PER_CENTSEC, y, y);
	}
 }

/*
//  Estᠣomentada porque este desarrollo env�los datos por SCIF y no los guarda en un fichero.
// Salva una fila completa de la tabla de trabajos en un fichero del tipo nombreTile1, nombreTile2, etc. donde Y es el numero de fila.
 void salvarFila(struct Nodo_trabajo ** tablaTrabajos, int y, int tabla_trabajos_MAX_X, int dummy_cierra_conexion){
 	int i, longitud;
 	char fileName[300];
 	FILE * file;
 	sprintf(fileName, "%s.%s%d", prefijo_fichero_temporal, nombreFichero, y);
 	if ((file = fopen(fileName, "w"))==NULL) {
    	printf("No se puede escribir en el fichero de fila temporal: %s\n", fileName);
    	exit(1);
  	}
 	for(i=0; i<tabla_trabajos_MAX_X; i++){
 		longitud=tablaTrabajos[i][y].fragmento_horizontal->longitud;
 		fwrite(&longitud, sizeof(int), 1, file);
 		fwrite(tablaTrabajos[i][y].fragmento_horizontal->contenido, sizeof(struct Nodo_valor), longitud, file);
// 		printf("Tamaño: %d [%d + %d]\n", sizeof(struct Nodo_valor), sizeof(long), sizeof(char));

 		longitud=tablaTrabajos[i][y].fragmento_vertical->longitud;
 		fwrite(&longitud, sizeof(int), 1, file);
 		fwrite(tablaTrabajos[i][y].fragmento_vertical->contenido, sizeof(struct Nodo_valor), longitud, file);

 		libera_fragmento(&(tablaTrabajos[i][y].fragmento_horizontal));
 		libera_fragmento(&(tablaTrabajos[i][y].fragmento_vertical));
 		tablaTrabajos[i][y].fragmento_horizontal = NULL;
 		tablaTrabajos[i][y].fragmento_vertical = NULL;
 	}
 	fclose(file);
	// printf("Row saved %d\n", y);

	if (TRAZA) {
		printf("%08ld,F,%d,%d\n", my_clock()/CLOCKS_PER_CENTSEC, y, y);
	}
 }
*/

#pragma GCC diagnostic ignored "-Wunused-result"
/** Recupera una fila completa de la tabla de trabajos a partir de un fichero del tipo nombreTile1, nombreTile2, etc. donde Y es el número de fila. */
  void recuperarFila(struct Nodo_trabajo ** tablaTrabajos, int y, int tabla_trabajos_MAX_X){
 	int i, longitud;
 	char fileName[300];
 	FILE * file;
 	sprintf(fileName, "%s.%s%d", prefijo_fichero_temporal, nombreFichero, y);
 	if ((file = fopen(fileName, "r"))==NULL) {
    	printf("No se puede leer del fichero de fila temporal: %s\n", fileName);
    	exit(1);
  	}
 	for(i=0; i<tabla_trabajos_MAX_X; i++){
 		tablaTrabajos[i][y].fragmento_horizontal=(struct Nodo_fragmento *) malloc_shared(sizeof(struct Nodo_fragmento));
 		tablaTrabajos[i][y].fragmento_vertical=(struct Nodo_fragmento *) malloc_shared(sizeof(struct Nodo_fragmento));

 		fread(&longitud, sizeof(int), 1, file);
 		tablaTrabajos[i][y].fragmento_horizontal->longitud=longitud;
 		tablaTrabajos[i][y].fragmento_horizontal->contenido=(struct Nodo_valor *) malloc_shared(sizeof(struct Nodo_valor)*longitud);
 		fread(tablaTrabajos[i][y].fragmento_horizontal->contenido, sizeof(struct Nodo_valor), longitud, file);

 		fread(&longitud, sizeof(int), 1, file);
 		tablaTrabajos[i][y].fragmento_vertical->longitud=longitud;
 		tablaTrabajos[i][y].fragmento_vertical->contenido=(struct Nodo_valor *) malloc_shared(sizeof(struct Nodo_valor)*longitud);
 		fread(tablaTrabajos[i][y].fragmento_vertical->contenido, sizeof(struct Nodo_valor), longitud, file);
 	}
 	fclose(file);
 }

#pragma GCC diagnostic warning "-Wunused-result"

 void liberarFila(struct Nodo_trabajo ** tablaTrabajos, int y, int tabla_trabajos_MAX_X){
 	int i;
 	for(i=0; i<tabla_trabajos_MAX_X; i++){
 		libera_fragmento(&(tablaTrabajos[i][y].fragmento_horizontal));
 		libera_fragmento(&(tablaTrabajos[i][y].fragmento_vertical));
 		tablaTrabajos[i][y].fragmento_horizontal = NULL;
 		tablaTrabajos[i][y].fragmento_vertical = NULL;
 	}
 }

/** Imprime el estado de la tabla de trabajos (de altura y anchura definidos, puede ser subtabla). */
 void imprime_estado_tabla_trabajos(struct Nodo_trabajo ** tabla_trabajos, int anchura, int altura){
 	printf("  ");
 	for (int i = 0; i < anchura; i++){
 		printf("   %02d: ",i);
 	}
 	printf("\n");
 	for (int j = 0; j < altura; j++){
 		printf("%02d: ",j);
 		for (int i = 0; i < anchura; i++){
 			if (tabla_trabajos[i][j].fragmento_horizontal!=NULL){
 				printf("[X]");
 			} else {
 				printf("[-]");
 			}
			if (tabla_trabajos[i][j].fragmento_vertical!=NULL){
 				printf("[X]");
 			} else {
 				printf("[-]");
 			}
 			printf(" ");
 		}
 		printf("\n");
 	}
 }

/** Escribe la tabla de trabajos a un fichero dado. */
void escribe_tabla_fichero(char* file, struct Nodo_trabajo ** tabla_trabajos, int anchura, int altura){
	FILE *fichero;
	if((fichero = fopen(file, "w"))==NULL) {
    	printf("No se puede crear el fichero para guardar la tabla de trabajos: %s\n", file);
    	exit(1);
  	}

 	for (int i = 0; i < anchura; i++){
// 		fprintf(fichero, "I=%d (%d)\n", i, anchura);
 		for (int j = 0; j < altura; j++){
// 			fprintf(fichero, "J=%d (%d)\n", j, altura);
 			if (tabla_trabajos[i][j].fragmento_horizontal!=NULL){
 				imprime_fragmento_fichero(&fichero, tabla_trabajos[i][j].fragmento_horizontal);
 			}
			if (tabla_trabajos[i][j].fragmento_vertical!=NULL){
 				imprime_fragmento_fichero(&fichero, tabla_trabajos[i][j].fragmento_vertical);
			}
 		}
 	}
	fclose(fichero); //Cerramos el fichero
}

void imprime_fragmento_fichero(FILE **fichero, Nodo_fragmento* fragmento){
	fprintf(*fichero, "[");
  	for (unsigned int i = 0; i<fragmento->longitud; i++){
  		struct Nodo_valor valor = fragmento->contenido[i];
  		fprintf(*fichero, "(%d,%d,%d),", valor.valor, valor.horizontal,valor.vertical);
  	}
  	fprintf(*fichero, "]\n");
}


/*
 * NWSW_LSA_512_controller.c
 *
 *  Created on: Mar 8, 2021
 *      Author: galvez
 */

#include <time.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include "MyLib.h"
#include "priority_queue.h"
#include "job_table.h"
#include "NWSW_LSA_512_worker.h"
#include "NWSW_LSA_512_controller.h"

/* Variables privadas
 */
pthread_t * workers_set = NULL;

/** Controlador del algoritmo paralelo Needleman-Wunsch. Éste será ÃÅ¡NICO y poseerá la tabla de trabajos en
 * memoria compartida, así como cualquier dato de uso general por el resto de procesos. */

/** Lanza un algoritmo Needleman-Wunsch paralelo, haciendo uso del proceso controlador del Needleman-Wunsch, que inicializará
 *  la matriz de resultados, colocará los trabajos y lanzará y coordinará a los núcleos, devolviendo las cadenas de alineamiento. */
void MC64NWSW_Controlador(char * prefijo, int global_bool, estado_algoritmo stage, char* horizontal, char* vertical, Nodo_resultado * resultado,
	short valor_insert, short valor_delete, int valor_matchreplace, short valor_gapextend, char* matriz_puntuaciones,
	int valor_k_horizontal, int valor_k_vertical, int num_tiles){

	/* Variables para cálculos de benchmarking */
	long inicio, intermedio, final;
	long inicio_backward;

	/* Guardamos las configuraciones del algoritmo en constantes globales */
	INSERT_COST = valor_insert;
	DELETE_COST = valor_delete;
	MATCHREPLACE_COST = valor_matchreplace;
	GAPEXTEND_COST = valor_gapextend;
	SIZE_SUBTRABAJO_H = valor_k_horizontal;
	SIZE_SUBTRABAJO_V = valor_k_vertical;
	NUMERO_TILES = num_tiles;

	/* Guardamos el tipo de ejecución. */
	global = global_bool;
	estado = stage;
	if (estado == BACKWARD)
		NUMERO_TILES = 1; // Me da igual el parámetro que me pasen, con uno me basta y sobra

	/* Antes de nada, creamos la tabla de trabajos, en memoria NO COMPARTIDA */
    int longHorizontal=strlen(horizontal);
	int longVertical=strlen(vertical);
	tabla_trabajos_MAX_X = total_subtrabajos(longHorizontal, SIZE_SUBTRABAJO_H);
	tabla_trabajos_MAX_Y = total_subtrabajos(longVertical, SIZE_SUBTRABAJO_V);

	tabla_trabajos = crea_tabla_trabajos(tabla_trabajos_MAX_X, tabla_trabajos_MAX_Y);

	/* Creamos las colas de procesos y de trabajos: tamaño máximo no puede ser menor que el mínimo
	 * (+1 porque así lo requiere la implementación o se reemplazará el primer valor). */
	cola_procesos = create_queue(max(NUMERO_TILES+1, MinQueueSize));
	//Máximo número de trabajos: diagonal mayor
	cola_trabajos = createQueuePrty(max(max(tabla_trabajos_MAX_X, tabla_trabajos_MAX_Y)+1, MinQueueSize));

	//Guardamos los valores de las referencias
	secuencia_horizontal = horizontal;
	secuencia_vertical = vertical;

	//Creamos la matriz para el alfabeto en función del fichero: sino, usa valores por defecto
	if (matriz_puntuaciones==NULL){
		init_matrix_default(&tabla_puntuacion, MATCHREPLACE_COST);
	} else {
		parse_matrix(matriz_puntuaciones, &tabla_puntuacion, MATCHREPLACE_COST);
	}

	/* NUEVA VERSIÓN: Entre fase y fase, se liberan todos los tiles por si se necesita ejecutar otro proceso.
	 * Así pues, ahora lanzaremos y relanzaremos los tiles entre fase y fase. */

	char * fichero_temporal_name;
	int size_fichero_temporal = strlen(prefijo) + 1 + strlen(nombreFichero) + 4;
	fichero_temporal_name = (char *)malloc_safe_agr248(sizeof(char) * size_fichero_temporal);
	sprintf(fichero_temporal_name, "%s.%sData", prefijo, nombreFichero);
    if (estado == FORWARD) {
/* FASE 1: CÁLCULO DE LAS SUBMATRICES */

		/* Inicializamos el máximo: 0 en [0,0] */
		maximo_valor = maximo_indice_x = maximo_indice_y = maximo_posicion_x = maximo_posicion_y = 0;

		//Inicializamos los primeros trabajos
		inicializar_trabajos(longVertical, longHorizontal);
		//Encolamos el primero, el [0,0]
    	encolar_trabajo(0,0);

	    //Comportamiento del clock(): devuelve ciclos que lleva corriendo la aplicación
	    inicio = my_clock();
	    if (TIEMPOS) {
	    	printf("INICIO: %.2f\n", (float)inicio/CLOCKS_PER_SEC);
	    }

		if (TRAZA) {
			printf("S-topped W-orking B-block finalized C-ontroller D-edicated\n");
			printf("%08ld\n", (my_clock()-inicio)/CLOCKS_PER_CENTSEC);
			printf("%d\n", tabla_trabajos_MAX_X);
			printf("%d\n", tabla_trabajos_MAX_Y);
			printf("%d\n", SIZE_SUBTRABAJO_H);
			printf("%d\n", SIZE_SUBTRABAJO_V);
			printf("%08ld,C,0,0\n", (my_clock()-inicio)/CLOCKS_PER_CENTSEC); //Núcleo dedicado el 0,0
			printf("%08ld,D,0,2\n", (my_clock()-inicio)/CLOCKS_PER_CENTSEC); //Núcleo dedicado el 0,2
		}

	    //Realiza el algoritmo en sí: el coordinador irá entregando y recibiendo trabajos
	    //hasta que se lleva a cabo el último trabajo, [tabla_trabajos_MAX_X-1, tabla_trabajos_MAX_Y-1]

	    //Siempre realizamos la primera fase, aunque fuese un único cuadro, se calculará dos veces todo para tener el máximo
	    //y no diferenciar casos entre global/local: esto es un algoritmo paralelo, mal parametrizado no será igual
	    inicializar_hilos(NUMERO_TILES);
	    coordinar_tiles();
	    finalizar_hilos();

	    intermedio = my_clock();
	    if (TIEMPOS) {
		    printf("FIN FASE 1 : %.3f [TOTAL: %.3f]\n", (float)intermedio/MILLIS_PER_SEC, (float)(intermedio-inicio)/MILLIS_PER_SEC);
	    }

		/* Inicializamos el máximo (valor, índice, posición) a partir del fichero */
		salva_maximos_fichero(fichero_temporal_name, maximo_valor, maximo_indice_x, maximo_indice_y, maximo_posicion_x, maximo_posicion_y, inicio, intermedio);
    }
    if (estado == BACKWARD) {
/* FASE 2: CÁLCULO DEL ALINEAMIENTO CON RETROCESO */
		inicio_backward = my_clock();
		/* Inicializamos el máximo (valor, índice, posición) a partir del fichero */
		recupera_maximos_fichero(fichero_temporal_name, &maximo_valor, &maximo_indice_x, &maximo_indice_y, &maximo_posicion_x, &maximo_posicion_y, &inicio, &intermedio);

	    //Antes de nada, encolamos el trabajo con el máximo, que será el primero para calcular el alineamiento
	    //En el caso de alineamiento global, será una de las "paredes" finales
	    //En el caso del alineamiento local, el alineamiento comenzará en el trabajo que tenga el máximo
    	encolar_trabajo(maximo_indice_x,maximo_indice_y);

    	//Calcula el alineamiento haciendo uso de los trabajadores y coordinando la dirección del algoritmo
    	inicializar_hilos(NUMERO_TILES); // Será = 1
    	controller_calcula_alineamiento(&(resultado->datos.alineamientos.alineamiento_query), &(resultado->datos.alineamientos.alineamiento_subject),
    		&(resultado->datos.alineamientos.align_s1_start), &(resultado->datos.alineamientos.align_s1_end),
    		&(resultado->datos.alineamientos.align_s2_start), &(resultado->datos.alineamientos.align_s2_end));
    	finalizar_hilos();

		// Borramos los ficheros temporales creados para salvar datos de la tabla de trabajos en SSD + fichero de datos
		char comando[100];
		sprintf(comando, "rm %s.%s*", prefijo, nombreFichero);

		#pragma GCC diagnostic ignored "-Wunused-result"
		system(comando);
		#pragma GCC diagnostic warning "-Wunused-result"
    }

    //Guardamos el resultado final de la puntuación
    resultado->datos.puntuacion_total = maximo_valor;

	// Limpiamos memoria: colas y tabla trabajos
	dispose_queue(cola_procesos);
	disposeQueuePrty(cola_trabajos);
	libera_tabla_trabajos(&tabla_trabajos, tabla_trabajos_MAX_X, tabla_trabajos_MAX_Y);

    // Guardamos tiempos
 	resultado->inicio = inicio;
 	resultado->intermedio = intermedio;
	if (estado == BACKWARD) {
		#pragma GCC diagnostic ignored "-Wmaybe-uninitialized"
	    final = my_clock() - inicio_backward + intermedio; // Hay que sincronizar final con intermedio, ya que al ser dos programas, el clock se resetea al reiniciarse >_<, de he resumarlo
		#pragma GCC diagnostic warning "-Wmaybe-uninitialized"
	    resultado->final = final;
	    if (TIEMPOS) {
	    	printf("FIN: %.3f [TOTAL: %.3f]\n", (float)final/MILLIS_PER_SEC, (float)(final-inicio)/MILLIS_PER_SEC);
	    }
	} else {
		printf("Inicio: %ld\nIntermedio: %ld\nDiff: %.3f\n", inicio, resultado->intermedio, (float)(resultado->intermedio - inicio)/MILLIS_PER_SEC);
	}
}

/** Lanza al resto de procesadores y les asigna el código binario de "NW_worker" (o "SW_worker" en caso de ser local). */
void inicializar_hilos(int num_tiles){

	// Verificamos que hay tiles suficientes: DE MOMENTO HE QUITADO EL USO DE SHARED_TILES
	int restantes = mylib_proc_remaining();
/*	if (num_tiles > restantes - SHARED_TILES) {
    	ilib_die("Error, sólo hay %d tiles disponibles, contando los %d reservados para el sistema, no puede lanzarse el algoritmo.\n", restantes, SHARED_TILES);
	} */
	if (num_tiles > restantes) {
    	mylib_die("Error, sólo hay %d tiles disponibles, no puede lanzarse el algoritmo.\n", restantes);
	}
	//printf("%d tiles disponibles\n", restantes);

	//Creamos la estructura para los parámetros: lanzaramos num_tiles procesos con el binario correspondiente
	mylibProcParam params;

	if (estado == FORWARD)
		params.t_params.estado = '1'; // Primera fase
	else
		params.t_params.estado = '2'; // Segunda fase


	if (global) {
		params.nameOfFunction = worker_main;
	} else {
		params.nameOfFunction = worker_main;
	}

	/* Arranque del procesamiento paralelo
	 */
	inicializar_mutexes(num_tiles);
	workers_set = mylib_proc_spawn(num_tiles, params);
	if (workers_set == NULL){
    	mylib_die("Error crítico, el algoritmo no ha podido lanzar los procesos paralelos.\n");
	}

	//Abrimos los canales de comunicación de los procesos para utilizar el canal común SINK
	define_canal_sumidero(); // Llamada de efecto nulo

	//Insertamos todos los procesos en la cola, que se llenará
	for (int i=0;i<num_tiles;i++){
		push(cola_procesos, i);
	}
}

 /** Bucle principal de la fase 1 que entrega y recibe trabajos entre los procesadores hasta que se reciba el último trabajo.
  *  Su comportamiento será el siguiente:
  *   - Mientras no recibamos el último fragmento (que sólo lo tendrá un procesador y será siempre el último), continuamos
  *   - Dentro del bucle, priorizamos la entrega de trabajos disponibles al resto de núcleos ociosos
  *   - En caso de que no haya trabajos o núcleos disponibles, quedamos a la espera de recibir nuevos trabajos finalizados
  * */
void coordinar_tiles(){
	int fin = 0;
	while (!fin){
		//Mientras haya trabajo disponible, seguimos enviando trabajos
//		printf("Cola procesos ->");
//		print_queue(cola_procesos);
//		printf("Cola trabajos -> ");
//		print_queue(cola_trabajos);
		while (!is_empty(cola_procesos) && !isEmptyQueuePrty(cola_trabajos)){
			enviar_trabajo(0, 0); //[0,0] para indicar que hay que calcular la submatriz completa guardando sólo los finales (fase 1)
		}
		//printf("Espero recibir trabajo hecho.\n");
		fin = recibir_trabajo();
		//printf("Recibo trabajo hecho.\n");
	}
	//printf("Acabo de coordinar tiles\n");
}

 /** Encola un trabajo de coordenadas [x,y] a la cola de trabajos pendientes. */
 void encolar_trabajo(int x, int y){
	//El trabajo será de la forma Y * maxX + x (ej: [2,0] -> 2, [4,3] para maxX=10 -> 34)
	int indice_trabajo = x + y * tabla_trabajos_MAX_X;
	pushQueuePrty(cola_trabajos, indice_trabajo);
//	printf("-- %08ld -- ", clock()/CLOCKS_PER_CENTSEC); printQueuePrty(cola_trabajos);
 }

/** Desencola un trabajo de la cola de trabajos pendientes. */
void desencolar_trabajo(int* x, int* y){
	//El trabajo será de la forma Y * maxX + x (ej: [2,0] -> 2, [4,3] para maxX=10 -> 34)
	int indice_trabajo = popQueuePrty(cola_trabajos);
	*y = indice_trabajo / tabla_trabajos_MAX_X; //división entera: coordenada Y
	*x = indice_trabajo % tabla_trabajos_MAX_X; //resto división entera: coordenada X
//	printf("-- %08ld -- ", clock()/CLOCKS_PER_CENTSEC); printQueuePrty(cola_trabajos);
}

/** Inicializa la tabla de trabajos. */
void inicializar_trabajos(int longVertical, int longHorizontal){
	/* Construir el juego inicial de trozos activos */
    /* Construcción de fragmentos verticales de la columna 0 (inicialización - phase.1) */
    for (int y=1; y<=longVertical; y += SIZE_SUBTRABAJO_V){
        /* Construir un trozo. Cuidado que el último trozo puede tener una longitud distinta a la normal */
        int longColumna=1+min(SIZE_SUBTRABAJO_V, longVertical-y+1);
        // Inicializamos el trozo vertical (horizontal = 0)
        Nodo_fragmento* trozoVertical = crea_fragmento_shared(longColumna);

		if (global) {
        	inicializa_fragmento(trozoVertical, 0, - INSERT_COST - GAPEXTEND_COST);
        	//La primera posición vale [0,0,0]
        	if (y==1) {
 				trozoVertical->contenido[0].valor = 0;
 				trozoVertical->contenido[0].horizontal = 0;
 				trozoVertical->contenido[0].vertical = 0;
        	}
        } else { //Smith-Waterman inicializa a 0 todo
        	inicializa_fragmento(trozoVertical, 0, 0);
    	}

        /* IMPORTANTE: Un trozo tiene un valor de y superior en 1 a donde realmente empieza.
         * Se supone que la posición y-1 es una diagonal necesaria sólo por las características propias del algoritmo */
        /* En los primeros trozos que se añaden no hay que controlar la eliminación de los que ya no sirvan,
         * y es por ello por lo que se hace directamente sin invocar a añadirTrozo().
         * TÃÂ©ngase en cuenta que para poder proceder así, hay que hacerlo antes de activar los Tiles. */
        int posFrag = total_subtrabajos(y-1, SIZE_SUBTRABAJO_V); //FragV[0,Y] donde Y es el trabajo de [0, pos/TOTAL]
        set_fragmento_columna(tabla_trabajos, 0, posFrag, trozoVertical, tabla_trabajos_MAX_X); //Añadimos los trabajos a la tabla
    }

	/* Construcción de fragmentos horizontales de la fila 0 del mismo modo (inicialización - phase.2) */
    for (int x=1; x<=longHorizontal; x += SIZE_SUBTRABAJO_H){
    	int longFila=1+min(SIZE_SUBTRABAJO_H, longHorizontal-x+1);
    	// Inicializamos el trozo horizontal (horizontal = 1)
        Nodo_fragmento* trozoHorizontal = crea_fragmento_shared(longFila);

        if (global) {
        	inicializa_fragmento(trozoHorizontal, 1, - DELETE_COST - GAPEXTEND_COST);
        	//La primera posición vale [0,0,0]
        	if (x==1) {
 				trozoHorizontal->contenido[0].valor = 0;
 				trozoHorizontal->contenido[0].horizontal = 0;
 				trozoHorizontal->contenido[0].vertical = 0;
        	}
        } else { //Smith-Waterman inicializa a 0 todo
        	inicializa_fragmento(trozoHorizontal, 1, 0);
    	}

        int posFrag = total_subtrabajos(x-1, SIZE_SUBTRABAJO_H); //FragH[X,0] donde X es el trabajo de [0, pos/TOTAL]
        set_fragmento_fila(tabla_trabajos, posFrag, 0, trozoHorizontal, tabla_trabajos_MAX_X); //Añadimos los trabajos a la tabla
    }
}

/** Define el canal sumidero (sink) y establece como emisores del mismo a todos los procesos del grupo grupo_trabajadores.
 * El sumidero se define en 3 fases:
 * ilib_rawchan_start_sink() define el receptor (proceso 0 de ILIB_GROUP_SIBLINGS);
 * ilib_rawchan_add_sink_sender() para definir a cada emisor (uso por cada uno);
 * ilib_rawchan_finish_sink() para finalizar la definición. */
void define_canal_sumidero(){
	/*
  int receptor = ilib_group_rank(ILIB_GROUP_SIBLINGS); //siempre será 0 pero por si las moscas
  ilibSink sink;

  // Creamos el canal: receptor, proceso 0 (controlador)
  ilib_rawchan_start_sink(ILIB_GROUP_SIBLINGS, receptor, CANAL_SINK, &sink);

  // Comenzamos a añadir todos los emisores, de 0 a NUMERO_TILES del grupo grupo_trabajadores
  for (int emisor = 0; emisor < NUMERO_TILES; emisor++)
    ilib_rawchan_add_sink_sender(grupo_trabajadores, emisor, CANAL_SINK, &sink);

  // Finalizamos la definición
  ilib_rawchan_finish_sink(&sink);

  //Abrimos el canal sumidero por el lado del receptor
  ilib_rawchan_open_receiver(CANAL_SINK, puerto_recepcion);
  */
}

/** Envía el próximo trabajo disponible al primer núcleo de la cola de trabajos.
 *  Los parámetros limite_x y limite_y sólo son utilizados en la fase de retroceso, donde nos puede interesar
 *  calcular sólo una parte del trabajo. Durante la primera parte tomarán valor [0,0] indicando que se calcula todo. */
void enviar_trabajo(int limite_x, int limite_y){
 	//Extraemos un trabajo, de índice [x,y]
 	int x, y;
 	desencolar_trabajo(&x, &y);

 	//Extraemos un proceso, de rango rango_proceso
 	int rango_proceso = pop(cola_procesos);

	pthread_mutex_lock(peer_to_peer_mutex+rango_proceso);
	 //printf("Bloqueo controlador para enviar\n");
		//Asignamos columna y fila del trabajo y coordenadas del trabajo
		Nodo_trabajo_completo * trabajo = &(peer_to_peer[rango_proceso].datos.trabajo);
		trabajo->fragmento_horizontal = get_fragmento_fila(tabla_trabajos, x, y, tabla_trabajos_MAX_X, tabla_trabajos_MAX_Y);
		trabajo->fragmento_vertical = get_fragmento_columna(tabla_trabajos, x, y, tabla_trabajos_MAX_X, tabla_trabajos_MAX_Y);
		trabajo->x = x;
		trabajo->y = y;
		trabajo->limite_x = limite_x;
		trabajo->limite_y = limite_y;
		trabajo->secuencia_horizontal = secuencia_horizontal;
		trabajo->secuencia_vertical = secuencia_vertical;
		trabajo->tabla_puntuacion = tabla_puntuacion;
		trabajo->INSERT_COST = INSERT_COST;
		trabajo->DELETE_COST = DELETE_COST;
		trabajo->GAPEXTEND_COST = GAPEXTEND_COST;
		trabajo->SIZE_SUBTRABAJO_H = SIZE_SUBTRABAJO_H;
		trabajo->SIZE_SUBTRABAJO_V = SIZE_SUBTRABAJO_V;
		trabajo->direccion = direccion_proxima;
		trabajo->valor_maximo = maximo_valor;

	//Enviamos el mensaje con la estructurda (grupo_trabajadores, rango desencolado, tag MENSAJE_TAG)
	//Es un envío bloqueante, por lo que mientras no haya núcleos en espera, esperará
	peer_to_peer[rango_proceso].datos_disponibles=1;
	pthread_cond_signal(peer_to_peer_datos_disponibles_cond+rango_proceso);
	//if (ilib_msg_send(grupo_trabajadores, rango_proceso, MENSAJE_TAG, &trabajo, sizeof(trabajo)) != ILIB_SUCCESS){
	//	ilib_die("Error crítico, el trabajador %d no ha podido recoger el trabajo.", ilib_group_rank(ILIB_GROUP_SIBLINGS));
	//}

	pthread_mutex_unlock(peer_to_peer_mutex+rango_proceso);
	 //printf("Desloqueo controlador para enviar\n");

	if (TRAZA) {
		int rango_para_traza = rango_proceso+1; //Ya que 0,0 es el controlador
		if (rango_para_traza > 15) //0,2 es el dedicado
			rango_para_traza++;
		int posicion_tile_x_traza = rango_para_traza % DIMENSION_TARJETA;
		int posicion_tile_y_traza = rango_para_traza / DIMENSION_TARJETA;
		// [work:%d,%d] con , x, y
		printf("%08ld,W,%d,%d\n", my_clock()/CLOCKS_PER_CENTSEC, posicion_tile_x_traza, posicion_tile_y_traza);
	}
 }

/** Envía un trabajo de índice x=-1, indicando así al núcleo que el algoritmo ha terminado y puede finalizar. */
void enviar_finalizacion(int rango){
	//Creamos un trabajo cuyo índice x sea -1

	pthread_mutex_lock(peer_to_peer_mutex+rango);

		//Asignamos columna y fila del trabajo y coordenadas del trabajo
		Nodo_trabajo_completo * trabajo = &(peer_to_peer[rango].datos.trabajo);
		trabajo->x = -1;
	peer_to_peer[rango].datos_disponibles=1;
	pthread_cond_signal(peer_to_peer_datos_disponibles_cond+rango);

	pthread_mutex_unlock(peer_to_peer_mutex+rango);
 }

 /** Envía la señal de finalización a todos los núcleos del sistema. */
 void finalizar_hilos(){
	for (int i=0;i<NUMERO_TILES;i++){
		printf("Enviando finalizacion a %d\n", i);
		enviar_finalizacion(i);
	}
	// Hacer un join y cerrar los mutexes
	mylib_wait_proc_vanishes(workers_set, NUMERO_TILES);
	workers_set=NULL;
	finalizar_mutexes(NUMERO_TILES);
	//ilib_finish(); // Finaliza el entorno iLib
 }

/** Recibe un trabajo ya finalizado por parte de alguno de los núcleos, o espera a que ÃÂ©stos envíen
 * dicho trabajo (bloqueante). Además, devuelve 0 mientras el último trabajo realizado no sea [x,y]!=[tabla_trabajos_MAX_X-1, tabla_trabajos_MAX_Y-1],
 * en otro caso, 1 (es el último trabajo, es necesario salvarlo y guardarlo para la fase 2). */
int recibir_trabajo(){
	int ultimo_trabajo = 0;
	//Extraemos un puntero que apunta al trabajo
	pthread_mutex_lock(&sumidero_mutex);
	// printf("Bloqueo controlador para recibir\n");
		while(!sumidero.datos_disponibles){
			pthread_cond_wait(&sumidero_datos_disponibles_cond, &sumidero_mutex);
		}
		Nodo_trabajo_completo* trabajo = sumidero.datos.trabajo;
		sumidero.datos_disponibles=0;
		pthread_cond_signal(&sumidero_libre_cond);
	pthread_mutex_unlock(&sumidero_mutex);
	// printf("Desloqueo controlador para recibir\n");

//	printf("Recibido trabajo [%d, %d] de %d\n", trabajo->x, trabajo->y, trabajo->indice_worker);

//	imprime_estado_tabla_trabajos(tabla_trabajos, tabla_trabajos_MAX_X, tabla_trabajos_MAX_Y);
	int x = trabajo->x, y = trabajo->y;

	// Guardamos el nuevo máximo si es mayor en SW
	if (!global && (trabajo->valor_maximo >= maximo_valor)) {
		maximo_indice_x = x;
		maximo_indice_y = y;
		maximo_posicion_x = trabajo->valor_maximo_x;
		maximo_posicion_y = trabajo->valor_maximo_y;
		maximo_valor = trabajo->valor_maximo;
	}

	// En NW, si estamos ante una última fila o columna, buscamos donde está el máximo para no obligar a empezar en el último trabajo, sino
	// desde la última fila/columna que dÃÂ© mayor máximo (permite ----- finales)
	if (global && (x == tabla_trabajos_MAX_X-1)) {
		int posicion_y;
		long valor = get_max_fragmento(trabajo->fragmento_vertical, &posicion_y);
		if (valor >= maximo_valor) {
			maximo_indice_x = x;
			maximo_indice_y = y;
			maximo_valor = valor;
			maximo_posicion_x = trabajo->fragmento_horizontal->longitud - 1; //ÃÅ¡ltima posición, X - 1
			maximo_posicion_y = posicion_y;
		}
	}

	if (global && (y == tabla_trabajos_MAX_Y-1)) {
		int posicion_x;
		long valor = get_max_fragmento(trabajo->fragmento_horizontal, &posicion_x);
		if (valor >= maximo_valor) {
			maximo_indice_x = x;
			maximo_indice_y = y;
			maximo_valor = valor;
			maximo_posicion_x = posicion_x;
			maximo_posicion_y = trabajo->fragmento_vertical->longitud - 1; //ÃÅ¡ltima posición, Y - 1
		}
	}

	if (TRAZA) {
		printf("%08ld,B,%d,%d\n", my_clock()/CLOCKS_PER_CENTSEC, x, y);

		int rango_para_traza = trabajo->indice_worker+1; //Ya que 0,0 es el controlador
		if (rango_para_traza > 2) //2,0 es el dedicado
			rango_para_traza++;
		int posicion_tile_x_traza = rango_para_traza % DIMENSION_TARJETA;
		int posicion_tile_y_traza = rango_para_traza / DIMENSION_TARJETA;
		printf("%08ld,S,%d,%d\n", my_clock()/CLOCKS_PER_CENTSEC, posicion_tile_x_traza, posicion_tile_y_traza);
	}

	//Una vez recibido, asignamos los trabajos a la tabla (controlando los límites):
	//última fila de trabajo [x,y] es fila inicial de trabajo [x+1,y]
	//última columna de trabajo [x,y] es columna inicial de trabajo [x,y+1]
	//Además, si alguno es el segundo fragmento de un trabajo que está ahora disponible, lo encolamos
	if (x<tabla_trabajos_MAX_X-1) {
		set_fragmento_columna(tabla_trabajos, x+1, y, trabajo->fragmento_vertical, tabla_trabajos_MAX_X);
		if (is_fragmento_fila(tabla_trabajos, x+1, y)!=NULL){
			//Encolamos, si es el último tambiÃÂ©n
			encolar_trabajo(x+1, y);
		}
	} else {
		libera_fragmento(&(trabajo->fragmento_vertical));
	}
	// printf("Controller working IV\n");
	if (y<tabla_trabajos_MAX_Y-1) {
		// printf("Setting fragment %d\n", y);
		set_fragmento_fila(tabla_trabajos, x, y+1, trabajo->fragmento_horizontal, tabla_trabajos_MAX_X);

			// printf("Valor de y: %d\n", y);
		if (is_fragmento_columna(tabla_trabajos, x, y+1)!=NULL){
			//Encolamos, si es el último tambiÃÂ©n
			// printf("Enqueuing a work\n");
			encolar_trabajo(x, y+1);
		}
		// printf("Controller working IV-I\n");
	} else {
		// printf("Freeing fragment\n");
		libera_fragmento(&(trabajo->fragmento_horizontal));
	}
	// printf("Controller working V\n");

	// Si es el último trabajo, salvamos la fila
	if ((x == tabla_trabajos_MAX_X-1) && (y == tabla_trabajos_MAX_Y-1)) {
 		pthread_t s;
 		struct args_salvarFila * args = (struct args_salvarFila *) malloc(sizeof(struct args_salvarFila));
 		args->tablaTrabajos = tabla_trabajos;
 		args->y = y;
 		args->tabla_trabajos_MAX_X = tabla_trabajos_MAX_X;
 		args->cierra_conexion = 1;
 		args->hilo_actual = NULL;
 		pthread_create(&s, NULL, pthread_salvarFila, (void *)args); // Este hilo liberará la memoria del último ejecutado.
 		pthread_join(s, NULL);
 		// No pthread_exit() is required.
 		//salvarFila(tabla_trabajos, y, tabla_trabajos_MAX_X, 1);
		ultimo_trabajo = 1;
	}

//	imprime_estado_tabla_trabajos(tabla_trabajos, tabla_trabajos_MAX_X, tabla_trabajos_MAX_Y);
//	printf("\n");
//	printf("Trabajo finalizado, trabajos pendientes: %d\n", getQueuePrtySize(cola_trabajos));
	//Encolamos el proceso, que ahora estará ocioso
	push(cola_procesos, trabajo->indice_worker);

	//Una vez utilizada, liberamos la memoria de la estructura
	free(trabajo);

	/* Condición de salida:
	* Si ÃÂ©ste es el último trabajo ([x,y]==[tabla_trabajos_MAX_X-1, tabla_trabajos_MAX_Y-1]), devolvemos 1, en el resto de casos, 0 */
	return ultimo_trabajo;
}

/** Bucle principal de la fase 2 que calcula las dos secuencias del alineamiento de atrás hacia delante coordinando a los núcleos.
  *  Su comportamiento será el siguiente:
  *   - Mientras no hayamos llegado al primer subtrabajo (que sólo lo tendrá un procesador y será siempre el último), continuamos
  *   (en caso a de al. local, el final vendrá dado por le subtrabajo al que lleguemos a un 0)
  *   - Dentro del bucle, damos el trabajo a un núcleo, que será el único activo
  *   - Y esperamos a la resolución del mismo
  *  Finalmente, invierte la cadena y devuelve el alineamiento.
  * */
void controller_calcula_alineamiento(char ** al_horizontal, char ** al_vertical, long * align_s1_start, long * align_s1_end, long * align_s2_start, long * align_s2_end){
	int fin = 0;
	int pos_final_x = maximo_posicion_x + 1; //El último trabajo no se calcula completo: sólo lo necesario
	int pos_final_y = maximo_posicion_y + 1; //+1 porque antes era posición.ahora tamaño fila/columna
	char* fragmento_aln_vert; //fragmento de alineamiento temporal vertical
	char* fragmento_aln_hor; //fragmento de alineamiento temporal horizontal

    //Creamos las cadenas que tendrá las secuencias alineadas finales DEL REVÉS (ya que vamos desde el final hacia adelante)
    //Tamaño máximo: cadena1+cadena2
    int tam_max = strlen(secuencia_horizontal) + strlen(secuencia_vertical);
	*al_vertical = (char *) malloc_safe_agr248(sizeof(char) * tam_max);
	*al_horizontal = (char *) malloc_safe_agr248(sizeof(char) * tam_max);
	int inicializadas = 0; //Para hacer un strcopy en lugar de strcat la primera vez

	// Ajuste del comienzo del alineamiento global:
	// Si no comenzamos en las últimas posiciones, el alineamiento comienza con ------ en un lado
	if (global) {
		int it=0;
		int pos_x_global = maximo_indice_x * SIZE_SUBTRABAJO_H + maximo_posicion_x;
		int pos_x_i = strlen(secuencia_horizontal)-1;
		int pos_y_global = maximo_indice_y * SIZE_SUBTRABAJO_V + maximo_posicion_y;
		int pos_y_j = strlen(secuencia_vertical)-1;

		// EN NINGÃÅ¡N CASO PUEDEN DARSE AMBAS CONDICIONES
		//---- en la segunda secuencia (ajustamos horizontal)
		while (pos_x_i >= pos_x_global) {
			(*al_horizontal)[it] = secuencia_horizontal[pos_x_i];
			(*al_vertical)[it] = '-';
			pos_x_i--;
			it++;
			inicializadas = 1;
		}
		//---- en la primera secuencia (ajustamos vertical)
		while (pos_y_j >= pos_y_global) {
			(*al_horizontal)[it] = '-';
			(*al_vertical)[it] = secuencia_vertical[pos_y_j];
			pos_y_j--;
			it++;
			inicializadas = 1;
		}
		if (inicializadas) {
			(*al_horizontal)[it] = FINCAD;
			(*al_vertical)[it] = FINCAD;
		}

		*align_s1_end = pos_x_global;
		*align_s2_end = pos_y_global;
	} else {
		// Es desde el maximo_posicion, contando el desplazamiento de los trabajos
		*align_s1_end = maximo_indice_x * SIZE_SUBTRABAJO_H + maximo_posicion_x;
		*align_s2_end = maximo_indice_y * SIZE_SUBTRABAJO_V + maximo_posicion_y;
	}

	direccion_proxima = DIRECCION_DIAGONAL; //Empezamos en diagonal siempre

	while (!fin){
		//Mientras haya trabajo disponible, seguimos enviando trabajos
//		printf("Cola procesos -> ");
//		print_queue(cola_procesos);
//		printf("Cola trabajos -> ");
//		print_queue(cola_trabajos);

		//TODO: EN LA ACTUAL IMPLEMENTACIÓN, SÓLO SE EJECUTARÁ UNA VEZ PUES SÓLO HAY UN PROCESO QUE CALCULA
		if (!is_empty(cola_procesos) && !isEmptyQueuePrty(cola_trabajos)){
			enviar_trabajo(pos_final_x, pos_final_y);
		}
		fin = recibir_alineamiento(&fragmento_aln_hor, &fragmento_aln_vert, &pos_final_x, &pos_final_y, align_s1_start, align_s2_start);

		if (inicializadas){
			strcat(*al_horizontal, fragmento_aln_hor);
			strcat(*al_vertical, fragmento_aln_vert);
		} else {
			strcpy(*al_horizontal, fragmento_aln_hor);
			strcpy(*al_vertical, fragmento_aln_vert);
		}

		// Liberamos las cadenas temporales de memoria compartida
		free(fragmento_aln_hor);
		free(fragmento_aln_vert);

		inicializadas = 1; //a partir de ahora, se concatenarán
	}

	//Invertimos la cadena del alineamiento y la devolvemos
    reverse(*al_horizontal);
    reverse(*al_vertical);
}

/** Recibe un fragmento de alineamiento ya calculado por parte de alguno de los núcleos, o espera a que ÃÂ©stos
 * envíen dicho trabajo (bloqueante). Además, control de fin:
 * Alineamiento global: devuelve 0 mientras [x,y]!=[0, 0], en otro caso, 1 (es el último trabajo).
 * Alineamiento local: devuelve 0 mientras limite_x != -1, en otro caso, 1 (se ha alcanzado el 0 de inicio alineamiento). */
int recibir_alineamiento (char ** frag_aln_horizontal, char ** frag_aln_vertical, int * pos_x_final, int * pos_y_final, long * fin_alineamiento_x, long * fin_alineamiento_y){
	//Extraemos un puntero que apunta al alineamiento
	pthread_mutex_lock(&sumidero_mutex);
		while(!sumidero.datos_disponibles){
			pthread_cond_wait(&sumidero_datos_disponibles_cond, &sumidero_mutex);
		}
		Nodo_alineamiento* alineamiento = sumidero.datos.alineamiento;
		sumidero.datos_disponibles=0;
		pthread_cond_signal(&sumidero_libre_cond);
	pthread_mutex_unlock(&sumidero_mutex);

	int x = alineamiento->x, y = alineamiento->y;

	*frag_aln_horizontal = alineamiento->fragmento_horizontal;
	*frag_aln_vertical = alineamiento->fragmento_vertical;

	*pos_x_final = alineamiento->limite_x;
	*pos_y_final = alineamiento->limite_y;

	direccion_proxima = alineamiento->direccion;

	//Metemos al proceso a la cola
	push(cola_procesos, alineamiento->indice_worker);

	// Global: Si estamos en el primer subtrabajo, hemos terminado, liberamos y retornamos
	// Local: Si se nos devuelve *pos_x_final == -1, hemos alcanzado el 0, inicio de alineamiento
	if ((x == 0 && y == 0) || (*pos_x_final == -1)) {
		*fin_alineamiento_x = x * SIZE_SUBTRABAJO_H + alineamiento->valor_minimo_x; // donde alineamiento->valor_minimo_x != -1
		*fin_alineamiento_y = y * SIZE_SUBTRABAJO_V + alineamiento->valor_minimo_y;

		free(alineamiento);
		return 1;
	}

	//Una vez recibido, invocamos al siguiente fragmento en función de los límites resultantes
	//si [lim_x,limy] == 0 -> izquierda superior (completo -> [x=0,y=0])
	//si lim_x == 0 -> izquierda hasta límite lim_y (sólo algunas filas)
	//si lim_y == 0 -> arriba hasta límite lim_x (sólo algunas columnas)
	if (*pos_x_final == 0 && *pos_y_final == 0){
		if (x>0 && y>0) {
			encolar_trabajo(x-1, y-1);
		} else if (x>0) { //y==0
			encolar_trabajo(x-1, y); //Rellenaremos huecos en todo el subtrabajo (límite vertical y=0)
			*pos_y_final = 1; // Para que no sea [0,0] y crea que es completo, pasamos que hay una única fila
		} else { //y>0 & x==0
			encolar_trabajo(x, y-1); //Rellenaremos huecos en todo el subtrabajo (límite horizontal x=0)
			*pos_x_final = 1; // Para que no sea [0,0] y crea que es completo, pasamos que hay una única columna
		}
	} else if (*pos_x_final == 0) {
		encolar_trabajo(x-1, y);
	} else if (*pos_y_final == 0) {
		encolar_trabajo(x, y-1);
	} else {
		printf("ERROR: ESTO NUNCA DEBERÁA DE OCURRIR PERO MIENTRAS DEPURO PUEDE QUE SÁ");
		exit(0);
	}

	//Una vez utilizada, liberamos la memoria de la estructura
	free(alineamiento);

	return 0; //ya que si era fin, ocurrió antes (1)
}

/** Escribe los datos de la tabla de trabajos en el fichero seleccionado. */
 void escribe_tabla_trabajos(char * fichero){
	escribe_tabla_fichero(fichero, tabla_trabajos, tabla_trabajos_MAX_X, tabla_trabajos_MAX_Y);
 }


/** Función auxiliar para el cálculo del mínimo entre dos números. */
int min(int a, int b){
	return (a < b) ? a : b;
}

/** Función auxiliar para el cálculo del máximo entre dos números. */
int max(int a, int b){
	return (a > b) ? a : b;
}

/** Función auxiliar que modifica una cadena y la convierte en su inversa. */
void reverse(char *t) {
	int i,j;
	for(i = 0, j = strlen(t)-1; i < j; i++, j--) {
		char c; c = t[i]; t[i] = t[j]; t[j] = c;
	}
}

/** Función auxiliar que calcula ceil(total/k) en enteros para sacar el número de trabajosc a partir
 *  del tamaño total y el subtamaño de cada uno. */
int total_subtrabajos(int total, int k) {
	int div = total / k;
	if (total % k > 0)
		div++;
	return div;
}

/** Escribe los valores temporales necesarios para la segunda fase. */
void salva_maximos_fichero(char* file_name, long maximo_valor, int maximo_indice_x, int maximo_indice_y, int maximo_posicion_x, int maximo_posicion_y, clock_t tiempo_inicio, clock_t tiempo_intermedio) {
	FILE *fichero;
	if((fichero = fopen(file_name, "w"))==NULL) {
    	printf("No se pueden salvar los valores temporales en el fichero: %s\n", file_name);
    	exit(1);
  	}
  	fprintf(fichero, "%ld\n%d\n%d\n%d\n%d\n%ld\n%ld", maximo_valor, maximo_indice_x, maximo_indice_y, maximo_posicion_x, maximo_posicion_y, tiempo_inicio, tiempo_intermedio);
	fclose(fichero); //Cerramos el fichero
}


#pragma GCC diagnostic ignored "-Wunused-result"
/** Lee los valores temporales necesarios para la segunda fase. */
void recupera_maximos_fichero(char * file_name, long *maximo_valor, int *maximo_indice_x, int *maximo_indice_y, int *maximo_posicion_x, int *maximo_posicion_y, clock_t *tiempo_inicio, clock_t *tiempo_intermedio) {
	FILE *fichero;
	if((fichero = fopen(file_name, "r"))==NULL) {
    	printf("No se puede abrir el fichero de valores temporales: %s\n", file_name);
    	exit(1);
  	}
  	fscanf(fichero, "%ld\n%d\n%d\n%d\n%d\n%ld\n%ld", maximo_valor, maximo_indice_x, maximo_indice_y, maximo_posicion_x, maximo_posicion_y, tiempo_inicio, tiempo_intermedio);
	fclose(fichero); //Cerramos el fichero
}
#pragma GCC diagnostic warning "-Wunused-result"

/* Manipulacion de exclusiones mutuas para comunicacion mediante memoria compartida
 */
 void inicializar_mutexes(int num_tiles){
 	// Inicializamos el "canal" sumidero
 	pthread_mutex_init(&sumidero_mutex, NULL);
 	sumidero.datos_disponibles = 0;
 	pthread_cond_init(&sumidero_datos_disponibles_cond, NULL);
 	pthread_cond_init(&sumidero_libre_cond, NULL);
 	// Inicializamos los canales punto a punto
 	peer_to_peer_mutex = (pthread_mutex_t *)malloc(sizeof(pthread_mutex_t)*num_tiles);
 	peer_to_peer = (peer_to_peer_struct *) malloc(sizeof(peer_to_peer_struct)*num_tiles);
 	peer_to_peer_datos_disponibles_cond = (pthread_cond_t *) malloc(sizeof(pthread_cond_t)*num_tiles);
 	for(int i=0; i<num_tiles; i++){
 		pthread_mutex_init(peer_to_peer_mutex+i, NULL);
 		peer_to_peer[i].datos_disponibles = 0;
 		pthread_cond_init(peer_to_peer_datos_disponibles_cond+i, NULL);
 	}
 }
 void finalizar_mutexes(int num_tiles){
 	pthread_mutex_destroy(&sumidero_mutex);
 	pthread_cond_destroy(&sumidero_datos_disponibles_cond);
 	pthread_cond_destroy(&sumidero_libre_cond);
 	for(int i=0; i<num_tiles; i++){
 		pthread_mutex_destroy(peer_to_peer_mutex+i);
 		pthread_cond_destroy(peer_to_peer_datos_disponibles_cond+i);
 	}
 	free(peer_to_peer_mutex);
 	peer_to_peer_mutex=NULL;
 	free(peer_to_peer);
 	peer_to_peer=NULL;
 	free(peer_to_peer_datos_disponibles_cond);
 	peer_to_peer_datos_disponibles_cond=NULL;
 }



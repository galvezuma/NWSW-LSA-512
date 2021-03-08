/*
 ============================================================================
 Name        : NWSW-LSA-512.c
 Author      : SGR
 Version     :
 Copyright   : Código abierto
 Description : Hello World in C, Ansi-style
 ============================================================================
 */

/*
 ============================================================================
 Name        : XPhi_NWSW_Final.c
 Author      : SGR
 ============================================================================
 */

/* Algoritmo Needleman-Wunsch paralelo, que utiliza un controlador sobre el resto de nÃºcleos. */

#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include "NWSW_LSA_512_controller.h"
#include "priority_queue.h"

// DEBUG CON CLUSTAL2-P2_TEST: COMUNICACIÓN POR DEVICES
//#include <fcntl.h>

void calcula_informacion_extendida(char * alineamiento_query, char * alineamiento_subject, long * align_num_ids, long * align_num_gaps, long * align_long);
void leer_fichero_fasta(char* file, char** cadena, char* cabecera);
void escribir_fichero_alineamiento(char* file, char* cadena1, char* cabecera1, char* cadena2, char* cabecera2);
void escribir_alineamiento(FILE** fichero, char* cadena, char* cabecera);
void escribir_fichero_puntuacion(char* file, long score);
void escribir_fichero_info(char* file, char* cabecera1, char* cabecera2, long INSERT_COST, long DELETE_COST, long MATCHREPLACE_COST, long GAPEXTEND_COST,
	int SIZE_SUBTRABAJO_H, int SIZE_SUBTRABAJO_V, int global, estado_algoritmo stage, Nodo_resultado * resultado, int informacion_extendida,
	long long_s1, long long_s2, long align_num_ids, long align_num_gaps, long aling_long);
void escribir_fichero_info_compacto(char* file, char* cabecera1, char* cabecera2, long INSERT_COST, long DELETE_COST, long MATCHREPLACE_COST, long GAPEXTEND_COST,
	int SIZE_SUBTRABAJO_H, int SIZE_SUBTRABAJO_V, int global, estado_algoritmo stage, Nodo_resultado * resultado, int informacion_extendida,
	long long_s1, long long_s2, long align_num_ids, long align_num_gaps, long aling_long);


/** A partir de dos ficheros con UNA ÚNICA SECUENCIA FASTA CADA UNO, extrae las secuencias y las lanza al NW paralelo.
 *  Argumentos:
 *  argv[0]: nombre del programa
 *  argv[-] (opción): --1pass (por defecto) o --2pass: Indica la fase del algoritmo a realizar y devolverá resultados diferentes en cada una (1: score, 2: alineamiento)
 *  argv[-] (opción): --global (por defecto) o --local: Lanza el algoritmo de alineamiento global Needleman-Wunsch (global) o el de alineamiento Smith-Waterman (local)
 *  argv[-] (opción): --extendedinfo (por defecto) o --basicinfo: Devuelve en el fichero de info información básica o más información (algo más costoso en la vuelta atrás)
 *  argv[1]: nombre del fichero FASTA con la secuencia query (original a comparar)
 *  argv[2]: nombre del fichero FASTA con la secuencia subject (con la que comparar)
 *  argv[3]: nombre del fichero de salida de texto con la secuencias FASTA alineadas
 *  argv[4]: nombre del fichero de salida de texto con la información del alineamiento
 *  -- OPCIONALES (en caso de no indicar, valores por defecto): --
 *  argv[5]: INSERT_COST: Coste de la operación de apertura de una inserción (hueco en la secuencia original [horizontal])
 *  argv[6]: DELETE_COST: Coste de la operación de apertura de una delección (hueco en la secuencia con la que comparar [vertical])
 *  argv[7]: MATCHREPLACE_COST: Puntuación por defecto en caso de que ambos nucleótidos coincidan, pero no estÃ©n en la matriz
 *  argv[8]: GAPEXTEND_COST: Coste de la extensión de un hueco una vez abierta
 *  argv[9]: fichero de matriz de puntuaciones a utilizar para el cálculo del alineamiento
 *  argv[10]: SIZE_SUBTRABAJO_H: Tamaño horizontal de un trabajo (K horizontal)
 *  rgv[11]: SIZE_SUBTRABAJO_V: Tamaño vertical de un trabajo (K vertical)
 *  argv[12]: NUMERO_TILES: Número de tiles a utilizar para la primera fase (<= remaining - RESERVED_TILES). Si no indicado, 59
 * */
int main(int argc, char *argv[]) {
	/* Secuencias y sus cabeceras de las secuencias para escribir al final */
	char* cadena_query;
    char cabecera_cadena_query[MAX_SIZEHEADER];
	char* cadena_subject;
	char cabecera_cadena_subject[MAX_SIZEHEADER];

	/* Resultado: alineamientos o puntuación final de la primera fase */
    Nodo_resultado resultado;

    /* Parámetros del algoritmo */
    long INSERT_COST = INSERT_COST_DEFAULT, DELETE_COST = DELETE_COST_DEFAULT;
    long MATCHREPLACE_COST = MATCHREPLACE_COST_DEFAULT, GAPEXTEND_COST = GAPEXTEND_COST_DEFAULT;
    int SIZE_SUBTRABAJO_H = SIZE_SUBTRABAJO_H_DEFAULT, SIZE_SUBTRABAJO_V = SIZE_SUBTRABAJO_V_DEFAULT;
    int NUMERO_TILES = 6; // NUMERO_TILES_DEFAULT;
    char* matriz_puntuaciones = NULL;
    estado_algoritmo stage = FORWARD; //Por defecto, primera pasada
    int global = 1; //Por defecto, el algoritmo NW global
    int informacion_extendida = 1; //Por defecto, se da información extendida

    /* EstadÃ­sticas de configuración extendida */
    long long_s1=0, long_s2=0;
    long align_num_ids=0, align_num_gaps=0, align_long=0;

    int ops = 0; //Si hay alguna opción al principio, desplazamos todo
    /* Antes de nada, extraemos si tiene opciones */

  	if ((argc>=2) && !strcmp(argv[1+ops], "--1pass")) {
   		stage = FORWARD; //Primera pasada
   		ops++;
   	} else if ((argc>=2) && !strcmp(argv[1+ops], "--2pass")) {
       	stage = BACKWARD; //Segunda pasada
       	ops++;
   	}

   	if ((argc>=2) && !strcmp(argv[1+ops], "--global")) {
   		global = 1; //Needleman-Wunsch
   		ops++;
   	} else if ((argc>=2) && !strcmp(argv[1+ops], "--local")) {
       	global = 0; //Smith-Waterman
       	ops++;
   	}

   	if ((argc>=2) && !strcmp(argv[1+ops], "--extendedinfo")) {
   		informacion_extendida = 1; //Información más extensa
   		ops++;
   	} else if ((argc>=2) && !strcmp(argv[1+ops], "--basicinfo")) {
       	informacion_extendida = 0; //Información básica
       	ops++;
   	}

    /* Primera fase: filtrado de parámetros */
	if (argc < ops+4) {
		printf("ERROR: Usage: MultiCore64-NW [--1pass || --2pass] [--global || --local] [--extendedinfo || --basicinfo] fasta_query fasta_subject output info_file [INSERT_COST] [DELETE_COST] [MATCHREPLACE_COST] [GAPEXTEND_COST] [matrix] [K HORIZONTAL SIZE] [K VERTICAL SIZE] [NUM TILES]\n");
		exit(-1);
	}
	// Parámetros opcionales:
	if (argc >= ops+13)
		NUMERO_TILES = atoi(argv[ops+12]);
	if (argc >= ops+12)
		SIZE_SUBTRABAJO_V = atoi(argv[ops+11]);
	if (argc >= ops+11)
		SIZE_SUBTRABAJO_H = atoi(argv[ops+10]);
	if (argc >= ops+10)
		matriz_puntuaciones = argv[ops+9];
	if (argc >= ops+9)
		GAPEXTEND_COST = atoi(argv[ops+8]);
	if (argc >= ops+8)
		MATCHREPLACE_COST = atoi(argv[ops+7]);
	if (argc >= ops+7)
		DELETE_COST = atoi(argv[ops+6]);
	if (argc >= ops+6)
		INSERT_COST = atoi(argv[ops+5]);

	// El prefijo temporal será el mismo nombre del fichero de salida, esto es, argv[ops+3]
	int size_prefijo = strlen(argv[ops+3]);
	prefijo_fichero_temporal = (char *)malloc_safe_agr248(sizeof(char) * size_prefijo);
	strcpy(prefijo_fichero_temporal, argv[ops+3]);

	leer_fichero_fasta(argv[ops+1], &cadena_query, cabecera_cadena_query);
	leer_fichero_fasta(argv[ops+2], &cadena_subject, cabecera_cadena_subject);


	//Inicializamos la librería
	mylib_init();

    //Lanzamos el algoritmo, que escribirá en los punteros las cadenas resultantes de alineamiento
	    MC64NWSW_Controlador(prefijo_fichero_temporal, global, stage, cadena_query, cadena_subject, &resultado,
    		INSERT_COST, DELETE_COST, MATCHREPLACE_COST, GAPEXTEND_COST, matriz_puntuaciones,
    		SIZE_SUBTRABAJO_H, SIZE_SUBTRABAJO_V, NUMERO_TILES);



    /* Escribimos los resultados: si es segunda fase, alineamientos, si es primera fase, puntuación total del alineamiento. */
	if (stage == FORWARD) {
    	escribir_fichero_puntuacion(argv[ops+3], resultado.datos.puntuacion_total);
    }
    if (stage == BACKWARD) {
   		escribir_fichero_alineamiento(argv[ops+3], resultado.datos.alineamientos.alineamiento_query, cabecera_cadena_query,
   			resultado.datos.alineamientos.alineamiento_subject, cabecera_cadena_subject);

   		if (informacion_extendida) {
   			long_s1 = strlen(cadena_query);
    		long_s2 = strlen(cadena_subject);
   			calcula_informacion_extendida(resultado.datos.alineamientos.alineamiento_query, resultado.datos.alineamientos.alineamiento_subject,
   				&align_num_ids, &align_num_gaps, &align_long);
   		}
    }

    /* Escribimos el fichero de informacion: si es para ClustalW el normal (necesita sacar luego el valor del maximum_score de ahÃ­). */
#ifdef CLUSTALW
   escribir_fichero_info(argv[ops+4], cabecera_cadena_query, cabecera_cadena_subject, INSERT_COST, DELETE_COST, MATCHREPLACE_COST, GAPEXTEND_COST,
		SIZE_SUBTRABAJO_H, SIZE_SUBTRABAJO_V, global, stage, &resultado, informacion_extendida, long_s1, long_s2, align_num_ids, align_num_gaps, align_long);
#else
   escribir_fichero_info_compacto(argv[ops+4], cabecera_cadena_query, cabecera_cadena_subject, INSERT_COST, DELETE_COST, MATCHREPLACE_COST, GAPEXTEND_COST,
		SIZE_SUBTRABAJO_H, SIZE_SUBTRABAJO_V, global, stage, &resultado, informacion_extendida, long_s1, long_s2, align_num_ids, align_num_gaps, align_long);
#endif

/*    // DEPURACIÃ“N
    printf("H_final [length=%d]: %s\n", strlen(alineamiento_query), alineamiento_query);
    printf("V_final [length=%d]: %s\n", strlen(alineamiento_subject), alineamiento_subject); */

	/* MENSAJES DE FINALIZACIÃ“N PARA DETECCIÃ“N EN SCRIPT */
	if (stage == FORWARD) { //DONE1 score infofile
		printf("DONE1 %ld %s\n", resultado.datos.puntuacion_total, argv[ops+4]);
	}
	if (stage == BACKWARD) { //DONE2 alignmentfile infofile
		printf("DONE2 %s %s\n", argv[ops+3], argv[ops+4]);
	}

	//Finalizamos el pthread
	mylib_finish();
}

/* Recorre el alineamiento y recoge parte de la información nueva adicional sobre el mismo (nÃºmero de IDS y nÃºmero de gaps) */
void calcula_informacion_extendida(char * alineamiento_query, char * alineamiento_subject, long * align_num_ids, long * align_num_gaps, long * align_long) {
	long align_num_ids_temp = 0;
	long align_num_gaps_temp = 0;
	unsigned int i;

	*align_long = strlen(alineamiento_query); // debe ser igual a alineamiento_subject
	for (i=0; i < *align_long; i++) {
		if (alineamiento_query[i] == alineamiento_subject[i])
			align_num_ids_temp++; // Ya que no puede haber -/-
		else if ((alineamiento_query[i] == '-') || (alineamiento_subject[i]  == '-'))
			align_num_gaps_temp++;
	}
	*align_num_ids = align_num_ids_temp;
	*align_num_gaps = align_num_gaps_temp;
}

/** Lee una secuencia de un fichero FASTA, extrayendo la secuencia y su cabecera. */
void leer_fichero_fasta(char* file, char** cadena_ptr, char* cabecera){
	FILE *fichero;
	if((fichero = fopen(file, "r"))==NULL) {
    	printf("No se puede abrir el fichero de secuencias FASTA: %s\n", file);
    	exit(1);
  	}

  	char temp[MAX_SIZELINE]; //Suponemos que cada lÃ­nea no será de más de MAX_SIZELINE caracteres
  	//Guardamos la primera lÃ­nea (cabecera FASTA)

	#pragma GCC diagnostic ignored "-Wunused-result"
  	fgets(cabecera, MAX_SIZEHEADER, fichero);
	#pragma GCC diagnostic warning "-Wunused-result"
  	if (cabecera[0]!='>') { //Si no comienza por >, no es cabecera FASTA, secuencia sin cabecera
  		cabecera[0]=FINCAD;
  	}

  	fseek(fichero, 0, SEEK_END); //Nos vamos al final para calcular el tamaÃ±o a reservar de la cadena
  	int size_primera_linea = strlen(cabecera);
  	long size = ftell(fichero) - size_primera_linea;
  	int posicion_cadena = 0;

  	//Creamos el espacio para la cadena (memoria compartida)
  	*cadena_ptr = (char *) malloc_shared (sizeof(char) * (size + 1));

  	fseek(fichero, size_primera_linea, SEEK_SET); //Continuamos leyendo
  	while (fgets(temp, MAX_SIZELINE, fichero)!=NULL) {
  		int ultimo = strlen(temp)-1;
  		if (temp[ultimo-1]=='\n' || temp[ultimo-1]=='\r'){ //Eliminamos salto de lÃ­nea (UNIX: \r\n)
  			ultimo = ultimo-1;
		} else if (temp[ultimo]=='\n' || temp[ultimo]=='\r'){ //Eliminamos salto de lÃ­nea (Windows: \n)
			// ultimo = ultimo;
  		}

  		int pos_temp = 0;
  		char caracter;
  		while (pos_temp < ultimo && (caracter = temp[pos_temp])!=FINCAD) {
  			(*cadena_ptr)[posicion_cadena] = toupper(caracter); // La letra SIEMPRE EN MAYÃšSCULAS para evitar conflictos minus-mayus al alinear
  			posicion_cadena++;
  			pos_temp++;
  		}
  	}
  	(*cadena_ptr)[posicion_cadena] = FINCAD;

	fclose(fichero); //Cerramos el fichero
}

/** Escribe un alineamiento por parejas FASTA en un fichero, a partir de las secuencias resultantes del
 *  alineamiento y cabeceras de ambas secuencias. */
void escribir_fichero_alineamiento(char* file, char* cadena1, char* cabecera1, char* cadena2, char* cabecera2){
	FILE *fichero;
	if((fichero = fopen(file, "w"))==NULL) {
    	printf("No se puede escribir en el fichero de alineamiento FASTA: %s\n", file);
    	exit(1);
  	}

  	escribir_alineamiento(&fichero, cadena1, cabecera1);
  	escribir_alineamiento(&fichero, cadena2, cabecera2);

	fclose(fichero); //Cerramos el fichero
}

/** Escribe un un alineamiento a partir de su secuencia resultante y cabecera en un fichero ya abierto. */
void escribir_alineamiento(FILE** fichero, char* cadena, char* cabecera){
	fprintf(*fichero, "%s", cabecera);
	int longitud = strlen(cadena);

	//Escribimos tandas de SIZELINE_ALIGNMENT como máximo en cada lÃ­nea
	int i = 0;
	while (i < longitud) {
		int j = 0;
		while (j<SIZELINE_ALIGNMENT && i < longitud) {
			fprintf(*fichero, "%c", cadena[i]);
			j++;
			i++;
		}
		fprintf(*fichero, "\n");
	}
	fprintf(*fichero, "\n"); //salto de lÃ­nea final
}

/** Escribe el valor long de la puntuación total pasada como parámetro como resultado de un alineamiento. */
void escribir_fichero_puntuacion(char* file, long score){
	FILE *fichero;
	if((fichero = fopen(file, "w"))==NULL) {
    	printf("No se puede escribir en el fichero con el valor de la puntuación total del alineamiento: %s\n", file);
    	exit(1);
  	}

  	fprintf(fichero, "%ld", score);

	fclose(fichero); //Cerramos el fichero
}

/** Escribe la información del alineamiento solicitado en un fichero. */
void escribir_fichero_info(char* file, char* cabecera1, char* cabecera2, long INSERT_COST, long DELETE_COST, long MATCHREPLACE_COST, long GAPEXTEND_COST,
	int SIZE_SUBTRABAJO_H, int SIZE_SUBTRABAJO_V, int global, estado_algoritmo stage, Nodo_resultado * resultado, int informacion_extendida,
	long long_s1, long long_s2, long align_num_ids, long align_num_gaps, long aling_long){

	FILE *fichero;
	if((fichero = fopen(file, "w"))==NULL) {
    	printf("No se puede escribir en el fichero de de información: %s\n", file);
    	exit(1);
  	}

	fprintf(fichero, "### ALIGNMENT RESULTS FROM MC64-NWSW ###\n");
	if (global) {
		fprintf(fichero, "GLOBAL ALIGNMENT\n\n");
	} else {
		fprintf(fichero, "LOCAL ALIGNMENT\n\n");
	}
    fprintf(fichero, "Sequence 1 header: %s\n", cabecera1);
    fprintf(fichero, "Sequence 2 header: %s\n", cabecera2);
    fprintf(fichero, "Insert operation cost: %ld\n", INSERT_COST);
	fprintf(fichero, "Delete operation cost: %ld\n", DELETE_COST);
	fprintf(fichero, "Match/Replace operation cost: %ld\n", MATCHREPLACE_COST);
	fprintf(fichero, "Gap extend operation cost: %ld\n", GAPEXTEND_COST);

	fprintf(fichero, "Horizontal K-size: %d\n", SIZE_SUBTRABAJO_H);
	fprintf(fichero, "Vertical K-size: %d\n", SIZE_SUBTRABAJO_V);

	if (stage == FORWARD) {
		fprintf(fichero, "\n  SCORING ONLY CALCULATION\n");
	} else if (stage == BACKWARD) {
		fprintf(fichero, "\n  FULL ALIGMENT CALCULATION\n");
	}
	fprintf(fichero, "Maximum score: %ld\n\n", resultado->datos.puntuacion_total);
	fprintf(fichero, "Initial time: %.3f seconds\n", (float)resultado->inicio/MILLIS_PER_SEC);
	fprintf(fichero, "First pass time: %.3f seconds\n", (float)(resultado->intermedio - resultado->inicio)/MILLIS_PER_SEC);
	if (stage == BACKWARD) {
		fprintf(fichero, "Total time: %.3f seconds\n", (float)(resultado->final - resultado->inicio)/MILLIS_PER_SEC);
	}

	if (stage == BACKWARD && informacion_extendida) {
		fprintf(fichero, "\n  EXTENDED STADISTICS\n");
		fprintf(fichero, "Alignment length: %ld\n", aling_long);
		fprintf(fichero, "    Identities: %ld (%.2f%%)\n", align_num_ids, 100 * (float)align_num_ids/aling_long);
		fprintf(fichero, "    Gaps: %ld (%.2f%%)\n", align_num_gaps, 100 * (float)align_num_gaps/aling_long);
		fprintf(fichero, "Sequence 1:\n");
		fprintf(fichero, "    Length: %ld\n", long_s1);
		fprintf(fichero, "    Alignment start position: %ld\n", resultado->datos.alineamientos.align_s1_start);
		fprintf(fichero, "    Alignment end position: %ld\n", resultado->datos.alineamientos.align_s1_end);
		fprintf(fichero, "Sequence 2:\n");
		fprintf(fichero, "    Length: %ld\n", long_s2);
		fprintf(fichero, "    Alignment start position: %ld\n", resultado->datos.alineamientos.align_s2_start);
		fprintf(fichero, "    Alignment end position: %ld\n", resultado->datos.alineamientos.align_s2_end);
	}

	fclose(fichero); //Cerramos el fichero
}

/** Escribe la información del alineamiento solicitado en un fichero en su versión compacta para separar. */
void escribir_fichero_info_compacto(char* file, char* cabecera1, char* cabecera2, long INSERT_COST, long DELETE_COST, long MATCHREPLACE_COST, long GAPEXTEND_COST,
	int SIZE_SUBTRABAJO_H, int SIZE_SUBTRABAJO_V, int global, estado_algoritmo stage, Nodo_resultado * resultado, int informacion_extendida,
	long long_s1, long long_s2, long align_num_ids, long align_num_gaps, long aling_long){

	char separador = '|'; // constante a la función

	FILE *fichero;
	if((fichero = fopen(file, "w"))==NULL) {
    	printf("No se puede escribir en el fichero de de información: %s\n", file);
    	exit(1);
  	}

	fprintf(fichero, "### ALIGNMENT RESULTS FROM MC64-NWSW: ");
	if (global) {
		fprintf(fichero, "GLOBAL ALIGNMENT ###\n");
	} else {
		fprintf(fichero, "LOCAL ALIGNMENT ###\n");
	}
	fprintf(fichero, "### COMPACTED INFO FORMAT ###\n");
	fprintf(fichero, "# seq1_header\n");
	fprintf(fichero, "# seq2_header\n");
	fprintf(fichero, "# INSERT_COST%cDELETE_COST%cMATCHREPLACE_COST%cGAPEXTEND_COST%cHORIZONTAL_K_SIZE%cVERTICAL_K_SIZE\n",
		separador, separador, separador, separador, separador);
	fprintf(fichero, "# stage%cmaximum_score%cinitial_time%c1pass_time%ctotal_time\n",
		separador, separador, separador, separador);
	fprintf(fichero, "# EXTENDED: alignment_length%cnum_identities%cnum_gaps%cs1_length%cs1_align_start%cs1_align_end%cs2_length%cs2_align_start%cs2_align_end\n",
		separador, separador, separador, separador, separador, separador, separador, separador);
	fprintf(fichero, "%s%s", cabecera1, cabecera2);
	fprintf(fichero, "%ld%c%ld%c%ld%c%ld%c%d%c%d\n",
		INSERT_COST, separador, DELETE_COST, separador, MATCHREPLACE_COST, separador, GAPEXTEND_COST, separador, SIZE_SUBTRABAJO_H, separador, SIZE_SUBTRABAJO_V);
	fprintf(fichero, "%s%c%ld%c%.2f%c%.3f",
		((stage == FORWARD) ? "1pass" : "2pass"), separador, resultado->datos.puntuacion_total, separador, (float)resultado->inicio, separador, (float)(resultado->final - resultado->inicio)/MILLIS_PER_SEC);
	if (stage == BACKWARD)
		fprintf(fichero, "%c%.2f\n", separador, (float)(resultado->final - resultado->inicio)/MILLIS_PER_SEC);
	else
		fprintf(fichero, "\n");
	if (informacion_extendida)
		fprintf(fichero, "EXTENDED: %ld%c%ld%c%ld%c%ld%c%ld%c%ld%c%ld%c%ld%c%ld\n",
			aling_long, separador, align_num_ids, separador, align_num_gaps, separador, long_s1, separador, resultado->datos.alineamientos.align_s1_start, separador, resultado->datos.alineamientos.align_s1_end,
			separador, long_s2, separador, resultado->datos.alineamientos.align_s2_start, separador, resultado->datos.alineamientos.align_s2_end);
	fclose(fichero); //Cerramos el fichero
}



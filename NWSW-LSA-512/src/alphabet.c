/*
 * alphabet.c
 *
 *  Created on: Mar 8, 2021
 *      Author: galvez
 */

/* Gestión del alfabeto de las secuencias a alinear y la matriz de sustitución. */
/* EL SISTEMA NO ADMITE MATRICES CON DIFERENCIA DE MAYÚSCULAS Y MINÚSCULAS, COSA QUE NO HA OCURRIDO HASTA AHORA, PERO NUNCA SE SABE. */

#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include "MyLib.h"
#include "alphabet.h"

/* De uso privado */
unsigned short * codifica_alfabeto(const char * alfabeto, int longAlfabeto);

/** Asigna al alfabeto y matriz los valores por defecto como constantes. */
void init_matrix_default(struct Matriz_puntuaciones** ptr_tabla_puntuacion, int MATCHREPLACE_COST) {
	struct Matriz_puntuaciones* tabla_puntuacion = (struct Matriz_puntuaciones*) malloc_shared(sizeof(struct Matriz_puntuaciones));

	// Copiamos los valores a dos nuevos arrays en memoria compartida para que accedan los workers
	unsigned int alfabeto_length = strlen(alfabeto_default);
	tabla_puntuacion->alfabeto_horizontal = (char*) malloc_shared(sizeof(char) * alfabeto_length);
	memcpy(tabla_puntuacion->alfabeto_horizontal, alfabeto_default, alfabeto_length);
	tabla_puntuacion->alfabeto_vertical = (char*) malloc_shared(sizeof(char) * alfabeto_length);
	memcpy(tabla_puntuacion->alfabeto_vertical, alfabeto_default, alfabeto_length);

	// Copiamos los valores a una nueva matriz en memoria compartida para que accedan los workers
	// Nota IMPORTANTE: +1 para poder asignar valores MATCHREPLACE a los nucleótidos que no aparecen en la matriz
	// haciendo uso de una nueva fila y columna, de valores MATCHREPLACE_COST
	unsigned int anchura = alfabeto_length+1;
	unsigned int altura = alfabeto_length+1;
	tabla_puntuacion->matriz = (int**) malloc_shared(sizeof(int*) * anchura);
    for (unsigned int i = 0; i < anchura; i++) {
       (tabla_puntuacion->matriz)[i] = (int*) malloc_shared(sizeof(int) * altura);
       // Asignamos valores
       for (unsigned int j = 0; j < altura; j++) {
       		if (i==anchura-1 || j==altura-1) { // Última fila o columna
       			(tabla_puntuacion->matriz)[i][j] = - MATCHREPLACE_COST;
       		} else {
           		(tabla_puntuacion->matriz)[i][j] = tabla_puntuacion_default[i][j];
       		}
       }
    }

    //Convertimos los alfabetos para su indexación indirecta en la matriz
    tabla_puntuacion->codificacion_horizontal = codifica_alfabeto(alfabeto_default, alfabeto_length);
    tabla_puntuacion->codificacion_vertical = codifica_alfabeto(alfabeto_default, alfabeto_length);

    tabla_puntuacion->id = 1;

    //Finalmente, asignamos la tabla
    *ptr_tabla_puntuacion = tabla_puntuacion;
}

/** Asigna un alfabeto y la matriz a partir de la ruta de un fichero de matrices estándar del NCBI. */
void parse_matrix(char* matrix_file, struct Matriz_puntuaciones** ptr_tabla_puntuacion, int MATCHREPLACE_COST) {
	FILE *fichero;
	if(matrix_file==NULL || (fichero = fopen(matrix_file, "r"))==NULL) {
		printf("Matriz %s couldn't be found, using default values\n", matrix_file);
    	init_matrix_default(ptr_tabla_puntuacion, MATCHREPLACE_COST);
    	return;
  	}

	struct Matriz_puntuaciones* tabla_puntuacion = (struct Matriz_puntuaciones*) malloc_shared(sizeof(struct Matriz_puntuaciones));

	// Creamos los arrays de tamaño máximo MAX_ALFABETO
	tabla_puntuacion->alfabeto_horizontal = (char*) malloc_shared(sizeof(char) * MAX_ALFABETO);
	tabla_puntuacion->alfabeto_vertical = (char*) malloc_shared(sizeof(char) * MAX_ALFABETO);

  	/* Nota: Algoritmo de lectura de matriz basado en el mismo en Java de BioJava */
	int rows = 0, cols = 0;

	// Primera lectura: sacamos los alfabetos y calculamos el tamaño de la matriz
  	char temp[MAX_SIZELINE]; //Suponemos que cada línea no será de más de MAX_SIZELINE caracteres
	while (fgets(temp, MAX_SIZELINE, fichero)!=NULL) {
		if (temp[0]!='#') { //Los comentarios son ignorados
			if (temp[0]==' ') { //Es el primer alfabeto
				unsigned int longitud = strlen(temp);
				for (unsigned int i=0; i<longitud; i++){
					if (temp[i]!=' ' && temp[i]!='\n' && temp[i]!='\r') { //Suponemos que es una letra
						tabla_puntuacion->alfabeto_horizontal[cols] = temp[i];
						cols++;
					}
				}
			} else if (temp[0]!='\n') { //Es la matriz, cogemos el primer caracter
				tabla_puntuacion->alfabeto_vertical[rows] = temp[0];
				rows++;
			}
		}
  	}
  	//Cerramos ya que nunca llegarán a MAX_ALFABETO
	tabla_puntuacion->alfabeto_horizontal[cols] = FINCAD;
	tabla_puntuacion->alfabeto_vertical[rows] = FINCAD;

    //Convertimos los alfabetos para su indexación indirecta en la matriz
    tabla_puntuacion->codificacion_horizontal = codifica_alfabeto(tabla_puntuacion->alfabeto_horizontal, cols);
    tabla_puntuacion->codificacion_vertical = codifica_alfabeto(tabla_puntuacion->alfabeto_vertical, rows);

	// Segunda lectura: creamos y rellenamos la matriz
	// Nota IMPORTANTE: +1 para poder asignar valores MATCHREPLACE a los nucleótidos que no aparecen en la matriz
	// haciendo uso de una nueva fila y columna, de valores MATCHREPLACE_COST
	unsigned int anchura = cols+1;
	unsigned int altura = rows+1;
	rows=0; cols=0;
	fseek(fichero, 0, SEEK_SET);

	tabla_puntuacion->matriz = (int**) malloc_shared(sizeof(int*) * anchura);
    for (unsigned int i = 0; i < anchura; i++) {
       (tabla_puntuacion->matriz)[i] = (int*) malloc_shared(sizeof(int) * altura);
       // Los valores ya los asignaremos luego
    }

    // Última fila y columna (no contemplada en el fichero de la matriz) con MATCHREPLACE_COST
    for (unsigned int i = 0; i < anchura; i++)
       (tabla_puntuacion->matriz)[i][altura-1] = MATCHREPLACE_COST;
    for (unsigned int j = 0; j < altura; j++)
       (tabla_puntuacion->matriz)[anchura-1][j] = MATCHREPLACE_COST;

	char* linea;
	char* nuevalinea;
	int valor;
	int fin_iteracion;
	while (fgets(temp, MAX_SIZELINE, fichero)!=NULL) {
		if (temp[0]!='#' && (temp[0]!=' ' && temp[0]!='\n' && temp[0]!='\r')) { //Comentarios y primera línea
			//Descartamos el primer elemento y el resto vamos añadiendo
			linea = temp;
			nuevalinea = linea;
			nuevalinea++;
			cols=0;
			//Mientras siga habiendo enteros (linea!=nuevalinea), añadimos a la matriz
			fin_iteracion = 0;
			while (!fin_iteracion) {
				linea=nuevalinea;
				valor = strtol(linea, &nuevalinea, 0);
				if (linea!=nuevalinea) {
					(tabla_puntuacion->matriz)[cols][rows]=valor;
					cols++;
				} else {
					fin_iteracion=1;
				}
			}
			rows++;
		}
  	}

    tabla_puntuacion->id = 1;

    //Finalmente, asignamos la tabla
    *ptr_tabla_puntuacion = tabla_puntuacion;
}

/** A partir de un alfabeto ya determinado, crea un alfabeto codificado para realizar
 *  indexaciones indirectas en la matriz de su tabla de puntuaciones:
 *   Ya que no aparecerán caracteres mayores de 128, este nuevo array convertirá un caracter
 *   al índice dentro de la matriz de alfabetos para su secuencia (horizontal o vertical),
 *   creando una relación biunívoca. */
unsigned short * codifica_alfabeto(const char * alfabeto, int longAlfabeto){
	unsigned short i, j;
	unsigned short * alfabetoConvertido;
	alfabetoConvertido = (unsigned short *)malloc_shared(MAX_ALFABETO * sizeof(unsigned short));

	// Inicialización, todos a fuera de rango
	for(i=0;i<MAX_ALFABETO; i++)
		alfabetoConvertido[i] = FUERA_DE_ALFABETO;

	// Recorrido parte buenos: guardamos los índices de la matriz en sus posiciones ASCII de la versión en maýusculas y en minúsculas
	for(i=0; i<longAlfabeto; i++) {
		alfabetoConvertido[tolower(alfabeto[i])] = i;
		alfabetoConvertido[toupper(alfabeto[i])] = i;
	}

	// Recorrido parte malos: los que no son "buenos": índice i (siguiente al número de símbolos)
	for(j=0; j<MAX_ALFABETO; j++)
	    if(alfabetoConvertido[j]==FUERA_DE_ALFABETO)
			alfabetoConvertido[j] = i;

	// Asignamos código de guión especial
//	alfabetoConvertido['-'] = CODIGO_GUION;
	// Asignamos código de guión a columna de "Otros": versión óptima en tiempo. Se podría controlar esta condición aparte pero sufriría el rendimiento para un caso no habitual: su valor será MATCHREPLACE
	alfabetoConvertido['-'] = i;

    return alfabetoConvertido;
}

/** Realiza una copia completa de la matriz de puntuaciones y sus alfabetos en MEMORIA LOCAL:
 *  Necesario para optimizar accesos en workers, evitando continuos accesos a memoria compartida. */
void matrix_copy(struct Matriz_puntuaciones** tabla_copia, struct Matriz_puntuaciones* tabla_original){
	struct Matriz_puntuaciones* tabla_puntuacion = (struct Matriz_puntuaciones*) malloc_safe_agr248(sizeof(struct Matriz_puntuaciones));

	// Alfabetos: hemos de copiar el NULL del final
	unsigned int horizontal_length = strlen(tabla_original->alfabeto_horizontal)+1;//+1 para el NULL final, evitará desbordamientos
	tabla_puntuacion->alfabeto_horizontal = (char*) malloc_safe_agr248(sizeof(char) * horizontal_length);
	memcpy(tabla_puntuacion->alfabeto_horizontal, tabla_original->alfabeto_horizontal, horizontal_length);
	unsigned int vertical_length = strlen(tabla_original->alfabeto_vertical)+1;//+1 para el NULL final, evitará desbordamientos
	tabla_puntuacion->alfabeto_vertical = (char*) malloc_safe_agr248(sizeof(char) * vertical_length);
	memcpy(tabla_puntuacion->alfabeto_vertical, tabla_original->alfabeto_vertical, vertical_length);

	// Alfabetos codificados
	tabla_puntuacion->codificacion_horizontal = (unsigned short *) malloc_safe_agr248(sizeof(unsigned short) * MAX_ALFABETO);
	memcpy(tabla_puntuacion->codificacion_horizontal, tabla_original->codificacion_horizontal, MAX_ALFABETO * sizeof(unsigned short));
	tabla_puntuacion->codificacion_vertical = (unsigned short *) malloc_safe_agr248(sizeof(unsigned short) * MAX_ALFABETO);
	memcpy(tabla_puntuacion->codificacion_vertical, tabla_original->codificacion_vertical, MAX_ALFABETO * sizeof(unsigned short));

	// Valores
	// Nota IMPORTANTE: +1 para poder asignar valores MATCHREPLACE a los nucleótidos que no aparecen en la matriz
	// haciendo uso de una nueva fila y columna, de valores MATCHREPLACE_COST
	// El +1 ya lo tiene de antes, luego no hay que añadirlo
	unsigned int anchura = horizontal_length;
	unsigned int altura = vertical_length;
	tabla_puntuacion->matriz = (int**) malloc_safe_agr248(sizeof(int*) * anchura);
    for (unsigned int i = 0; i < anchura; i++) {
       (tabla_puntuacion->matriz)[i] = (int*) malloc_safe_agr248(sizeof(int) * altura);
       // Asignamos valores
       for (unsigned int j = 0; j < anchura; j++) {
           (tabla_puntuacion->matriz)[i][j] = tabla_original->matriz[i][j];
       }
    }

	tabla_puntuacion->id = tabla_original->id;

	//Finalmente, asignamos la tabla
	*tabla_copia = tabla_puntuacion;
}

/** Libera la memoria de una matriz copiada previamente. */
void matrix_free(struct Matriz_puntuaciones** ptr_tabla){
	struct Matriz_puntuaciones* tabla = *ptr_tabla;
	// Liberamos la matriz: Anchura + 1 por el extra MATCHREPLACE
	unsigned int anchura = strlen(tabla->alfabeto_horizontal) + 1;

	for (unsigned int x = 0; x < anchura; x++)
    	free(tabla->matriz[x]);
    free(tabla->matriz);

    // Liberamos los alfabetos y codificaciones
	free(tabla->alfabeto_horizontal);
	free(tabla->alfabeto_vertical);
	free(tabla->codificacion_horizontal);
	free(tabla->codificacion_vertical);

	// Liberamos el propio puntero a la matriz
	free(*ptr_tabla);

	*ptr_tabla = NULL;
}

/** Imprime la matriz de puntuaciones. */
 void imprime_matriz_puntuaciones(struct Matriz_puntuaciones* tabla_puntuacion){
	printf("Matrix ID: %d", tabla_puntuacion->id);
	// Nota IMPORTANTE: +1 para poder asignar valores MATCHREPLACE a los nucleótidos que no aparecen en la matriz
	// haciendo uso de una nueva fila y columna, de valores MATCHREPLACE_COST
 	unsigned int altura = strlen(tabla_puntuacion->alfabeto_vertical)+1;
 	unsigned int anchura = strlen(tabla_puntuacion->alfabeto_horizontal)+1;
 	for (unsigned int i = 0; i < anchura; i++){
 		if (i == anchura-1) {
 			printf("\tM-R");
 		} else {
 			printf("\t%c", tabla_puntuacion->alfabeto_horizontal[i]);
 		}
 	}
 	for (unsigned int j = 0; j < altura; j++){
 		if (j == altura-1) {
 			printf("\nM-R");
 		} else {
 			printf("\n%c", tabla_puntuacion->alfabeto_vertical[j]);
 	    }
 		for (unsigned int i = 0; i < anchura; i++){
 			printf("\t%d", tabla_puntuacion->matriz[i][j]);
 		}
 	}
 	printf("\n");
 }



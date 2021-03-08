/*
 * alphabet.h
 *
 *  Created on: Mar 8, 2021
 *      Author: galvez
 */

#ifndef ALPHABET_H_
#define ALPHABET_H_

#include "definitions.h"

#define CODIGO_GUION 127
#define FUERA_DE_ALFABETO 126
#define MAX_ALFABETO 128 //Número máximo de caracteres en un alfabeto: tamaño del array

/* Alfabeto completo con el que trabajará al algoritmo (depende de la matriz). */

/* VALORES CONSTANTES a usar por defecto */
static const char alfabeto_default[6] = "ACGTN"; //Alfabeto por defecto: A, C, G, T, N
static const int tabla_puntuacion_default[5][5]= {{+1, -1, -1, -1, -2},{-1, +1, -1, -1, -2},{-1, -1, +1, -1, -2},{-1, -1, -1, +1, -2},{-2, -2, -2, -2, -1}};
//Matriz: sub NUC4.2 (EDNAMAT):
//const int tablaPuntuacion_default[5][5]= {{+5, -4, -4, -4, -2},{-4, +5, -4, -4, -2},{-4, -4, +5, -4, -2},{-4, -4, -4, +5, -2},{-2, -2, -2, -2, -1}};

extern void init_matrix_default(struct Matriz_puntuaciones** ptr_tabla_puntuacion, int MATCHREPLACE_COST);

extern void parse_matrix(char* matrix_file, struct Matriz_puntuaciones** ptr_tabla_puntuacion, int MATCHREPLACE_COST);

extern void matrix_copy(struct Matriz_puntuaciones** tabla_copia, struct Matriz_puntuaciones* tabla_original);

extern void imprime_matriz_puntuaciones(struct Matriz_puntuaciones* tabla_puntuacion);

extern void matrix_free(struct Matriz_puntuaciones** ptr_tabla);

#endif /* ALPHABET_H_ */

/*
 * NWSW-LSA-512.h
 *
 *  Created on: Mar 9, 2021
 *      Author: galvez
 */

#ifndef NWSW_LSA_512_H_
#define NWSW_LSA_512_H_

#include <stdint.h>
#include <stdio.h>

/*****************************/
/* CLI PARAMETERS MANAGEMENT */
/*****************************/

#define VERSION 1.2 // Current program version
#define MAX_FILENAME_SIZE 512 // Maximum length of the full path of a file

/* User options */
/* parameters that can be specified by CLI */
#define __HELP "--help"
#define __1PASS "--1pass"
#define __2PASS "--2pass"
#define __NJ "--nj"
#define __UPGMA "--upgma"
#define __GLOBAL "--global"
#define __LOCAL "--local"
//#define __EXTENDED "--extended"
//#define __BASIC "--basic"
#define __CISC "--cisc"
#define __128 "--128"
#define __256 "--256"
#define __512 "--512"
#define __INSERT "--insert"
#define __DELETE "--delete"
#define __MATCHREPLACE "--matchreplace"
#define __THREADS "--threads"
#define __PARALLEL "--parallel"
#define __GAPEXTEND "--gapextend"
#define __MATRIX "--matrix"
#define __VERBOSE "--verbose"

/* Enumerations for parameters with exclusive values */
enum Pass { ONLY_SCORE, FULL_ALIGNMENT }; // Calculus of score or full alignment
enum Tree { NEIGHBOUR_JOIN, UPGMA, NONE_TREE };
enum Algorithm {NEEDLEMAN_WUNSCH, SMITH_WATERMAN }; // Alignment global or local
enum Vectorization { CISC, SSE3, AVX2, AVX512 };
//enum Info {BASIC, EXTENDED}; // It should be produced Basic or Extended information file as output

/* Default values for the alignment algorithm. */
#define DEFAULT_PASS ONLY_SCORE
#define DEFAULT_TREE NONE_TREE
#define DEFAULT_ALGORITHM NEEDLEMAN_WUNSCH
#define DEFAULT_VECTORIZATION SSE3
//#define DEFAULT_INFO EXTENDED
#define DEFAULT_INSERT_COST 4 // Cost of opening a gap in the Query sequence [horizontal]
#define DEFAULT_DELETE_COST 4 // Cost of opening a gap in the Subject sequence [vertical]
#define DEFAULT_MATCHREPLACE_COST 1 // Default cost when a nucleotide to compare is not in the score matrix
#define DEFAULT_GAPEXTEND_COST 1 // Cost of extending a gap once it has been opened
#define DEFAULT_THREADS getNumberOfCores() // Number of cores given by the system
#define DEFAULT_PARALLEL 1 // Number of pairwise alignments to be calculated in parallel when multifasta is provided
#define DEFAULT_VERBOSE 0 // Default verbose value is FALSE

/* Structure to contain user parameters */
struct UserParameters {
	enum Pass pass;
	enum Tree tree;
	enum Algorithm algorithm;
	enum Vectorization vectorization;
//	enum Info info;
	uint32_t insert;
	uint32_t delete;
	uint32_t matchReplace __attribute__((aligned (4)));
	uint32_t gapExtend __attribute__((aligned (4)));
	char matrixFilename[MAX_FILENAME_SIZE];
	uint16_t threads;
	uint16_t parallel;
	unsigned char verbose;
	/* Single pairwise filenames */
	char queryFilename[MAX_FILENAME_SIZE];				// No default. Mandatory
	char subjectFilename[MAX_FILENAME_SIZE];			// No default. Mandatory
//	char infoFilename[MAX_FILENAME_SIZE];				// No default. Removed
	char alignFilename[MAX_FILENAME_SIZE];
	/* Multiple pairwise filenames */
	char multifastaFilename[MAX_FILENAME_SIZE];			// No default. Mandatory
	char newickFilename[MAX_FILENAME_SIZE];				// No default. Mandatory
};

struct UserParameters defaultUserParameters();
void displayUserParameters(struct UserParameters *, FILE *);
int startsWith(char *text, char *subquery);
char * afterEquals(char *str, char *prefix);
uint32_t parseInt(char * str, unsigned char *ok);
int isEmptyOrNull(char *str);

#define enumPassToString(value) ((value) == ONLY_SCORE)? __1PASS : __2PASS
#define enumTreeToString(value) ((value) == NEIGHBOUR_JOIN)? __NJ : (((value) == UPGMA)? __UPGMA : "None")
#define enumAlgorithmToString(value) ((value) == NEEDLEMAN_WUNSCH)? __GLOBAL : __LOCAL
#define enumVectorizationToString(value) ((value) == CISC)? __CISC : ((value) == SSE3)? __128 : ((value) == AVX2)? __256 : __512
//#define enumInfoToString(value) ((value) == BASIC)? __BASIC : __EXTENDED

#endif /* NWSW_LSA_512_H_ */

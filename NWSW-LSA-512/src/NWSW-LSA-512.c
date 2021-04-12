/*
 ============================================================================
 Name        : NWSW-LSA-512.c
 Author      : SGR
 Version     :
 Copyright   : Open source
 Description : Global and local alignment in C, Ansi-style
 ============================================================================
 */

/* Algorithm Needleman-Wunsch and Smith-Waterman parallel with vectorization, that uses shared memory. */

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <ctype.h>
#include <unistd.h>
#include <limits.h>

#include "NWSW-LSA-512.h"
#include "Utilities.h"
#include "GlobalData.h"
#include "Sequence.h"
#include "JobTable.h"
#include "SingleAlignment.h"
#include "Worker.h"

/** DISPLAYS USAGE */
void printHelp() {
	printf("nwsw_lsa v%.1f. A program to perform single or multiple pairwise alignments\n", VERSION);
	printf("nwsw_lsa is a program to perform single or multiple pairwise alignments.\n");
	printf("Usage \t (single pairwise):  nwsw_lsa [options] query.fa subject.fa [align.fa]\n");
	printf("\t (multiple pairwise): nwsw_lsa {%s | %s} [options] sequences.fa newick.tree\n", __NJ, __UPGMA);
	printf("*** Options:\n");
	printf("[%s]            \t Displays this helping information and exits.\n", __HELP);
	printf("[%s | %s] \t Calculates only the score (%s) or the full alignment (%s).\n\t\t\t If %s then a file align.fa must be indicated. Default %s.\n", __1PASS, __2PASS, __1PASS, __2PASS, __2PASS, enumPassToString(DEFAULT_PASS));
	printf("[%s | %s] \t Performs many to many pairwise alignments and produces a Newick tree.\n\t\t\t This overrides [%s | %s] option.\n", __NJ, __UPGMA, __1PASS, __2PASS);
	printf("[%s | %s]\t Calculates NW (%s) or SW (%s). Both with affine gaps. Default %s\n", __GLOBAL, __LOCAL, __GLOBAL, __LOCAL, enumAlgorithmToString(DEFAULT_ALGORITHM));
//	printf("[%s | %s]\t Provides a file with basic info on execution or with extended information (%s). Default %s\n", __EXTENDED, __BASIC, __EXTENDED, enumInfoToString(DEFAULT_INFO));
	printf("[%s=value]    \t Cost value for open gap in the query sequence. By default %d.\n", __INSERT, DEFAULT_INSERT_COST);
	printf("[%s=value]    \t Cost value for open gap in the subject sequence. By default %d.\n", __DELETE, DEFAULT_DELETE_COST);
	printf("[%s=value] \t Cost value when any nucleotide to compare is not in the score matrix. By default %d.\n", __MATCHREPLACE, DEFAULT_MATCHREPLACE_COST);
	printf("[%s=value]    \t Cost value to extend a gap in any sequence. By default %d.\n", __GAPEXTEND, DEFAULT_GAPEXTEND_COST);
	printf("[%s=file]     \t Filename with the score matrix to use. By default NUC4.4.\n", __MATRIX);
	printf("[%s=value]   \t Number of threads to use in the Forward stage of the alignment.\n", __THREADS);
	printf("[%s=value]   \t Number of alignments to be calculated in parallel. Default %d.\n", __PARALLEL, DEFAULT_PARALLEL);
	printf("[%s]         \t Displays messages on algorithm execution.\n", __VERBOSE);
	printf("Single Pairwise mode (no %s or %s are specified):\n", __NJ, __UPGMA);
	printf("\t*** Input files:\n");
	printf("\tquery.fa   \t Filename with Query sequence in fasta format.\n");
	printf("\tsubject.fa \t Filename with Subject sequence in fasta format.\n");
	printf("\t*** Output files:\n");
//	printf("align.info \t Filename where to put basic or extended information on the alignment process.\n");
	printf("\t[align.fa] \t Filename where to put the final alignment. Mandatory when %s is used.\n", __2PASS);
	printf("Multiple Pairwise mode (%s or %s are specified):\n", __NJ, __UPGMA);
	printf("\t*** Input files:\n");
	printf("\tsequences.fa   \t Filename with many sequences in multifasta format.\n");
	printf("\t*** Output files:\n");
	printf("\tnewick.tree \t Filename where to put the final Newick tree.\n");
}

int singlePairwise(struct UserParameters *ptrUserParams);
int multiplePairwise(struct UserParameters *ptrUserParams);

/** MAIN PROGRAM: OBTAINS PARAMETERS AND CALCULATES AN ALIGNMENT OR MULTIPLE ALIGNMENTS */
//
int main(int argc, char *argv[]) {
	// This sentence makes the output to console in the same order than expected. Required by Eclipse
	setvbuf(stdout, NULL, _IONBF, 0);
	struct UserParameters userParams = defaultUserParameters();
	int pos=1;
	unsigned char goodValue;
	unsigned char ok = 1;
	int positionalParameter = 0;
	char filenames[3][MAX_FILENAME_SIZE];
	for(int i=0; i<3; i++) strcpy(filenames[i], "");
	while (pos < argc) {
		if (! strcmp(__HELP, argv[pos]) ) {
			printHelp(); exit(EXIT_SUCCESS);
		} else if(! strcmp(__1PASS, argv[pos]) ) {
			userParams.pass = ONLY_SCORE; pos++;
		} else if(! strcmp(__2PASS, argv[pos]) ) {
			userParams.pass = FULL_ALIGNMENT; pos++;
		} else if(! strcmp(__NJ, argv[pos]) ) {
			userParams.tree = NEIGHBOUR_JOIN; pos++;
		} else if(! strcmp(__UPGMA, argv[pos]) ) {
			userParams.tree = UPGMA; pos++;
		} else if(! strcmp(__GLOBAL, argv[pos]) ) {
			userParams.algorithm = NEEDLEMAN_WUNSCH; pos++;
		} else if(! strcmp(__LOCAL, argv[pos]) ) {
			userParams.algorithm = SMITH_WATERMAN; pos++;
//		} else if(! strcmp(__EXTENDED, argv[pos]) ) {
//			userParams.info = EXTENDED; pos++;
//		} else if(! strcmp(__BASIC, argv[pos]) ) {
//			userParams.info = BASIC; pos++;
		} else if(! strcmp(__VERBOSE, argv[pos]) ) {
			userParams.verbose = 1; pos++;
		} else if(startsWith(argv[pos], __INSERT) ) {
			userParams.insert = parseInt(afterEquals(argv[pos], __INSERT), &goodValue);
			if (!goodValue) { fprintf(stderr, "Unknown value in: %s\n", argv[pos]); ok = 0; }
			pos++;
		} else if(startsWith(argv[pos], __MATCHREPLACE) ) {
			userParams.matchReplace = parseInt(afterEquals(argv[pos], __MATCHREPLACE), &goodValue);
			if (!goodValue) { fprintf(stderr, "Unknown value in: %s\n", argv[pos]); ok = 0; }
			pos++;
		} else if(startsWith(argv[pos], __DELETE) ) {
			userParams.delete = parseInt(afterEquals(argv[pos], __DELETE), &goodValue);
			if (!goodValue) { fprintf(stderr, "Unknown value in: %s\n", argv[pos]); ok = 0; }
			pos++;
		} else if(startsWith(argv[pos], __GAPEXTEND) ) {
			userParams.gapExtend = parseInt(afterEquals(argv[pos], __GAPEXTEND), &goodValue);
			if (!goodValue) { fprintf(stderr, "Unknown value in: %s\n", argv[pos]); ok = 0; }
			pos++;
		} else if(startsWith(argv[pos], __THREADS) ) {
			userParams.threads = parseInt(afterEquals(argv[pos], __THREADS), &goodValue);
			if (!goodValue || userParams.threads < 1) { fprintf(stderr, "Invalid value in: %s\n", argv[pos]); ok = 0; }
			pos++;
		} else if(startsWith(argv[pos], __PARALLEL) ) {
			userParams.parallel = parseInt(afterEquals(argv[pos], __PARALLEL), &goodValue);
			if (!goodValue || userParams.parallel < 1) { fprintf(stderr, "Invalid value in: %s\n", argv[pos]); ok = 0; }
			pos++;
		} else if(startsWith(argv[pos], __MATRIX) ) {
			char *ptrFilename = afterEquals(argv[pos], __MATRIX);
			if (ptrFilename != NULL)
				strcpy(userParams.matrixFilename, ptrFilename);
			if (access(userParams.matrixFilename, R_OK)) {
				fprintf(stderr, "Unable to read file '%s' in: %s\n", userParams.matrixFilename, argv[pos]); ok = 0;
			}
			pos++;
		} else if(startsWith(argv[pos], "-") ){
			fprintf(stderr, "Unknown parameter: %s\n", argv[pos++]); ok = 0;
		} else {
			// We assume it is a positional parameter
			if (positionalParameter < 3) {
				strcpy(filenames[positionalParameter], argv[pos]);
			} else {
				fprintf(stderr, "Too many positional parameters: %s\n", argv[pos]); ok=0;
			}
			positionalParameter++; pos++;
		}
	}
	// Sets filenames depending on --nj and --upgma
	if (ok && userParams.tree == NONE_TREE) { // PAIRWISE ALIGNMENT
		if (userParams.parallel != DEFAULT_PARALLEL)
			fprintf(stderr, "Ignoring %s parameter (%d)\n", __PARALLEL, userParams.parallel);
		if (strcmp(filenames[0], "")) { // It is the Query sequence filename
			strcpy(userParams.queryFilename, filenames[0]);
			if (access(userParams.queryFilename, R_OK)) {
				fprintf(stderr, "Unable to read file '%s' in: %s\n", userParams.queryFilename, argv[pos]); ok = 0;
			}
		} else {
			fprintf(stderr, "Missing Query sequence filename\n"); ok = 0;
		}
		if (strcmp(filenames[1], "")) { // It is the Subject sequence filename
			strcpy(userParams.subjectFilename, filenames[1]);
			if (access(userParams.subjectFilename, R_OK)) {
				fprintf(stderr, "Unable to read file '%s' in: %s\n", userParams.subjectFilename, argv[pos]); ok = 0;
			}
		} else {
			fprintf(stderr, "Missing Subject sequence filename\n"); ok = 0;
		}
//			} else if (positionalParameter == 2) { // It is the Information output filename
//				strcpy(userParams.infoFilename, argv[pos]);
		 // It is the Alignment output filename
		strcpy(userParams.alignFilename, filenames[2]);
		// Checks consistency between --2pass and resulting alignment filename
		if (userParams.pass == FULL_ALIGNMENT && isEmptyOrNull(userParams.alignFilename) ) {
			fprintf(stderr, "The filename to store the alignment must be specified\n"); ok = 0;
		} else if (userParams.pass == ONLY_SCORE && ! isEmptyOrNull(userParams.alignFilename) ) {
			fprintf(stderr, "The filename to store the alignment must not be specified when calculating the score only\n"); ok = 0;
		}
	} else if (ok) {
		if (strcmp(filenames[0], "")) { // It is the Multifasta sequence filename
			strcpy(userParams.multifastaFilename, filenames[0]);
			if (access(userParams.multifastaFilename, R_OK)) {
				fprintf(stderr, "Unable to read file '%s' in: %s\n", userParams.multifastaFilename, argv[pos]); ok = 0;
			}
		} else {
			fprintf(stderr, "Missing Multifasta sequence filename\n"); ok = 0;
		}
		if (strcmp(filenames[1], "")) { // It is the Newick output  filename
			strcpy(userParams.newickFilename, filenames[1]);
		} else {
			fprintf(stderr, "Missing Newick output filename\n"); ok = 0;
		}
	}
	if (! ok)
		fatalError0("Please, use the option --help to display correct usage.\n");

	//
	/* Everything goes OK */
	if (userParams.tree == NONE_TREE)
		return singlePairwise(&userParams);
	else {
		userParams.pass = ONLY_SCORE;
		return multiplePairwise(&userParams);
	}
}


/* Returns a structure with default parameters. Returns it by copy */
struct UserParameters defaultUserParameters() {
	struct UserParameters userParams;
	userParams.pass 		= DEFAULT_PASS;
	userParams.tree			= DEFAULT_TREE;
	userParams.algorithm 	= DEFAULT_ALGORITHM;
//	userParams.info 		= DEFAULT_INFO;
	userParams.insert 		= DEFAULT_INSERT_COST;
	userParams.delete  		= DEFAULT_DELETE_COST;
	userParams.matchReplace = DEFAULT_MATCHREPLACE_COST;
	userParams.gapExtend	= DEFAULT_GAPEXTEND_COST;
	strcpy(userParams.matrixFilename, "");
	userParams.threads 		= DEFAULT_THREADS;
	userParams.parallel		= DEFAULT_PARALLEL;
	userParams.verbose		= DEFAULT_VERBOSE;
	//userParams.queryFilename			// No default. Mandatory
	//userParams.subjectFilename			// No default. Mandatory
	//userParams.infoFilename				// No default. Removed
	strcpy(userParams.alignFilename, "");
	return userParams;
}

/* Displays the content of a structure with the User's parameters */
void displayUserParameters(struct UserParameters * ptrUserParams, FILE * out) {
	fprintf(out, "User parameters are:\n");
	if (ptrUserParams->tree == NONE_TREE)
		fprintf(out, "Pass: %s\n", enumPassToString(ptrUserParams->pass));
	else
		fprintf(out, "Guide tree: %s\n", enumTreeToString(ptrUserParams->tree));
	fprintf(out, "Algorithm: %s\n", enumAlgorithmToString(ptrUserParams->algorithm));
//	fprintf(out, "Output info: %s\n", enumInfoToString(ptrUserParams->info));
	fprintf(out, "Insert: %d\n", ptrUserParams->insert);
	fprintf(out, "Delete: %d\n", ptrUserParams->delete);
	fprintf(out, "Match-Replace: %d\n", ptrUserParams->matchReplace);
	fprintf(out, "Gap extend: %d\n", ptrUserParams->gapExtend);
	fprintf(out, "Threads: %d\n", ptrUserParams->threads);
	if (ptrUserParams->tree != NONE_TREE)
		fprintf(out, "Parallel: %d\n", ptrUserParams->parallel);
	fprintf(out, "Verbose: %s\n", ptrUserParams->verbose? "True" : "False");
	fprintf(out, "Score matrix: %s\n", strcmp(ptrUserParams->matrixFilename, "")? ptrUserParams->matrixFilename : "Internal NUC4.4");
	fprintf(out, "Query filename: %s\n", ptrUserParams->queryFilename);
	fprintf(out, "Subject filename: %s\n", ptrUserParams->subjectFilename);
//	fprintf(out, "Information output filename: %s\n", ptrUserParams->infoFilename);
	fprintf(out, "Alignment output filename: %s\n", ptrUserParams->alignFilename);
}

/* Returns TRUE if str starts with the characters of prefix. Otherwise FALSE */
int startsWith(char *str, char *prefix) {
	return strncmp(prefix, str, strlen(prefix)) == 0;
}

/* Returns a pointer to the next character of the '=' sign located at position strlen(prefix).
 * If that position does not hold a '=' the NULL is returned
 */
char * afterEquals(char *str, char *prefix) {
	if (str[strlen(prefix)] != '=') return NULL;
	return str+strlen(prefix)+1;
}

/* Parses a text composed only by digits into an unsigned integer and returns it.
 * If the string is empty or contains invalid digits returns FALSE into ok.
 * Otherwise returns TRUE into ok.
 */
uint32_t parseInt(char * str, unsigned char *ok) {
	// If the string is empty: Error
	if (str == NULL || strlen(str) == 0) { *ok=0; return 0; }
	int ret = 0;
	for (int i = 0; str[i] != '\0'; ++i)
		if (str[i]<'0' || str[i]>'9') { *ok=0; return 0; } // If it isn't a digit: Error
		else ret = ret * 10 + str[i] - '0';
	*ok=1;
	return ret;
}

/* Returns TRUE if str is NULL or "" */
int isEmptyOrNull(char *str) {
	return (str == NULL) || str[0] == 0;
}



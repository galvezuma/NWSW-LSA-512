/*
 ============================================================================
 Name        : NWSW-LSA-512.c
 Author      : SGR
 Version     :
 Copyright   : Código abierto
 Description : Global and local alignment in C, Ansi-style
 ============================================================================
 */

/* Algoritmo Needleman-Wunsch y Smith-Waterman paralelo con vectorización, que utiliza memoria compartida. */

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <ctype.h>
#include <unistd.h>

#include "NWSW-LSA-512.h"
#include "Utilities.h"
#include "GlobalData.h"

/** DISPLAYS USAGE */
void printHelp() {
	printf("nwsw_lsa v%.1f. Usage: nwsw_lsa [options] query.fa subject.fa align.info [align.fa]\n", VERSION);
	printf("*** Options:\n");
	printf("[%s]            \t Displays this helping information and exits.\n", __HELP);
	printf("[%s | %s] \t Calculates only the score (%s) or the full alignment (%s).\n\t\t\t If %s then a file align.fa must be indicated. Default %s\n", __1PASS, __2PASS, __1PASS, __2PASS, __2PASS, enumPassToString(DEFAULT_PASS));
	printf("[%s | %s]\t Calculates NW (%s) or SW (%s). Both with affine gaps. Default %s\n", __GLOBAL, __LOCAL, __GLOBAL, __LOCAL, enumAlgorithmToString(DEFAULT_ALGORITHM));
	printf("[%s | %s]\t Provides a file with basic info on execution or with extended information (%s). Default %s\n", __EXTENDED, __BASIC, __EXTENDED, enumInfoToString(DEFAULT_INFO));
	printf("[%s=value]    \t Cost value for open gap in the query sequence. By default %d.\n", __INSERT, DEFAULT_INSERT_COST);
	printf("[%s=value]    \t Cost value for open gap in the subject sequence. By default %d.\n", __DELETE, DEFAULT_DELETE_COST);
	printf("[%s=value]     \t Cost value when any nucleotide to compare is not in the score matrix. By default %d.\n", __MATCHREPLACE, DEFAULT_MATCHREPLACE_COST);
	printf("[%s=value]    \t Cost value to extend a gap in any sequence. By default %d.\n", __GAPEXTEND, DEFAULT_GAPEXTEND_COST);
	printf("[%s=file]     \t Filename with the score matrix to use. By default NUC4.4.\n", __MATRIX);
	printf("[%s=value]   \t Number of threads to use in the Forward stage of the alignment.\n", __THREADS);
	printf("[%s]         \t Displays messages on algorithm execution.\n", __VERBOSE);
	printf("*** Input files:\n");
	printf("query.fa   \t Filename with Query sequence in fasta format.\n");
	printf("subject.fa \t Filename with Subject sequence in fasta format.\n");
	printf("*** Output files\n");
	printf("align.info \t Filename where to put basic or extended information on the alignment process.\n");
	printf("[align.fa] \t Filename where to put the final alignment. Mandatory when --2pass is used.\n");
}

/** MAIN PROGRAM: OBTAINS PARAMETERS AND CALCULATES AN ALIGNMENT  */
int main(int argc, char *argv[]) {
	// This sentence makes the output to console in the same order than expected. Required by Eclipse
	setvbuf(stdout, NULL, _IONBF, 0);
	struct UserParameters userParams = defaultUserParameters();
	int pos=1;
	unsigned char goodValue;
	unsigned char ok = 1;
	int positionalParameter = 0;
	while (pos < argc) {
		if (! strcmp(__HELP, argv[pos]) ) {
			printHelp(); exit(EXIT_SUCCESS);
		} else if(! strcmp(__1PASS, argv[pos]) ) {
			userParams.pass = ONLY_SCORE; pos++;
		} else if(! strcmp(__2PASS, argv[pos]) ) {
			userParams.pass = FULL_ALIGNMENT; pos++;
		} else if(! strcmp(__GLOBAL, argv[pos]) ) {
			userParams.algorithm = NEEDELMAN_WUNSCH; pos++;
		} else if(! strcmp(__LOCAL, argv[pos]) ) {
			userParams.algorithm = SMITH_WATERMAN; pos++;
		} else if(! strcmp(__EXTENDED, argv[pos]) ) {
			userParams.info = EXTENDED; pos++;
		} else if(! strcmp(__BASIC, argv[pos]) ) {
			userParams.info = BASIC; pos++;
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
			if (!goodValue) { fprintf(stderr, "Unknown value in: %s\n", argv[pos]); ok = 0; }
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
			if (positionalParameter == 0) { // It is the Query sequence filename
				strcpy(userParams.queryFilename, argv[pos]);
				if (access(userParams.queryFilename, R_OK)) {
					fprintf(stderr, "Unable to read file '%s' in: %s\n", userParams.queryFilename, argv[pos]); ok = 0;
				}
			} else if (positionalParameter == 1) { // It is the Subject sequence filename
				strcpy(userParams.subjectFilename, argv[pos]);
				if (access(userParams.subjectFilename, R_OK)) {
					fprintf(stderr, "Unable to read file '%s' in: %s\n", userParams.subjectFilename, argv[pos]); ok = 0;
				}
			} else if (positionalParameter == 2) { // It is the Information output filename
				strcpy(userParams.infoFilename, argv[pos]);
			} else if (positionalParameter == 3) { // It is the Alignment output filename
				strcpy(userParams.alignFilename, argv[pos]);
			} else {
				fprintf(stderr, "Too many positional parameters: %s\n", argv[pos]); ok=0;
			}
			positionalParameter++; pos++;
		}
	}
	// Additional Checks of parameters
	// 1.- Consistency between --2pass and resulting alignment filename
	if (ok) {
		if (userParams.pass == FULL_ALIGNMENT && isEmptyOrNull(userParams.alignFilename) ) {
			fprintf(stderr, "The filename to store the alignment must be specified\n"); ok = 0;
		} else if (userParams.pass == ONLY_SCORE && ! isEmptyOrNull(userParams.alignFilename) ) {
			fprintf(stderr, "The filename to store the alignment must not be specified when calculating the score only\n"); ok = 0;
		}
		if (positionalParameter < 1) { fprintf(stderr, "Missing Query sequence filename\n"); ok = 0; }
		if (positionalParameter < 2) { fprintf(stderr, "Missing Subject sequence filename\n"); ok = 0; }
		if (positionalParameter < 3) { fprintf(stderr, "Missing Information output filename\n"); ok = 0; }
	}
	// If some parameter is wrong, we exit
	if (! ok) {
		fatalError0("Please, use the option --help to display correct usage.\n");
	}
	//
	// If everything goes OK then we display the parameters and continue
	if (userParams.verbose) {
		fprintf(stdout, "nwsw_lsa v%.1f. Universidad de Málaga. Spain.\n", VERSION);
		displayUserParameters(&userParams, stdout);
	}

	/* MAIN PROCESSING */
	struct GlobalData globalData;
	// Reads the Query and Subject sequence Fasta files
	readFastaFile(userParams.queryFilename, &(globalData.query));
	readFastaFile(userParams.subjectFilename, &(globalData.subject));
	if (userParams.verbose) {
		fprintf(stdout, "****** INPUT SEQUENCE ******\n");
		fprintf(stdout, "Query (%ld nt): %s\n", strlen(globalData.query.data), globalData.query.name);
		fprintf(stdout, "Subject (%ld nt) %s\n", strlen(globalData.subject.data), globalData.subject.name);
	}

	return(EXIT_SUCCESS);
}

/* Returns a structure with default parameters. Returns it by copy */
struct UserParameters defaultUserParameters() {
	struct UserParameters userParams;
	userParams.pass 		= DEFAULT_PASS;
	userParams.algorithm 	= DEFAULT_ALGORITHM;
	userParams.info 		= DEFAULT_INFO;
	userParams.insert 		= DEFAULT_INSERT_COST;
	userParams.delete  		= DEFAULT_DELETE_COST;
	userParams.matchReplace = DEFAULT_MATCHREPLACE_COST;
	userParams.gapExtend	= DEFAULT_GAPEXTEND_COST;
	strcpy(userParams.matrixFilename, "");
	userParams.threads 		= DEFAULT_THREADS;
	userParams.verbose		= DEFAULT_VERBOSE;
	//userParams.queryFilename			// No default. Mandatory
	//userParams.subjectFilename			// No default. Mandatory
	//userParams.infoFilename				// No default. Mandatory
	strcpy(userParams.alignFilename, "");
	return userParams;
}

/* Displays the content of a structure with the User's parameters */
void displayUserParameters(struct UserParameters * ptrUserParams, FILE * out) {
	fprintf(out, "User parameters are:\n");
	fprintf(out, "Pass: %s\n", enumPassToString(ptrUserParams->pass));
	fprintf(out, "Algorithm: %s\n", enumAlgorithmToString(ptrUserParams->algorithm));
	fprintf(out, "Output info: %s\n", enumInfoToString(ptrUserParams->info));
	fprintf(out, "Insert: %d\n", ptrUserParams->insert);
	fprintf(out, "Delete: %d\n", ptrUserParams->delete);
	fprintf(out, "Match-Replace: %d\n", ptrUserParams->matchReplace);
	fprintf(out, "Gap extend: %d\n", ptrUserParams->gapExtend);
	fprintf(out, "Threads: %d\n", ptrUserParams->threads);
	fprintf(out, "Verbose: %s\n", ptrUserParams->verbose? "True" : "False");
	fprintf(out, "Score matrix: %s\n", strcmp(ptrUserParams->matrixFilename, "")? ptrUserParams->matrixFilename : "Internal NUC4.4");
	fprintf(out, "Query filename: %s\n", ptrUserParams->queryFilename);
	fprintf(out, "Subject filename: %s\n", ptrUserParams->subjectFilename);
	fprintf(out, "Information output filename: %s\n", ptrUserParams->infoFilename);
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

/** Reads a Sequence from a Fasta file. Takes the header and the content. */
void readFastaFile(char* filename, struct Sequence *seq) {
	char* fastaHeader = seq->name;
	char** fastaContent = &(seq->data);
	char * ok;
	FILE *file;
	if((file = fopen(filename, "r")) == NULL) {
    	fatalError1("Unable to open the input Fasta file: %s\n", filename); // This HALTS
  	}

  	// The first line is the header (fastaHeader FASTA)
  	ok = fgets(fastaHeader, MAX_SIZEHEADER, file);
  	if ((! ok) || (fastaHeader[0]!='>')) { // If it is empty o does not starts with '>', then it is not Fasta
    	fatalError1("Fasta file %s should begin with >.\n", filename); // This HALTS
  	}
  	int sizeFirstLine = strlen(fastaHeader);
  	// Remove the CR and LF at the end of the header
  	while(fastaHeader[strlen(fastaHeader)-1] == '\n' || fastaHeader[strlen(fastaHeader)-1] == '\r')
  		fastaHeader[strlen(fastaHeader)-1] = '\0';

  	// Let's allocate space for the Fasta content (the sequence in nucleotides)
  	// Actually, we allocate more space than required because we store also de CR and LF
  	fseek(file, 0, SEEK_END); // Let's go to the end of the file to estimate the maximum size of the Fasta sequence
  	long size = ftell(file) - sizeFirstLine;
  	*fastaContent = (char *) internalMalloc(sizeof(char) * (size + 1));

  	// We will iterate through the lines of the file
  	int offset = 0; //Initially the allocated memory is empty
  	fseek(file, sizeFirstLine, SEEK_SET); // We start reading after the header of the Fasta file
  	while (fgets((*fastaContent)+offset, MAX_SIZELINE, file) != NULL) { // Data is directly read over the allocated memory
  		int tempOffset = offset;
  		offset += strlen((*fastaContent)+offset);
  		while ((tempOffset < offset) && ( (*fastaContent)[offset-1]=='\n' || (*fastaContent)[offset-1]=='\r' ) ){ // CR and LF are removed
  			offset--;
		}
  	}
  	(*fastaContent)[offset] = '\0'; // The Fasta string must be finished (delimited by \0)

	fclose(file); // The file should be closed
}



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
#include <limits.h>

#include "NWSW-LSA-512.h"
#include "Utilities.h"
#include "Fasta.h"
#include "GlobalData.h"
#include "Sequence.h"
#include "ScoreMatrix.h"
#include "JobTable.h"
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
	printf("[%s=value]     \t Cost value when any nucleotide to compare is not in the score matrix. By default %d.\n", __MATCHREPLACE, DEFAULT_MATCHREPLACE_COST);
	printf("[%s=value]    \t Cost value to extend a gap in any sequence. By default %d.\n", __GAPEXTEND, DEFAULT_GAPEXTEND_COST);
	printf("[%s=file]     \t Filename with the score matrix to use. By default NUC4.4.\n", __MATRIX);
	printf("[%s=value]   \t Number of threads to use in the Forward stage of the alignment.\n", __THREADS);
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
int executePairwise(struct GlobalData *gd);

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


int multiplePairwise(struct UserParameters *ptrUserParams) {
	/* LOAD FILES */
	// The next structure is used throughout the functions of the program
	struct GlobalData globalData;
	// Loads globally the User parameters
	copyUserParameters(&globalData, ptrUserParams);

	// Reads the Score Matrix file or loads the default one
	if (! strcmp(ptrUserParams->matrixFilename, "")) {
		loadDefaultMatrix(&globalData);
	} else {
		loadUserMatrix(ptrUserParams->matrixFilename, &globalData);
	}
	// Checks that in multi-pairwise the horizontal and vertical alphabets are used interchangeably and, therefore,
	// they should be the same
	if (strcmp(globalData.scoreMatrix.horizontalAlphabet, globalData.scoreMatrix.verticalAlphabet))
		fatalError2("In multi-pairwise alignments, both horizontal and vertical alphabets in matrix should be the same: %s versus %s",
					globalData.scoreMatrix.horizontalAlphabet,
					globalData.scoreMatrix.verticalAlphabet);

	// Reads a file with many fasta sequences
	struct NodeListSequence *ptrNodesFasta = readMultipleFastaFile(ptrUserParams->multifastaFilename);
	if (ptrNodesFasta == NULL || ptrNodesFasta->ptrNext == NULL)
		fatalError1("Multifasta file %s must have two sequences at least.\n", ptrUserParams->multifastaFilename);

	// VERBOSE
	if (ptrUserParams->verbose) {
		fprintf(stdout, "nwsw_lsa v%.1f. Universidad de Málaga. Spain.\n", VERSION);
		displayUserParameters(ptrUserParams, stdout);
		// Display sequences in multifasta file
		fprintf(stdout, "****** INPUT SEQUENCES ******\n");
		int i = 0;
		for (struct NodeListSequence *ptr = ptrNodesFasta; ptr != NULL; ptr = ptr->ptrNext)
			fprintf(stdout, "Sequence %d (%ld nt): %s\n", i++, strlen(ptr->sequence.data), ptr->sequence.name);
		fprintf(stdout, "****** SCORE MATRIX ******\n");
		// Display matrix content
		displayScoreMatrix(stdout, &(globalData.scoreMatrix));
	}

	/* INITIALIZATION */
	// Checks the letters in the sequences and put them in Uppercase
	// Check is made only with horizontal codification
	int i = 0;
	for (struct NodeListSequence *ptr = ptrNodesFasta; ptr != NULL; ptr = ptr->ptrNext, i++) {
		int numInvalidLettersInSequence;
		// Checks the Query
		numInvalidLettersInSequence = toUpperCodeAndCheck(	&ptr->sequence,
													globalData.scoreMatrix.horizontalAlphabet,
													globalData.scoreMatrix.horizontalCodification
												 );
		if (numInvalidLettersInSequence != 0) fprintf(stderr, "Sequence %d has %d letters not found in the alphabet.\n", i, numInvalidLettersInSequence);
	}
	const int numSeq = i;

	// PREPARE STRUCTURE TO STORE SCORES
	int scores[numSeq-1][numSeq-1];
	/*****************************************/
	/* MAIN CALL TO MULTI-PAIRWISE ALIGNMENT */
	/*****************************************/
	i = 0;
	for (struct NodeListSequence *ptrI = ptrNodesFasta; ptrI->ptrNext != NULL; ptrI = ptrI->ptrNext, i++) {
		int j = 0;
		for (struct NodeListSequence *ptrJ = ptrI->ptrNext; ptrJ != NULL; ptrJ = ptrJ->ptrNext, j++) {
			globalData.query = ptrI->sequence;
			globalData.subject = ptrJ->sequence;
			/* MAIN CALL TO EXECUTE PAIRWISE ALIGNMENT */
			int score = executePairwise(&globalData);
			scores[i][j] = scores[j][i] = score;
			fprintf(stdout, "The score %dx%d is: %d\n", i, j, score);
		}
	}
	/* RELEASE MEMORY */
	freeNodeListSequenceStruct(&ptrNodesFasta);
	freeScoreMatrixStruct(&(globalData.scoreMatrix));
	return(EXIT_SUCCESS);

}

/********************************************************/
/* Main Function to calculate Single Pairwise alignment */
/********************************************************/
//
int singlePairwise(struct UserParameters *ptrUserParams) {
	/* LOAD FILES */
	// The next structure is used throughout the functions of the program
	struct GlobalData globalData;
	// Loads globally the User parameters
	copyUserParameters(&globalData, ptrUserParams);

	// Reads the Score Matrix file or loads the default one
	if (! strcmp(ptrUserParams->matrixFilename, "")) {
		loadDefaultMatrix(&globalData);
	} else {
		loadUserMatrix(ptrUserParams->matrixFilename, &globalData);
	}

	// Reads the Query and Subject sequence Fasta files
	readSingleFastaFile(ptrUserParams->queryFilename, &(globalData.query));
	readSingleFastaFile(ptrUserParams->subjectFilename, &(globalData.subject));

	// VERBOSE
	if (ptrUserParams->verbose) {
		fprintf(stdout, "nwsw_lsa v%.1f. Universidad de Málaga. Spain.\n", VERSION);
		displayUserParameters(ptrUserParams, stdout);
		fprintf(stdout, "****** INPUT SEQUENCES ******\n");
		fprintf(stdout, "Query (%ld nt): %s\n", strlen(globalData.query.data), globalData.query.name);
		fprintf(stdout, "Subject (%ld nt) %s\n", strlen(globalData.subject.data), globalData.subject.name);
		fprintf(stdout, "****** SCORE MATRIX ******\n");
		displayScoreMatrix(stdout, &(globalData.scoreMatrix));
	}

	/* INITIALIZATION */
	// Checks the letters in the sequences and put them in Uppercase
	int invalidLettersQuery, invalidLettersSubject;
	// Checks the Query
	invalidLettersQuery = toUpperCodeAndCheck(	&(globalData.query),
												globalData.scoreMatrix.horizontalAlphabet,
												globalData.scoreMatrix.horizontalCodification
											 );
	if (invalidLettersQuery != 0) fprintf(stderr, "Query has %d letters not found in the alphabet.\n", invalidLettersQuery);

	// Checks the subject
	invalidLettersSubject = toUpperCodeAndCheck(&(globalData.subject),
												globalData.scoreMatrix.verticalAlphabet,
												globalData.scoreMatrix.verticalCodification
											   );
	if (invalidLettersSubject != 0) fprintf(stderr, "Subject has %d letters not found in the alphabet.\n", invalidLettersSubject);

	/***********************************/
	/* MAIN CALL TO PAIRWISE ALIGNMENT */
	/***********************************/
	int score = executePairwise(&globalData);
	fprintf(stdout, "The score is: %d\n", score);
	/***********************************/
	/***********************************/

	/* RELEASE MEMORY */
	freeScoreMatrixStruct(&(globalData.scoreMatrix));
	freeSequenceStruct(&(globalData.query));
	freeSequenceStruct(&(globalData.subject));
	return(EXIT_SUCCESS);
}

int executePairwise(struct GlobalData *gd) {
	/* ESTIMATES THE FRAGMENTS SIZE AND CREATES THE TABLE OF JOBS */
	createJobTable(gd);
	/* INITIALIZES THE STACK OF JOBS */
	gd->stackOfJobs = createStack(gd, max(gd->jobTable.numFragments_X, gd->jobTable.numFragments_Y));
	push(&gd->stackOfJobs, getJob(&gd->jobTable, 0, 0)); // The job at position 0,0 is the only available by now
	/* CHECKS SUPPORT FOR VECTORIZATION */
	checkSupport(gd);


	/* MAIN PROCESSING */
//	displayJob(getJob(&globalData.jobTable, 0, 0));
//	displayJob(getJob(&globalData.jobTable, 1, 0));
//	printf("%d, %d\n", globalData.jobTable.numFragments_X, globalData.jobTable.numFragments_Y);
		long initTime = myClock();
		gd->bestScore = -INFINITE;
		// PREPARE SYNCHRONIZATION
		int error;
		gd->jobTableFulfilled = 0;
		if (gd->verbose) {
			error = pthread_mutex_init(&gd->verboseStdOut_mutex, NULL);
			if (error) fatalError0("Unable to initialize the verbose mutex.\n");
		}
		error = pthread_mutex_init(&gd->globalDataAccess_mutex, NULL);
		if (error) fatalError0("Unable to initialize the main mutex.\n");
		error = pthread_cond_init(&gd->jobAvailable_condition, NULL);
		if (error) fatalError0("Unable to initialize the main synchronization condition.\n");
		// CREATE AND LAUNCH THREADS
		if (gd->verbose) printf("Creating %d threads.\n", gd->threads);
		pthread_t * threads= (pthread_t *) internalMalloc(sizeof(pthread_t) * gd->threads);
		for(int i=0; i<gd->threads; i++){
			error = pthread_create(&(threads[i]), NULL, threadWorker, (void *)gd);
			if (error) fatalError1("Unable to create the thread num. %d\n", i);
		}
		// WAIT FOR THREADS
		for(int i=0; i<gd->threads; i++){
			error = pthread_join(threads[i], NULL);
			if (error) fatalError1("Unable to join to the thread num. %d\n", i);
		}
		// DESTROY THREADS AND SYNCHRONIZATION ELEMENTS
		error = pthread_cond_destroy(&gd->jobAvailable_condition);
		if (error) fatalError0("Unable to destroy the main synchronization condition.\n");
		error = pthread_mutex_destroy(&gd->globalDataAccess_mutex);
		if (error) fatalError0("Unable to destroy the main mutex.\n");
		if (gd->verbose) {
			error = pthread_mutex_destroy(&gd->verboseStdOut_mutex);
			if (error) fatalError0("Unable to destroy the verbose mutex.\n");
		}
		internalFree((void **) &threads);
		// VERBOSE DISPLAY TIME
		long endTime = myClock();
		if (gd->verbose)
			fprintf(stdout, "Time execution: %.3f\n", (float)(endTime-initTime)/1000);

	/* RELEASE MEMORY */
	freeStack(&gd->stackOfJobs);
	freeJobTableStruct(&gd->jobTable, gd->pass);
	return gd->bestScore;
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



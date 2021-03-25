/*
 * ScoreMatrix.c
 *
 *  Created on: Mar 11, 2021
 *      Author: galvez
 */

#include <stdio.h>
#include <ctype.h>
#include <string.h>

#include "Utilities.h"
#include "GlobalData.h"
#include "ScoreMatrix.h"

// Loads the matrix from an already open FILE
void loadMatrixFromFILE(FILE *file, struct GlobalData *gd);
// Converts a string (alphabet) into an array of indexes
void recodeAlphabet(uint8_t *codedAlphabet, const char *alphabet);

/* Loads the basic NUC4.4 score matrix */
// Creates a temporal file with the data and loads it as an User score matrix though a common internal function
void loadDefaultMatrix( struct GlobalData *gd) {
	// Opens a temporal FILE
	FILE* tmp = tmpfile();
	if (tmp == NULL)
		fatalError0("Unable to create temporal file.\n");
	// Writes the NUC4.4
	fprintf(tmp, "# This matrix was created by Todd Lowe   12/10/92            \n");
	fprintf(tmp, "#                                                            \n");
	fprintf(tmp, "# Uses ambiguous nucleotide codes, probabilities rounded to  \n");
	fprintf(tmp, "#  nearest integer                                           \n");
	fprintf(tmp, "#                                                            \n");
	fprintf(tmp, "# Lowest score = -4, Highest score = 5                       \n");
	fprintf(tmp, "#                                                            \n");
	fprintf(tmp, "    A   T   G   C   S   W   R   Y   K   M   B   V   H   D   N\n");
	fprintf(tmp, "A   5  -4  -4  -4  -4   1   1  -4  -4   1  -4  -1  -1  -1  -2\n");
	fprintf(tmp, "T  -4   5  -4  -4  -4   1  -4   1   1  -4  -1  -4  -1  -1  -2\n");
	fprintf(tmp, "G  -4  -4   5  -4   1  -4   1  -4   1  -4  -1  -1  -4  -1  -2\n");
	fprintf(tmp, "C  -4  -4  -4   5   1  -4  -4   1  -4   1  -1  -1  -1  -4  -2\n");
	fprintf(tmp, "S  -4  -4   1   1  -1  -4  -2  -2  -2  -2  -1  -1  -3  -3  -1\n");
	fprintf(tmp, "W   1   1  -4  -4  -4  -1  -2  -2  -2  -2  -3  -3  -1  -1  -1\n");
	fprintf(tmp, "R   1  -4   1  -4  -2  -2  -1  -4  -2  -2  -3  -1  -3  -1  -1\n");
	fprintf(tmp, "Y  -4   1  -4   1  -2  -2  -4  -1  -2  -2  -1  -3  -1  -3  -1\n");
	fprintf(tmp, "K  -4   1   1  -4  -2  -2  -2  -2  -1  -4  -1  -3  -3  -1  -1\n");
	fprintf(tmp, "M   1  -4  -4   1  -2  -2  -2  -2  -4  -1  -3  -1  -1  -3  -1\n");
	fprintf(tmp, "B  -4  -1  -1  -1  -1  -3  -3  -1  -1  -3  -1  -2  -2  -2  -1\n");
	fprintf(tmp, "V  -1  -4  -1  -1  -1  -3  -1  -3  -3  -1  -2  -1  -2  -2  -1\n");
	fprintf(tmp, "H  -1  -1  -4  -1  -3  -1  -3  -1  -3  -1  -2  -2  -1  -2  -1\n");
	fprintf(tmp, "D  -1  -1  -1  -4  -3  -1  -1  -3  -1  -3  -2  -2  -2  -1  -1\n");
	fprintf(tmp, "N  -2  -2  -2  -2  -1  -1  -1  -1  -1  -1  -1  -1  -1  -1  -1\n");

    rewind(tmp); // Positions at the beginning.
	loadMatrixFromFILE(tmp, gd); // Calls common function to load from a FILE
    fclose(tmp); // Though it is assumed that the created temporary file will automatically be deleted after the termination of program.
}

/* Loads an User matrix given its filename */
void loadUserMatrix(const char *filename, struct GlobalData *gd) {
	// Tries to open the User file
	FILE* file = fopen(filename, "rt");
	if (file == NULL)
		fatalError0("Unable to create temporal file.\n");
	loadMatrixFromFILE(file, gd); // Calls common function to load from a FILE
	fclose(file);

}

// Loads the matrix from an already open FILE
void loadMatrixFromFILE(FILE *file, struct GlobalData *gd) {

	struct ScoreMatrix* scoreMatrix = & gd->scoreMatrix;

  	/* Note: Algorithm to read a matrix based on the same one in BioJava */
	unsigned int rows = 0, cols = 0;

	// FIRST pass: we obtain the alphabets and calculate the size of the matrix
  	char line[MAX_SIZELINE]; // We assume each line is, at most, MAX_SIZELINE length in bytes
	while (fgets(line, MAX_SIZELINE, file)!=NULL) {
		if (line[0]=='#' || line[0]=='\n') continue; // Comments and blank lines are ignored
		if (line[0]==' ') { // Horizontal alphabet. All the letters in the same line
			for (size_t i=0; i<strlen(line); i++) {
				if (line[i]!=' ' && line[i]!='\n' && line[i]!='\r')  // We assume it is a valid alphabet character
					scoreMatrix->horizontalAlphabet[cols++] = line[i];
			}
		} else { // Vertical alphabet. A letter per line.
			scoreMatrix->verticalAlphabet[rows++] = line[0]; // We assume it is the matrix with a character alphabet at position 0
		}
  	}
  	// We close the string securely (their length will never reach MAX_LENGTH_ALPHABET
	scoreMatrix->horizontalAlphabet[cols] = '\0';
	scoreMatrix->verticalAlphabet[rows] = '\0';

    // Alphabets are coded to to point directly to the matrix
    recodeAlphabet(scoreMatrix->horizontalCodification, scoreMatrix->horizontalAlphabet);
    recodeAlphabet(scoreMatrix->verticalCodification, scoreMatrix->verticalAlphabet);

	// SECOND pass: the matrix is created and fulfilled
	// IMPORTANT note: The matrix has an additional row and column to store the MATCHREPLACE values associated
    // to the nucleotides (letters) that do not appear in the alphabets.
    // This is why width and height have +1 over cols and rows respectively.
	unsigned int width = cols+1;
	unsigned int height = rows+1;
	fseek(file, 0, SEEK_SET); // It's the same than rewind(file);

	// Firstly, the empty matrix is created
	scoreMatrix->matrix = (int**) internalMalloc(sizeof(int*) * width);
    for (unsigned int i = 0; i < width; i++) {
       (scoreMatrix->matrix)[i] = (int*) internalMalloc(sizeof(int) * height);
    }

    // Secondly, the last row and the last column (not included in the original FILE) are fulfilled with MATCHREPLACE_COST
    for (unsigned int i = 0; i < width; i++)
       (scoreMatrix->matrix)[i][height-1] = -gd->matchReplace;
    for (unsigned int j = 0; j < height; j++)
       (scoreMatrix->matrix)[width-1][j] = -gd->matchReplace;

	char* ptr;
	char* newPtr;
	int value;

	rows = 0;
	while (fgets(line, MAX_SIZELINE, file)!=NULL) {
		if (line[0]=='#' || line[0]==' ' || line[0]=='\n' || line[0]=='\r') continue; // Skip comments, blank lines and horizontal alphabet
		// First letter is discarded. The rest must be numbers to add to the matrix
		ptr = newPtr = line+1;
		cols=0;
		// While integers are available (ptr !=newPtr), we add them to the matrix
		value = strtol(ptr, &newPtr, 0);
		while (ptr != newPtr) {
			(scoreMatrix->matrix)[cols][rows]=value;
			cols++;
			ptr=newPtr;
			value = strtol(ptr, &newPtr, 0);
		}
		rows++;
  	}

    // scoreMatrix->id = 1; // Just in case of adaptive matrix
}

/* Input is an alphabet already loaded.
 * Creates a coded alphabet as a direct index to the matrix of the ScoreMatrix:
 *   Because there is no character greater than 127, this new array converts a character into an index
 *   to the array (alphabet) where such a character is associated to its sequence (horizontal or vertical).
 *   This way a biunivocal (one-to-one) relationship is created.
 */
void recodeAlphabet(uint8_t *codedAlphabet, const char *alphabet) {
	uint8_t i;

	// Initializes all positions as Out of Range
	for(i=0;i<MAX_LENGTH_ALPHABET; i++)
		codedAlphabet[i] = OUT_OF_ALPHABET;

	/* Iterates though the VALID letters of the Alphabet. */
	// The indexes of the array are stored in their ASCII positions. This can be made for both upper and lower case
	for(i=0; i<strlen(alphabet); i++) {
		codedAlphabet[tolower(alphabet[i])] = i;
		codedAlphabet[toupper(alphabet[i])] = i;
	} // Value of 'i' after the loop will be used later

	/* Iterates though the INVALID letters (those that are not in the Alphabet). */
	// The indexes of the invalid letters point beyond the last valid letter of the alphabet
	for(uint8_t j=0; j<MAX_LENGTH_ALPHABET; j++)
	    if(codedAlphabet[j]==OUT_OF_ALPHABET)
			codedAlphabet[j] = i; // 'i' contains the index of the next position beyond the last valid letter in the alphabet
}

void displayScoreMatrix(FILE *out, struct ScoreMatrix *scoreMatrix) {
	// The scoreMatrix has one more row and column to address invalid letters
 	unsigned int width = strlen(scoreMatrix->horizontalAlphabet)+1;
 	unsigned int height = strlen(scoreMatrix->verticalAlphabet)+1;

 	// Displays the first line with the horizontal alphabet
 	for (unsigned int i = 0; i < width; i++)
 		fprintf(out, "\t%s", (i == width-1)? "M-R" : toString(scoreMatrix->horizontalAlphabet[i]));
 	fprintf(out, "\n");

 	// Displays next lines, one by one.
 	for (unsigned int j = 0; j < height; j++) {
 		// Firstly the alphabet letter ...
 		fprintf(out, "%s", (j == height-1)? "M-R" : toString(scoreMatrix->verticalAlphabet[j]));
 		// ... and then the scores
 		for (unsigned int i = 0; i < width; i++)
 			fprintf(out, "\t%d", scoreMatrix->matrix[i][j]);
 	 	fprintf(out, "\n");
 	}
}

// Frees the memory used by a struct ScoreMatrix
void freeScoreMatrixStruct(struct ScoreMatrix *scoreMatrix) {
	// The scoreMatrix has one more row and column to address invalid letters
 	unsigned int width = strlen(scoreMatrix->horizontalAlphabet)+1;
    for (unsigned int i = 0; i < width; i++) {
       internalFree((void **) &((scoreMatrix->matrix)[i]));
    }
    internalFree((void **) &(scoreMatrix->matrix));
}

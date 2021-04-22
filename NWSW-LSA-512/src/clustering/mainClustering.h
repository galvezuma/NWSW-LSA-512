/*
 * mainNJ.h
 *
 *  Created on: Apr 7, 2021
 *      Author: galvez
 */

#ifndef CLUSTERING_MAINCLUSTERING_H_
#define CLUSTERING_MAINCLUSTERING_H_

int mainNJ(int numSequences, double *ptrArray, char **arraySpeciesName, int verbose, char *filename);
int mainUPGMA(int numSequences, double *ptrArray, char **arraySpeciesName, int verbose, char *filename);

#endif /* CLUSTERING_MAINCLUSTERING_H_ */

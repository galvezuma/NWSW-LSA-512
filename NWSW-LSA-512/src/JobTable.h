/*
 * JobTable.h
 *
 *  Created on: Mar 12, 2021
 *      Author: galvez
 */

#ifndef JOBTABLE_H_
#define JOBTABLE_H_

#include "NWSW-LSA-512.h"
struct GlobalData;

/* Default values for fragments size. */
#define DEFAULT_FRAGMENTSIZE_X 8*500
#define DEFAULT_FRAGMENTSIZE_Y 8*500

struct Node {
	int t; // End in horizontal gap
	int u; // End in vertical gap
	int s; // Best value
};

struct Job {
	unsigned x, y;
	uint32_t realSize_X;
	uint32_t realSize_Y;
	struct Node * volatile ptrRow;
	struct Node * volatile ptrColumn;
};

struct JobTable {
	uint32_t fragmentSize_X;
	uint32_t fragmentSize_Y;
	int numFragments_X; // Num of X fragments of the logical Job table. No matter how many pass (--1pass or --2pass)
	int numFragments_Y; // Num of Y fragments of the logical Job table. No matter how many pass (--1pass or --2pass)
	struct Node **rows;
	struct Node **columns;
	struct Job *jobs;
};

void createJobTable(struct GlobalData *gd);
void freeJobTableStruct(struct JobTable * jobTable, enum Pass pass);
struct Job *getJob(struct JobTable *jobTable, int posX, int posY);

/* DEBUG FUNCTIONS */
void displayJobFromJobTable(struct JobTable *jobTable, int x, int y);
void displayJob(struct Job *job);
void displayNode(struct Node *node);

#endif /* JOBTABLE_H_ */

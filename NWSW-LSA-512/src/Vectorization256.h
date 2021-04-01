/*
 * Vectorization256.h
 *
 *  Created on: Mar 25, 2021
 *      Author: galvez
 */

#ifndef VECTORIZATION256_H_
#define VECTORIZATION256_H_

#include "NWSW-LSA-512.h"
#include "GlobalData.h"
#include "JobTable.h"

int processJob_256(struct GlobalData *gd, struct Job *job, struct Node *retFragX, struct Node *retFragY);

#endif /* VECTORIZATION256_H_ */

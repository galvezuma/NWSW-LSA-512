/*
 * Vectorization512.h
 *
 *  Created on: Mar 25, 2021
 *      Author: galvez
 */

#ifndef VECTORIZATION512_H_
#define VECTORIZATION512_H_

#include "NWSW-LSA-512.h"
#include "GlobalData.h"
#include "JobTable.h"

int processJob_512(struct GlobalData *gd, struct Job *job, struct Node *retFragX, struct Node *retFragY);

#endif /* VECTORIZATION512_H_ */

/*
 * Vectorization128.h
 *
 *  Created on: Mar 29, 2021
 *      Author: galvez
 */

#ifndef VECTORIZATION128_H_
#define VECTORIZATION128_H_

#include "NWSW-LSA-512.h"
#include "GlobalData.h"
#include "JobTable.h"

int processJob_128(struct GlobalData *gd, struct Job *job, struct Node *retFragX, struct Node *retFragY);

#endif /* VECTORIZATION128_H_ */

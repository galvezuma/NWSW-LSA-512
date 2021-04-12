/*
 * VectorizationKNC.h
 *
 *  Created on: Apr 5, 2021
 *      Author: galvez
 */

#ifndef VECTORIZATIONKNC_H_
#define VECTORIZATIONKNC_H_

#include "NWSW-LSA-512.h"
#include "GlobalData.h"
#include "JobTable.h"

int processJob_KNC(struct GlobalData *gd, struct Job *job, struct Node *retFragX, struct Node *retFragY);

#endif /* VECTORIZATIONKNC_H_ */

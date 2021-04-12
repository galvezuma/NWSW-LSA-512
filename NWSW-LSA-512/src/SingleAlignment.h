/*
 * Alignment.h
 *
 *  Created on: Apr 7, 2021
 *      Author: galvez
 */

#ifndef SINGLEALIGNMENT_H_
#define SINGLEALIGNMENT_H_

#include "NWSW-LSA-512.h"
#include "GlobalData.h"

int singlePairwise(struct UserParameters *ptrUserParams);
int executePairwise(struct GlobalData *gd);

#endif /* SINGLEALIGNMENT_H_ */

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
#include "BackwardsPairwise.h"

struct PairwiseAlignmentResult {
	struct FastaPairwiseAlignment fastaPairwiseAlignment;
	int score;
	double normalizedScore;
};

int singlePairwise(struct UserParameters *ptrUserParams);
struct PairwiseAlignmentResult executePairwise(struct GlobalData *gd);

#endif /* SINGLEALIGNMENT_H_ */

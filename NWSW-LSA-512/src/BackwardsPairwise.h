/*
 * BackwardsPairwise.h
 *
 *  Created on: Apr 12, 2021
 *      Author: galvez
 */

#ifndef BACKWARDSPAIRWISE_H_
#define BACKWARDSPAIRWISE_H_

#include "Sequence.h"
#include "GlobalData.h"

struct FastaPairwiseAlignment {
	struct Sequence *sequence[2];
	char *alignment[2];
};

struct FastaPairwiseAlignment getFastaAlignment(struct GlobalData * globalData);

#endif /* BACKWARDSPAIRWISE_H_ */

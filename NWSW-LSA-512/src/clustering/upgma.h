/*
 * upgma.h
 *
 *  Created on: Apr 22, 2021
 *      Author: galvez
 */

#ifndef CLUSTERING_UPGMA_H_
#define CLUSTERING_UPGMA_H_

/* THIS FOLLOWS THE SAME STRUCTURE THAN THE NJ */

/* Find the nearest clusters in dmat and write their index in c1 and c2. */
void upgma_find_nearest_clusters(const dist_matrix *dmat, uint32_t *c1, uint32_t *c2);

/* Create a new dist_matrix joining the two specified clusters c1, c2 in a new one called new_name. */
dist_matrix *upgma_join_clusters(const dist_matrix *dmat, const char *new_name, uint32_t c1, uint32_t c2);

btree_storage *upgma_tree_init(const dist_matrix *dmat, btree_node **leafs);
void upgma_tree_add_node(const dist_matrix *dmat, btree_storage *storage, btree_node **partial_trees, const char *name, uint32_t c1, uint32_t c2);


#endif /* CLUSTERING_UPGMA_H_ */

/*
 * mainNJ.c
 *
 *  Created on: Apr 7, 2021
 *      Author: galvez
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <assert.h>
#include <inttypes.h>
#include <math.h>

#include "../clustering/dist_matrix.h"
#include "../clustering/neighbour_joining.h"
#include "../clustering/upgma.h"
#include "../clustering/phylogenetic_tree.h"
#include "../clustering/utilities.h"


dist_matrix *loadDistMatrix(int species_count, double *ptrArray, char **arraySpeciesName);

int mainNJ(int numSequences, double *ptrArray, char **arraySpeciesName, int verbose, char *filename) {

        dist_matrix *dmat = loadDistMatrix(numSequences, ptrArray, arraySpeciesName);

        if (!dmat) {
            // Error, stop here
            return EXIT_FAILURE;
        }

        double u[dmat->species_count];

        /* Ensure that cluster_name has enough space for the longest possible name */
        char cluster_name[2 + 3 * sizeof(dmat->species_count) + 1];
        uint32_t cluster_id = 1;

        btree_node *partial_trees[dmat->species_count];
        btree_storage *tree_storage = nj_tree_init(dmat, partial_trees);

        while (dmat->species_count >= 2) {
            if (verbose) {
                dist_matrix_print(dmat);
                printf("\n");
                btree_print_trees(partial_trees, dmat->species_count);
                printf("\n");
            }

            /* Compute the average distance of each clusters from the others */
            dist_matrix_compute_avg_distances(dmat, u);

            /* Find the pair of nearest clusters */
            uint32_t c1, c2;
            nj_find_nearest_clusters(dmat, u, &c1, &c2);

            /* Generate a name for the new cluster */
            unsigned long result = snprintf(cluster_name, sizeof(cluster_name), "C_%" PRIu32, cluster_id);
            assert(result > 0 && result < sizeof(cluster_name));

            if (verbose) {
                fprintf(stderr, "Joining clusters '%s' and '%s' in '%s'.\n", dmat->species_names[c1], dmat->species_names[c2], cluster_name);
            }

            /* Add a node for the new cluster to the array of partial trees */
            nj_tree_add_node(dmat, u, tree_storage, partial_trees, cluster_name, c1, c2);

            /* Create a new dist_matrix joining the specified clusters */
            dist_matrix *joined = nj_join_clusters(dmat, cluster_name, c1, c2);

            if (joined == NULL) {
                /* Error, stop here */
                break;
            }

            /* Release the old distance matrix */
            dist_matrix_free(dmat);
            dmat = joined;

            cluster_id++;
        }

        btree_node *phyl_tree = partial_trees[0];

        if (verbose) {
            btree_print_tree(phyl_tree);
            printf("\n\n");
        }
        exportToNewick(filename, phyl_tree);
        dist_matrix_free(dmat);
        btree_storage_free(tree_storage);

    return EXIT_SUCCESS;
}

// This is very similar to mainNJ but modifying the way to unify clusters.
//
int mainUPGMA(int numSequences, double *ptrArray, char **arraySpeciesName, int verbose, char *filename) {
	// This loads the square matrix given by the main program (MultiAlignment.c).
    dist_matrix *dmat = loadDistMatrix(numSequences, ptrArray, arraySpeciesName);
    if (!dmat) return EXIT_FAILURE; // This should never happens because we do not open any file

    /* Ensure that cluster_name has enough space for the longest possible name */
    char cluster_name[2 + 3 * sizeof(dmat->species_count) + 1];
    uint32_t cluster_id = 1;

    btree_node *partial_trees[dmat->species_count];
    btree_storage *tree_storage = nj_tree_init(dmat, partial_trees);

    while (dmat->species_count >= 2) {
        if (verbose) {
            dist_matrix_print(dmat); printf("\n");
            btree_print_trees(partial_trees, dmat->species_count); printf("\n");
        }

        /* Find the pair of nearest clusters */
        uint32_t c1, c2;
        upgma_find_nearest_clusters(dmat, &c1, &c2);

        /* Generate a name for the new cluster */
        unsigned long result = snprintf(cluster_name, sizeof(cluster_name), "C_%" PRIu32, cluster_id);
        assert(result > 0 && result < sizeof(cluster_name));

        if (verbose) fprintf(stderr, "Joining clusters '%s' and '%s' in '%s'.\n", dmat->species_names[c1], dmat->species_names[c2], cluster_name);

        /* Add a node for the new cluster to the array of partial trees */
        upgma_tree_add_node(dmat, tree_storage, partial_trees, cluster_name, c1, c2);

        /* Create a new dist_matrix joining the specified clusters */
        dist_matrix *joined = upgma_join_clusters(dmat, cluster_name, c1, c2);

        if (joined == NULL) break; /* Error, stop here */

        /* Release the old distance matrix */
        dist_matrix_free(dmat);
        dmat = joined;

        cluster_id++;
    }

    btree_node *phyl_tree = partial_trees[0];

    if (verbose) {
        btree_print_tree(phyl_tree);
        printf("\n\n");
    }
    exportToNewick(filename, phyl_tree);
    dist_matrix_free(dmat);
    btree_storage_free(tree_storage);

return EXIT_SUCCESS;
}


dist_matrix *loadDistMatrix(int species_count, double *ptrArray, char **arraySpeciesName) {

    dist_matrix *dmat = dist_matrix_init(species_count);

    if (!dmat) {
        perror("NJ: Unable to create distance matrix");
        return NULL;
    }

    for (uint32_t i = 0; i < species_count; i++) {
        /* species name: up to 512 alphabetic or whitespace characters */
        char species_name[513];
        strcpy(species_name, arraySpeciesName[i]);

        trim_trailing_space(species_name);

        dist_matrix_set_species_name(dmat, i, species_name);
        dmat->cluster_sizes[i] = 1;

        for (uint32_t j = 0; j < i; j++) {
            double *element = dist_matrix_element(dmat, i, j);
            *element = ptrArray[(i-1) * (species_count - 1) + j];
        }
    }

    return dmat;
}

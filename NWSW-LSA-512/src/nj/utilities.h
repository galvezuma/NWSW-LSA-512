/*
 * utilities.h
 *
 *  Created on: Apr 7, 2021
 *      Author: galvez
 */

#ifndef NJ_UTILITIES_H_
#define NJ_UTILITIES_H_

#include <stddef.h>

#define member_size(type, member) sizeof(((type *)0)->member)

char *neigh_strdup(const char *s);

size_t trim_trailing_space(char *s);
size_t filename_copy(const char *path, char *dest, size_t size);


#endif /* NJ_UTILITIES_H_ */

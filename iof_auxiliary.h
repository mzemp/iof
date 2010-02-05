/*
** auxiliary.c
**
** Various auxiliary functions and definitions
**
** Written by Marcel Zemp
*/

#ifndef IOF_AUXILIARY_H
#define IOF_AUXILIARY_H

/*
** Physical constants
*/ 

#define k_BOLTZMANN 1.3806503e-23
#define PROTON_MASS 1.6726216e-27

/*
** Function for flipping bytes depending on endianness
*/

void reorder(void *, size_t, size_t);

#endif /* IOF_AUXILIARY_H */

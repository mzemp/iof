/*
** auxiliary.h
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
** Structures
*/

typedef struct cosmological_parameters {

    double OmegaM0;
    double OmegaDM0; 
    double OmegaB0;
    double OmegaL0;
    double OmegaK0;
    double OmegaR0;
    } COSMOLOGICAL_PARAMETERS;

typedef struct unit_system {

    double LBox;
    double GNewton;
    double Hubble0;
    double rhocrit0;
    } UNIT_SYSTEM;

/*
** Function for flipping bytes depending on endianness
*/

void reorder(void *, size_t, size_t);

/*
** Function for correcting positions in periodic boundary conditions cubes
*/

double correct_position(double, double, double);

/*
** Function for putting positions back in box
*/

double put_in_box(double, double, double);

#endif /* IOF_AUXILIARY_H */

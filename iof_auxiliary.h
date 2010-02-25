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
** Structures
*/

typedef struct physical_constants {

    double GNewton;
    double GNewton_Cosmology;
    double k_Boltzmann;
    double m_Proton;
    double m_Neutron;
    double m_Electron;
    double Mo;
    double Lo;
    double rho_crit;
    double rho_crit_Cosmology;
    } PHYSICAL_CONSTANTS;

typedef struct conversion_factors {

    double km_per_s_2_kpc_per_Gyr;
    double kpc_per_Gyr_2_km_per_s;
    } CONVERSION_FACTORS;

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

typedef struct coordinate_transformation {

    /*
    ** Unit scale factors
    */
    double L_usf;
    double T_usf;
    double V_usf;
    double M_usf;
    /*
    ** Coordinate system scale factors
    */
    double L_cssf;
    double V_cssf;
    /*
    ** Coordinate system shifts
    */
    double L_css[3];
    double V_css[3];
    } COORDINATE_TRANSFORMATION;

/*
** Constant structures
*/

const PHYSICAL_CONSTANTS PhysicalConstants;
const CONVERSION_FACTORS ConversionFactors;

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

/*
** Setting default values for coordinate transforamtion (i.e. no transformation)
*/

void set_default_values_coordinate_transformation(COORDINATE_TRANSFORMATION *);

/*
** Function for calculating units transformation
*/

void calculate_units_transformation(UNIT_SYSTEM, UNIT_SYSTEM, COORDINATE_TRANSFORMATION *);

#endif /* IOF_AUXILIARY_H */

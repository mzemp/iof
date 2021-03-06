/*
** iof_auxiliary.h
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

	double c_light;
	double au;
	double yr;
	double G_Newton;
	double k_Boltzmann;
	double uamu;
	double m_proton;
	double m_neutron;
	double m_electron;
	double eV;
	double mu_sun;
	double Mo;
	double Lo;
	double pc;
	double G_Newton_cosmology;
	double rho_crit;
	double rho_crit_cosmology;
	} PHYSICAL_CONSTANTS;

typedef struct conversion_factors {

	double km_per_s_2_kpc_per_Gyr;
	double kpc_per_Gyr_2_km_per_s;
	} CONVERSION_FACTORS;

typedef struct cosmological_parameters {

	double OmegaM0;
	double OmegaD0; 
	double OmegaB0;
	double OmegaL0;
	double OmegaK0;
	double OmegaR0;
	double h0_100;
	} COSMOLOGICAL_PARAMETERS;

typedef struct unit_system {

	double LBox;
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

/*
** Function for calculating E function used throughout cosmology
*/

double Ecosmo(double a, COSMOLOGICAL_PARAMETERS cp);

/*
** Function for calculating unit vectors for spherical coordinates
*/

void calculate_unit_vectors_spherical(double *, double *, double *, double *);

/*
** Function for calculating unit vectors for cylindrical coordinates
*/

void calculate_unit_vectors_cylindrical(double *, double *, double *, double *, double *);

#endif /* IOF_AUXILIARY_H */

/*
** auxiliary.c
**
** Various auxiliary functions and definitons
**
** Written by Marcel Zemp
*/

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <malloc.h>
#include <assert.h>
#include "iof_auxiliary.h"

/*
** Constant structures
*/

const PHYSICAL_CONSTANTS PhysicalConstants = {

    6.67428e-11,         // GNewton [m^3 s^[-2} kg^{-1}
    4.4985021530e-6,     // GNewton_Cosmology [kpc^3 Gyr^{-2} Mo^{-1}]
    1.3806504e-23,       // k_Boltzmann [J K^{-1}]
    1.672621637e-27,     // m_Proton [kg]
    1.67492729e-27,      // m_Neutron [kg]
    9.10938215e-31,      // m_Electron [kg]
    1.9884e30,           // Mo [kg]
    3.845e26,            // Lo [J s^{-1}]
    1.8783e-26,          // rho_crit [h_0^2 kg m^{-3}]
    2.7753662719e2       // rho_crit_Cosmology [h_0^2 Mo kpc^{-3}]
    };

const CONVERSION_FACTORS ConversionFactors = {

    1.02271216511,       // km_per_s_2_kpc_per_Gyr
    0.97779222216        // kpc_per_Gyr_2_km_per_s
    };

/*
** Function for flipping bytes depending on endianness
*/

void reorder(void *pointer, size_t size, size_t n){
    
    unsigned char *buffer, tmp;
    int i, j;
    
    buffer = pointer;
    
    for(i = 0; i < n; i++){
        for (j = 0; j < (size/2); j++) {
	    tmp = buffer[i*size+j];
	    buffer[i*size+j] = buffer[(i+1)*size-1-j];
	    buffer[(i+1)*size-1-j] = tmp;
	    }
	}
    }

/*
** Function for correcting positions in periodic boundary conditions cubes
*/

double correct_position(double c, double r, double l) {

    double d;

    d = r-c;
    if (d > 0.5*l) return r - l;
    else if (d < -0.5*l) return r + l;
    else return r;
    }

/*
** Function for putting positions back in box
*/

double put_in_box(double r, double lmin, double lmax) {

    double l;

    l = lmax - lmin;
    if (r < lmin) return r + l;
    else if (r > lmax) return r - l;
    else return r;
    }

/*
** Setting default values for coordinate transforamtion (i.e. no transformation)
*/

void set_default_values_coordinate_transformation(COORDINATE_TRANSFORMATION *ct) {

    int i;

    ct->L_usf = 1;
    ct->T_usf = 1;
    ct->V_usf = 1;
    ct->M_usf = 1;
    ct->L_cssf = 1;
    ct->V_cssf = 1;
    for (i = 0; i < 3; i++) {
	ct->L_css[i] = 0;
	ct->V_css[i] = 0;
	}
    }

/*
** Function for calculating units transformation
*/

void calculate_units_transformation(UNIT_SYSTEM fromus, UNIT_SYSTEM tous, COORDINATE_TRANSFORMATION *ct) {
    
    ct->L_usf = tous.LBox/fromus.LBox;
    ct->T_usf = fromus.Hubble0/tous.Hubble0;
    ct->V_usf = ct->L_usf/ct->T_usf;
    ct->M_usf = tous.rhocrit0/fromus.rhocrit0*pow(ct->L_usf,3);
    }

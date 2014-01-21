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

	299792458.,			 /* c_light [m s^{-1}] */
	149597870700.,		 /* au [m] */
	31557600.,			 /* yr [s] */
	6.67384e-11,		 /* G_Newton [m^3 s^{-2} kg^{-1}] */
	1.3806488e-23,		 /* k_Boltzmann [J K^{-1}] */
	1.660538921e-27,	 /* uamu [kg] */
	1.672621777e-27,	 /* m_proton [kg] */
	1.674927351e-27,	 /* m_neutron [kg] */
	9.10938291e-31,		 /* m_electron [kg] */
	1.602176565e-19,	 /* eV [J] */
	1.32712440018e20,	 /* mu_sun [m^3 s^{-2}] */
	1.98855e30,			 /* Mo [kg] */
	3.845e26,			 /* Lo [W] */
	3.08567758147e+16,	 /* pc [m] */
	4.4985021521e-6,	 /* G_Newton_cosmology [kpc^3 Gyr^{-2} Mo^{-1}] */
	1.8784e-26,			 /* rho_crit [h_0^2 kg m^{-3}] */
	2.7753662721e2		 /* rho_crit_cosmology [h_0^2 Mo kpc^{-3}] */
	};

const CONVERSION_FACTORS ConversionFactors = {

	1.02271216505,		 /* km_per_s_2_kpc_per_Gyr */
	0.977792221673		 /* kpc_per_Gyr_2_km_per_s */
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

/*
** Function for calculating E function used throughout cosmology
*/

double Ecosmo(double a, COSMOLOGICAL_PARAMETERS cp) {

	return sqrt(cp.OmegaM0*pow(a,-3) + cp.OmegaL0 + cp.OmegaK0*pow(a,-2) + cp.OmegaR0*pow(a,-4));
	}

/*
** Function for calculating unit vectors for spherical coordinates
*/

void calculate_unit_vectors_spherical(double pos[3], double erad[3], double ephi[3], double etheta[3]) {

	double dist;
	double cosphi, sinphi, costheta, sintheta;

	/*
	** Calculate cosphi & sinphi
	*/
	dist = sqrt(pos[0]*pos[0]+pos[1]*pos[1]);
	cosphi = pos[0]/dist;
	sinphi = pos[1]/dist;
	if ((pos[0] == 0) && (pos[1] == 0)) {
		cosphi = 1;
		sinphi = 0;
		}
	/*
	** Calculate costheta & sintheta
	*/
	dist = sqrt(pos[0]*pos[0]+pos[1]*pos[1]+pos[2]*pos[2]);
	costheta = pos[2]/dist;
	sintheta = sqrt(1-costheta*costheta);
	if (dist == 0) {
		cosphi = 0;
		sinphi = 0;
		costheta = 0;
		sintheta = 0;
		}
	/*
	** Calculate unit vectors
	*/
	erad[0] = sintheta*cosphi;
	erad[1] = sintheta*sinphi;
	erad[2] = costheta;
	ephi[0] = -sinphi;
	ephi[1] = cosphi;
	ephi[2] = 0;
	etheta[0] = -costheta*cosphi;
	etheta[1] = -costheta*sinphi;
	etheta[2] = sintheta;
	}

/*
** Function for calculating unit vectors for cylindrical coordinates
*/

void calculate_unit_vectors_cylindrical(double pos[3], double ezin[3], double erad[3], double ephi[3], double ez[3]) {

	int i;
	double dist, distv;
	double vec[3];

	/*
	** Normalise z component
	*/
	dist = 1.0/sqrt(ezin[0]*ezin[0]+ezin[1]*ezin[1]+ezin[2]*ezin[2]);
	for (i = 0; i < 3; i++) ez[i] = ezin[i]*dist;
	/*
	** Calculate radial component
	*/
	dist = pos[0]*ez[0]+pos[1]*ez[1]+pos[2]*ez[2];
	if (fabs(dist/sqrt(pos[0]*pos[0]+pos[1]*pos[1]+pos[2]*pos[2])) > 1-1e-15) {
		/*
		** Alingned with z component 
		** => asign random vector perpendicular to z component
		** and make sure that the radial component is never negative!
		*/
		vec[0] = rand();
		vec[1] = rand();
		vec[2] = rand();
		distv = vec[0]*ez[0]+vec[1]*ez[1]+vec[2]*ez[2];
		/* fprintf(stderr,"Close case\n"); */
		for (i = 0; i < 3; i++) erad[i] = vec[i]-distv*ez[i];
		if (pos[0]*erad[0]+pos[1]*erad[1]+pos[2]*erad[2] < 0) {
			/* fprintf(stderr,"Negative case\n"); */
			for (i = 0; i < 3; i++) erad[i] *= -1;
			}
		}
	else {
		for (i = 0; i < 3; i++) erad[i] = pos[i]-dist*ez[i];
		}
	/*
	** Normalise radial component
	*/
	dist = 1.0/sqrt(erad[0]*erad[0]+erad[1]*erad[1]+erad[2]*erad[2]);
	for (i = 0; i < 3; i++) erad[i] = erad[i]*dist;
	/*
	** Calculate phi component as a cross product e_z x e_rad
	*/
	ephi[0] = ez[1]*erad[2]-ez[2]*erad[1];
	ephi[1] = ez[2]*erad[0]-ez[0]*erad[2];
	ephi[2] = ez[0]*erad[1]-ez[1]*erad[0];
	}

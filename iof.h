/*
** iof.h
**
** written by Marcel Zemp
*/

#ifndef IOF_H
#define IOF_H

#ifdef __cplusplus
extern "C" {
#endif

#include <rpc/types.h>
#include <rpc/xdr.h>

#define kB 1.3806503e-23
#define PROTONMASS 1.6726216e-27

/*
** Tipsy definitions
*/

typedef struct tipsy_header {

    double time;
    int ntotal;
    int ndim;
    int ngas;
    int ndark;
    int nstar;
    } TIPSY_HEADER;

typedef struct gas_particle {

    float mass;
    float pos[3];
    float vel[3];
    float rho;
    float temp;
    float hsmooth;
    float metals;
    float phi;
    } GAS_PARTICLE;

typedef struct dark_particle {

    float mass;
    float pos[3];
    float vel[3];
    float eps;
    float phi;
    } DARK_PARTICLE;

typedef struct star_particle {

    float mass;
    float pos[3];
    float vel[3];
    float metals;
    float tform;
    float eps;
    float phi;
    } STAR_PARTICLE;

typedef struct gas_particle_dpp {

    float mass;
    double pos[3];
    float vel[3];
    float rho;
    float temp;
    float hsmooth;
    float metals;
    float phi;
    } GAS_PARTICLE_DPP;

typedef struct dark_particle_dpp {

    float mass;
    double pos[3];
    float vel[3];
    float eps;
    float phi;
    } DARK_PARTICLE_DPP;

typedef struct star_particle_dpp {

    float mass;
    double pos[3];
    float vel[3];
    float metals;
    float tform;
    float eps;
    float phi;
    } STAR_PARTICLE_DPP;

typedef struct tipsy_structure {

    TIPSY_HEADER *th;
    GAS_PARTICLE *gp;
    DARK_PARTICLE *dp;
    STAR_PARTICLE *sp;
    } TIPSY_STRUCTURE;

typedef struct tipsy_structure_dpp {

    TIPSY_HEADER *th;
    GAS_PARTICLE_DPP *gpdpp;
    DARK_PARTICLE_DPP *dpdpp;
    STAR_PARTICLE_DPP *spdpp;
    } TIPSY_STRUCTURE_DPP;

/*
** Gadget definitions from allvars.h
*/

typedef struct io_header {
    
    int npart[6];                        /*!< number of particles of each type in this file */
    double mass[6];                      /*!< mass of particles of each type. If 0, then the masses are explicitly
					   stored in the mass-block of the snapshot file, otherwise they are omitted */
    double time;                         /*!< time of snapshot file */
    double redshift;                     /*!< redshift of snapshot file */
    int flag_sfr;                        /*!< flags whether the simulation was including star formation */
    int flag_feedback;                   /*!< flags whether feedback was included (obsolete) */
    unsigned int npartTotal[6];          /*!< total number of particles of each type in this snapshot. This can be
					   different from npart if one is dealing with a multi-file snapshot. */
    int flag_cooling;                    /*!< flags whether cooling was included  */
    int num_files;                       /*!< number of files in multi-file snapshot */
    double BoxSize;                      /*!< box-size of simulation in case periodic boundaries were used */
    double Omega0;                       /*!< matter density in units of critical density */
    double OmegaLambda;                  /*!< cosmological constant parameter */
    double HubbleParam;                  /*!< Hubble parameter in units of 100 km/sec/Mpc */
    int flag_stellarage;                 /*!< flags whether the file contains formation times of star particles */
    int flag_metals;                     /*!< flags whether the file contains metallicity values for gas and star particles */
    unsigned int npartTotalHighWord[6];  /*!< High word of the total number of particles of each type */
    int  flag_entropy_instead_u;         /*!< flags that IC-file contains entropy instead of u */
    char fill[60];                       /*!< fills to 256 Bytes */
    } GADGET_HEADER;

/*
** ART header from io.h
*/

typedef struct {

    float aunin;
    float auni0;
    float amplt;
    float astep;
    int   istep;
    float partw;
    float tintg;
    float ekin;
    float ekin1;
    float ekin2;
    float au0;
    float aeu0;
    int   Nrow;
    int   Ngrid;
    int   Nspecies;
    int   Nseed;
    float OmM0;
    float OmL0;
    float h100;
    float Wp5;
    float OmK0;
    float OmB0;  
    float mass[10];
    unsigned int num[10];
    float zero1;
    float DelDC;
    float abox;   /* Scale factor in the box */
    float Hbox;   /* Hubble constant in the box */
    float zero2;
    float fill[75];
    } ART_HEADER;

/*
** Array definitions
*/

typedef struct array_header {
    
    int N[4];
    } ARRAY_HEADER;

typedef struct array_particle {
    
    int *ia;
    float *fa;
    double *da;
    } ARRAY_PARTICLE;

/*
** Sorting structures needed by read_gadget_binary
*/

typedef struct qsort_gp {

    int index;
    GAS_PARTICLE gp;
    } QSORT_GP;

typedef struct qsort_dp {

    int index;
    DARK_PARTICLE dp;
    } QSORT_DP;

typedef struct qsort_sp {

    int index;
    STAR_PARTICLE sp;
    } QSORT_SP;

/*
** Function definitions
*/

void copy_gp_to_gpdpp(const GAS_PARTICLE (*), GAS_PARTICLE_DPP (*));
void copy_dp_to_dpdpp(const DARK_PARTICLE (*), DARK_PARTICLE_DPP (*));
void copy_sp_to_spdpp(const STAR_PARTICLE (*), STAR_PARTICLE_DPP (*));
void copy_gpdpp_to_gp(const GAS_PARTICLE_DPP (*), GAS_PARTICLE (*));
void copy_dpdpp_to_dp(const DARK_PARTICLE_DPP (*), DARK_PARTICLE (*));
void copy_spdpp_to_sp(const STAR_PARTICLE_DPP (*), STAR_PARTICLE (*));

void read_tipsy_binary_header(FILE (*), TIPSY_HEADER (*));
void read_tipsy_binary_gas(FILE (*), GAS_PARTICLE (*));
void read_tipsy_binary_dark(FILE (*), DARK_PARTICLE (*));
void read_tipsy_binary_star(FILE (*), STAR_PARTICLE (*));
void read_tipsy_binary_gas_dpp(FILE (*), GAS_PARTICLE_DPP (*));
void read_tipsy_binary_dark_dpp(FILE (*), DARK_PARTICLE_DPP (*));
void read_tipsy_binary_star_dpp(FILE (*), STAR_PARTICLE_DPP (*));
void read_tipsy_binary(FILE (*), TIPSY_STRUCTURE (*));
void read_tipsy_binary_dpp(FILE (*), TIPSY_STRUCTURE_DPP (*));

void write_tipsy_binary_header(FILE (*), const TIPSY_HEADER (*));
void write_tipsy_binary_gas(FILE (*), const GAS_PARTICLE (*));
void write_tipsy_binary_dark(FILE (*), const DARK_PARTICLE (*));
void write_tipsy_binary_star(FILE (*), const STAR_PARTICLE (*));
void write_tipsy_binary_gas_dpp(FILE (*), const GAS_PARTICLE_DPP (*));
void write_tipsy_binary_dark_dpp(FILE (*), const DARK_PARTICLE_DPP (*));
void write_tipsy_binary_star_dpp(FILE (*), const STAR_PARTICLE_DPP (*));
void write_tipsy_binary(FILE (*), const TIPSY_STRUCTURE (*));
void write_tipsy_binary_dpp(FILE (*), const TIPSY_STRUCTURE_DPP (*));

void read_tipsy_standard_header(XDR (*), TIPSY_HEADER (*));
void read_tipsy_standard_gas(XDR (*), GAS_PARTICLE (*));
void read_tipsy_standard_dark(XDR (*), DARK_PARTICLE (*));
void read_tipsy_standard_star(XDR (*), STAR_PARTICLE (*));
void read_tipsy_standard_gas_dpp(XDR (*), GAS_PARTICLE_DPP (*));
void read_tipsy_standard_dark_dpp(XDR (*), DARK_PARTICLE_DPP (*));
void read_tipsy_standard_star_dpp(XDR (*), STAR_PARTICLE_DPP (*));
void read_tipsy_standard(FILE (*), TIPSY_STRUCTURE (*));
void read_tipsy_standard_dpp(FILE (*), TIPSY_STRUCTURE_DPP (*));

void write_tipsy_standard_header(XDR (*), TIPSY_HEADER (*));
void write_tipsy_standard_gas(XDR (*), GAS_PARTICLE (*));
void write_tipsy_standard_dark(XDR (*), DARK_PARTICLE (*));
void write_tipsy_standard_star(XDR (*), STAR_PARTICLE (*));
void write_tipsy_standard_gas_dpp(XDR (*), GAS_PARTICLE_DPP (*));
void write_tipsy_standard_dark_dpp(XDR (*), DARK_PARTICLE_DPP (*));
void write_tipsy_standard_star_dpp(XDR (*), STAR_PARTICLE_DPP (*));
void write_tipsy_standard(FILE (*), const TIPSY_STRUCTURE (*));
void write_tipsy_standard_dpp(FILE (*), const TIPSY_STRUCTURE_DPP (*));

void read_tipsy_ascii(FILE (*), TIPSY_STRUCTURE (*));
void read_tipsy_ascii_dpp(FILE (*), TIPSY_STRUCTURE_DPP (*));

void write_tipsy_ascii(FILE (*), const TIPSY_STRUCTURE (*));
void write_tipsy_ascii_dpp(FILE (*), const TIPSY_STRUCTURE_DPP (*));

void read_gadget_binary(FILE (*), TIPSY_STRUCTURE (*), double, double, double, double, double, double, double);
void write_gadget_binary(FILE (*), const TIPSY_STRUCTURE (*), double, double, double, double, double, double, double);

void allocate_array_particle(const ARRAY_HEADER (*), ARRAY_PARTICLE (*));
    
void read_array_header(XDR (*), ARRAY_HEADER (*));
void read_array_particle(XDR (*), const ARRAY_HEADER (*), ARRAY_PARTICLE (*));

void write_array_header(XDR (*), ARRAY_HEADER (*));
void write_array_particle(XDR (*), const ARRAY_HEADER (*), ARRAY_PARTICLE (*));

void reorder(void (*), size_t, size_t);

#ifdef __cplusplus
}
#endif

#endif /* IOF_H */

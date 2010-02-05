/*
** gadget.h
**
** Various reading and writing functions for gadget format
**
** Written by Marcel Zemp
*/

#ifndef IOF_GADGET_H
#define IOF_GADGET_H

/*
** Structures
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

typedef struct qsort_tgp {

    int index;
    TIPSY_GAS_PARTICLE tgp;
    } QSORT_TGP;

typedef struct qsort_tdp {

    int index;
    TIPSY_DARK_PARTICLE tdp;
    } QSORT_TDP;

typedef struct qsort_tsp {

    int index;
    TIPSY_STAR_PARTICLE tsp;
    } QSORT_TSP;

/*
** Reading and writing functions
*/

void read_gadget_nb(FILE *, TIPSY_STRUCTURE *, double, double, double, double, double, double, double);
void write_gadget_nb(FILE *, const TIPSY_STRUCTURE *, double, double, double, double, double, double, double);

#endif /* IOF_GADGET_H */

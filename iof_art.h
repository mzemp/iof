/*
** iof_art.h
**
** Various reading and writing functions for art format
**
** written by Marcel Zemp
*/

#ifndef IOF_ART_H
#define IOF_ART_H

#include <sfc.h>

/*
** Definitions
*/

#define ART_MAX_NUMBER_GAS_LEVELS 30
#define ART_MAX_NUMBER_DARK_LEVELS 10
#define ART_BANNER_LENGTH 45
#define ART_MAX_NUMBER_STAR_PROPERTIES 5
#define ART_MAX_NUMBER_GAS_BLOCKS 2

const double art_cell_delta[8][3];

/*
** Structures
*/

typedef struct art_header {

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

typedef struct art_coordinates {

    double r[3];
    double v[3];
    } ART_COORDINATES;

typedef struct art_star_properties {

    double mass;
    double initialmass;
    double tform;
    double metallicitySNII;
    double metallicitySNIa;
    } ART_STAR_PROPERTIES;

typedef struct art_data {

    /*
    ** from general header
    */
    char Banner[ART_BANNER_LENGTH];
    ART_HEADER ah;
    /*
    ** from star properties file
    */
    double totalstellarmass, totalstellarinitialmass;
    /*
    ** from gas file
    */
    double refinementvolumemin[3], refinementvolumemax[3];
    double starformationvolumemin[3], starformationvolumemax[3];
    /*
    ** derived stuff to get data better organised
    */
    int massfromdata;
    int doswap;
    int gascontained, darkcontained, starcontained;
    int Lmingas, Lmaxgas, Nlevelgas;
    int Lmindark, Lmaxdark, Nleveldark;
    int Ndim;
    int Nparticleperrecord;
    int Nrecord;
    int Ndarklevel[ART_MAX_NUMBER_DARK_LEVELS];
    int Nstarproperties;
    int Nhydroproperties, Notherproperties;
    int Nrtchemspecies, Nchemspecies, Nrtdiskvars;
    long Ngas, Ndark, Nstar;
    long Ncell[ART_MAX_NUMBER_GAS_LEVELS], Ncellrefined[ART_MAX_NUMBER_GAS_LEVELS];
    double shift;
    double rootcelllength;
    double toplevelmassdark, toplevelsoftdark, refinementstepdark;
    double massdark[ART_MAX_NUMBER_DARK_LEVELS];
    double softdark[ART_MAX_NUMBER_DARK_LEVELS];
    /*
    ** ART preprocessor flags
    */
    int GRAVITY, HYDRO, STARFORM, ADVECT_SPECIES, ENRICH, ENRICH_SNIa, RADIATIVE_TRANSFER, ELECTRON_ION_NONEQUILIBRIUM;
    /*
    ** SFC info
    */
    SFC_INFO sfci;
    /*
    ** ART files
    */
    char HeaderFileName[256], CoordinatesDataFileName[256], StarPropertiesFileName[256], GasFileName[256];
    FILE *HeaderFile, *CoordinatesDataFile, *StarPropertiesFile[ART_MAX_NUMBER_STAR_PROPERTIES], *GasFile[ART_MAX_NUMBER_GAS_BLOCKS];
    } ART_DATA;

/*
** Reading and writing functions
*/

void read_art_nb_general_header(ART_DATA *);
void read_art_nb_star_header(ART_DATA *, int index);
void read_art_nb_gas_header(ART_DATA *, int index);

#endif /* IOF_ART_H */

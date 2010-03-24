/*
** iof_art.h
**
** Various handling functions for art format
**
** Written by Marcel Zemp
*/

#ifndef IOF_ART_H
#define IOF_ART_H

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
    float magic1;
    float DelDC;
    float abox;   /* Scale factor in the box */
    float Hbox;   /* Hubble constant in the box */
    float magic2;
    float fill[75];
    } ART_HEADER;

typedef struct art_gas_properties {

    double gas_density;
    double gas_energy;
    double momentum[3];
    double pressure;
    double gamma;
    double internal_energy;
    double electron_internal_energy;
    double potential;
    double potential_hydro;
    double HI_density;
    double HII_density;
    double HeI_density;
    double HeII_density;
    double HeIII_density;
    double H2_density;
    double metal_density_SNII;
    double metal_density_SNIa;
    } ART_GAS_PROPERTIES;

typedef struct art_star_properties {

    double mass;
    double initial_mass;
    double t_form;
    double metallicity_SNII;
    double metallicity_SNIa;
    } ART_STAR_PROPERTIES;

typedef struct art_coordinates {

    double r[3];
    double v[3];
    } ART_COORDINATES;

typedef struct art_sfc_info {

    int nDim;
    int num_grid;
    int sfc_order;
    int nBitsPerDim;
    int nBits;
    int max_sfc_index;
    } ART_SFC_INFO;

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
    long int Ngas, Ndark, Nstar;
    long int Ncell[ART_MAX_NUMBER_GAS_LEVELS], Ncellrefined[ART_MAX_NUMBER_GAS_LEVELS];
    double shift;
    double rootcelllength;
    double toplevelmassdark, toplevelsoftdark, refinementstepdark;
    double massdark[ART_MAX_NUMBER_DARK_LEVELS];
    double softdark[ART_MAX_NUMBER_DARK_LEVELS];
    /*
    ** ART preprocessor flags
    */
    int GRAVITY, HYDRO, ADVECT_SPECIES, STARFORM, ENRICH, ENRICH_SNIa, RADIATIVE_TRANSFER, ELECTRON_ION_NONEQUILIBRIUM;
    /*
    ** SFC info
    */
    ART_SFC_INFO asfci;
    /*
    ** ART files
    */
    char HeaderFileName[256], CoordinatesDataFileName[256], StarPropertiesFileName[256], GasFileName[256];
    FILE *HeaderFile, *CoordinatesDataFile, *StarPropertiesFile[ART_MAX_NUMBER_STAR_PROPERTIES], *GasFile[ART_MAX_NUMBER_GAS_BLOCKS];
    } ART_DATA;

/*
** Functions
*/

void set_default_values_art_data(ART_DATA *);
void prepare_art_data(ART_DATA *);

void read_art_nb_general_header(ART_DATA *);
void read_art_nb_gas_header(ART_DATA *, int index);
void read_art_nb_gas_header_level(ART_DATA *, int, int **);
void read_art_nb_star_header(ART_DATA *, int index);

void read_art_nb_coordinates_record(ART_DATA, ART_COORDINATES (*));
void read_art_nb_gas_properties(ART_DATA, ART_GAS_PROPERTIES *); 
void read_art_nb_star_properties(ART_DATA, ART_STAR_PROPERTIES *);

void move_art_nb_gas_filepositions_level_begin(ART_DATA, int);
void move_art_nb_gas_filepositions_level_end(ART_DATA, int);
void move_art_nb_star_filepositions_begin(ART_DATA);
void move_art_nb_star_filepositions_end(ART_DATA);

#endif /* IOF_ART_H */

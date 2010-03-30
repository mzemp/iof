/*
** iof_art.c
**
** Various handling functions for art format
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
#include "iof_art.h"

/*
** Definitions
*/

const double art_cell_delta[8][3] = {
	{ -0.5, -0.5, -0.5 }, {  0.5, -0.5, -0.5 }, { -0.5,  0.5, -0.5 }, {  0.5,  0.5, -0.5 }, 
	{ -0.5, -0.5,  0.5 }, {  0.5, -0.5,  0.5 }, { -0.5,  0.5,  0.5 }, {  0.5,  0.5,  0.5 }
    };

/*
** Functions
*/

void set_default_values_art_data(ART_DATA *ad) {

    int i;

    ad->doswap = 0;
    ad->gascontained = 0;
    ad->darkcontained = 0;
    ad->starcontained = 0;
    ad->Nparticleperrecord = 0;
    ad->Nrecord = 0;
    ad->Ndim = 3;
    ad->Ngas = 0;
    ad->Ndark = 0;
    ad->Nstar = 0;
    ad->refinementstepdark = 2;
    ad->toplevelmassdark = -1;
    ad->toplevelsoftdark = -1;
    ad->Lmingas = 0;
    ad->Lmaxgas = 0;
    ad->Nlevelgas = 0;
    ad->Lmindark = 0;
    ad->Lmaxdark = 0;
    ad->Nleveldark = 0;
    ad->shift = 1;
    ad->rootcelllength = 1;
    for (i = 0; i < ART_MAX_NUMBER_DARK_LEVELS; i++) {
	ad->Ndarklevel[i] = 0;
	ad->massdark[i] = 0;
	ad->softdark[i] = 0;
	}
    for (i = 0; i < ART_MAX_NUMBER_GAS_LEVELS; i++) {
	ad->Ncell[i] = 0;
	ad->Ncellrefined[i] = 0;
	}
    ad->Nhydroproperties = 0;
    ad->Notherproperties = 0;
    ad->Nstarproperties = 0;
    ad->Nrtchemspecies = 0;
    ad->Nchemspecies = 0;
    ad->Nrtdiskvars = 0;
    ad->GRAVITY = 1;
    ad->HYDRO = 1;
    ad->ADVECT_SPECIES = 1;
    ad->STARFORM = 1;
    ad->ENRICH = 1;
    ad->ENRICH_SNIa = 1;
    ad->RADIATIVE_TRANSFER = 1;
    ad->ELECTRON_ION_NONEQUILIBRIUM = 0;
    ad->asfci.nDim = 0;
    ad->asfci.num_grid = 0;
    ad->asfci.sfc_order = -1;
    ad->asfci.nBitsPerDim = 0;
    ad->asfci.nBits = 0;
    ad->asfci.max_sfc_index = 0;
    ad->HeaderFile = NULL;
    ad->CoordinatesDataFile = NULL;
    for (i = 0; i < ART_MAX_NUMBER_STAR_PROPERTIES; i++) {
	ad->StarPropertiesFile[i] = NULL;
	}
    for (i = 0; i < ART_MAX_NUMBER_GAS_BLOCKS; i++) {
	ad->GasFile[i] = NULL;
	}
    }

void prepare_art_data(ART_DATA *ad) {

    int i;

    /*
    ** Derive number of properties from flags
    */
    ad->Nhydroproperties = 0;
    ad->Notherproperties = 0;
    ad->Nstarproperties = 0;
    if (ad->gascontained) {
	if (ad->HYDRO) {
	    if (ad->ADVECT_SPECIES) {
		if (ad->RADIATIVE_TRANSFER) ad->Nrtchemspecies = 6;
		else ad->Nrtchemspecies = 0;
		if (ad->ENRICH) {
		    if (ad->ENRICH_SNIa) ad->Nchemspecies = ad->Nrtchemspecies + 2;
		    else ad->Nchemspecies = ad->Nrtchemspecies + 1;
		    }
		else ad->Nchemspecies = ad->Nrtchemspecies;
		}
	    else ad->Nchemspecies = 0;
	    if (ad->ELECTRON_ION_NONEQUILIBRIUM) ad->Nhydroproperties = 6 + ad->Ndim + ad->Nchemspecies;
	    else ad->Nhydroproperties = 5 + ad->Ndim + ad->Nchemspecies;
	    }
	if (ad->RADIATIVE_TRANSFER) ad->Nrtdiskvars = 6;
	else ad->Nrtdiskvars = 0;
	if (ad->GRAVITY) ad->Notherproperties++;
	if (ad->HYDRO) ad->Notherproperties++;
	if (ad->RADIATIVE_TRANSFER) ad->Notherproperties += ad->Nrtdiskvars;
	}
    if (ad->starcontained) {
	ad->Nstarproperties = 3;
	if (ad->ENRICH) {
	    if (ad->ENRICH_SNIa) ad->Nstarproperties = ad->Nstarproperties + 2;
	    else ad->Nstarproperties++;
	    }
	}
    /*
    ** Read all headers
    */
    ad->HeaderFile = fopen(ad->HeaderFileName,"r");
    assert(ad->HeaderFile != NULL);
    read_art_nb_general_header(ad);
    fclose(ad->HeaderFile);
    if (ad->gascontained) {
	ad->GasFile[0] = fopen(ad->GasFileName,"r");
	assert(ad->GasFile[0] != NULL);
	read_art_nb_gas_header(ad,0);
	if (ad->GRAVITY || ad->RADIATIVE_TRANSFER) {
	    ad->GasFile[1] = fopen(ad->GasFileName,"r");
	    assert(ad->GasFile[1] != NULL);
	    read_art_nb_gas_header(ad,1);
	    }
	}
    if (ad->darkcontained || ad->starcontained) {
	ad->CoordinatesDataFile = fopen(ad->CoordinatesDataFileName,"r");
	assert(ad->CoordinatesDataFile != NULL);
	}
    if (ad->starcontained) {
	for (i = 0; i < ad->Nstarproperties; i++) {
	    ad->StarPropertiesFile[i] = fopen(ad->StarPropertiesFileName,"r");
	    assert(ad->StarPropertiesFile[i] != NULL);
	    read_art_nb_star_header(ad,i);
	    }
	}
    }

void read_art_nb_general_header(ART_DATA *ad) {

    int i, index;
    int header, trailer;

    ad->doswap = 0;
    assert(fread(&header,sizeof(int),1,ad->HeaderFile) == 1);
    if (header != sizeof(ad->Banner)+sizeof(ART_HEADER)) {
	ad->doswap = 1;
	reorder(&header,sizeof(int),1);
	}
    assert(header == sizeof(ad->Banner)+sizeof(ART_HEADER));
    assert(fread(ad->Banner,sizeof(char),sizeof(ad->Banner),ad->HeaderFile) == sizeof(ad->Banner));
    assert(fread(&ad->ah,sizeof(ART_HEADER),1,ad->HeaderFile) == 1);
    if (ad->doswap) reorder(&ad->ah,4,sizeof(ART_HEADER)/4);
    assert(fread(&trailer,sizeof(int),1,ad->HeaderFile) == 1);
    if (ad->doswap) reorder(&trailer,sizeof(int),1);
    assert(header == trailer);
    /*
    ** Set some derived quantities
    */
    ad->asfci.nDim = ad->Ndim;
    ad->asfci.num_grid = ad->ah.Ngrid;
    ad->Ngas = 0;
    if (ad->darkcontained) {
	if (ad->starcontained == 1) {
	    ad->Lmaxdark = ad->ah.Nspecies-2;
	    ad->Ndark = ad->ah.num[ad->Lmaxdark];
	    ad->Nstar = ad->ah.num[ad->ah.Nspecies-1] - ad->Ndark;
	    }
	else {
	    ad->Lmaxdark = ad->ah.Nspecies-1;
	    ad->Ndark = ad->ah.num[ad->Lmaxdark];
	    ad->Nstar = 0;
	    }
	ad->Lmindark = 0;
	ad->Nleveldark = ad->Lmaxdark-ad->Lmindark+1;
	ad->Nparticleperrecord = ad->ah.Nrow*ad->ah.Nrow;
	ad->Nrecord = (ad->Ndark+ad->Nstar+ad->Nparticleperrecord-1)/ad->Nparticleperrecord;
	for (i = ad->Lmindark; i <= ad->Lmaxdark; i++) {
	    index = ad->Lmaxdark-i;
	    if (index == 0) {
		ad->Ndarklevel[i] = ad->ah.num[index];
		}
	    else {
		ad->Ndarklevel[i] = ad->ah.num[index]-ad->ah.num[index-1];
		}
	    }
	}
    }

void read_art_nb_gas_header(ART_DATA *ad, int index) {

    int i;
    int header, trailer;
    int idummy;
    float fdummy;
    double ddummy;
    char cdummy[256];

    /*
    ** jobname
    */
    assert(fread(&header,sizeof(int),1,ad->GasFile[index]) == 1);
    if (ad->doswap) reorder(&header,sizeof(int),1);
    assert(fread(&cdummy,sizeof(char),256,ad->GasFile[index]) == 256);
    assert(fread(&trailer,sizeof(int),1,ad->GasFile[index]) == 1);
    if (ad->doswap) reorder(&trailer,sizeof(int),1);
    assert(header == trailer);
    /* 
    ** istep, t, dt, adum, ainit 
    */
    assert(fread(&header,sizeof(int),1,ad->GasFile[index]) == 1);
    if (ad->doswap) reorder(&header,sizeof(int),1);
    assert(fread(&idummy,sizeof(int),1,ad->GasFile[index]) == 1);
    assert(fread(&ddummy,sizeof(double),1,ad->GasFile[index]) == 1);
    assert(fread(&ddummy,sizeof(double),1,ad->GasFile[index]) == 1);
    assert(fread(&fdummy,sizeof(float),1,ad->GasFile[index]) == 1);
    assert(fread(&fdummy,sizeof(float),1,ad->GasFile[index]) == 1);
    assert(fread(&trailer,sizeof(int),1,ad->GasFile[index]) == 1);
    if (ad->doswap) reorder(&trailer,sizeof(int),1);
    assert(header == trailer);
    /* 
    ** boxh, Om0, Oml0, Omb0, hubble
    */
    assert(fread(&header,sizeof(int),1,ad->GasFile[index]) == 1);
    if (ad->doswap) reorder(&header,sizeof(int),1);
    assert(fread(&fdummy,sizeof(float),1,ad->GasFile[index]) == 1);
    assert(fread(&fdummy,sizeof(float),1,ad->GasFile[index]) == 1);
    assert(fread(&fdummy,sizeof(float),1,ad->GasFile[index]) == 1);
    assert(fread(&fdummy,sizeof(float),1,ad->GasFile[index]) == 1);
    assert(fread(&fdummy,sizeof(float),1,ad->GasFile[index]) == 1);
    assert(fread(&trailer,sizeof(int),1,ad->GasFile[index]) == 1);
    if (ad->doswap) reorder(&trailer,sizeof(int),1);
    assert(header == trailer);
    /*
    ** Nextra (some old crap - should be zero!)
    */
    assert(fread(&header,sizeof(int),1,ad->GasFile[index]) == 1);
    if (ad->doswap) reorder(&header,sizeof(int),1);
    assert(fread(&idummy,sizeof(int),1,ad->GasFile[index]) == 1);
    if (ad->doswap) reorder(&idummy,sizeof(int),1);
    assert(idummy == 0);
    assert(fread(&trailer,sizeof(int),1,ad->GasFile[index]) == 1);
    if (ad->doswap) reorder(&trailer,sizeof(int),1);
    assert(header == trailer);
    /*
    ** extra
    */
    assert(fread(&header,sizeof(int),1,ad->GasFile[index]) == 1);
    if (ad->doswap) reorder(&header,sizeof(int),1);
    assert(fread(&trailer,sizeof(int),1,ad->GasFile[index]) == 1);
    if (ad->doswap) reorder(&trailer,sizeof(int),1);
    assert(header == trailer);
    /*
    ** lextra
    */
    assert(fread(&header,sizeof(int),1,ad->GasFile[index]) == 1);
    if (ad->doswap) reorder(&header,sizeof(int),1);
    assert(fread(&trailer,sizeof(int),1,ad->GasFile[index]) == 1);
    if (ad->doswap) reorder(&trailer,sizeof(int),1);
    assert(header == trailer);
    /* 
    ** minlevel, maxlevel
    */
    assert(fread(&header,sizeof(int),1,ad->GasFile[index]) == 1);
    if (ad->doswap) reorder(&header,sizeof(int),1);
    assert(fread(&idummy,sizeof(int),1,ad->GasFile[index]) == 1);
    if (ad->doswap) reorder(&idummy,sizeof(int),1);
    ad->Lmingas = idummy;
    assert(fread(&idummy,sizeof(int),1,ad->GasFile[index]) == 1);
    if (ad->doswap) reorder(&idummy,sizeof(int),1);
    ad->Lmaxgas = idummy;
    assert(fread(&trailer,sizeof(int),1,ad->GasFile[index]) == 1);
    if (ad->doswap) reorder(&trailer,sizeof(int),1);
    assert(header == trailer);
    ad->Nlevelgas = ad->Lmaxgas-ad->Lmingas+1;
    assert(ad->Lmingas >= 0);
    assert(ad->Lmaxgas >= ad->Lmingas);
    /*
    ** tl
    */
    assert(fread(&header,sizeof(int),1,ad->GasFile[index]) == 1);
    if (ad->doswap) reorder(&header,sizeof(int),1);
    for (i = 0; i < ad->Nlevelgas; i++) {
	assert(fread(&ddummy,sizeof(double),1,ad->GasFile[index]) == 1);
	}
    assert(fread(&trailer,sizeof(int),1,ad->GasFile[index]) == 1);
    if (ad->doswap) reorder(&trailer,sizeof(int),1);
    assert(header == trailer);
    /* 
    ** dtl
    */
    assert(fread(&header,sizeof(int),1,ad->GasFile[index]) == 1);
    if (ad->doswap) reorder(&header,sizeof(int),1);
    for (i = 0; i < ad->Nlevelgas; i++) {
	assert(fread(&ddummy,sizeof(double),1,ad->GasFile[index]) == 1);
	}
    assert(fread(&trailer,sizeof(int),1,ad->GasFile[index]) == 1);
    if (ad->doswap) reorder(&trailer,sizeof(int),1);
    assert(header == trailer);
    /*
    ** tl_old
    */
    assert(fread(&header,sizeof(int),1,ad->GasFile[index]) == 1);
    if (ad->doswap) reorder(&header,sizeof(int),1);
    for (i = 0; i < ad->Nlevelgas; i++) {
	assert(fread(&ddummy,sizeof(double),1,ad->GasFile[index]) == 1);
	}
    assert(fread(&trailer,sizeof(int),1,ad->GasFile[index]) == 1);
    if (ad->doswap) reorder(&trailer,sizeof(int),1);
    assert(header == trailer);
    /* 
    ** dtl_old
    */
    assert(fread(&header,sizeof(int),1,ad->GasFile[index]) == 1);
    if (ad->doswap) reorder(&header,sizeof(int),1);
    for (i = 0; i < ad->Nlevelgas; i++) {
	assert(fread(&ddummy,sizeof(double),1,ad->GasFile[index]) == 1);
	}
    assert(fread(&trailer,sizeof(int),1,ad->GasFile[index]) == 1);
    if (ad->doswap) reorder(&trailer,sizeof(int),1);
    assert(header == trailer);
    /* 
    ** iSO (sweep direction for flux solver)
    */
    if (ad->HYDRO) {
	assert(fread(&header,sizeof(int),1,ad->GasFile[index]) == 1);
	if (ad->doswap) reorder(&header,sizeof(int),1);
	for (i = 0; i < ad->Nlevelgas; i++) {
	    assert(fread(&idummy,sizeof(int),1,ad->GasFile[index]) == 1);
	    }
	assert(fread(&trailer,sizeof(int),1,ad->GasFile[index]) == 1);
	if (ad->doswap) reorder(&trailer,sizeof(int),1);
	assert(header == trailer);
	}
    /*
    ** space filling curve (SFC) order
    */
    assert(fread(&header,sizeof(int),1,ad->GasFile[index]) == 1);
    if (ad->doswap) reorder(&header,sizeof(int),1);
    assert(fread(&idummy,sizeof(int),1,ad->GasFile[index]) == 1);
    if (ad->doswap) reorder(&idummy,sizeof(int),1);
    ad->asfci.sfc_order = idummy;
    assert(fread(&trailer,sizeof(int),1,ad->GasFile[index]) == 1);
    if (ad->doswap) reorder(&trailer,sizeof(int),1);
    assert(header == trailer);
    /* 
    ** refinement volume 
    */
    assert(fread(&header,sizeof(int),1,ad->GasFile[index]) == 1);
    if (ad->doswap) reorder(&header,sizeof(int),1);
    for (i = 0; i < ad->Ndim; i++) {
	assert(fread(&fdummy,sizeof(float),1,ad->GasFile[index]) == 1);
	if (ad->doswap) reorder(&fdummy,sizeof(float),1);
	ad->refinementvolumemin[i] = fdummy;
	}
    for (i = 0; i < ad->Ndim; i++) {
	assert(fread(&fdummy,sizeof(float),1,ad->GasFile[index]) == 1);
	if (ad->doswap) reorder(&fdummy,sizeof(float),1);
	ad->refinementvolumemax[i] = fdummy;
	}
    assert(fread(&trailer,sizeof(int),1,ad->GasFile[index]) == 1);
    if (ad->doswap) reorder(&trailer,sizeof(int),1);
    assert(header == trailer);
    /* 
    ** star formation volume 
    */
    if (ad->STARFORM) {
	assert(fread(&header,sizeof(int),1,ad->GasFile[index]) == 1);
	if (ad->doswap) reorder(&header,sizeof(int),1);
	for (i = 0; i < ad->Ndim; i++) {
	    assert(fread(&fdummy,sizeof(float),1,ad->GasFile[index]) == 1);
	    if (ad->doswap) reorder(&fdummy,sizeof(float),1);
	    ad->starformationvolumemin[i] = fdummy;
	    }
	for (i = 0; i < ad->Ndim; i++) {
	    assert(fread(&fdummy,sizeof(float),1,ad->GasFile[index]) == 1);
	    if (ad->doswap) reorder(&fdummy,sizeof(float),1);
	    ad->starformationvolumemax[i] = fdummy;
	    }
	assert(fread(&trailer,sizeof(int),1,ad->GasFile[index]) == 1);
	if (ad->doswap) reorder(&trailer,sizeof(int),1);
	assert(header == trailer);
	}
    else {
	for (i = 0; i < ad->Ndim; i++) {
	    ad->starformationvolumemin[i] = 0;
	    ad->starformationvolumemax[i] = 0;
	    }
	}
    }

void read_art_nb_gas_header_level(ART_DATA *ad, int level, int **cellrefinedin) {

    int header, trailer;
    int i;
    int idummy;
    int *cellrefined = NULL;
    long int lidummy;

    assert(ad->Lmingas == 0);
    cellrefined = *cellrefinedin;
    /*
    ** Number of cells
    */
    assert(fread(&header,sizeof(int),1,ad->GasFile[0]) == 1);
    if (ad->doswap) reorder(&header,sizeof(int),1);
    if (level == ad->Lmingas) {
	assert(fread(&idummy,sizeof(int),1,ad->GasFile[0]) == 1);
	if (ad->doswap) reorder(&idummy,sizeof(int),1);
	ad->Ncell[0] = idummy;
	assert(ad->Ncell[0] == ad->ah.Ngrid*ad->ah.Ngrid*ad->ah.Ngrid);
	}
    else {
	assert(fread(&lidummy,sizeof(long int),1,ad->GasFile[0]) == 1);
	if (ad->doswap) reorder(&lidummy,sizeof(long int),1);
	ad->Ncell[level] = lidummy;
	}
    assert(fread(&trailer,sizeof(int),1,ad->GasFile[0]) == 1);
    if (ad->doswap) reorder(&trailer,sizeof(int),1);
    assert(header == trailer);
    /*
    ** Cellrefined array
    */
    cellrefined = realloc(cellrefined,ad->Ncell[level]*sizeof(int));
    assert(cellrefined != NULL);
    assert(fread(&header,sizeof(int),1,ad->GasFile[0]) == 1);
    if (ad->doswap) reorder(&header,sizeof(int),1);
    assert(fread(cellrefined,sizeof(int),ad->Ncell[level],ad->GasFile[0]) == ad->Ncell[level]);
    if (ad->doswap) reorder(cellrefined,sizeof(int),ad->Ncell[level]);
    assert(fread(&trailer,sizeof(int),1,ad->GasFile[0]) == 1);
    if (ad->doswap) reorder(&trailer,sizeof(int),1);
    assert(header == trailer);
    if (level > ad->Lmingas) assert(ad->Ncell[level] == ad->Ncellrefined[level-1]*8);
    if (level == ad->Lmingas) {
	ad->Ngas = 0;
	for (i = 0; i < ad->Ncell[level]; i++) {
	    ad->Ngas += cellrefined[i];
	    if (cellrefined[i] > 1) {
		ad->Ncellrefined[level]++;
		cellrefined[i] = 1;
		}
	    else if (cellrefined[i] == 1) {
		cellrefined[i] = 0;
		}
	    }
	}
    else {
	for (i = 0; i < ad->Ncell[level]; i++) {
	    if (cellrefined[i] == 1) ad->Ncellrefined[level]++;
	    }
	}
    /*
    ** Get second file ready if needed
    */
    if (ad->GRAVITY || ad->RADIATIVE_TRANSFER) {
	if (level == ad->Lmingas) lidummy = 4*sizeof(int) + sizeof(int) + ad->Ncell[level]*sizeof(int);
	else lidummy = 4*sizeof(int) + sizeof(long int) + ad->Ncell[level]*sizeof(int);
	assert(fseek(ad->GasFile[1],lidummy,SEEK_CUR) == 0);
	}
    *cellrefinedin = cellrefined;
    }

void read_art_nb_star_header(ART_DATA *ad, int index) {

    int header, trailer;
    int idummy;
    double ddummy;

    /*
    ** st, sa
    */
    assert(fread(&header,sizeof(int),1,ad->StarPropertiesFile[index]) == 1);
    if (ad->doswap) reorder(&header,sizeof(int),1);
    assert(fread(&ddummy,sizeof(double),1,ad->StarPropertiesFile[index]) == 1);
    assert(fread(&ddummy,sizeof(double),1,ad->StarPropertiesFile[index]) == 1);
    assert(fread(&trailer,sizeof(int),1,ad->StarPropertiesFile[index]) == 1);
    if (ad->doswap) reorder(&trailer,sizeof(int),1);
    assert(header == trailer);
    /*
    ** Nstar
    */
    assert(fread(&header,sizeof(int),1,ad->StarPropertiesFile[index]) == 1);
    if (ad->doswap) reorder(&header,sizeof(int),1);
    assert(fread(&idummy,sizeof(int),1,ad->StarPropertiesFile[index]) == 1);
    if (ad->doswap) reorder(&idummy,sizeof(int),1);
    assert(idummy == ad->Nstar);
    assert(fread(&trailer,sizeof(int),1,ad->StarPropertiesFile[index]) == 1);
    if (ad->doswap) reorder(&trailer,sizeof(int),1);
    assert(header == trailer);
    /*
    ** total stellar mass, total stellar initial mass
    */
    assert(fread(&header,sizeof(int),1,ad->StarPropertiesFile[index]) == 1);
    if (ad->doswap) reorder(&header,sizeof(int),1);
    assert(fread(&ddummy,sizeof(double),1,ad->StarPropertiesFile[index]) == 1);
    if (ad->doswap) reorder(&ddummy,sizeof(double),1);
    ad->totalstellarmass = ddummy;
    assert(fread(&ddummy,sizeof(double),1,ad->StarPropertiesFile[index]) == 1);
    if (ad->doswap) reorder(&ddummy,sizeof(double),1);
    ad->totalstellarinitialmass = ddummy;
    assert(fread(&trailer,sizeof(int),1,ad->StarPropertiesFile[index]) == 1);
    if (ad->doswap) reorder(&trailer,sizeof(int),1);
    assert(header == trailer);
    }

void read_art_nb_coordinates_record(ART_DATA ad, ART_COORDINATES *coordinates) {

    int i;
    float fdummy;

    assert(ad.CoordinatesDataFile != NULL);
    for (i = 0; i < ad.Nparticleperrecord; i++) {
	assert(fread(&fdummy,sizeof(float),1,ad.CoordinatesDataFile) == 1);
	if (ad.doswap) reorder(&fdummy,sizeof(float),1);
	coordinates[i].r[0] = fdummy;
	}
    for (i = 0; i < ad.Nparticleperrecord; i++) {
	assert(fread(&fdummy,sizeof(float),1,ad.CoordinatesDataFile) == 1);
	if (ad.doswap) reorder(&fdummy,sizeof(float),1);
	coordinates[i].r[1] = fdummy;
	}
    for (i = 0; i < ad.Nparticleperrecord; i++) {
	assert(fread(&fdummy,sizeof(float),1,ad.CoordinatesDataFile) == 1);
	if (ad.doswap) reorder(&fdummy,sizeof(float),1);
	coordinates[i].r[2] = fdummy;
	}
    for (i = 0; i < ad.Nparticleperrecord; i++) {
	assert(fread(&fdummy,sizeof(float),1,ad.CoordinatesDataFile) == 1);
	if (ad.doswap) reorder(&fdummy,sizeof(float),1);
	coordinates[i].v[0] = fdummy;
	}
    for (i = 0; i < ad.Nparticleperrecord; i++) {
	assert(fread(&fdummy,sizeof(float),1,ad.CoordinatesDataFile) == 1);
	if (ad.doswap) reorder(&fdummy,sizeof(float),1);
	coordinates[i].v[1] = fdummy;
	}
    for (i = 0; i < ad.Nparticleperrecord; i++) {
	assert(fread(&fdummy,sizeof(float),1,ad.CoordinatesDataFile) == 1);
	if (ad.doswap) reorder(&fdummy,sizeof(float),1);
	coordinates[i].v[2] = fdummy;
	}
    }

void read_art_nb_gas_properties(ART_DATA ad, ART_GAS_PROPERTIES *agp) {

    int index;
    float *cellhydroproperties, *cellotherproperties;

    cellhydroproperties = malloc(ad.Nhydroproperties*sizeof(float));
    assert(cellhydroproperties != NULL);
    cellotherproperties = malloc(ad.Notherproperties*sizeof(float));
    assert(cellotherproperties != NULL);

    assert(fread(&cellhydroproperties,sizeof(float),ad.Nhydroproperties,ad.GasFile[0]) == ad.Nhydroproperties);
    if (ad.doswap) reorder(&cellhydroproperties,sizeof(float),ad.Nhydroproperties);
    if (ad.GRAVITY || ad.RADIATIVE_TRANSFER) {
	assert(fread(&cellotherproperties,sizeof(float),ad.Notherproperties,ad.GasFile[1]) == ad.Notherproperties);
	if (ad.doswap) reorder(&cellotherproperties,sizeof(float),ad.Notherproperties);
	}

    agp->gas_density = cellhydroproperties[0];
    agp->gas_energy = cellhydroproperties[1];
    agp->momentum[0] = cellhydroproperties[2];
    agp->momentum[1] = cellhydroproperties[3];
    agp->momentum[2] = cellhydroproperties[4];
    agp->pressure = cellhydroproperties[5];
    agp->gamma = cellhydroproperties[6];
    agp->internal_energy = cellhydroproperties[7];
    agp->electron_internal_energy = 0;
    if (ad.ELECTRON_ION_NONEQUILIBRIUM) agp->electron_internal_energy = cellhydroproperties[8];
    agp->HI_density = 0;
    agp->HII_density = 0;
    agp->HeI_density = 0;
    agp->HeII_density = 0;
    agp->HeIII_density = 0;
    agp->H2_density = 0;
    if (ad.Nrtchemspecies > 0) {
	assert(ad.Nrtchemspecies == 6);
	if (ad.ELECTRON_ION_NONEQUILIBRIUM) index = 9;
	else index = 8;
	/*
	** What is written out for species j is X_j * gas_density where X_j = N_j / N_Nucleons
	** => need to multiply with number of nucleons of species j to get the mass density
	*/
	agp->HI_density    = cellhydroproperties[index];
	agp->HII_density   = cellhydroproperties[index+1];
	agp->HeI_density   = cellhydroproperties[index+2]*4;
	agp->HeII_density  = cellhydroproperties[index+3]*4;
	agp->HeIII_density = cellhydroproperties[index+4]*4;
	agp->H2_density    = cellhydroproperties[index+5]*2;
	}
    agp->metal_density_SNII = 0;
    agp->metal_density_SNIa = 0;
    if (ad.ENRICH) {
	if (ad.ENRICH_SNIa) {
	    index = ad.Nhydroproperties-2;
	    agp->metal_density_SNII = cellhydroproperties[index];
	    agp->metal_density_SNIa = cellhydroproperties[index+1];
	    }
	else {
	    index = ad.Nhydroproperties-1;
	    agp->metal_density_SNII = cellhydroproperties[index];
	    }
	}
    agp->potential = 0;
    agp->potential_hydro = 0;
    if (ad.GRAVITY || ad.RADIATIVE_TRANSFER) {
	if (ad.GRAVITY) agp->potential = cellotherproperties[0];
	if (ad.HYDRO) agp->potential_hydro = cellotherproperties[1];
	}
    free(cellhydroproperties);
    free(cellotherproperties);
    }

void read_art_nb_star_properties(ART_DATA ad, ART_STAR_PROPERTIES *asp) {

    int i;
    float fdummy;

    asp->mass = 0;
    asp->initial_mass = 0;
    asp->t_form = 0;
    asp->metallicity_SNII = 0;
    asp->metallicity_SNIa = 0;
    for (i = 0; i < ad.Nstarproperties; i++) {
	assert(fread(&fdummy,sizeof(float),1,ad.StarPropertiesFile[i]) == 1);
	if (ad.doswap) reorder(&fdummy,sizeof(float),1);
	if (i == 0) asp->mass = fdummy;
	else if (i == 1) asp->initial_mass = fdummy;
	else if (i == 2) asp->t_form = fdummy;
	else if (i == 3) asp->metallicity_SNII = fdummy;
	else if (i == 4) asp->metallicity_SNIa = fdummy;
	}
    }

void move_art_nb_gas_filepositions_level_begin(ART_DATA ad, int level) {

    int header, trailer;
		
    fread(&header,sizeof(int),1,ad.GasFile[0]);
    if (ad.doswap) reorder(&header,sizeof(int),1);
    
    fread(&header,sizeof(int),1,ad.GasFile[1]);
    if (ad.doswap) reorder(&header,sizeof(int),1);
    assert(fseek(ad.GasFile[1],ad.Ncell[level]*ad.Nhydroproperties*sizeof(float),SEEK_CUR) == 0);
    fread(&trailer,sizeof(int),1,ad.GasFile[1]);
    if (ad.doswap) reorder(&trailer,sizeof(int),1);
    assert(header == trailer);
    
    fread(&header,sizeof(int),1,ad.GasFile[1]);
    if (ad.doswap) reorder(&header,sizeof(int),1);
    }

void move_art_nb_gas_filepositions_level_end(ART_DATA ad, int level) {

    int header, trailer;
    
    assert(fread(&trailer,sizeof(int),1,ad.GasFile[0]) == 1);
    if (ad.doswap) reorder(&trailer,sizeof(int),1);
    
    fread(&header,sizeof(int),1,ad.GasFile[0]);
    if (ad.doswap) reorder(&header,sizeof(int),1);
    assert(fseek(ad.GasFile[0],ad.Ncell[level]*ad.Notherproperties*sizeof(float),SEEK_CUR) == 0);
    assert(fread(&trailer,sizeof(int),1,ad.GasFile[0]) == 1);
    if (ad.doswap) reorder(&trailer,sizeof(int),1);
    assert(header == trailer);
    
    assert(fread(&trailer,sizeof(int),1,ad.GasFile[1]) == 1);
    if (ad.doswap) reorder(&trailer,sizeof(int),1);

    if (level == ad.Lmaxgas) {
	fclose(ad.GasFile[0]);
	if (ad.GRAVITY || ad.RADIATIVE_TRANSFER) {
	    fclose(ad.GasFile[1]);
	    }
	}
    }

void move_art_nb_star_filepositions_begin(ART_DATA ad) {

    int i;
    int header;
    long int seekamount;

    for (i = 0; i < ad.Nstarproperties; i++) {
	assert(fread(&header,sizeof(int),1,ad.StarPropertiesFile[i]) == 1);
	if (ad.doswap) reorder(&header,sizeof(int),1);
	assert(header == ad.Nstar*sizeof(float));
	seekamount = i*(ad.Nstar*sizeof(float) + 2*sizeof(int));
	assert(fseek(ad.StarPropertiesFile[i],seekamount,SEEK_CUR) == 0);
	}
    }

void move_art_nb_star_filepositions_end(ART_DATA ad) {

    int i;
    int trailer;

    for (i = 0; i < ad.Nstarproperties; i++) {
	assert(fread(&trailer,sizeof(int),1,ad.StarPropertiesFile[i]) == 1);
	if (ad.doswap) reorder(&trailer,sizeof(int),1);
	assert(trailer == ad.Nstar*sizeof(float));
	fclose(ad.StarPropertiesFile[i]);
	}
    }

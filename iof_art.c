/*
** iof_art.c
**
** Various reading and writing functions for art format
**
** written by Marcel Zemp
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
** Reading and writing functions
*/

void read_art_nb_general_header(ART_DATA *ad) {

    int i, L;
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
	if (ad->toplevelmassdark == -1) {
	    ad->massfromdata = 1;
	    ad->toplevelmassdark = ad->ah.mass[ad->Lmaxdark];
	    }
	assert(ad->toplevelsoftdark >= 0);
	assert(ad->toplevelmassdark >= 0);
	for (i = 0; i <= ad->Lmaxdark; i++) {
	    L = ad->Lmaxdark-i;
	    if (i == 0) {
		ad->Ndarklevel[L] = ad->ah.num[i];
		}
	    else {
		ad->Ndarklevel[L] = ad->ah.num[i]-ad->ah.num[i-1];
		}
	    if (ad->massfromdata == 1) {
		ad->massdark[L] = ad->ah.mass[i];
		}
	    }
	}
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
    ad->sfci.sfc_order = idummy;
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
    /*
    ** Ncell[0]
    */
    assert(fread(&header,sizeof(int),1,ad->GasFile[index]) == 1);
    if (ad->doswap) reorder(&header,sizeof(int),1);
    assert(fread(&idummy,sizeof(int),1,ad->GasFile[index]) == 1);
    if (ad->doswap) reorder(&idummy,sizeof(int),1);
    ad->Ncell[0] = idummy;
    assert(ad->Ncell[0] == ad->ah.Ngrid*ad->ah.Ngrid*ad->ah.Ngrid);
    assert(fread(&trailer,sizeof(int),1,ad->GasFile[index]) == 1);
    if (ad->doswap) reorder(&trailer,sizeof(int),1);
    assert(header == trailer);
    }

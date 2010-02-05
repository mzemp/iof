/*
** gadget.c
**
** Various reading and writing functions for gadget format
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
#include "iof_tipsy.h"
#include "iof_gadget.h"

/*
** Comparison functions
*/

int comp_tgp(const void *a, const void *b) {

    QSORT_TGP *aa = (QSORT_TGP *) a;
    QSORT_TGP *bb = (QSORT_TGP *) b;

    return (aa->index < bb->index)?-1:1;
    }

int comp_tdp(const void *a, const void *b) {

    QSORT_TDP *aa = (QSORT_TDP *) a;
    QSORT_TDP *bb = (QSORT_TDP *) b;

    return (aa->index < bb->index)?-1:1;
    }

int comp_tsp(const void *a, const void *b) {

    QSORT_TSP *aa = (QSORT_TSP *) a;
    QSORT_TSP *bb = (QSORT_TSP *) b;

    return (aa->index < bb->index)?-1:1;
    }

/*
** Reading and writing functions
*/

void read_gadget_nb(FILE *fp, TIPSY_STRUCTURE *ts, double a, double dx, double dy, double dz, double dof, double mmw, double uvf) {

    int i, j, dummy1, dummy2, Npwm, Ntotgad;
    int *gasindex, *darkindex, *starindex;
    float temp;
    double inverse_sqrt_a, deltapos[3];
    GADGET_HEADER gh;
    TIPSY_HEADER *th;
    TIPSY_GAS_PARTICLE *tgp;
    TIPSY_DARK_PARTICLE *tdp;
    TIPSY_STAR_PARTICLE *tsp;
    QSORT_TGP *qtgp;
    QSORT_TDP *qtdp;
    QSORT_TSP *qtsp;

    th = ts->th;
    tgp = ts->tgp;
    tdp = ts->tdp;
    tsp = ts->tsp;
    deltapos[0] = dx;
    deltapos[1] = dy;
    deltapos[2] = dz;
    /*
    ** Read in header
    */
    th = realloc(th,sizeof(TIPSY_HEADER));
    assert(th != NULL);
    assert(fread(&dummy1,sizeof(int),1,fp) == 1);
    assert(fread(&gh,sizeof(GADGET_HEADER),1,fp) == 1);
    assert(fread(&dummy2,sizeof(int),1,fp) == 1);
    assert(dummy1 == dummy2 && dummy2 == sizeof(GADGET_HEADER));
    if (gh.num_files > 1) {
	fprintf(stderr,"Sorry, this snapshot contains more than one file!\n");
	fprintf(stderr,"This is not supported by this program.\n");
	exit(0);
	}
    th->time = gh.time;
    th->ntotal = gh.npart[0] + gh.npart[1] + gh.npart[4];
    th->ndim = 3;
    th->ngas = gh.npart[0];
    th->ndark = gh.npart[1];
    th->nstar = gh.npart[4];
    Npwm = 0;
    for (i = 0; i < 6; i++) {
	if (gh.mass[i] == 0) {
	    Npwm += gh.npart[i];
	    }
	}
    Ntotgad = 0;
    for (i = 0; i < 6; i++) {
	Ntotgad += gh.npart[i];
	}
    /*
    ** Determine scale factor
    */
    if (a == 1) {
	inverse_sqrt_a = 1;
	}
    else if (a == -1) {
	inverse_sqrt_a = 1/sqrt(th->time);
	}
    else {
	inverse_sqrt_a = 1/sqrt(a);
	}
    /*
    ** Allocate memory
    */
    if(th->ngas != 0) {
	tgp = realloc(tgp,th->ngas*sizeof(TIPSY_GAS_PARTICLE));
	assert(tgp != NULL);
	}
    else {
	tgp = NULL;
	}
    if(th->ndark != 0) {
	tdp = realloc(tdp,th->ndark*sizeof(TIPSY_DARK_PARTICLE));
	assert(tdp != NULL);
	}
    else {
	tdp = NULL;
	}
    if(th->nstar != 0) {
	tsp = realloc(tsp,th->nstar*sizeof(TIPSY_STAR_PARTICLE));
	assert(tsp != NULL);
	}
    else {
	tsp = NULL;
	}
    /*
    ** Read in positions
    */
    assert(fread(&dummy1,sizeof(int),1,fp) == 1); 
    for(i = 0; i < th->ngas; i++) {
	for(j = 0; j < 3; j++) {
	    assert(fread(&temp,sizeof(float),1,fp) == 1);
	    tgp[i].pos[j] = temp + deltapos[j];
	    }
	}
    for(i = 0; i < th->ndark; i++) {
	for(j = 0; j < 3; j++) {
	    assert(fread(&temp,sizeof(float),1,fp) == 1);
	    tdp[i].pos[j] = temp + deltapos[j];
	    }
	}
    fseek(fp,3*gh.npart[2]*sizeof(float),SEEK_CUR);
    fseek(fp,3*gh.npart[3]*sizeof(float),SEEK_CUR);
    for(i = 0; i < th->nstar; i++) {
	for(j = 0; j < 3; j++) {
	    assert(fread(&temp,sizeof(float),1,fp) == 1);
	    tsp[i].pos[j] = temp + deltapos[j];
	    }
	}
    fseek(fp,3*gh.npart[5]*sizeof(float),SEEK_CUR);
    assert(fread(&dummy2,sizeof(int),1,fp) == 1);
    assert(dummy1 == dummy2 && dummy2 == 3*Ntotgad*sizeof(float));
    /*
    ** Read in velocities
    */
    assert(fread(&dummy1,sizeof(int),1,fp) == 1); 
    for(i = 0; i < th->ngas; i++) {
	for(j = 0; j < 3; j++) {
	    assert(fread(&tgp[i].vel[j],sizeof(float),1,fp) == 1);
	    tgp[i].vel[j] *= inverse_sqrt_a;
	    }
	}
    for(i = 0; i < th->ndark; i++) {
	for(j = 0; j < 3; j++) {
	    assert(fread(&tdp[i].vel[j],sizeof(float),1,fp) == 1);
	    tdp[i].vel[j] *= inverse_sqrt_a;
	    }
	}
    fseek(fp,3*gh.npart[2]*sizeof(float),SEEK_CUR);
    fseek(fp,3*gh.npart[3]*sizeof(float),SEEK_CUR);
    for(i = 0; i < th->nstar; i++) {
	for(j = 0; j < 3; j++) {
	    assert(fread(&tsp[i].vel[j],sizeof(float),1,fp) == 1);
	    tsp[i].vel[j] *= inverse_sqrt_a;
	    }
	}
    fseek(fp,3*gh.npart[5]*sizeof(float),SEEK_CUR);
    assert(fread(&dummy2,sizeof(int),1,fp) == 1);
    assert(dummy1 == dummy2 && dummy2 == 3*Ntotgad*sizeof(float));
    /*
    ** Read in indices
    */
    gasindex = malloc(th->ngas*sizeof(int));
    darkindex = malloc(th->ndark*sizeof(int));
    starindex = malloc(th->nstar*sizeof(int));
    assert(gasindex != NULL);
    assert(darkindex != NULL);
    assert(starindex != NULL);
    assert(fread(&dummy1,sizeof(int),1,fp) == 1); 
    for(i = 0; i < th->ngas; i++) {
	assert(fread(&gasindex[i],sizeof(int),1,fp) == 1);
	}
    for(i = 0; i < th->ndark; i++) {
	assert(fread(&darkindex[i],sizeof(int),1,fp) == 1);
	}
    fseek(fp,gh.npart[2]*sizeof(int),SEEK_CUR);
    fseek(fp,gh.npart[3]*sizeof(int),SEEK_CUR);
    for(i = 0; i < th->nstar; i++) {
	assert(fread(&starindex[i],sizeof(int),1,fp) == 1);
	}
    fseek(fp,gh.npart[5]*sizeof(int),SEEK_CUR);
    assert(fread(&dummy2,sizeof(int),1,fp) == 1);
    assert(dummy1 == dummy2 && dummy2 == Ntotgad*sizeof(int));
    /*
    ** Read in masses
    */
    if (Npwm > 0) {
	assert(fread(&dummy1,sizeof(int),1,fp) == 1);
	if (gh.mass[0] == 0) {
	    for(i = 0; i < th->ngas; i++) {
		assert(fread(&tgp[i].mass,sizeof(float),1,fp) == 1);
		}
	    }
	else {
	    for(i = 0; i < th->ngas; i++) {
		tgp[i].mass = gh.mass[0];
		}
	    }
	if (gh.mass[1] == 0) {
	    for(i = 0; i < th->ndark; i++) {
		assert(fread(&tdp[i].mass,sizeof(float),1,fp) == 1);
		}
	    }
	else {
	    for(i = 0; i < th->ndark; i++) {
		tdp[i].mass = gh.mass[1];
		}
	    }
	if (gh.mass[2] == 0) {
	    fseek(fp,gh.npart[2]*sizeof(float),SEEK_CUR);
	    }
	if (gh.mass[3] == 0) {
	    fseek(fp,gh.npart[3]*sizeof(float),SEEK_CUR);
	    }
	if (gh.mass[4] == 0) {
	    for(i = 0; i < th->nstar; i++) {
		assert(fread(&tsp[i].mass,sizeof(float),1,fp) == 1);
		}
	    }
	else {
	    for(i = 0; i < th->nstar; i++) {
		tsp[i].mass = gh.mass[4];
		}
	    }
	if (gh.mass[5] == 0) {
	    fseek(fp,gh.npart[5]*sizeof(float),SEEK_CUR);
	    }
	assert(fread(&dummy2,sizeof(int),1,fp) == 1);
	assert(dummy1 == dummy2 && dummy2 == Npwm*sizeof(float));
	}
    else {
	for(i = 0; i < th->ngas; i++) {
	    tgp[i].mass = gh.mass[0];
	    }
	for(i = 0; i < th->ndark; i++) {
	    tdp[i].mass = gh.mass[1];
	    }
	for(i = 0; i < th->nstar; i++) {
	    tsp[i].mass = gh.mass[4];
	    }
	}
    /*
    ** Read in gas temperatures
    */
    if (th->ngas > 0) {
	assert(fread(&dummy1,sizeof(int),1,fp) == 1);
	for(i = 0; i < th->ngas; i++) {
	    assert(fread(&temp,sizeof(float),1,fp) == 1);
	    tgp[i].temp = temp*(uvf*uvf)*(2*mmw*PROTON_MASS)/(dof*k_BOLTZMANN);
	    }
	assert(fread(&dummy2,sizeof(int),1,fp) == 1);
	assert(dummy1 == dummy2 && dummy2 == th->ngas*sizeof(float));
	}
    /*
    ** Set values for parameters that are normally not in the gadget file
    */
    for(i = 0; i < th->ngas; i++) {
	tgp[i].rho = 0;
	tgp[i].hsmooth = 0;
	tgp[i].metals = 0;
	tgp[i].phi = 0;
	}
    for(i = 0; i < th->ndark; i++) {
	tdp[i].eps = 0;
	tdp[i].phi = 0;
	}
    for(i = 0; i < th->nstar; i++) {
	tsp[i].metals = 0;
	tsp[i].tform = 0;
	tsp[i].eps = 0;
	tsp[i].phi = 0;
	}
    /*
    ** Read in possible additional stuff and overwrite
    */
    fread(&dummy1,sizeof(int),1,fp);
    if (feof(fp) == 0) {
	if (th->ngas > 0) {
	    /*
	    ** Read in gas densities
	    */
	    for(i = 0; i < th->ngas; i++) {
		assert(fread(&tgp[i].rho,sizeof(float),1,fp) == 1);
		}
	    assert(fread(&dummy2,sizeof(int),1,fp) == 1);
	    assert(dummy1 == dummy2 && dummy2 == th->ngas*sizeof(float));
	    /*
	    ** Read in softening lengths
	    */
	    fread(&dummy1,sizeof(int),1,fp);
	    if (feof(fp) == 0) {
		for(i = 0; i < th->ngas; i++) {
		    assert(fread(&tgp[i].hsmooth,sizeof(float),1,fp) == 1);
		    }
		assert(fread(&dummy2,sizeof(int),1,fp) == 1);
		assert(dummy1 == dummy2 && dummy2 == th->ngas*sizeof(float));
		fread(&dummy1,sizeof(int),1,fp);
		}
	    }
	/*
	** Read in potenitals
	*/
	if (feof(fp) == 0) {
	    for(i = 0; i < th->ngas; i++) {
		assert(fread(&tgp[i].phi,sizeof(float),1,fp) == 1);
		}
	    for(i = 0; i < th->ndark; i++) {
		assert(fread(&tdp[i].phi,sizeof(float),1,fp) == 1);
		}
	    for(i = 0; i < th->nstar; i++) {
		assert(fread(&tsp[i].phi,sizeof(float),1,fp) == 1);
		}
	    assert(fread(&dummy2,sizeof(int),1,fp) == 1);
	    assert(dummy1 == dummy2 && dummy2 == th->ntotal*sizeof(float));
	    }
	}
    /*
    ** Sort particles
    */
    qtgp = malloc(th->ngas*sizeof(QSORT_TGP));
    assert(qtgp != NULL);
    for (i = 0; i < th->ngas; i++) {
	qtgp[i].index = gasindex[i];
	qtgp[i].tgp = tgp[i];
	}
    qsort(qtgp,th->ngas,sizeof(QSORT_TGP),comp_tgp);
    for (i = 0; i < th->ngas; i++) {
	tgp[i] = qtgp[i].tgp;
	}
    qtdp = malloc(th->ndark*sizeof(QSORT_TDP));
    assert(qtdp != NULL);
    for (i = 0; i < th->ndark; i++) {
	qtdp[i].index = darkindex[i];
	qtdp[i].tdp = tdp[i];
	}
    qsort(qtdp,th->ndark,sizeof(QSORT_TDP),comp_tdp);
    for (i = 0; i < th->ndark; i++) {
	tdp[i] = qtdp[i].tdp;
	}
    qtsp = malloc(th->nstar*sizeof(QSORT_TSP));
    assert(qtsp != NULL);
    for (i = 0; i < th->nstar; i++) {
	qtsp[i].index = starindex[i];
	qtsp[i].tsp = tsp[i];
	}
    qsort(qtsp,th->nstar,sizeof(QSORT_TSP),comp_tsp);
    for (i = 0; i < th->nstar; i++) {
	tsp[i] = qtsp[i].tsp;
	}
    free(gasindex);
    free(darkindex);
    free(starindex);
    free(qtgp);
    free(qtdp);
    free(qtsp);
    /*
    ** Return pointers
    */
    ts->th = th;
    ts->tgp = tgp;
    ts->tdp = tdp;
    ts->tsp = tsp;
    }

void write_gadget_nb(FILE *fp, const TIPSY_STRUCTURE *ts, double a, double dx, double dy, double dz, double dof, double mmw, double uvf) {

    int i, j, dummy, index;
    float temp;
    double sqrt_a, deltapos[3];
    GADGET_HEADER gh;
    TIPSY_HEADER *th;
    TIPSY_GAS_PARTICLE *tgp;
    TIPSY_DARK_PARTICLE *tdp;
    TIPSY_STAR_PARTICLE *tsp;

    th = ts->th;
    tgp = ts->tgp;
    tdp = ts->tdp;
    tsp = ts->tsp;
    deltapos[0] = dx;
    deltapos[1] = dy;
    deltapos[2] = dz;
    /*
    ** Determine scale factor
    */
    if (a == 1) {
	sqrt_a = 1;
	}
    else if (a == -1) {
	sqrt_a = sqrt(th->time);
	}
    else {
	sqrt_a = sqrt(a);
	}
    /*
    ** Initialise header
    */
    gh.npart[0] = th->ngas;
    gh.npart[1] = th->ndark;
    gh.npart[2] = 0;
    gh.npart[3] = 0;
    gh.npart[4] = th->nstar;
    gh.npart[5] = 0;
    for (i = 0; i < 6; i++) {
	gh.mass[i] = 0;
	}
    gh.time = th->time;
    gh.redshift = 0;
    gh.flag_sfr = 0;
    gh.flag_feedback = 0;
    gh.npartTotal[0] = th->ngas;
    gh.npartTotal[1] = th->ndark;
    gh.npartTotal[2] = 0;
    gh.npartTotal[3] = 0;
    gh.npartTotal[4] = th->nstar;
    gh.npartTotal[5] = 0;
    gh.flag_cooling = 0;
    gh.num_files = 1;
    gh.BoxSize = 0;
    gh.Omega0 = 0;
    gh.OmegaLambda = 0;
    gh.HubbleParam = 0;
    gh.flag_stellarage = 0;
    gh.flag_metals = 0;
    for (i = 0; i < 6; i++) {
	gh.npartTotalHighWord[i] = 0;
	}
    gh.flag_entropy_instead_u = 0;
    /*
    ** Write out header
    */
    dummy = sizeof(GADGET_HEADER);
    assert(fwrite(&dummy,sizeof(int),1,fp) == 1);
    assert(fwrite(&gh,sizeof(GADGET_HEADER),1,fp) == 1);
    assert(fwrite(&dummy,sizeof(int),1,fp) == 1);
    /*
    ** Write out positions
    */
    dummy = 3*th->ntotal*sizeof(float);
    assert(fwrite(&dummy,sizeof(int),1,fp) == 1); 
    for(i = 0; i < th->ngas; i++) {
	for(j = 0; j < 3; j++) {
	    temp = tgp[i].pos[j] + deltapos[j];
	    assert(fwrite(&temp,sizeof(float),1,fp) == 1);
	    }
	}
    for(i = 0; i < th->ndark; i++) {
	for(j = 0; j < 3; j++) {
	    temp = tdp[i].pos[j] + deltapos[j];
	    assert(fwrite(&temp,sizeof(float),1,fp) == 1);
	    }
	}
    for(i = 0; i < th->nstar; i++) {
	for(j = 0; j < 3; j++) {
	    temp = tsp[i].pos[j] + deltapos[j];
	    assert(fwrite(&temp,sizeof(float),1,fp) == 1);
	    }
	}
    assert(fwrite(&dummy,sizeof(int),1,fp) == 1);
    /*
    ** Write out velocities
    */
    dummy = 3*th->ntotal*sizeof(float);
    assert(fwrite(&dummy,sizeof(int),1,fp) == 1);
    for(i = 0; i < th->ngas; i++) {
	for(j = 0; j < 3; j++) {
	    temp = tgp[i].vel[j]*sqrt_a;
	    assert(fwrite(&temp,sizeof(float),1,fp) == 1);
	    }
	}
    for(i = 0; i < th->ndark; i++) {
	for(j = 0; j < 3; j++) {
	    temp = tdp[i].vel[j]*sqrt_a;
	    assert(fwrite(&temp,sizeof(float),1,fp) == 1);
	    }
	}
    for(i = 0; i < th->nstar; i++) {
	for(j = 0; j < 3; j++) {
	    temp = tsp[i].vel[j]*sqrt_a;
	    assert(fwrite(&temp,sizeof(float),1,fp) == 1);
	    }
	}
    assert(fwrite(&dummy,sizeof(int),1,fp) == 1);
    /*
    ** Write out indices
    */
    dummy = th->ntotal*sizeof(int);
    assert(fwrite(&dummy,sizeof(int),1,fp) == 1);
    index = 0;
    for(i = 0; i < th->ngas; i++) {
	index++;
	assert(fwrite(&index,sizeof(int),1,fp) == 1);
	}
    for(i = 0; i < th->ndark; i++) {
	index++;
	assert(fwrite(&index,sizeof(int),1,fp) == 1);
	}
    for(i = 0; i < th->nstar; i++) {
	index++;
	assert(fwrite(&index,sizeof(int),1,fp) == 1);
	}
    assert(fwrite(&dummy,sizeof(int),1,fp) == 1);
    /*
    ** Write out masses
    */
    dummy = th->ntotal*sizeof(float);
    assert(fwrite(&dummy,sizeof(int),1,fp) == 1);
    for(i = 0; i < th->ngas; i++) {
	assert(fwrite(&tgp[i].mass,sizeof(float),1,fp) == 1);
	}
    for(i = 0; i < th->ndark; i++) {
	assert(fwrite(&tdp[i].mass,sizeof(float),1,fp) == 1);
	}
    for(i = 0; i < th->nstar; i++) {
	assert(fwrite(&tsp[i].mass,sizeof(float),1,fp) == 1);
	}
    assert(fwrite(&dummy,sizeof(int),1,fp) == 1);
    /*
    ** Write out gas specific stuff
    */
    if (th->ngas > 0) {
	dummy = th->ngas*sizeof(float);
	assert(fwrite(&dummy,sizeof(int),1,fp) == 1);
	for(i = 0; i < th->ngas; i++) {
	    temp = tgp[i].temp*(dof*k_BOLTZMANN)/(2*mmw*PROTON_MASS)/(uvf*uvf);
	    assert(fwrite(&temp,sizeof(float),1,fp) == 1);
	    }
	assert(fwrite(&dummy,sizeof(int),1,fp) == 1);
	dummy = th->ngas*sizeof(float);
	assert(fwrite(&dummy,sizeof(int),1,fp) == 1);
	for(i = 0; i < th->ngas; i++) {
	    assert(fwrite(&tgp[i].rho,sizeof(float),1,fp) == 1);
	    }
	assert(fwrite(&dummy,sizeof(int),1,fp) == 1);
	dummy = th->ngas*sizeof(float);
	assert(fwrite(&dummy,sizeof(int),1,fp) == 1);
	for(i = 0; i < th->ngas; i++) {
	    assert(fwrite(&tgp[i].hsmooth,sizeof(float),1,fp) == 1);
	    }
	assert(fwrite(&dummy,sizeof(int),1,fp) == 1);
	}
    /*
    ** Write out potenitals
    */
    dummy = th->ntotal*sizeof(float);
    assert(fwrite(&dummy,sizeof(int),1,fp) == 1);
    for(i = 0; i < th->ngas; i++) {
	assert(fwrite(&tgp[i].phi,sizeof(float),1,fp) == 1);
	}
    for(i = 0; i < th->ndark; i++) {
	assert(fwrite(&tdp[i].phi,sizeof(float),1,fp) == 1);
	}
    for(i = 0; i < th->nstar; i++) {
	assert(fwrite(&tsp[i].phi,sizeof(float),1,fp) == 1);
	}
    assert(fwrite(&dummy,sizeof(int),1,fp) == 1);
    }

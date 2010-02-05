/*
** tipsy.c
**
** Various reading and writing functions for tipsy format
**
** Written by Marcel Zemp
*/

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <malloc.h>
#include <assert.h>
#include "iof_tipsy.h"

/*
** Copy functions
*/

void copy_tgp_to_tgpdpp(const TIPSY_GAS_PARTICLE *tgp, TIPSY_GAS_PARTICLE_DPP *tgpdpp) {

    int j;

    tgpdpp->mass = tgp->mass;
    for (j = 0; j < 3; j++) {
	tgpdpp->pos[j] = tgp->pos[j];
	}
    for (j = 0; j < 3; j++) {
	tgpdpp->vel[j] = tgp->vel[j];
	}
    tgpdpp->rho = tgp->rho;
    tgpdpp->temp = tgp->temp;
    tgpdpp->hsmooth = tgp->hsmooth;
    tgpdpp->metals = tgp->metals;
    tgpdpp->phi = tgp->phi;
    }

void copy_tdp_to_tdpdpp(const TIPSY_DARK_PARTICLE *tdp, TIPSY_DARK_PARTICLE_DPP *tdpdpp) {

    int j;

    tdpdpp->mass = tdp->mass;
    for (j = 0; j < 3; j++) {
	tdpdpp->pos[j] = tdp->pos[j];
	}
    for (j = 0; j < 3; j++) {
	tdpdpp->vel[j] = tdp->vel[j];
	}
    tdpdpp->eps = tdp->eps;
    tdpdpp->phi = tdp->phi;
    }

void copy_tsp_to_tspdpp(const TIPSY_STAR_PARTICLE *tsp, TIPSY_STAR_PARTICLE_DPP *tspdpp) {

    int j;

    tspdpp->mass = tsp->mass;
    for (j = 0; j < 3; j++) {
	tspdpp->pos[j] = tsp->pos[j];
	}
    for (j = 0; j < 3; j++) {
	tspdpp->vel[j] = tsp->vel[j];
	}
    tspdpp->metals = tsp->metals;
    tspdpp->tform = tsp->tform;
    tspdpp->eps = tsp->eps;
    tspdpp->phi = tsp->phi;
    }

void copy_tgpdpp_to_tgp(const TIPSY_GAS_PARTICLE_DPP *tgpdpp, TIPSY_GAS_PARTICLE *tgp) {

    int j;

    tgp->mass = tgpdpp->mass;
    for (j = 0; j < 3; j++) {
	tgp->pos[j] = tgpdpp->pos[j];
	}
    for (j = 0; j < 3; j++) {
	tgp->vel[j] = tgpdpp->vel[j];
	}
    tgp->rho = tgpdpp->rho;
    tgp->temp = tgpdpp->temp;
    tgp->hsmooth = tgpdpp->hsmooth;
    tgp->metals = tgpdpp->metals;
    tgp->phi = tgpdpp->phi;
    }

void copy_tdpdpp_to_tdp(const TIPSY_DARK_PARTICLE_DPP *tdpdpp, TIPSY_DARK_PARTICLE *tdp) {

    int j;

    tdp->mass = tdpdpp->mass;
    for (j = 0; j < 3; j++) {
	tdp->pos[j] = tdpdpp->pos[j];
	}
    for (j = 0; j < 3; j++) {
	tdp->vel[j] = tdpdpp->vel[j];
	}
    tdp->eps = tdpdpp->eps;
    tdp->phi = tdpdpp->phi;
    }

void copy_tspdpp_to_tsp(const TIPSY_STAR_PARTICLE_DPP *tspdpp, TIPSY_STAR_PARTICLE *tsp) {

    int j;

    tsp->mass = tspdpp->mass;
    for (j = 0; j < 3; j++) {
	tsp->pos[j] = tspdpp->pos[j];
	}
    for (j = 0; j < 3; j++) {
	tsp->vel[j] = tspdpp->vel[j];
	}
    tsp->metals = tspdpp->metals;
    tsp->tform = tspdpp->tform;
    tsp->eps = tspdpp->eps;
    tsp->phi = tspdpp->phi;
    }

/*
** Reading and writing functions
*/

void read_tipsy_nb_header(FILE *fp, TIPSY_HEADER *th) {

    assert(fread(th,sizeof(TIPSY_HEADER),1,fp) == 1);
    }

void read_tipsy_nb_gas(FILE *fp, TIPSY_GAS_PARTICLE *tgp) {

    assert(fread(tgp,sizeof(TIPSY_GAS_PARTICLE),1,fp) == 1);
    }

void read_tipsy_nb_dark(FILE *fp, TIPSY_DARK_PARTICLE *tdp) {

    assert(fread(tdp,sizeof(TIPSY_DARK_PARTICLE),1,fp) == 1);
    }

void read_tipsy_nb_star(FILE *fp, TIPSY_STAR_PARTICLE *tsp) {

    assert(fread(tsp,sizeof(TIPSY_STAR_PARTICLE),1,fp) == 1);
    }

void read_tipsy_nb_gas_dpp(FILE *fp, TIPSY_GAS_PARTICLE_DPP *tgpdpp) {

    assert(fread(tgpdpp,sizeof(TIPSY_GAS_PARTICLE_DPP),1,fp) == 1);
    }

void read_tipsy_nb_dark_dpp(FILE *fp, TIPSY_DARK_PARTICLE_DPP *tdpdpp) {

    assert(fread(tdpdpp,sizeof(TIPSY_DARK_PARTICLE_DPP),1,fp) == 1);
    }

void read_tipsy_nb_star_dpp(FILE *fp, TIPSY_STAR_PARTICLE_DPP *tspdpp) {

    assert(fread(tspdpp,sizeof(TIPSY_STAR_PARTICLE_DPP),1,fp) == 1);
    }

void read_tipsy_nb(FILE *fp, TIPSY_STRUCTURE *ts) {

    TIPSY_HEADER *th;
    TIPSY_GAS_PARTICLE *tgp;
    TIPSY_DARK_PARTICLE *tdp;
    TIPSY_STAR_PARTICLE *tsp;

    th = ts->th;
    tgp = ts->tgp;
    tdp = ts->tdp;
    tsp = ts->tsp;
    /*
    ** Read in header
    */
    th = realloc(th,sizeof(TIPSY_HEADER));
    assert(th != NULL);
    assert(fread(th,sizeof(TIPSY_HEADER),1,fp) == 1);
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
    ** Read in particles
    */
    assert(fread(tgp,sizeof(TIPSY_GAS_PARTICLE),th->ngas,fp) == th->ngas);
    assert(fread(tdp,sizeof(TIPSY_DARK_PARTICLE),th->ndark,fp) == th->ndark);
    assert(fread(tsp,sizeof(TIPSY_STAR_PARTICLE),th->nstar,fp) == th->nstar);
    /*
    ** Return pointers
    */
    ts->th = th;
    ts->tgp = tgp;
    ts->tdp = tdp;
    ts->tsp = tsp;
    }

void read_tipsy_nb_dpp(FILE *fp, TIPSY_STRUCTURE_DPP *tsdpp) {

    TIPSY_HEADER *th;
    TIPSY_GAS_PARTICLE_DPP *tgpdpp;
    TIPSY_DARK_PARTICLE_DPP *tdpdpp;
    TIPSY_STAR_PARTICLE_DPP *tspdpp;

    th = tsdpp->th;
    tgpdpp = tsdpp->tgpdpp;
    tdpdpp = tsdpp->tdpdpp;
    tspdpp = tsdpp->tspdpp;
    /*
    ** Read in header
    */
    th = realloc(th,sizeof(TIPSY_HEADER));
    assert(th != NULL);
    assert(fread(th,sizeof(TIPSY_HEADER),1,fp) == 1);
    /*
    ** Allocate memory
    */
    if(th->ngas != 0) {
	tgpdpp = realloc(tgpdpp,th->ngas*sizeof(TIPSY_GAS_PARTICLE_DPP));
	assert(tgpdpp != NULL);
	}
    else {
	tgpdpp = NULL;
	}
    if(th->ndark != 0) {
	tdpdpp = realloc(tdpdpp,th->ndark*sizeof(TIPSY_DARK_PARTICLE_DPP));
	assert(tdpdpp != NULL);
	}
    else {
	tdpdpp = NULL;
	}
    if(th->nstar != 0) {
	tspdpp = realloc(tspdpp,th->nstar*sizeof(TIPSY_STAR_PARTICLE_DPP));
	assert(tspdpp != NULL);
	}
    else {
	tspdpp = NULL;
	}
    /*
    ** Read in particles
    */
    assert(fread(tgpdpp,sizeof(TIPSY_GAS_PARTICLE_DPP),th->ngas,fp) == th->ngas);
    assert(fread(tdpdpp,sizeof(TIPSY_DARK_PARTICLE_DPP),th->ndark,fp) == th->ndark);
    assert(fread(tspdpp,sizeof(TIPSY_STAR_PARTICLE_DPP),th->nstar,fp) == th->nstar);
    /*
    ** Return pointers
    */
    tsdpp->th = th;
    tsdpp->tgpdpp = tgpdpp;
    tsdpp->tdpdpp = tdpdpp;
    tsdpp->tspdpp = tspdpp;
    }

void write_tipsy_nb_header(FILE *fp, const TIPSY_HEADER *th) {

    assert(fwrite(th,sizeof(TIPSY_HEADER),1,fp) == 1);
    }

void write_tipsy_nb_gas(FILE *fp, const TIPSY_GAS_PARTICLE *tgp) {

    assert(fwrite(tgp,sizeof(TIPSY_GAS_PARTICLE),1,fp) == 1);
    }

void write_tipsy_nb_dark(FILE *fp, const TIPSY_DARK_PARTICLE *tdp) {

    assert(fwrite(tdp,sizeof(TIPSY_DARK_PARTICLE),1,fp) == 1);
    }

void write_tipsy_nb_star(FILE *fp, const TIPSY_STAR_PARTICLE *tsp) {

    assert(fwrite(tsp,sizeof(TIPSY_STAR_PARTICLE),1,fp) == 1);
    }

void write_tipsy_nb_gas_dpp(FILE *fp, const TIPSY_GAS_PARTICLE_DPP *tgpdpp) {

    assert(fwrite(tgpdpp,sizeof(TIPSY_GAS_PARTICLE_DPP),1,fp) == 1);
    }

void write_tipsy_nb_dark_dpp(FILE *fp, const TIPSY_DARK_PARTICLE_DPP *tdpdpp) {

    assert(fwrite(tdpdpp,sizeof(TIPSY_DARK_PARTICLE_DPP),1,fp) == 1);
    }

void write_tipsy_nb_star_dpp(FILE *fp, const TIPSY_STAR_PARTICLE_DPP *tspdpp) {

    assert(fwrite(tspdpp,sizeof(TIPSY_STAR_PARTICLE_DPP),1,fp) == 1);
    }

void write_tipsy_nb(FILE *fp, const TIPSY_STRUCTURE *ts) {

    TIPSY_HEADER *th;
    TIPSY_GAS_PARTICLE *tgp;
    TIPSY_DARK_PARTICLE *tdp;
    TIPSY_STAR_PARTICLE *tsp;

    th = ts->th;
    tgp = ts->tgp;
    tdp = ts->tdp;
    tsp = ts->tsp;
    /*
    ** Write out header
    */
    assert(fwrite(th,sizeof(TIPSY_HEADER),1,fp) == 1);
    /*
    ** Write out particles
    */
    assert(fwrite(tgp,sizeof(TIPSY_GAS_PARTICLE),th->ngas,fp) == th->ngas);
    assert(fwrite(tdp,sizeof(TIPSY_DARK_PARTICLE),th->ndark,fp) == th->ndark);
    assert(fwrite(tsp,sizeof(TIPSY_STAR_PARTICLE),th->nstar,fp) == th->nstar);
    }

void write_tipsy_nb_dpp(FILE *fp, const TIPSY_STRUCTURE_DPP *tsdpp) {

    TIPSY_HEADER *th;
    TIPSY_GAS_PARTICLE_DPP *tgpdpp;
    TIPSY_DARK_PARTICLE_DPP *tdpdpp;
    TIPSY_STAR_PARTICLE_DPP *tspdpp;

    th = tsdpp->th;
    tgpdpp = tsdpp->tgpdpp;
    tdpdpp = tsdpp->tdpdpp;
    tspdpp = tsdpp->tspdpp;
    /*
    ** Write out header
    */
    assert(fwrite(th,sizeof(TIPSY_HEADER),1,fp) == 1);
    /*
    ** Write out particles
    */
    assert(fwrite(tgpdpp,sizeof(TIPSY_GAS_PARTICLE_DPP),th->ngas,fp) == th->ngas);
    assert(fwrite(tdpdpp,sizeof(TIPSY_DARK_PARTICLE_DPP),th->ndark,fp) == th->ndark);
    assert(fwrite(tspdpp,sizeof(TIPSY_STAR_PARTICLE_DPP),th->nstar,fp) == th->nstar);
    }

void read_tipsy_xdr_header(XDR *xdrs, TIPSY_HEADER *th) {

    int pad;

    assert(xdr_double(xdrs,&th->time) == 1);
    assert(xdr_int(xdrs,&th->ntotal) == 1);
    assert(xdr_int(xdrs,&th->ndim) == 1);
    assert(xdr_int(xdrs,&th->ngas) == 1);
    assert(xdr_int(xdrs,&th->ndark) == 1);
    assert(xdr_int(xdrs,&th->nstar) == 1);
    assert(xdr_int(xdrs,&pad) == 1);
    }

void read_tipsy_xdr_gas(XDR *xdrs, TIPSY_GAS_PARTICLE *tgp) {

    int j;

    assert(xdr_float(xdrs,&tgp->mass) == 1);
    for (j = 0; j < 3; j++) {
	assert(xdr_float(xdrs,&tgp->pos[j]) == 1);
	}
    for (j = 0; j < 3; j++) {
	assert(xdr_float(xdrs,&tgp->vel[j]) == 1);
	}
    assert(xdr_float(xdrs,&tgp->rho) == 1);
    assert(xdr_float(xdrs,&tgp->temp) == 1);
    assert(xdr_float(xdrs,&tgp->hsmooth) == 1);
    assert(xdr_float(xdrs,&tgp->metals) == 1);
    assert(xdr_float(xdrs,&tgp->phi) == 1);
    }

void read_tipsy_xdr_dark(XDR *xdrs, TIPSY_DARK_PARTICLE *tdp) {

    int j;

    assert(xdr_float(xdrs,&tdp->mass) == 1);
    for (j = 0; j < 3; j++) {
	assert(xdr_float(xdrs,&tdp->pos[j]) == 1);
	}
    for (j = 0; j < 3; j++) {
	assert(xdr_float(xdrs,&tdp->vel[j]) == 1);
	}
    assert(xdr_float(xdrs,&tdp->eps) == 1);
    assert(xdr_float(xdrs,&tdp->phi) == 1);
    }

void read_tipsy_xdr_star(XDR *xdrs, TIPSY_STAR_PARTICLE *tsp) {

    int j;

    assert(xdr_float(xdrs,&tsp->mass) == 1);
    for (j = 0; j < 3; j++) {
	assert(xdr_float(xdrs,&tsp->pos[j]) == 1);
	}
    for (j = 0; j < 3; j++) {
	assert(xdr_float(xdrs,&tsp->vel[j]) == 1);
	}
    assert(xdr_float(xdrs,&tsp->metals) == 1);
    assert(xdr_float(xdrs,&tsp->tform) == 1);
    assert(xdr_float(xdrs,&tsp->eps) == 1);
    assert(xdr_float(xdrs,&tsp->phi) == 1);
    }

void read_tipsy_xdr_gas_dpp(XDR *xdrs, TIPSY_GAS_PARTICLE_DPP *tgpdpp) {

    int j;

    assert(xdr_float(xdrs,&tgpdpp->mass) == 1);
    for (j = 0; j < 3; j++) {
	assert(xdr_double(xdrs,&tgpdpp->pos[j]) == 1);
	}
    for (j = 0; j < 3; j++) {
	assert(xdr_float(xdrs,&tgpdpp->vel[j]) == 1);
	}
    assert(xdr_float(xdrs,&tgpdpp->rho) == 1);
    assert(xdr_float(xdrs,&tgpdpp->temp) == 1);
    assert(xdr_float(xdrs,&tgpdpp->hsmooth) == 1);
    assert(xdr_float(xdrs,&tgpdpp->metals) == 1);
    assert(xdr_float(xdrs,&tgpdpp->phi) == 1);
    }

void read_tipsy_xdr_dark_dpp(XDR *xdrs, TIPSY_DARK_PARTICLE_DPP *tdpdpp) {

    int j;

    assert(xdr_float(xdrs,&tdpdpp->mass) == 1);
    for (j = 0; j < 3; j++) {
	assert(xdr_double(xdrs,&tdpdpp->pos[j]) == 1);
	}
    for (j = 0; j < 3; j++) {
	assert(xdr_float(xdrs,&tdpdpp->vel[j]) == 1);
	}
    assert(xdr_float(xdrs,&tdpdpp->eps) == 1);
    assert(xdr_float(xdrs,&tdpdpp->phi) == 1);
    }

void read_tipsy_xdr_star_dpp(XDR *xdrs, TIPSY_STAR_PARTICLE_DPP *tspdpp) {

    int j;

    assert(xdr_float(xdrs,&tspdpp->mass) == 1);
    for (j = 0; j < 3; j++) {
	assert(xdr_double(xdrs,&tspdpp->pos[j]) == 1);
	}
    for (j = 0; j < 3; j++) {
	assert(xdr_float(xdrs,&tspdpp->vel[j]) == 1);
	}
    assert(xdr_float(xdrs,&tspdpp->metals) == 1);
    assert(xdr_float(xdrs,&tspdpp->tform) == 1);
    assert(xdr_float(xdrs,&tspdpp->eps) == 1);
    assert(xdr_float(xdrs,&tspdpp->phi) == 1);
    }

void read_tipsy_xdr(FILE *fp, TIPSY_STRUCTURE *ts) {

    int i;

    XDR xdrs;
    TIPSY_HEADER *th;
    TIPSY_GAS_PARTICLE *tgp;
    TIPSY_DARK_PARTICLE *tdp;
    TIPSY_STAR_PARTICLE *tsp;

    th = ts->th;
    tgp = ts->tgp;
    tdp = ts->tdp;
    tsp = ts->tsp;
    /*
    ** Read in header
    */
    th = realloc(th,sizeof(TIPSY_HEADER));
    assert(th != NULL);
    xdrstdio_create(&xdrs,fp,XDR_DECODE);
    read_tipsy_xdr_header(&xdrs,th);
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
    ** Read in gas particles
    */
    for (i = 0; i < th->ngas; i++) {
	read_tipsy_xdr_gas(&xdrs,&tgp[i]);
	}
    /*
    ** Read in dark matter particles
    */
    for (i = 0; i < th->ndark; i++) {
	read_tipsy_xdr_dark(&xdrs,&tdp[i]);
	}
    /*
    ** Read in star particles
    */
    for (i = 0; i < th->nstar; i++) {
	read_tipsy_xdr_star(&xdrs,&tsp[i]);
	}
    xdr_destroy(&xdrs);
    /*
    ** Return pointers
    */
    ts->th = th;
    ts->tgp = tgp;
    ts->tdp = tdp;
    ts->tsp = tsp;
    }

void read_tipsy_xdr_dpp(FILE *fp, TIPSY_STRUCTURE_DPP *tsdpp) {

    int i;

    XDR xdrs;
    TIPSY_HEADER *th;
    TIPSY_GAS_PARTICLE_DPP *tgpdpp;
    TIPSY_DARK_PARTICLE_DPP *tdpdpp;
    TIPSY_STAR_PARTICLE_DPP *tspdpp;

    th = tsdpp->th;
    tgpdpp = tsdpp->tgpdpp;
    tdpdpp = tsdpp->tdpdpp;
    tspdpp = tsdpp->tspdpp;
    /*
    ** Read in header
    */
    th = realloc(th,sizeof(TIPSY_HEADER));
    assert(th != NULL);
    xdrstdio_create(&xdrs,fp,XDR_DECODE);
    read_tipsy_xdr_header(&xdrs,th);
    /*
    ** Allocate memory
    */
    if(th->ngas != 0) {
	tgpdpp = realloc(tgpdpp,th->ngas*sizeof(TIPSY_GAS_PARTICLE_DPP));
	assert(tgpdpp != NULL);
	}
    else {
	tgpdpp = NULL;
	}
    if(th->ndark != 0) {
	tdpdpp = realloc(tdpdpp,th->ndark*sizeof(TIPSY_DARK_PARTICLE_DPP));
	assert(tdpdpp != NULL);
	}
    else {
	tdpdpp = NULL;
	}
    if(th->nstar != 0) {
	tspdpp = realloc(tspdpp,th->nstar*sizeof(TIPSY_STAR_PARTICLE_DPP));
	assert(tspdpp != NULL);
	}
    else {
	tspdpp = NULL;
	}
    /*
    ** Read in gas particles
    */
    for (i = 0; i < th->ngas; i++) {
	read_tipsy_xdr_gas_dpp(&xdrs,&tgpdpp[i]);
	}
    /*
    ** Read in dark matter particles
    */
    for (i = 0; i < th->ndark; i++) {
	read_tipsy_xdr_dark_dpp(&xdrs,&tdpdpp[i]);
	}
    /*
    ** Read in star particles
    */
    for (i = 0; i < th->nstar; i++) {
	read_tipsy_xdr_star_dpp(&xdrs,&tspdpp[i]);
	}
    xdr_destroy(&xdrs);
    /*
    ** Return pointers
    */
    tsdpp->th = th;
    tsdpp->tgpdpp = tgpdpp;
    tsdpp->tdpdpp = tdpdpp;
    tsdpp->tspdpp = tspdpp;
    }

void write_tipsy_xdr_header(XDR *xdrs, TIPSY_HEADER *th) {

    int pad = 0;

    assert(xdr_double(xdrs,&th->time) == 1);
    assert(xdr_int(xdrs,&th->ntotal) == 1);
    assert(xdr_int(xdrs,&th->ndim) == 1);
    assert(xdr_int(xdrs,&th->ngas) == 1);
    assert(xdr_int(xdrs,&th->ndark) == 1);
    assert(xdr_int(xdrs,&th->nstar) == 1);
    assert(xdr_int(xdrs,&pad) == 1);
    }

void write_tipsy_xdr_gas(XDR *xdrs, TIPSY_GAS_PARTICLE *tgp) {

    int j;

    assert(xdr_float(xdrs,&tgp->mass) == 1);
    for (j = 0; j < 3; j++) {
	assert(xdr_float(xdrs,&tgp->pos[j]) == 1);
	}
    for (j = 0; j < 3; j++) {
	assert(xdr_float(xdrs,&tgp->vel[j]) == 1);
	}
    assert(xdr_float(xdrs,&tgp->rho) == 1);
    assert(xdr_float(xdrs,&tgp->temp) == 1);
    assert(xdr_float(xdrs,&tgp->hsmooth) == 1);
    assert(xdr_float(xdrs,&tgp->metals) == 1);
    assert(xdr_float(xdrs,&tgp->phi) == 1);
    }

void write_tipsy_xdr_dark(XDR *xdrs, TIPSY_DARK_PARTICLE *tdp) {

    int j;

    assert(xdr_float(xdrs,&tdp->mass) == 1);
    for (j = 0; j < 3; j++) {
	assert(xdr_float(xdrs,&tdp->pos[j]) == 1);
	}
    for (j = 0; j < 3; j++) {
	assert(xdr_float(xdrs,&tdp->vel[j]) == 1);
	}
    assert(xdr_float(xdrs,&tdp->eps) == 1);
    assert(xdr_float(xdrs,&tdp->phi) == 1);
    }

void write_tipsy_xdr_star(XDR *xdrs, TIPSY_STAR_PARTICLE *tsp) {

    int j;

    assert(xdr_float(xdrs,&tsp->mass) == 1);
    for (j = 0; j < 3; j++) {
	assert(xdr_float(xdrs,&tsp->pos[j]) == 1);
	}
    for (j = 0; j < 3; j++) {
	assert(xdr_float(xdrs,&tsp->vel[j]) == 1);
	}
    assert(xdr_float(xdrs,&tsp->metals) == 1);
    assert(xdr_float(xdrs,&tsp->tform) == 1);
    assert(xdr_float(xdrs,&tsp->eps) == 1);
    assert(xdr_float(xdrs,&tsp->phi) == 1);
    }

void write_tipsy_xdr_gas_dpp(XDR *xdrs, TIPSY_GAS_PARTICLE_DPP *tgpdpp) {

    int j;

    assert(xdr_float(xdrs,&tgpdpp->mass) == 1);
    for (j = 0; j < 3; j++) {
	assert(xdr_double(xdrs,&tgpdpp->pos[j]) == 1);
	}
    for (j = 0; j < 3; j++) {
	assert(xdr_float(xdrs,&tgpdpp->vel[j]) == 1);
	}
    assert(xdr_float(xdrs,&tgpdpp->rho) == 1);
    assert(xdr_float(xdrs,&tgpdpp->temp) == 1);
    assert(xdr_float(xdrs,&tgpdpp->hsmooth) == 1);
    assert(xdr_float(xdrs,&tgpdpp->metals) == 1);
    assert(xdr_float(xdrs,&tgpdpp->phi) == 1);
    }

void write_tipsy_xdr_dark_dpp(XDR *xdrs, TIPSY_DARK_PARTICLE_DPP *tdpdpp) {

    int j;

    assert(xdr_float(xdrs,&tdpdpp->mass) == 1);
    for (j = 0; j < 3; j++) {
	assert(xdr_double(xdrs,&tdpdpp->pos[j]) == 1);
	}
    for (j = 0; j < 3; j++) {
	assert(xdr_float(xdrs,&tdpdpp->vel[j]) == 1);
	}
    assert(xdr_float(xdrs,&tdpdpp->eps) == 1);
    assert(xdr_float(xdrs,&tdpdpp->phi) == 1);
    }

void write_tipsy_xdr_star_dpp(XDR *xdrs, TIPSY_STAR_PARTICLE_DPP *tspdpp) {

    int j;

    assert(xdr_float(xdrs,&tspdpp->mass) == 1);
    for (j = 0; j < 3; j++) {
	assert(xdr_double(xdrs,&tspdpp->pos[j]) == 1);
	}
    for (j = 0; j < 3; j++) {
	assert(xdr_float(xdrs,&tspdpp->vel[j]) == 1);
	}
    assert(xdr_float(xdrs,&tspdpp->metals) == 1);
    assert(xdr_float(xdrs,&tspdpp->tform) == 1);
    assert(xdr_float(xdrs,&tspdpp->eps) == 1);
    assert(xdr_float(xdrs,&tspdpp->phi) == 1);
    }

void write_tipsy_xdr(FILE *fp, const TIPSY_STRUCTURE *ts) {

    int i;

    XDR xdrs;
    TIPSY_HEADER *th;
    TIPSY_GAS_PARTICLE *tgp;
    TIPSY_DARK_PARTICLE *tdp;
    TIPSY_STAR_PARTICLE *tsp;

    th = ts->th;
    tgp = ts->tgp;
    tdp = ts->tdp;
    tsp = ts->tsp;
    /*
    ** Write out header
    */
    xdrstdio_create(&xdrs,fp,XDR_ENCODE);
    write_tipsy_xdr_header(&xdrs,th);
    /*
    ** Write gas particles
    */
    for (i = 0; i < th->ngas; i++) {
	write_tipsy_xdr_gas(&xdrs,&tgp[i]);
	}
    /*
    ** Write dark matter particles
    */
    for (i = 0; i < th->ndark; i++) {
	write_tipsy_xdr_dark(&xdrs,&tdp[i]);
	}
    /*
    ** Write star particles
    */
    for (i = 0; i < th->nstar; i++) {
	write_tipsy_xdr_star(&xdrs,&tsp[i]);
	}
    xdr_destroy(&xdrs);
    }

void write_tipsy_xdr_dpp(FILE *fp, const TIPSY_STRUCTURE_DPP *tsdpp) {

    int i;

    XDR xdrs;
    TIPSY_HEADER *th;
    TIPSY_GAS_PARTICLE_DPP *tgpdpp;
    TIPSY_DARK_PARTICLE_DPP *tdpdpp;
    TIPSY_STAR_PARTICLE_DPP *tspdpp;

    th = tsdpp->th;
    tgpdpp = tsdpp->tgpdpp;
    tdpdpp = tsdpp->tdpdpp;
    tspdpp = tsdpp->tspdpp;
    /*
    ** Write out header
    */
    xdrstdio_create(&xdrs,fp,XDR_ENCODE);
    write_tipsy_xdr_header(&xdrs,th);
    /*
    ** Write gas particles
    */
    for (i = 0; i < th->ngas; i++) {
	write_tipsy_xdr_gas_dpp(&xdrs,&tgpdpp[i]);
	}
    /*
    ** Write dark matter particles
    */
    for (i = 0; i < th->ndark; i++) {
	write_tipsy_xdr_dark_dpp(&xdrs,&tdpdpp[i]);
	}
    /*
    ** Write star particles
    */
    for (i = 0; i < th->nstar; i++) {
	write_tipsy_xdr_star_dpp(&xdrs,&tspdpp[i]);
	}
    xdr_destroy(&xdrs);
    }

/*
** Tipsy ascii functions
*/

void read_tipsy_ascii(FILE *fp, TIPSY_STRUCTURE *ts) {

    int i, j;
    TIPSY_HEADER *th;
    TIPSY_GAS_PARTICLE *tgp;
    TIPSY_DARK_PARTICLE *tdp;
    TIPSY_STAR_PARTICLE *tsp;
    
    th = ts->th;
    tgp = ts->tgp;
    tdp = ts->tdp;
    tsp = ts->tsp;
    /*
    ** Read in header
    */
    th = realloc(th,sizeof(TIPSY_HEADER));
    assert(th != NULL);
    assert(fscanf(fp,"%d %d %d",&th->ntotal,&th->ngas,&th->nstar) == 3);
    assert(fscanf(fp,"%d",&th->ndim) == 1);
    assert(fscanf(fp,"%lf",&th->time) == 1);
    th->ndark = th->ntotal - th->ngas - th->nstar;
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
    ** Read in masses
    */
    for (i = 0; i < th->ngas; i++) {
	assert(fscanf(fp,"%f",&tgp[i].mass) == 1);
	}
    for (i = 0; i < th->ndark; i++) {
	assert(fscanf(fp,"%f",&tdp[i].mass) == 1);
	}
    for (i = 0; i < th->nstar; i++) {
	assert(fscanf(fp,"%f",&tsp[i].mass) == 1);
	}
    /*
    ** Read in positions
    */
    for (j = 0; j < 3; j++) {
	for (i = 0; i < th->ngas; i++) {
	    assert(fscanf(fp,"%f",&tgp[i].pos[j]) == 1);
	    }
	for (i = 0; i < th->ndark; i++) {
	    assert(fscanf(fp,"%f",&tdp[i].pos[j]) == 1);
	    }
	for (i = 0; i < th->nstar; i++) {
	    assert(fscanf(fp,"%f",&tsp[i].pos[j]) == 1);
	    }
	}
    /*
    ** Read in velocities
    */
    for (j = 0; j < 3; j++) {
	for (i = 0; i < th->ngas; i++) {
	    assert(fscanf(fp,"%f",&tgp[i].vel[j]) == 1);
	    }
	for (i = 0; i < th->ndark; i++) {
	    assert(fscanf(fp,"%f",&tdp[i].vel[j]) == 1);
	    }
	for (i = 0; i < th->nstar; i++) {
	    assert(fscanf(fp,"%f",&tsp[i].vel[j]) == 1);
	    }
	}
    /*
    ** Read in softenings of dark matter and star particles
    */
    for (i = 0; i < th->ndark; i++) {
	assert(fscanf(fp,"%f",&tdp[i].eps) == 1);
	}
    for (i = 0; i < th->nstar; i++) {
	assert(fscanf(fp,"%f",&tsp[i].eps) == 1);
	}
    /*
    ** Read in gas specific stuff
    */
    for (i = 0; i < th->ngas; i++) {
	assert(fscanf(fp,"%f",&tgp[i].rho) == 1);
	}
    for (i = 0; i < th->ngas; i++) {
	assert(fscanf(fp,"%f",&tgp[i].temp) == 1);
	}
    for (i = 0; i < th->ngas; i++) {
	assert(fscanf(fp,"%f",&tgp[i].hsmooth) == 1);
	}
    for (i = 0; i < th->ngas; i++) {
	assert(fscanf(fp,"%f",&tgp[i].metals) == 1);
	}
    /*
    ** Read in star specific stuff
    */
    for (i = 0; i < th->nstar; i++) {
	assert(fscanf(fp,"%f",&tsp[i].metals) == 1);
	}
    for (i = 0; i < th->nstar; i++) {
	assert(fscanf(fp,"%f",&tsp[i].tform) == 1);
	}
    /*
    ** Read in potentials
    */
    for (i = 0; i < th->ngas; i++) {
	assert(fscanf(fp,"%f",&tgp[i].phi) == 1);
	}
    for (i = 0; i < th->ndark; i++) {
	assert(fscanf(fp,"%f",&tdp[i].phi) == 1);
	}
    for (i = 0; i < th->nstar; i++) {
	assert(fscanf(fp,"%f",&tsp[i].phi) == 1);
	}
    /*
    ** Return pointers
    */
    ts->th = th;
    ts->tgp = tgp;
    ts->tdp = tdp;
    ts->tsp = tsp;
    }

void read_tipsy_ascii_dpp(FILE *fp, TIPSY_STRUCTURE_DPP *tsdpp) {

    int i, j;
    TIPSY_HEADER *th;
    TIPSY_GAS_PARTICLE_DPP *tgpdpp;
    TIPSY_DARK_PARTICLE_DPP *tdpdpp;
    TIPSY_STAR_PARTICLE_DPP *tspdpp;
    
    th = tsdpp->th;
    tgpdpp = tsdpp->tgpdpp;
    tdpdpp = tsdpp->tdpdpp;
    tspdpp = tsdpp->tspdpp;
    /*
    ** Read in header
    */
    th = realloc(th,sizeof(TIPSY_HEADER));
    assert(th != NULL);
    assert(fscanf(fp,"%d %d %d",&th->ntotal,&th->ngas,&th->nstar) == 3);
    assert(fscanf(fp,"%d",&th->ndim) == 1);
    assert(fscanf(fp,"%lf",&th->time) == 1);
    th->ndark = th->ntotal - th->ngas - th->nstar;
    /*
    ** Allocate memory
    */
    if(th->ngas != 0) {
	tgpdpp = realloc(tgpdpp,th->ngas*sizeof(TIPSY_GAS_PARTICLE_DPP));
	assert(tgpdpp != NULL);
	}
    else {
	tgpdpp = NULL;
	}
    if(th->ndark != 0) {
	tdpdpp = realloc(tdpdpp,th->ndark*sizeof(TIPSY_DARK_PARTICLE_DPP));
	assert(tdpdpp != NULL);
	}
    else {
	tdpdpp = NULL;
	}
    if(th->nstar != 0) {
	tspdpp = realloc(tspdpp,th->nstar*sizeof(TIPSY_STAR_PARTICLE_DPP));
	assert(tspdpp != NULL);
	}
    else {
	tspdpp = NULL;
	}
    /*
    ** Read in masses
    */
    for (i = 0; i < th->ngas; i++) {
	assert(fscanf(fp,"%f",&tgpdpp[i].mass) == 1);
	}
    for (i = 0; i < th->ndark; i++) {
	assert(fscanf(fp,"%f",&tdpdpp[i].mass) == 1);
	}
    for (i = 0; i < th->nstar; i++) {
	assert(fscanf(fp,"%f",&tspdpp[i].mass) == 1);
	}
    /*
    ** Read in positions
    */
    for (j = 0; j < 3; j++) {
	for (i = 0; i < th->ngas; i++) {
	    assert(fscanf(fp,"%lf",&tgpdpp[i].pos[j]) == 1);
	    }
	for (i = 0; i < th->ndark; i++) {
	    assert(fscanf(fp,"%lf",&tdpdpp[i].pos[j]) == 1);
	    }
	for (i = 0; i < th->nstar; i++) {
	    assert(fscanf(fp,"%lf",&tspdpp[i].pos[j]) == 1);
	    }
	}
    /*
    ** Read in velocities
    */
    for (j = 0; j < 3; j++) {
	for (i = 0; i < th->ngas; i++) {
	    assert(fscanf(fp,"%f",&tgpdpp[i].vel[j]) == 1);
	    }
	for (i = 0; i < th->ndark; i++) {
	    assert(fscanf(fp,"%f",&tdpdpp[i].vel[j]) == 1);
	    }
	for (i = 0; i < th->nstar; i++) {
	    assert(fscanf(fp,"%f",&tspdpp[i].vel[j]) == 1);
	    }
	}
    /*
    ** Read in softenings of dark matter and star particles
    */
    for (i = 0; i < th->ndark; i++) {
	assert(fscanf(fp,"%f",&tdpdpp[i].eps) == 1);
	}
    for (i = 0; i < th->nstar; i++) {
	assert(fscanf(fp,"%f",&tspdpp[i].eps) == 1);
	}
    /*
    ** Read in gas specific stuff
    */
    for (i = 0; i < th->ngas; i++) {
	assert(fscanf(fp,"%f",&tgpdpp[i].rho) == 1);
	}
    for (i = 0; i < th->ngas; i++) {
	assert(fscanf(fp,"%f",&tgpdpp[i].temp) == 1);
	}
    for (i = 0; i < th->ngas; i++) {
	assert(fscanf(fp,"%f",&tgpdpp[i].hsmooth) == 1);
	}
    for (i = 0; i < th->ngas; i++) {
	assert(fscanf(fp,"%f",&tgpdpp[i].metals) == 1);
	}
    /*
    ** Read in star specific stuff
    */
    for (i = 0; i < th->nstar; i++) {
	assert(fscanf(fp,"%f",&tspdpp[i].metals) == 1);
	}
    for (i = 0; i < th->nstar; i++) {
	assert(fscanf(fp,"%f",&tspdpp[i].tform) == 1);
	}
    /*
    ** Read in potentials
    */
    for (i = 0; i < th->ngas; i++) {
	assert(fscanf(fp,"%f",&tgpdpp[i].phi) == 1);
	}
    for (i = 0; i < th->ndark; i++) {
	assert(fscanf(fp,"%f",&tdpdpp[i].phi) == 1);
	}
    for (i = 0; i < th->nstar; i++) {
	assert(fscanf(fp,"%f",&tspdpp[i].phi) == 1);
	}
    /*
    ** Return pointers
    */
    tsdpp->th = th;
    tsdpp->tgpdpp = tgpdpp;
    tsdpp->tdpdpp = tdpdpp;
    tsdpp->tspdpp = tspdpp;
    }

void write_tipsy_ascii(FILE *fp, const TIPSY_STRUCTURE *ts) {

    int i, j;
    TIPSY_HEADER *th;
    TIPSY_GAS_PARTICLE *tgp;
    TIPSY_DARK_PARTICLE *tdp;
    TIPSY_STAR_PARTICLE *tsp;
    
    th = ts->th;
    tgp = ts->tgp;
    tdp = ts->tdp;
    tsp = ts->tsp;
    /*
    ** Write out header
    */
    assert(fprintf(fp,"%d %d %d\n",th->ntotal,th->ngas,th->nstar) > 0);
    assert(fprintf(fp,"%d\n",3) > 0);
    assert(fprintf(fp,"%.6e\n",th->time) > 0);
    /*
    ** Write out masses
    */
    for (i = 0; i < th->ngas; i++) {
	assert(fprintf(fp,"%.6e\n",tgp[i].mass) > 0);
	}
    for (i = 0; i < th->ndark; i++) {
	assert(fprintf(fp,"%.6e\n",tdp[i].mass) > 0);
	}
    for (i = 0; i < th->nstar; i++) {
	assert(fprintf(fp,"%.6e\n",tsp[i].mass) > 0);
	}
    /*
    ** Write out positions
    */
    for (j = 0; j < 3; j++) {
	for (i = 0; i < th->ngas; i++) {
	    assert(fprintf(fp,"%.6e\n",tgp[i].pos[j]) > 0);
	    }
	for (i = 0; i < th->ndark; i++) {
	    assert(fprintf(fp,"%.6e\n",tdp[i].pos[j]) > 0);
	    }
	for (i = 0; i < th->nstar; i++) {
	    assert(fprintf(fp,"%.6e\n",tsp[i].pos[j]) > 0);
	    }
	}
    /*
    ** Write out velocities
    */
    for (j = 0; j < 3; j++) {
	for (i = 0; i < th->ngas; i++) {
	    assert(fprintf(fp,"%.6e\n",tgp[i].vel[j]) > 0);
	    }
	for (i = 0; i < th->ndark; i++) {
	    assert(fprintf(fp,"%.6e\n",tdp[i].vel[j]) > 0);
	    }
	for (i = 0; i < th->nstar; i++) {
	    assert(fprintf(fp,"%.6e\n",tsp[i].vel[j]) > 0);
	    }
	}
     /*
    ** Write out softenings for dark matter and star particles
    */
    for (i = 0; i < th->ndark; i++) {
	assert(fprintf(fp,"%.6e\n",tdp[i].eps) > 0);
	}
    for (i = 0; i < th->nstar; i++) {
	assert(fprintf(fp,"%.6e\n",tsp[i].eps) > 0);
	}
    /*
    ** Write out gas specific stuff
    */
    for (i = 0; i < th->ngas; i++) {
	assert(fprintf(fp,"%.6e\n",tgp[i].rho) > 0);
	}
    for (i = 0; i < th->ngas; i++) {
	assert(fprintf(fp,"%.6e\n",tgp[i].temp) > 0);
	}
    for (i = 0; i < th->ngas; i++) {
	assert(fprintf(fp,"%.6e\n",tgp[i].hsmooth) > 0);
	}
    for (i = 0; i < th->ngas; i++) {
	assert(fprintf(fp,"%.6e\n",tgp[i].metals) > 0);
	}
    /*
    ** Write out star specific stuff
    */
    for (i = 0; i < th->nstar; i++) {
	assert(fprintf(fp,"%.6e\n",tsp[i].metals) > 0);
	}
    for (i = 0; i < th->nstar; i++) {
	assert(fprintf(fp,"%.6e\n",tsp[i].tform) > 0);
	}
    /*
    ** Write out potentials
    */
    for (i = 0; i < th->ngas; i++) {
	assert(fprintf(fp,"%.6e\n",tgp[i].phi) > 0);
	}
    for (i = 0; i < th->ndark; i++) {
	assert(fprintf(fp,"%.6e\n",tdp[i].phi) > 0);
	}
    for (i = 0; i < th->nstar; i++) {
	assert(fprintf(fp,"%.6e\n",tsp[i].phi) > 0);
	}
    }

void write_tipsy_ascii_dpp(FILE *fp, const TIPSY_STRUCTURE_DPP *tsdpp) {

    int i, j;
    TIPSY_HEADER *th;
    TIPSY_GAS_PARTICLE_DPP *tgpdpp;
    TIPSY_DARK_PARTICLE_DPP *tdpdpp;
    TIPSY_STAR_PARTICLE_DPP *tspdpp;
    
    th = tsdpp->th;
    tgpdpp = tsdpp->tgpdpp;
    tdpdpp = tsdpp->tdpdpp;
    tspdpp = tsdpp->tspdpp;
    /*
    ** Write out header
    */
    assert(fprintf(fp,"%d %d %d\n",th->ntotal,th->ngas,th->nstar) > 0);
    assert(fprintf(fp,"%d\n",3) > 0);
    assert(fprintf(fp,"%.6e\n",th->time) > 0);
    /*
    ** Write out masses
    */
    for (i = 0; i < th->ngas; i++) {
	assert(fprintf(fp,"%.6e\n",tgpdpp[i].mass) > 0);
	}
    for (i = 0; i < th->ndark; i++) {
	assert(fprintf(fp,"%.6e\n",tdpdpp[i].mass) > 0);
	}
    for (i = 0; i < th->nstar; i++) {
	assert(fprintf(fp,"%.6e\n",tspdpp[i].mass) > 0);
	}
    /*
    ** Write out positions
    */
    for (j = 0; j < 3; j++) {
	for (i = 0; i < th->ngas; i++) {
	    assert(fprintf(fp,"%.14e\n",tgpdpp[i].pos[j]) > 0);
	    }
	for (i = 0; i < th->ndark; i++) {
	    assert(fprintf(fp,"%.14e\n",tdpdpp[i].pos[j]) > 0);
	    }
	for (i = 0; i < th->nstar; i++) {
	    assert(fprintf(fp,"%.14e\n",tspdpp[i].pos[j]) > 0);
	    }
	}
    /*
    ** Write out velocities
    */
    for (j = 0; j < 3; j++) {
	for (i = 0; i < th->ngas; i++) {
	    assert(fprintf(fp,"%.6e\n",tgpdpp[i].vel[j]) > 0);
	    }
	for (i = 0; i < th->ndark; i++) {
	    assert(fprintf(fp,"%.6e\n",tdpdpp[i].vel[j]) > 0);
	    }
	for (i = 0; i < th->nstar; i++) {
	    assert(fprintf(fp,"%.6e\n",tspdpp[i].vel[j]) > 0);
	    }
	}
     /*
    ** Write out softenings for dark matter and star particles
    */
    for (i = 0; i < th->ndark; i++) {
	assert(fprintf(fp,"%.6e\n",tdpdpp[i].eps) > 0);
	}
    for (i = 0; i < th->nstar; i++) {
	assert(fprintf(fp,"%.6e\n",tspdpp[i].eps) > 0);
	}
    /*
    ** Write out gas specific stuff
    */
    for (i = 0; i < th->ngas; i++) {
	assert(fprintf(fp,"%.6e\n",tgpdpp[i].rho) > 0);
	}
    for (i = 0; i < th->ngas; i++) {
	assert(fprintf(fp,"%.6e\n",tgpdpp[i].temp) > 0);
	}
    for (i = 0; i < th->ngas; i++) {
	assert(fprintf(fp,"%.6e\n",tgpdpp[i].hsmooth) > 0);
	}
    for (i = 0; i < th->ngas; i++) {
	assert(fprintf(fp,"%.6e\n",tgpdpp[i].metals) > 0);
	}
    /*
    ** Write out star specific stuff
    */
    for (i = 0; i < th->nstar; i++) {
	assert(fprintf(fp,"%.6e\n",tspdpp[i].metals) > 0);
	}
    for (i = 0; i < th->nstar; i++) {
	assert(fprintf(fp,"%.6e\n",tspdpp[i].tform) > 0);
	}
    /*
    ** Write out potentials
    */
    for (i = 0; i < th->ngas; i++) {
	assert(fprintf(fp,"%.6e\n",tgpdpp[i].phi) > 0);
	}
    for (i = 0; i < th->ndark; i++) {
	assert(fprintf(fp,"%.6e\n",tdpdpp[i].phi) > 0);
	}
    for (i = 0; i < th->nstar; i++) {
	assert(fprintf(fp,"%.6e\n",tspdpp[i].phi) > 0);
	}
    }

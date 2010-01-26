/*
** iof.c
**
** Various reading & writing functions
**
** written by Marcel Zemp
*/

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <malloc.h>
#include <assert.h>
#include "iof.h"

/*
** Copy functions
*/

void copy_gp_to_gpdpp(const GAS_PARTICLE *gp, GAS_PARTICLE_DPP *gpdpp) {

    int j;

    gpdpp->mass = gp->mass;
    for (j = 0; j < 3; j++) {
	gpdpp->pos[j] = gp->pos[j];
	}
    for (j = 0; j < 3; j++) {
	gpdpp->vel[j] = gp->vel[j];
	}
    gpdpp->rho = gp->rho;
    gpdpp->temp = gp->temp;
    gpdpp->hsmooth = gp->hsmooth;
    gpdpp->metals = gp->metals;
    gpdpp->phi = gp->phi;
    }

void copy_dp_to_dpdpp(const DARK_PARTICLE *dp, DARK_PARTICLE_DPP *dpdpp) {

    int j;

    dpdpp->mass = dp->mass;
    for (j = 0; j < 3; j++) {
	dpdpp->pos[j] = dp->pos[j];
	}
    for (j = 0; j < 3; j++) {
	dpdpp->vel[j] = dp->vel[j];
	}
    dpdpp->eps = dp->eps;
    dpdpp->phi = dp->phi;
    }

void copy_sp_to_spdpp(const STAR_PARTICLE *sp, STAR_PARTICLE_DPP *spdpp) {

    int j;

    spdpp->mass = sp->mass;
    for (j = 0; j < 3; j++) {
	spdpp->pos[j] = sp->pos[j];
	}
    for (j = 0; j < 3; j++) {
	spdpp->vel[j] = sp->vel[j];
	}
    spdpp->metals = sp->metals;
    spdpp->tform = sp->tform;
    spdpp->eps = sp->eps;
    spdpp->phi = sp->phi;
    }

void copy_gpdpp_to_gp(const GAS_PARTICLE_DPP *gpdpp, GAS_PARTICLE *gp) {

    int j;

    gp->mass = gpdpp->mass;
    for (j = 0; j < 3; j++) {
	gp->pos[j] = gpdpp->pos[j];
	}
    for (j = 0; j < 3; j++) {
	gp->vel[j] = gpdpp->vel[j];
	}
    gp->rho = gpdpp->rho;
    gp->temp = gpdpp->temp;
    gp->hsmooth = gpdpp->hsmooth;
    gp->metals = gpdpp->metals;
    gp->phi = gpdpp->phi;
    }

void copy_dpdpp_to_dp(const DARK_PARTICLE_DPP *dpdpp, DARK_PARTICLE *dp) {

    int j;

    dp->mass = dpdpp->mass;
    for (j = 0; j < 3; j++) {
	dp->pos[j] = dpdpp->pos[j];
	}
    for (j = 0; j < 3; j++) {
	dp->vel[j] = dpdpp->vel[j];
	}
    dp->eps = dpdpp->eps;
    dp->phi = dpdpp->phi;
    }

void copy_spdpp_to_sp(const STAR_PARTICLE_DPP *spdpp, STAR_PARTICLE *sp) {

    int j;

    sp->mass = spdpp->mass;
    for (j = 0; j < 3; j++) {
	sp->pos[j] = spdpp->pos[j];
	}
    for (j = 0; j < 3; j++) {
	sp->vel[j] = spdpp->vel[j];
	}
    sp->metals = spdpp->metals;
    sp->tform = spdpp->tform;
    sp->eps = spdpp->eps;
    sp->phi = spdpp->phi;
    }

/*
** Tipsy binary functions
*/

void read_tipsy_binary_header(FILE *fp, TIPSY_HEADER *th) {

    assert(fread(th,sizeof(TIPSY_HEADER),1,fp) == 1);
    }

void read_tipsy_binary_gas(FILE *fp, GAS_PARTICLE *gp) {

    assert(fread(gp,sizeof(GAS_PARTICLE),1,fp) == 1);
    }

void read_tipsy_binary_dark(FILE *fp, DARK_PARTICLE *dp) {

    assert(fread(dp,sizeof(DARK_PARTICLE),1,fp) == 1);
    }

void read_tipsy_binary_star(FILE *fp, STAR_PARTICLE *sp) {

    assert(fread(sp,sizeof(STAR_PARTICLE),1,fp) == 1);
    }

void read_tipsy_binary_gas_dpp(FILE *fp, GAS_PARTICLE_DPP *gpdpp) {

    assert(fread(gpdpp,sizeof(GAS_PARTICLE_DPP),1,fp) == 1);
    }

void read_tipsy_binary_dark_dpp(FILE *fp, DARK_PARTICLE_DPP *dpdpp) {

    assert(fread(dpdpp,sizeof(DARK_PARTICLE_DPP),1,fp) == 1);
    }

void read_tipsy_binary_star_dpp(FILE *fp, STAR_PARTICLE_DPP *spdpp) {

    assert(fread(spdpp,sizeof(STAR_PARTICLE_DPP),1,fp) == 1);
    }

void read_tipsy_binary(FILE *fp, TIPSY_STRUCTURE *ts) {

    TIPSY_HEADER *th;
    GAS_PARTICLE *gp;
    DARK_PARTICLE *dp;
    STAR_PARTICLE *sp;

    th = ts->th;
    gp = ts->gp;
    dp = ts->dp;
    sp = ts->sp;
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
	gp = realloc(gp,th->ngas*sizeof(GAS_PARTICLE));
	assert(gp != NULL);
	}
    else {
	gp = NULL;
	}
    if(th->ndark != 0) {
	dp = realloc(dp,th->ndark*sizeof(DARK_PARTICLE));
	assert(dp != NULL);
	}
    else {
	dp = NULL;
	}
    if(th->nstar != 0) {
	sp = realloc(sp,th->nstar*sizeof(STAR_PARTICLE));
	assert(sp != NULL);
	}
    else {
	sp = NULL;
	}
    /*
    ** Read in particles
    */
    assert(fread(gp,sizeof(GAS_PARTICLE),th->ngas,fp) == th->ngas);
    assert(fread(dp,sizeof(DARK_PARTICLE),th->ndark,fp) == th->ndark);
    assert(fread(sp,sizeof(STAR_PARTICLE),th->nstar,fp) == th->nstar);
    /*
    ** Return pointers
    */
    ts->th = th;
    ts->gp = gp;
    ts->dp = dp;
    ts->sp = sp;
    }

void read_tipsy_binary_dpp(FILE *fp, TIPSY_STRUCTURE_DPP *tsdpp) {

    TIPSY_HEADER *th;
    GAS_PARTICLE_DPP *gpdpp;
    DARK_PARTICLE_DPP *dpdpp;
    STAR_PARTICLE_DPP *spdpp;

    th = tsdpp->th;
    gpdpp = tsdpp->gpdpp;
    dpdpp = tsdpp->dpdpp;
    spdpp = tsdpp->spdpp;
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
	gpdpp = realloc(gpdpp,th->ngas*sizeof(GAS_PARTICLE_DPP));
	assert(gpdpp != NULL);
	}
    else {
	gpdpp = NULL;
	}
    if(th->ndark != 0) {
	dpdpp = realloc(dpdpp,th->ndark*sizeof(DARK_PARTICLE_DPP));
	assert(dpdpp != NULL);
	}
    else {
	dpdpp = NULL;
	}
    if(th->nstar != 0) {
	spdpp = realloc(spdpp,th->nstar*sizeof(STAR_PARTICLE_DPP));
	assert(spdpp != NULL);
	}
    else {
	spdpp = NULL;
	}
    /*
    ** Read in particles
    */
    assert(fread(gpdpp,sizeof(GAS_PARTICLE_DPP),th->ngas,fp) == th->ngas);
    assert(fread(dpdpp,sizeof(DARK_PARTICLE_DPP),th->ndark,fp) == th->ndark);
    assert(fread(spdpp,sizeof(STAR_PARTICLE_DPP),th->nstar,fp) == th->nstar);
    /*
    ** Return pointers
    */
    tsdpp->th = th;
    tsdpp->gpdpp = gpdpp;
    tsdpp->dpdpp = dpdpp;
    tsdpp->spdpp = spdpp;
    }

void write_tipsy_binary_header(FILE *fp, const TIPSY_HEADER *th) {

    assert(fwrite(th,sizeof(TIPSY_HEADER),1,fp) == 1);
    }

void write_tipsy_binary_gas(FILE *fp, const GAS_PARTICLE *gp) {

    assert(fwrite(gp,sizeof(GAS_PARTICLE),1,fp) == 1);
    }

void write_tipsy_binary_dark(FILE *fp, const DARK_PARTICLE *dp) {

    assert(fwrite(dp,sizeof(DARK_PARTICLE),1,fp) == 1);
    }

void write_tipsy_binary_star(FILE *fp, const STAR_PARTICLE *sp) {

    assert(fwrite(sp,sizeof(STAR_PARTICLE),1,fp) == 1);
    }

void write_tipsy_binary_gas_dpp(FILE *fp, const GAS_PARTICLE_DPP *gpdpp) {

    assert(fwrite(gpdpp,sizeof(GAS_PARTICLE_DPP),1,fp) == 1);
    }

void write_tipsy_binary_dark_dpp(FILE *fp, const DARK_PARTICLE_DPP *dpdpp) {

    assert(fwrite(dpdpp,sizeof(DARK_PARTICLE_DPP),1,fp) == 1);
    }

void write_tipsy_binary_star_dpp(FILE *fp, const STAR_PARTICLE_DPP *spdpp) {

    assert(fwrite(spdpp,sizeof(STAR_PARTICLE_DPP),1,fp) == 1);
    }

void write_tipsy_binary(FILE *fp, const TIPSY_STRUCTURE *ts) {

    TIPSY_HEADER *th;
    GAS_PARTICLE *gp;
    DARK_PARTICLE *dp;
    STAR_PARTICLE *sp;

    th = ts->th;
    gp = ts->gp;
    dp = ts->dp;
    sp = ts->sp;
    /*
    ** Write out header
    */
    assert(fwrite(th,sizeof(TIPSY_HEADER),1,fp) == 1);
    /*
    ** Write out particles
    */
    assert(fwrite(gp,sizeof(GAS_PARTICLE),th->ngas,fp) == th->ngas);
    assert(fwrite(dp,sizeof(DARK_PARTICLE),th->ndark,fp) == th->ndark);
    assert(fwrite(sp,sizeof(STAR_PARTICLE),th->nstar,fp) == th->nstar);
    }

void write_tipsy_binary_dpp(FILE *fp, const TIPSY_STRUCTURE_DPP *tsdpp) {

    TIPSY_HEADER *th;
    GAS_PARTICLE_DPP *gpdpp;
    DARK_PARTICLE_DPP *dpdpp;
    STAR_PARTICLE_DPP *spdpp;

    th = tsdpp->th;
    gpdpp = tsdpp->gpdpp;
    dpdpp = tsdpp->dpdpp;
    spdpp = tsdpp->spdpp;
    /*
    ** Write out header
    */
    assert(fwrite(th,sizeof(TIPSY_HEADER),1,fp) == 1);
    /*
    ** Write out particles
    */
    assert(fwrite(gpdpp,sizeof(GAS_PARTICLE_DPP),th->ngas,fp) == th->ngas);
    assert(fwrite(dpdpp,sizeof(DARK_PARTICLE_DPP),th->ndark,fp) == th->ndark);
    assert(fwrite(spdpp,sizeof(STAR_PARTICLE_DPP),th->nstar,fp) == th->nstar);
    }

/*
** Tipsy standard binary functions
*/

void read_tipsy_standard_header(XDR *xdrs, TIPSY_HEADER *th) {

    int pad;

    assert(xdr_double(xdrs,&th->time) == 1);
    assert(xdr_int(xdrs,&th->ntotal) == 1);
    assert(xdr_int(xdrs,&th->ndim) == 1);
    assert(xdr_int(xdrs,&th->ngas) == 1);
    assert(xdr_int(xdrs,&th->ndark) == 1);
    assert(xdr_int(xdrs,&th->nstar) == 1);
    assert(xdr_int(xdrs,&pad) == 1);
    }

void read_tipsy_standard_gas(XDR *xdrs, GAS_PARTICLE *gp) {

    int j;

    assert(xdr_float(xdrs,&gp->mass) == 1);
    for (j = 0; j < 3; j++) {
	assert(xdr_float(xdrs,&gp->pos[j]) == 1);
	}
    for (j = 0; j < 3; j++) {
	assert(xdr_float(xdrs,&gp->vel[j]) == 1);
	}
    assert(xdr_float(xdrs,&gp->rho) == 1);
    assert(xdr_float(xdrs,&gp->temp) == 1);
    assert(xdr_float(xdrs,&gp->hsmooth) == 1);
    assert(xdr_float(xdrs,&gp->metals) == 1);
    assert(xdr_float(xdrs,&gp->phi) == 1);
    }

void read_tipsy_standard_dark(XDR *xdrs, DARK_PARTICLE *dp) {

    int j;

    assert(xdr_float(xdrs,&dp->mass) == 1);
    for (j = 0; j < 3; j++) {
	assert(xdr_float(xdrs,&dp->pos[j]) == 1);
	}
    for (j = 0; j < 3; j++) {
	assert(xdr_float(xdrs,&dp->vel[j]) == 1);
	}
    assert(xdr_float(xdrs,&dp->eps) == 1);
    assert(xdr_float(xdrs,&dp->phi) == 1);
    }

void read_tipsy_standard_star(XDR *xdrs, STAR_PARTICLE *sp) {

    int j;

    assert(xdr_float(xdrs,&sp->mass) == 1);
    for (j = 0; j < 3; j++) {
	assert(xdr_float(xdrs,&sp->pos[j]) == 1);
	}
    for (j = 0; j < 3; j++) {
	assert(xdr_float(xdrs,&sp->vel[j]) == 1);
	}
    assert(xdr_float(xdrs,&sp->metals) == 1);
    assert(xdr_float(xdrs,&sp->tform) == 1);
    assert(xdr_float(xdrs,&sp->eps) == 1);
    assert(xdr_float(xdrs,&sp->phi) == 1);
    }

void read_tipsy_standard_gas_dpp(XDR *xdrs, GAS_PARTICLE_DPP *gpdpp) {

    int j;

    assert(xdr_float(xdrs,&gpdpp->mass) == 1);
    for (j = 0; j < 3; j++) {
	assert(xdr_double(xdrs,&gpdpp->pos[j]) == 1);
	}
    for (j = 0; j < 3; j++) {
	assert(xdr_float(xdrs,&gpdpp->vel[j]) == 1);
	}
    assert(xdr_float(xdrs,&gpdpp->rho) == 1);
    assert(xdr_float(xdrs,&gpdpp->temp) == 1);
    assert(xdr_float(xdrs,&gpdpp->hsmooth) == 1);
    assert(xdr_float(xdrs,&gpdpp->metals) == 1);
    assert(xdr_float(xdrs,&gpdpp->phi) == 1);
    }

void read_tipsy_standard_dark_dpp(XDR *xdrs, DARK_PARTICLE_DPP *dpdpp) {

    int j;

    assert(xdr_float(xdrs,&dpdpp->mass) == 1);
    for (j = 0; j < 3; j++) {
	assert(xdr_double(xdrs,&dpdpp->pos[j]) == 1);
	}
    for (j = 0; j < 3; j++) {
	assert(xdr_float(xdrs,&dpdpp->vel[j]) == 1);
	}
    assert(xdr_float(xdrs,&dpdpp->eps) == 1);
    assert(xdr_float(xdrs,&dpdpp->phi) == 1);
    }

void read_tipsy_standard_star_dpp(XDR *xdrs, STAR_PARTICLE_DPP *spdpp) {

    int j;

    assert(xdr_float(xdrs,&spdpp->mass) == 1);
    for (j = 0; j < 3; j++) {
	assert(xdr_double(xdrs,&spdpp->pos[j]) == 1);
	}
    for (j = 0; j < 3; j++) {
	assert(xdr_float(xdrs,&spdpp->vel[j]) == 1);
	}
    assert(xdr_float(xdrs,&spdpp->metals) == 1);
    assert(xdr_float(xdrs,&spdpp->tform) == 1);
    assert(xdr_float(xdrs,&spdpp->eps) == 1);
    assert(xdr_float(xdrs,&spdpp->phi) == 1);
    }

void read_tipsy_standard(FILE *fp, TIPSY_STRUCTURE *ts) {

    int i;

    XDR xdrs;
    TIPSY_HEADER *th;
    GAS_PARTICLE *gp;
    DARK_PARTICLE *dp;
    STAR_PARTICLE *sp;

    th = ts->th;
    gp = ts->gp;
    dp = ts->dp;
    sp = ts->sp;
    /*
    ** Read in header
    */
    th = realloc(th,sizeof(TIPSY_HEADER));
    assert(th != NULL);
    xdrstdio_create(&xdrs,fp,XDR_DECODE);
    read_tipsy_standard_header(&xdrs,th);
    /*
    ** Allocate memory
    */
    if(th->ngas != 0) {
	gp = realloc(gp,th->ngas*sizeof(GAS_PARTICLE));
	assert(gp != NULL);
	}
    else {
	gp = NULL;
	}
    if(th->ndark != 0) {
	dp = realloc(dp,th->ndark*sizeof(DARK_PARTICLE));
	assert(dp != NULL);
	}
    else {
	dp = NULL;
	}
    if(th->nstar != 0) {
	sp = realloc(sp,th->nstar*sizeof(STAR_PARTICLE));
	assert(sp != NULL);
	}
    else {
	sp = NULL;
	}
    /*
    ** Read in gas particles
    */
    for (i = 0; i < th->ngas; i++) {
	read_tipsy_standard_gas(&xdrs,&gp[i]);
	}
    /*
    ** Read in dark matter particles
    */
    for (i = 0; i < th->ndark; i++) {
	read_tipsy_standard_dark(&xdrs,&dp[i]);
	}
    /*
    ** Read in star particles
    */
    for (i = 0; i < th->nstar; i++) {
	read_tipsy_standard_star(&xdrs,&sp[i]);
	}
    xdr_destroy(&xdrs);
    /*
    ** Return pointers
    */
    ts->th = th;
    ts->gp = gp;
    ts->dp = dp;
    ts->sp = sp;
    }

void read_tipsy_standard_dpp(FILE *fp, TIPSY_STRUCTURE_DPP *tsdpp) {

    int i;

    XDR xdrs;
    TIPSY_HEADER *th;
    GAS_PARTICLE_DPP *gpdpp;
    DARK_PARTICLE_DPP *dpdpp;
    STAR_PARTICLE_DPP *spdpp;

    th = tsdpp->th;
    gpdpp = tsdpp->gpdpp;
    dpdpp = tsdpp->dpdpp;
    spdpp = tsdpp->spdpp;
    /*
    ** Read in header
    */
    th = realloc(th,sizeof(TIPSY_HEADER));
    assert(th != NULL);
    xdrstdio_create(&xdrs,fp,XDR_DECODE);
    read_tipsy_standard_header(&xdrs,th);
    /*
    ** Allocate memory
    */
    if(th->ngas != 0) {
	gpdpp = realloc(gpdpp,th->ngas*sizeof(GAS_PARTICLE_DPP));
	assert(gpdpp != NULL);
	}
    else {
	gpdpp = NULL;
	}
    if(th->ndark != 0) {
	dpdpp = realloc(dpdpp,th->ndark*sizeof(DARK_PARTICLE_DPP));
	assert(dpdpp != NULL);
	}
    else {
	dpdpp = NULL;
	}
    if(th->nstar != 0) {
	spdpp = realloc(spdpp,th->nstar*sizeof(STAR_PARTICLE_DPP));
	assert(spdpp != NULL);
	}
    else {
	spdpp = NULL;
	}
    /*
    ** Read in gas particles
    */
    for (i = 0; i < th->ngas; i++) {
	read_tipsy_standard_gas_dpp(&xdrs,&gpdpp[i]);
	}
    /*
    ** Read in dark matter particles
    */
    for (i = 0; i < th->ndark; i++) {
	read_tipsy_standard_dark_dpp(&xdrs,&dpdpp[i]);
	}
    /*
    ** Read in star particles
    */
    for (i = 0; i < th->nstar; i++) {
	read_tipsy_standard_star_dpp(&xdrs,&spdpp[i]);
	}
    xdr_destroy(&xdrs);
    /*
    ** Return pointers
    */
    tsdpp->th = th;
    tsdpp->gpdpp = gpdpp;
    tsdpp->dpdpp = dpdpp;
    tsdpp->spdpp = spdpp;
    }

void write_tipsy_standard_header(XDR *xdrs, TIPSY_HEADER *th) {

    int pad = 0;

    assert(xdr_double(xdrs,&th->time) == 1);
    assert(xdr_int(xdrs,&th->ntotal) == 1);
    assert(xdr_int(xdrs,&th->ndim) == 1);
    assert(xdr_int(xdrs,&th->ngas) == 1);
    assert(xdr_int(xdrs,&th->ndark) == 1);
    assert(xdr_int(xdrs,&th->nstar) == 1);
    assert(xdr_int(xdrs,&pad) == 1);
    }

void write_tipsy_standard_gas(XDR *xdrs, GAS_PARTICLE *gp) {

    int j;

    assert(xdr_float(xdrs,&gp->mass) == 1);
    for (j = 0; j < 3; j++) {
	assert(xdr_float(xdrs,&gp->pos[j]) == 1);
	}
    for (j = 0; j < 3; j++) {
	assert(xdr_float(xdrs,&gp->vel[j]) == 1);
	}
    assert(xdr_float(xdrs,&gp->rho) == 1);
    assert(xdr_float(xdrs,&gp->temp) == 1);
    assert(xdr_float(xdrs,&gp->hsmooth) == 1);
    assert(xdr_float(xdrs,&gp->metals) == 1);
    assert(xdr_float(xdrs,&gp->phi) == 1);
    }

void write_tipsy_standard_dark(XDR *xdrs, DARK_PARTICLE *dp) {

    int j;

    assert(xdr_float(xdrs,&dp->mass) == 1);
    for (j = 0; j < 3; j++) {
	assert(xdr_float(xdrs,&dp->pos[j]) == 1);
	}
    for (j = 0; j < 3; j++) {
	assert(xdr_float(xdrs,&dp->vel[j]) == 1);
	}
    assert(xdr_float(xdrs,&dp->eps) == 1);
    assert(xdr_float(xdrs,&dp->phi) == 1);
    }

void write_tipsy_standard_star(XDR *xdrs, STAR_PARTICLE *sp) {

    int j;

    assert(xdr_float(xdrs,&sp->mass) == 1);
    for (j = 0; j < 3; j++) {
	assert(xdr_float(xdrs,&sp->pos[j]) == 1);
	}
    for (j = 0; j < 3; j++) {
	assert(xdr_float(xdrs,&sp->vel[j]) == 1);
	}
    assert(xdr_float(xdrs,&sp->metals) == 1);
    assert(xdr_float(xdrs,&sp->tform) == 1);
    assert(xdr_float(xdrs,&sp->eps) == 1);
    assert(xdr_float(xdrs,&sp->phi) == 1);
    }

void write_tipsy_standard_gas_dpp(XDR *xdrs, GAS_PARTICLE_DPP *gpdpp) {

    int j;

    assert(xdr_float(xdrs,&gpdpp->mass) == 1);
    for (j = 0; j < 3; j++) {
	assert(xdr_double(xdrs,&gpdpp->pos[j]) == 1);
	}
    for (j = 0; j < 3; j++) {
	assert(xdr_float(xdrs,&gpdpp->vel[j]) == 1);
	}
    assert(xdr_float(xdrs,&gpdpp->rho) == 1);
    assert(xdr_float(xdrs,&gpdpp->temp) == 1);
    assert(xdr_float(xdrs,&gpdpp->hsmooth) == 1);
    assert(xdr_float(xdrs,&gpdpp->metals) == 1);
    assert(xdr_float(xdrs,&gpdpp->phi) == 1);
    }

void write_tipsy_standard_dark_dpp(XDR *xdrs, DARK_PARTICLE_DPP *dpdpp) {

    int j;

    assert(xdr_float(xdrs,&dpdpp->mass) == 1);
    for (j = 0; j < 3; j++) {
	assert(xdr_double(xdrs,&dpdpp->pos[j]) == 1);
	}
    for (j = 0; j < 3; j++) {
	assert(xdr_float(xdrs,&dpdpp->vel[j]) == 1);
	}
    assert(xdr_float(xdrs,&dpdpp->eps) == 1);
    assert(xdr_float(xdrs,&dpdpp->phi) == 1);
    }

void write_tipsy_standard_star_dpp(XDR *xdrs, STAR_PARTICLE_DPP *spdpp) {

    int j;

    assert(xdr_float(xdrs,&spdpp->mass) == 1);
    for (j = 0; j < 3; j++) {
	assert(xdr_double(xdrs,&spdpp->pos[j]) == 1);
	}
    for (j = 0; j < 3; j++) {
	assert(xdr_float(xdrs,&spdpp->vel[j]) == 1);
	}
    assert(xdr_float(xdrs,&spdpp->metals) == 1);
    assert(xdr_float(xdrs,&spdpp->tform) == 1);
    assert(xdr_float(xdrs,&spdpp->eps) == 1);
    assert(xdr_float(xdrs,&spdpp->phi) == 1);
    }

void write_tipsy_standard(FILE *fp, const TIPSY_STRUCTURE *ts) {

    int i;

    XDR xdrs;
    TIPSY_HEADER *th;
    GAS_PARTICLE *gp;
    DARK_PARTICLE *dp;
    STAR_PARTICLE *sp;

    th = ts->th;
    gp = ts->gp;
    dp = ts->dp;
    sp = ts->sp;
    /*
    ** Write out header
    */
    xdrstdio_create(&xdrs,fp,XDR_ENCODE);
    write_tipsy_standard_header(&xdrs,th);
    /*
    ** Write gas particles
    */
    for (i = 0; i < th->ngas; i++) {
	write_tipsy_standard_gas(&xdrs,&gp[i]);
	}
    /*
    ** Write dark matter particles
    */
    for (i = 0; i < th->ndark; i++) {
	write_tipsy_standard_dark(&xdrs,&dp[i]);
	}
    /*
    ** Write star particles
    */
    for (i = 0; i < th->nstar; i++) {
	write_tipsy_standard_star(&xdrs,&sp[i]);
	}
    xdr_destroy(&xdrs);
    }

void write_tipsy_standard_dpp(FILE *fp, const TIPSY_STRUCTURE_DPP *tsdpp) {

    int i;

    XDR xdrs;
    TIPSY_HEADER *th;
    GAS_PARTICLE_DPP *gpdpp;
    DARK_PARTICLE_DPP *dpdpp;
    STAR_PARTICLE_DPP *spdpp;

    th = tsdpp->th;
    gpdpp = tsdpp->gpdpp;
    dpdpp = tsdpp->dpdpp;
    spdpp = tsdpp->spdpp;
    /*
    ** Write out header
    */
    xdrstdio_create(&xdrs,fp,XDR_ENCODE);
    write_tipsy_standard_header(&xdrs,th);
    /*
    ** Write gas particles
    */
    for (i = 0; i < th->ngas; i++) {
	write_tipsy_standard_gas_dpp(&xdrs,&gpdpp[i]);
	}
    /*
    ** Write dark matter particles
    */
    for (i = 0; i < th->ndark; i++) {
	write_tipsy_standard_dark_dpp(&xdrs,&dpdpp[i]);
	}
    /*
    ** Write star particles
    */
    for (i = 0; i < th->nstar; i++) {
	write_tipsy_standard_star_dpp(&xdrs,&spdpp[i]);
	}
    xdr_destroy(&xdrs);
    }

/*
** Tipsy ascii functions
*/

void read_tipsy_ascii(FILE *fp, TIPSY_STRUCTURE *ts) {

    int i, j;
    TIPSY_HEADER *th;
    GAS_PARTICLE *gp;
    DARK_PARTICLE *dp;
    STAR_PARTICLE *sp;
    
    th = ts->th;
    gp = ts->gp;
    dp = ts->dp;
    sp = ts->sp;
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
	gp = realloc(gp,th->ngas*sizeof(GAS_PARTICLE));
	assert(gp != NULL);
	}
    else {
	gp = NULL;
	}
    if(th->ndark != 0) {
	dp = realloc(dp,th->ndark*sizeof(DARK_PARTICLE));
	assert(dp != NULL);
	}
    else {
	dp = NULL;
	}
    if(th->nstar != 0) {
	sp = realloc(sp,th->nstar*sizeof(STAR_PARTICLE));
	assert(sp != NULL);
	}
    else {
	sp = NULL;
	}
    /*
    ** Read in masses
    */
    for (i = 0; i < th->ngas; i++) {
	assert(fscanf(fp,"%f",&gp[i].mass) == 1);
	}
    for (i = 0; i < th->ndark; i++) {
	assert(fscanf(fp,"%f",&dp[i].mass) == 1);
	}
    for (i = 0; i < th->nstar; i++) {
	assert(fscanf(fp,"%f",&sp[i].mass) == 1);
	}
    /*
    ** Read in positions
    */
    for (j = 0; j < 3; j++) {
	for (i = 0; i < th->ngas; i++) {
	    assert(fscanf(fp,"%f",&gp[i].pos[j]) == 1);
	    }
	for (i = 0; i < th->ndark; i++) {
	    assert(fscanf(fp,"%f",&dp[i].pos[j]) == 1);
	    }
	for (i = 0; i < th->nstar; i++) {
	    assert(fscanf(fp,"%f",&sp[i].pos[j]) == 1);
	    }
	}
    /*
    ** Read in velocities
    */
    for (j = 0; j < 3; j++) {
	for (i = 0; i < th->ngas; i++) {
	    assert(fscanf(fp,"%f",&gp[i].vel[j]) == 1);
	    }
	for (i = 0; i < th->ndark; i++) {
	    assert(fscanf(fp,"%f",&dp[i].vel[j]) == 1);
	    }
	for (i = 0; i < th->nstar; i++) {
	    assert(fscanf(fp,"%f",&sp[i].vel[j]) == 1);
	    }
	}
    /*
    ** Read in softenings of dark matter and star particles
    */
    for (i = 0; i < th->ndark; i++) {
	assert(fscanf(fp,"%f",&dp[i].eps) == 1);
	}
    for (i = 0; i < th->nstar; i++) {
	assert(fscanf(fp,"%f",&sp[i].eps) == 1);
	}
    /*
    ** Read in gas specific stuff
    */
    for (i = 0; i < th->ngas; i++) {
	assert(fscanf(fp,"%f",&gp[i].rho) == 1);
	}
    for (i = 0; i < th->ngas; i++) {
	assert(fscanf(fp,"%f",&gp[i].temp) == 1);
	}
    for (i = 0; i < th->ngas; i++) {
	assert(fscanf(fp,"%f",&gp[i].hsmooth) == 1);
	}
    for (i = 0; i < th->ngas; i++) {
	assert(fscanf(fp,"%f",&gp[i].metals) == 1);
	}
    /*
    ** Read in star specific stuff
    */
    for (i = 0; i < th->nstar; i++) {
	assert(fscanf(fp,"%f",&sp[i].metals) == 1);
	}
    for (i = 0; i < th->nstar; i++) {
	assert(fscanf(fp,"%f",&sp[i].tform) == 1);
	}
    /*
    ** Read in potentials
    */
    for (i = 0; i < th->ngas; i++) {
	assert(fscanf(fp,"%f",&gp[i].phi) == 1);
	}
    for (i = 0; i < th->ndark; i++) {
	assert(fscanf(fp,"%f",&dp[i].phi) == 1);
	}
    for (i = 0; i < th->nstar; i++) {
	assert(fscanf(fp,"%f",&sp[i].phi) == 1);
	}
    /*
    ** Return pointers
    */
    ts->th = th;
    ts->gp = gp;
    ts->dp = dp;
    ts->sp = sp;
    }

void read_tipsy_ascii_dpp(FILE *fp, TIPSY_STRUCTURE_DPP *tsdpp) {

    int i, j;
    TIPSY_HEADER *th;
    GAS_PARTICLE_DPP *gpdpp;
    DARK_PARTICLE_DPP *dpdpp;
    STAR_PARTICLE_DPP *spdpp;
    
    th = tsdpp->th;
    gpdpp = tsdpp->gpdpp;
    dpdpp = tsdpp->dpdpp;
    spdpp = tsdpp->spdpp;
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
	gpdpp = realloc(gpdpp,th->ngas*sizeof(GAS_PARTICLE_DPP));
	assert(gpdpp != NULL);
	}
    else {
	gpdpp = NULL;
	}
    if(th->ndark != 0) {
	dpdpp = realloc(dpdpp,th->ndark*sizeof(DARK_PARTICLE_DPP));
	assert(dpdpp != NULL);
	}
    else {
	dpdpp = NULL;
	}
    if(th->nstar != 0) {
	spdpp = realloc(spdpp,th->nstar*sizeof(STAR_PARTICLE_DPP));
	assert(spdpp != NULL);
	}
    else {
	spdpp = NULL;
	}
    /*
    ** Read in masses
    */
    for (i = 0; i < th->ngas; i++) {
	assert(fscanf(fp,"%f",&gpdpp[i].mass) == 1);
	}
    for (i = 0; i < th->ndark; i++) {
	assert(fscanf(fp,"%f",&dpdpp[i].mass) == 1);
	}
    for (i = 0; i < th->nstar; i++) {
	assert(fscanf(fp,"%f",&spdpp[i].mass) == 1);
	}
    /*
    ** Read in positions
    */
    for (j = 0; j < 3; j++) {
	for (i = 0; i < th->ngas; i++) {
	    assert(fscanf(fp,"%lf",&gpdpp[i].pos[j]) == 1);
	    }
	for (i = 0; i < th->ndark; i++) {
	    assert(fscanf(fp,"%lf",&dpdpp[i].pos[j]) == 1);
	    }
	for (i = 0; i < th->nstar; i++) {
	    assert(fscanf(fp,"%lf",&spdpp[i].pos[j]) == 1);
	    }
	}
    /*
    ** Read in velocities
    */
    for (j = 0; j < 3; j++) {
	for (i = 0; i < th->ngas; i++) {
	    assert(fscanf(fp,"%f",&gpdpp[i].vel[j]) == 1);
	    }
	for (i = 0; i < th->ndark; i++) {
	    assert(fscanf(fp,"%f",&dpdpp[i].vel[j]) == 1);
	    }
	for (i = 0; i < th->nstar; i++) {
	    assert(fscanf(fp,"%f",&spdpp[i].vel[j]) == 1);
	    }
	}
    /*
    ** Read in softenings of dark matter and star particles
    */
    for (i = 0; i < th->ndark; i++) {
	assert(fscanf(fp,"%f",&dpdpp[i].eps) == 1);
	}
    for (i = 0; i < th->nstar; i++) {
	assert(fscanf(fp,"%f",&spdpp[i].eps) == 1);
	}
    /*
    ** Read in gas specific stuff
    */
    for (i = 0; i < th->ngas; i++) {
	assert(fscanf(fp,"%f",&gpdpp[i].rho) == 1);
	}
    for (i = 0; i < th->ngas; i++) {
	assert(fscanf(fp,"%f",&gpdpp[i].temp) == 1);
	}
    for (i = 0; i < th->ngas; i++) {
	assert(fscanf(fp,"%f",&gpdpp[i].hsmooth) == 1);
	}
    for (i = 0; i < th->ngas; i++) {
	assert(fscanf(fp,"%f",&gpdpp[i].metals) == 1);
	}
    /*
    ** Read in star specific stuff
    */
    for (i = 0; i < th->nstar; i++) {
	assert(fscanf(fp,"%f",&spdpp[i].metals) == 1);
	}
    for (i = 0; i < th->nstar; i++) {
	assert(fscanf(fp,"%f",&spdpp[i].tform) == 1);
	}
    /*
    ** Read in potentials
    */
    for (i = 0; i < th->ngas; i++) {
	assert(fscanf(fp,"%f",&gpdpp[i].phi) == 1);
	}
    for (i = 0; i < th->ndark; i++) {
	assert(fscanf(fp,"%f",&dpdpp[i].phi) == 1);
	}
    for (i = 0; i < th->nstar; i++) {
	assert(fscanf(fp,"%f",&spdpp[i].phi) == 1);
	}
    /*
    ** Return pointers
    */
    tsdpp->th = th;
    tsdpp->gpdpp = gpdpp;
    tsdpp->dpdpp = dpdpp;
    tsdpp->spdpp = spdpp;
    }

void write_tipsy_ascii(FILE *fp, const TIPSY_STRUCTURE *ts) {

    int i, j;
    TIPSY_HEADER *th;
    GAS_PARTICLE *gp;
    DARK_PARTICLE *dp;
    STAR_PARTICLE *sp;
    
    th = ts->th;
    gp = ts->gp;
    dp = ts->dp;
    sp = ts->sp;
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
	assert(fprintf(fp,"%.6e\n",gp[i].mass) > 0);
	}
    for (i = 0; i < th->ndark; i++) {
	assert(fprintf(fp,"%.6e\n",dp[i].mass) > 0);
	}
    for (i = 0; i < th->nstar; i++) {
	assert(fprintf(fp,"%.6e\n",sp[i].mass) > 0);
	}
    /*
    ** Write out positions
    */
    for (j = 0; j < 3; j++) {
	for (i = 0; i < th->ngas; i++) {
	    assert(fprintf(fp,"%.6e\n",gp[i].pos[j]) > 0);
	    }
	for (i = 0; i < th->ndark; i++) {
	    assert(fprintf(fp,"%.6e\n",dp[i].pos[j]) > 0);
	    }
	for (i = 0; i < th->nstar; i++) {
	    assert(fprintf(fp,"%.6e\n",sp[i].pos[j]) > 0);
	    }
	}
    /*
    ** Write out velocities
    */
    for (j = 0; j < 3; j++) {
	for (i = 0; i < th->ngas; i++) {
	    assert(fprintf(fp,"%.6e\n",gp[i].vel[j]) > 0);
	    }
	for (i = 0; i < th->ndark; i++) {
	    assert(fprintf(fp,"%.6e\n",dp[i].vel[j]) > 0);
	    }
	for (i = 0; i < th->nstar; i++) {
	    assert(fprintf(fp,"%.6e\n",sp[i].vel[j]) > 0);
	    }
	}
     /*
    ** Write out softenings for dark matter and star particles
    */
    for (i = 0; i < th->ndark; i++) {
	assert(fprintf(fp,"%.6e\n",dp[i].eps) > 0);
	}
    for (i = 0; i < th->nstar; i++) {
	assert(fprintf(fp,"%.6e\n",sp[i].eps) > 0);
	}
    /*
    ** Write out gas specific stuff
    */
    for (i = 0; i < th->ngas; i++) {
	assert(fprintf(fp,"%.6e\n",gp[i].rho) > 0);
	}
    for (i = 0; i < th->ngas; i++) {
	assert(fprintf(fp,"%.6e\n",gp[i].temp) > 0);
	}
    for (i = 0; i < th->ngas; i++) {
	assert(fprintf(fp,"%.6e\n",gp[i].hsmooth) > 0);
	}
    for (i = 0; i < th->ngas; i++) {
	assert(fprintf(fp,"%.6e\n",gp[i].metals) > 0);
	}
    /*
    ** Write out star specific stuff
    */
    for (i = 0; i < th->nstar; i++) {
	assert(fprintf(fp,"%.6e\n",sp[i].metals) > 0);
	}
    for (i = 0; i < th->nstar; i++) {
	assert(fprintf(fp,"%.6e\n",sp[i].tform) > 0);
	}
    /*
    ** Write out potentials
    */
    for (i = 0; i < th->ngas; i++) {
	assert(fprintf(fp,"%.6e\n",gp[i].phi) > 0);
	}
    for (i = 0; i < th->ndark; i++) {
	assert(fprintf(fp,"%.6e\n",dp[i].phi) > 0);
	}
    for (i = 0; i < th->nstar; i++) {
	assert(fprintf(fp,"%.6e\n",sp[i].phi) > 0);
	}
    }

void write_tipsy_ascii_dpp(FILE *fp, const TIPSY_STRUCTURE_DPP *tsdpp) {

    int i, j;
    TIPSY_HEADER *th;
    GAS_PARTICLE_DPP *gpdpp;
    DARK_PARTICLE_DPP *dpdpp;
    STAR_PARTICLE_DPP *spdpp;
    
    th = tsdpp->th;
    gpdpp = tsdpp->gpdpp;
    dpdpp = tsdpp->dpdpp;
    spdpp = tsdpp->spdpp;
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
	assert(fprintf(fp,"%.6e\n",gpdpp[i].mass) > 0);
	}
    for (i = 0; i < th->ndark; i++) {
	assert(fprintf(fp,"%.6e\n",dpdpp[i].mass) > 0);
	}
    for (i = 0; i < th->nstar; i++) {
	assert(fprintf(fp,"%.6e\n",spdpp[i].mass) > 0);
	}
    /*
    ** Write out positions
    */
    for (j = 0; j < 3; j++) {
	for (i = 0; i < th->ngas; i++) {
	    assert(fprintf(fp,"%.14e\n",gpdpp[i].pos[j]) > 0);
	    }
	for (i = 0; i < th->ndark; i++) {
	    assert(fprintf(fp,"%.14e\n",dpdpp[i].pos[j]) > 0);
	    }
	for (i = 0; i < th->nstar; i++) {
	    assert(fprintf(fp,"%.14e\n",spdpp[i].pos[j]) > 0);
	    }
	}
    /*
    ** Write out velocities
    */
    for (j = 0; j < 3; j++) {
	for (i = 0; i < th->ngas; i++) {
	    assert(fprintf(fp,"%.6e\n",gpdpp[i].vel[j]) > 0);
	    }
	for (i = 0; i < th->ndark; i++) {
	    assert(fprintf(fp,"%.6e\n",dpdpp[i].vel[j]) > 0);
	    }
	for (i = 0; i < th->nstar; i++) {
	    assert(fprintf(fp,"%.6e\n",spdpp[i].vel[j]) > 0);
	    }
	}
     /*
    ** Write out softenings for dark matter and star particles
    */
    for (i = 0; i < th->ndark; i++) {
	assert(fprintf(fp,"%.6e\n",dpdpp[i].eps) > 0);
	}
    for (i = 0; i < th->nstar; i++) {
	assert(fprintf(fp,"%.6e\n",spdpp[i].eps) > 0);
	}
    /*
    ** Write out gas specific stuff
    */
    for (i = 0; i < th->ngas; i++) {
	assert(fprintf(fp,"%.6e\n",gpdpp[i].rho) > 0);
	}
    for (i = 0; i < th->ngas; i++) {
	assert(fprintf(fp,"%.6e\n",gpdpp[i].temp) > 0);
	}
    for (i = 0; i < th->ngas; i++) {
	assert(fprintf(fp,"%.6e\n",gpdpp[i].hsmooth) > 0);
	}
    for (i = 0; i < th->ngas; i++) {
	assert(fprintf(fp,"%.6e\n",gpdpp[i].metals) > 0);
	}
    /*
    ** Write out star specific stuff
    */
    for (i = 0; i < th->nstar; i++) {
	assert(fprintf(fp,"%.6e\n",spdpp[i].metals) > 0);
	}
    for (i = 0; i < th->nstar; i++) {
	assert(fprintf(fp,"%.6e\n",spdpp[i].tform) > 0);
	}
    /*
    ** Write out potentials
    */
    for (i = 0; i < th->ngas; i++) {
	assert(fprintf(fp,"%.6e\n",gpdpp[i].phi) > 0);
	}
    for (i = 0; i < th->ndark; i++) {
	assert(fprintf(fp,"%.6e\n",dpdpp[i].phi) > 0);
	}
    for (i = 0; i < th->nstar; i++) {
	assert(fprintf(fp,"%.6e\n",spdpp[i].phi) > 0);
	}
    }

/*
** Gadget functions
*/

int comp_gp(const void *a, const void *b) {

    QSORT_GP *aa = (QSORT_GP *) a;
    QSORT_GP *bb = (QSORT_GP *) b;

    return (aa->index < bb->index)?-1:1;
    }

int comp_dp(const void *a, const void *b) {

    QSORT_DP *aa = (QSORT_DP *) a;
    QSORT_DP *bb = (QSORT_DP *) b;

    return (aa->index < bb->index)?-1:1;
    }

int comp_sp(const void *a, const void *b) {

    QSORT_SP *aa = (QSORT_SP *) a;
    QSORT_SP *bb = (QSORT_SP *) b;

    return (aa->index < bb->index)?-1:1;
    }

void read_gadget_binary(FILE *fp, TIPSY_STRUCTURE *ts, double a, double dx, double dy, double dz, double dof, double mmw, double uvf) {

    int i, j, dummy1, dummy2, Npwm, Ntotgad;
    int *gasindex, *darkindex, *starindex;
    float temp;
    double inverse_sqrt_a, deltapos[3];
    GADGET_HEADER gh;
    TIPSY_HEADER *th;
    GAS_PARTICLE *gp;
    DARK_PARTICLE *dp;
    STAR_PARTICLE *sp;
    QSORT_GP *qgp;
    QSORT_DP *qdp;
    QSORT_SP *qsp;

    th = ts->th;
    gp = ts->gp;
    dp = ts->dp;
    sp = ts->sp;
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
	gp = realloc(gp,th->ngas*sizeof(GAS_PARTICLE));
	assert(gp != NULL);
	}
    else {
	gp = NULL;
	}
    if(th->ndark != 0) {
	dp = realloc(dp,th->ndark*sizeof(DARK_PARTICLE));
	assert(dp != NULL);
	}
    else {
	dp = NULL;
	}
    if(th->nstar != 0) {
	sp = realloc(sp,th->nstar*sizeof(STAR_PARTICLE));
	assert(sp != NULL);
	}
    else {
	sp = NULL;
	}
    /*
    ** Read in positions
    */
    assert(fread(&dummy1,sizeof(int),1,fp) == 1); 
    for(i = 0; i < th->ngas; i++) {
	for(j = 0; j < 3; j++) {
	    assert(fread(&temp,sizeof(float),1,fp) == 1);
	    gp[i].pos[j] = temp + deltapos[j];
	    }
	}
    for(i = 0; i < th->ndark; i++) {
	for(j = 0; j < 3; j++) {
	    assert(fread(&temp,sizeof(float),1,fp) == 1);
	    dp[i].pos[j] = temp + deltapos[j];
	    }
	}
    fseek(fp,3*gh.npart[2]*sizeof(float),SEEK_CUR);
    fseek(fp,3*gh.npart[3]*sizeof(float),SEEK_CUR);
    for(i = 0; i < th->nstar; i++) {
	for(j = 0; j < 3; j++) {
	    assert(fread(&temp,sizeof(float),1,fp) == 1);
	    sp[i].pos[j] = temp + deltapos[j];
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
	    assert(fread(&gp[i].vel[j],sizeof(float),1,fp) == 1);
	    gp[i].vel[j] *= inverse_sqrt_a;
	    }
	}
    for(i = 0; i < th->ndark; i++) {
	for(j = 0; j < 3; j++) {
	    assert(fread(&dp[i].vel[j],sizeof(float),1,fp) == 1);
	    dp[i].vel[j] *= inverse_sqrt_a;
	    }
	}
    fseek(fp,3*gh.npart[2]*sizeof(float),SEEK_CUR);
    fseek(fp,3*gh.npart[3]*sizeof(float),SEEK_CUR);
    for(i = 0; i < th->nstar; i++) {
	for(j = 0; j < 3; j++) {
	    assert(fread(&sp[i].vel[j],sizeof(float),1,fp) == 1);
	    sp[i].vel[j] *= inverse_sqrt_a;
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
		assert(fread(&gp[i].mass,sizeof(float),1,fp) == 1);
		}
	    }
	else {
	    for(i = 0; i < th->ngas; i++) {
		gp[i].mass = gh.mass[0];
		}
	    }
	if (gh.mass[1] == 0) {
	    for(i = 0; i < th->ndark; i++) {
		assert(fread(&dp[i].mass,sizeof(float),1,fp) == 1);
		}
	    }
	else {
	    for(i = 0; i < th->ndark; i++) {
		dp[i].mass = gh.mass[1];
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
		assert(fread(&sp[i].mass,sizeof(float),1,fp) == 1);
		}
	    }
	else {
	    for(i = 0; i < th->nstar; i++) {
		sp[i].mass = gh.mass[4];
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
	    gp[i].mass = gh.mass[0];
	    }
	for(i = 0; i < th->ndark; i++) {
	    dp[i].mass = gh.mass[1];
	    }
	for(i = 0; i < th->nstar; i++) {
	    sp[i].mass = gh.mass[4];
	    }
	}
    /*
    ** Read in gas temperatures
    */
    if (th->ngas > 0) {
	assert(fread(&dummy1,sizeof(int),1,fp) == 1);
	for(i = 0; i < th->ngas; i++) {
	    assert(fread(&temp,sizeof(float),1,fp) == 1);
	    gp[i].temp = temp*(uvf*uvf)*(2*mmw*PROTONMASS)/(dof*kB);
	    }
	assert(fread(&dummy2,sizeof(int),1,fp) == 1);
	assert(dummy1 == dummy2 && dummy2 == th->ngas*sizeof(float));
	}
    /*
    ** Set values for parameters that are normally not in the gadget file
    */
    for(i = 0; i < th->ngas; i++) {
	gp[i].rho = 0;
	gp[i].hsmooth = 0;
	gp[i].metals = 0;
	gp[i].phi = 0;
	}
    for(i = 0; i < th->ndark; i++) {
	dp[i].eps = 0;
	dp[i].phi = 0;
	}
    for(i = 0; i < th->nstar; i++) {
	sp[i].metals = 0;
	sp[i].tform = 0;
	sp[i].eps = 0;
	sp[i].phi = 0;
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
		assert(fread(&gp[i].rho,sizeof(float),1,fp) == 1);
		}
	    assert(fread(&dummy2,sizeof(int),1,fp) == 1);
	    assert(dummy1 == dummy2 && dummy2 == th->ngas*sizeof(float));
	    /*
	    ** Read in softening lengths
	    */
	    fread(&dummy1,sizeof(int),1,fp);
	    if (feof(fp) == 0) {
		for(i = 0; i < th->ngas; i++) {
		    assert(fread(&gp[i].hsmooth,sizeof(float),1,fp) == 1);
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
		assert(fread(&gp[i].phi,sizeof(float),1,fp) == 1);
		}
	    for(i = 0; i < th->ndark; i++) {
		assert(fread(&dp[i].phi,sizeof(float),1,fp) == 1);
		}
	    for(i = 0; i < th->nstar; i++) {
		assert(fread(&sp[i].phi,sizeof(float),1,fp) == 1);
		}
	    assert(fread(&dummy2,sizeof(int),1,fp) == 1);
	    assert(dummy1 == dummy2 && dummy2 == th->ntotal*sizeof(float));
	    }
	}
    /*
    ** Sort particles
    */
    qgp = malloc(th->ngas*sizeof(QSORT_GP));
    assert(qgp != NULL);
    for (i = 0; i < th->ngas; i++) {
	qgp[i].index = gasindex[i];
	qgp[i].gp = gp[i];
	}
    qsort(qgp,th->ngas,sizeof(QSORT_GP),comp_gp);
    for (i = 0; i < th->ngas; i++) {
	gp[i] = qgp[i].gp;
	}
    qdp = malloc(th->ndark*sizeof(QSORT_DP));
    assert(qdp != NULL);
    for (i = 0; i < th->ndark; i++) {
	qdp[i].index = darkindex[i];
	qdp[i].dp = dp[i];
	}
    qsort(qdp,th->ndark,sizeof(QSORT_DP),comp_dp);
    for (i = 0; i < th->ndark; i++) {
	dp[i] = qdp[i].dp;
	}
    qsp = malloc(th->nstar*sizeof(QSORT_SP));
    assert(qsp != NULL);
    for (i = 0; i < th->nstar; i++) {
	qsp[i].index = starindex[i];
	qsp[i].sp = sp[i];
	}
    qsort(qsp,th->nstar,sizeof(QSORT_SP),comp_sp);
    for (i = 0; i < th->nstar; i++) {
	sp[i] = qsp[i].sp;
	}
    free(gasindex);
    free(darkindex);
    free(starindex);
    free(qgp);
    free(qdp);
    free(qsp);
    /*
    ** Return pointers
    */
    ts->th = th;
    ts->gp = gp;
    ts->dp = dp;
    ts->sp = sp;
    }

void write_gadget_binary(FILE *fp, const TIPSY_STRUCTURE *ts, double a, double dx, double dy, double dz, double dof, double mmw, double uvf) {

    int i, j, dummy, index;
    float temp;
    double sqrt_a, deltapos[3];
    GADGET_HEADER gh;
    TIPSY_HEADER *th;
    GAS_PARTICLE *gp;
    DARK_PARTICLE *dp;
    STAR_PARTICLE *sp;

    th = ts->th;
    gp = ts->gp;
    dp = ts->dp;
    sp = ts->sp;
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
	    temp = gp[i].pos[j] + deltapos[j];
	    assert(fwrite(&temp,sizeof(float),1,fp) == 1);
	    }
	}
    for(i = 0; i < th->ndark; i++) {
	for(j = 0; j < 3; j++) {
	    temp = dp[i].pos[j] + deltapos[j];
	    assert(fwrite(&temp,sizeof(float),1,fp) == 1);
	    }
	}
    for(i = 0; i < th->nstar; i++) {
	for(j = 0; j < 3; j++) {
	    temp = sp[i].pos[j] + deltapos[j];
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
	    temp = gp[i].vel[j]*sqrt_a;
	    assert(fwrite(&temp,sizeof(float),1,fp) == 1);
	    }
	}
    for(i = 0; i < th->ndark; i++) {
	for(j = 0; j < 3; j++) {
	    temp = dp[i].vel[j]*sqrt_a;
	    assert(fwrite(&temp,sizeof(float),1,fp) == 1);
	    }
	}
    for(i = 0; i < th->nstar; i++) {
	for(j = 0; j < 3; j++) {
	    temp = sp[i].vel[j]*sqrt_a;
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
	assert(fwrite(&gp[i].mass,sizeof(float),1,fp) == 1);
	}
    for(i = 0; i < th->ndark; i++) {
	assert(fwrite(&dp[i].mass,sizeof(float),1,fp) == 1);
	}
    for(i = 0; i < th->nstar; i++) {
	assert(fwrite(&sp[i].mass,sizeof(float),1,fp) == 1);
	}
    assert(fwrite(&dummy,sizeof(int),1,fp) == 1);
    /*
    ** Write out gas specific stuff
    */
    if (th->ngas > 0) {
	dummy = th->ngas*sizeof(float);
	assert(fwrite(&dummy,sizeof(int),1,fp) == 1);
	for(i = 0; i < th->ngas; i++) {
	    temp = gp[i].temp*(dof*kB)/(2*mmw*PROTONMASS)/(uvf*uvf);
	    assert(fwrite(&temp,sizeof(float),1,fp) == 1);
	    }
	assert(fwrite(&dummy,sizeof(int),1,fp) == 1);
	dummy = th->ngas*sizeof(float);
	assert(fwrite(&dummy,sizeof(int),1,fp) == 1);
	for(i = 0; i < th->ngas; i++) {
	    assert(fwrite(&gp[i].rho,sizeof(float),1,fp) == 1);
	    }
	assert(fwrite(&dummy,sizeof(int),1,fp) == 1);
	dummy = th->ngas*sizeof(float);
	assert(fwrite(&dummy,sizeof(int),1,fp) == 1);
	for(i = 0; i < th->ngas; i++) {
	    assert(fwrite(&gp[i].hsmooth,sizeof(float),1,fp) == 1);
	    }
	assert(fwrite(&dummy,sizeof(int),1,fp) == 1);
	}
    /*
    ** Write out potenitals
    */
    dummy = th->ntotal*sizeof(float);
    assert(fwrite(&dummy,sizeof(int),1,fp) == 1);
    for(i = 0; i < th->ngas; i++) {
	assert(fwrite(&gp[i].phi,sizeof(float),1,fp) == 1);
	}
    for(i = 0; i < th->ndark; i++) {
	assert(fwrite(&dp[i].phi,sizeof(float),1,fp) == 1);
	}
    for(i = 0; i < th->nstar; i++) {
	assert(fwrite(&sp[i].phi,sizeof(float),1,fp) == 1);
	}
    assert(fwrite(&dummy,sizeof(int),1,fp) == 1);
    }

/*
** Array allocation, read and write functions
*/

void allocate_array_particle(const ARRAY_HEADER *ah, ARRAY_PARTICLE *ap) {

    if (ah->N[1] != 0) {
	ap->ia = malloc(ah->N[1]*sizeof(int));
        assert(ap->ia != NULL);
        }
    else {
	ap->ia = NULL;
	}
    if (ah->N[2] != 0) {
	ap->fa = malloc(ah->N[2]*sizeof(float));
        assert(ap->fa != NULL);
        }
    else {
	ap->fa = NULL;
	}
    if (ah->N[3] != 0) {
	ap->da = malloc(ah->N[3]*sizeof(double));
        assert(ap->da != NULL);
        }
    else {
	ap->da = NULL;
	}
    }

void read_array_header(XDR *xdrs, ARRAY_HEADER *ah) {

    int i;

    for (i = 0; i < 4; i++) {
	assert(xdr_int(xdrs,&ah->N[i]) == 1);
	}
    }

void read_array_particle(XDR *xdrs, const ARRAY_HEADER *ah, ARRAY_PARTICLE *ap) {

    int i;

    for (i = 0; i < ah->N[1]; i++) {
	assert(xdr_int(xdrs,&ap->ia[i]) == 1);
	}
    for (i = 0; i < ah->N[2]; i++) {
	assert(xdr_float(xdrs,&ap->fa[i]) == 1);
	}
    for (i = 0; i < ah->N[3]; i++) {
	assert(xdr_double(xdrs,&ap->da[i]) == 1);
	}
    }

void write_array_header(XDR *xdrs, ARRAY_HEADER *ah) {

    int i;
    
    for (i = 0; i < 4; i++) {
	assert(xdr_int(xdrs,&ah->N[i]) == 1);
	}
    }

void write_array_particle(XDR *xdrs, const ARRAY_HEADER *ah, ARRAY_PARTICLE *ap) {

    int i;

    for (i = 0; i < ah->N[1]; i++) {
	assert(xdr_int(xdrs,&ap->ia[i]) == 1);
	}
    for (i = 0; i < ah->N[2]; i++) {
	assert(xdr_float(xdrs,&ap->fa[i]) == 1);
	}
    for (i = 0; i < ah->N[3]; i++) {
	assert(xdr_double(xdrs,&ap->da[i]) == 1);
	}
    }

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

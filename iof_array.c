/*
** array.c
**
** Various handling functions for arrays
**
** Written by Marcel Zemp
*/

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <malloc.h>
#include <assert.h>
#include "iof_array.h"

/*
** Functions
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

void read_array_xdr_header(XDR *xdrs, ARRAY_HEADER *ah) {

	int i;

	for (i = 0; i < 4; i++) {
		assert(xdr_int(xdrs,&ah->N[i]) == 1);
		}
	}

void read_array_xdr_particle(XDR *xdrs, const ARRAY_HEADER *ah, ARRAY_PARTICLE *ap) {

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

void write_array_xdr_header(XDR *xdrs, ARRAY_HEADER *ah) {

	int i;
	
	for (i = 0; i < 4; i++) {
		assert(xdr_int(xdrs,&ah->N[i]) == 1);
		}
	}

void write_array_xdr_particle(XDR *xdrs, const ARRAY_HEADER *ah, ARRAY_PARTICLE *ap) {

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

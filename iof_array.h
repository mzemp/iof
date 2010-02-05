/*
** array.c
**
** Various reading and writing functions for arrays
**
** written by Marcel Zemp
*/

#ifndef IOF_ARRAY_H
#define IOF_ARRAY_H

#include <rpc/types.h>
#include <rpc/xdr.h>

/*
** Structures
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
** Reading and writing functions
*/

void allocate_array_particle(const ARRAY_HEADER *, ARRAY_PARTICLE *);
    
void read_array_xdr_header(XDR *, ARRAY_HEADER *);
void read_array_xdr_particle(XDR *, const ARRAY_HEADER *, ARRAY_PARTICLE *);

void write_array_xdr_header(XDR *, ARRAY_HEADER *);
void write_array_xdr_particle(XDR *, const ARRAY_HEADER *, ARRAY_PARTICLE *);

#endif /* IOF_ARRAY_H */

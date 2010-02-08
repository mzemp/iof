/*
** auxiliary.c
**
** Various auxiliary functions and definitons
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

/*
** Function for correcting positions in periodic boundary conditions cubes
*/

double correct_position(double c, double r, double l) {

    double d;

    d = r-c;
    if (d > 0.5*l) return r - l;
    else if (d < -0.5*l) return r + l;
    else return r;
    }

/*
** Function for putting positions back in box
*/

double put_in_box(double r, double lmin, double lmax) {

    double l;

    l = lmax - lmin;
    if (r < lmin) return r + l;
    else if (r > lmax) return r - l;
    else return r;
    }

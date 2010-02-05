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

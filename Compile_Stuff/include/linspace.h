#if !defined(LINSPACE_H)
#define LINSPACE_H

#include <stdlib.h>

double * linspace(float start, float end, int n);


// Allocates and fills an array with 'n' linearly spaced doubles from 'start' to 'end' (inclusive). Returns pointer to array, NULL on error or invalid parameter

double * linspace(float start, float end, int n) { 

    double *ret = NULL;
    if ( n >= 2 ) { // Check if n is valid. If less than 2, that isn't an array, its a scalar

        double *arr = malloc(n * sizeof(*arr));
        if ( arr != NULL ) {

            double step = (end - start) / (n - 1); // Calculate the step size between elements

            for ( int i = 0 ; i < n ; i++ )
                arr[i] = start + (step * i);

            ret = arr;
        }
    }
    return ret;
}

#endif
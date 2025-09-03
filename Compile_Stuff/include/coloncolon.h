#if !defined(COLONCOLON_H)
#define COLONCOLON_H

#include <stdlib.h>
#include <math.h>

double * coloncolon(float start, float step, float end);


// Allocates and fills an array with 'step'-spaced doubles from 'start' to less than or equal to 'end', where the last
// element is NAN. Returns pointer to array, NULL on error or invalid parameter

double * coloncolon(float start, float step, float end) {

    double *ret = NULL;

    if ( step != 0 ) {

        int n = (int)((end - start) / step) + 1 + 1 ;    // Calculate number of elements. The last +1 is for the terminate element

            double *arr = malloc(n * sizeof(*arr));
            if ( arr != NULL ) {

                for ( int i = 0 ; i < n - 1 ; i++ )
                    arr[i] = start + (i * step);
                    
                arr[n - 1] = nan(""); // Terminate the array with the NAN
                ret = arr;
            }
    }    
    return ret;
}

#endif
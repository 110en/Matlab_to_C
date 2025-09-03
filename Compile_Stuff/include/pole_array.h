#if !defined(POLE_ARRAY_H)
#define POLE_ARRAY_H

#include <stdlib.h>
#include <math.h>

double max_array(const double *arr);
double min_array(const double *arr);


// Returns the maximum value of an array of doubles that is terminated by NaN.
// If the array is NULL, returns 0.
double max_array(const double *arr) {

    double ret = 0;

    if ( arr != NULL ) {

        double max = arr[0];
        int i = 0; 

        while ( !isnan(arr[i]) ) {

            if ( arr[i] > max )
                max = arr[i];
                
            i++;
        }
        ret = max;
    }
    return ret;
}


// Returns the minimum value of an array of doubles that is terminated by NaN.
// If the array is NULL, returns 0.
double min_array(const double *arr) {

    double ret = 0;

    if ( arr != NULL ) {

        double min = arr[0];
        int i = 0;

        while ( !isnan(arr[i]) ) {

            if ( arr[i] < min )
                min = arr[i];
                
            i++;
        }
        ret = min;
    }
    return ret;
}

#endif
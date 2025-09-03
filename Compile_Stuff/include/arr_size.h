#if !defined(ARR_SIZE_H)
#define ARR_SIZE_H

#include <stdlib.h>
#include <math.h>

int arr_size(double *arr);


// Finds the size of the array 'arr' of doubles that is terminated by a NaN (size doesn't include the NaN)
int arr_size(double *arr) {

    int size = 0;

    if ( arr != NULL ) {

        while ( !isnan(arr[size]) )
            size++;
    }
    return size;
}

#endif

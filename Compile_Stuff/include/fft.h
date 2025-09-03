#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <complex.h>
#include <omp.h>
#include "fftw3.h"


#if !defined(FFT_H)
#define FFT_H

void ftx(const fftw_complex *s, fftw_complex *fs, int M, int N);
void fty(const fftw_complex *s, fftw_complex *fs, int M, int N, int P);
void ftz(const fftw_complex *s, fftw_complex *fs, int M, int N, int P);
void iftx(const fftw_complex *fs, fftw_complex *s, int M, int N);
void ifty(const fftw_complex *fs, fftw_complex *s, int M, int N, int P);
static void ifftshift_1d(fftw_complex *data, int M);
static void fftshift_1d(fftw_complex *data, int M);


// Forward FFT along first axis (rows) for 2D array (M rows, N cols). Can also take in 3D arrays where N = # columns * # pages.
// Assumes that arrays are in column major order
void ftx(const fftw_complex *s, fftw_complex *fs, int M, int N) {

    // Create blank arrays to use for the FFTW plan
    fftw_complex *blank_in = fftw_malloc(M * sizeof(*blank_in));
    fftw_complex *blank_out = fftw_malloc(M * sizeof(*blank_out));
    if ( blank_in == NULL || blank_out == NULL ) {
        perror("Error allocating memory for blank arrays in ftx\n");
        fflush(stdout);
        exit(EXIT_FAILURE);
    }
    
    // Create the "plan" for the FFT, that makes an optimized FFT code for current machine.
    fftw_plan plan = fftw_plan_dft_1d(M, blank_in, blank_out, FFTW_FORWARD, FFTW_ESTIMATE);

    #pragma omp parallel
    {

        // Allocate buffers
        fftw_complex *arr = fftw_malloc(M * sizeof(*arr));
        fftw_complex *arr_out = fftw_malloc(M * sizeof(*arr_out));
        if ( arr != NULL && arr_out != NULL ) {            

                // Loop through each column (and technically also pages) to ensure all rows are processed
                #pragma omp for
                for ( int j = 0 ; j < N ; j++ ) {

                    // Extract a row from parameter to FFT input array
                    memcpy(arr, &s[j*M], M * sizeof(*s));

                    // ifftshift
                    ifftshift_1d(arr, M);

                    // FFT
                    fftw_execute_dft(plan, arr, arr_out);

                    // fftshift
                    fftshift_1d(arr_out, M);

                    // Store result
                    for ( int i = 0 ; i < M ; i++ )
                        fs[i + (j * M)] = arr_out[i];
                }
        }
        // Free memory no longer needed
        fftw_free(arr);
        fftw_free(arr_out);
    }
    // Free memory no longer needed
    fftw_destroy_plan(plan);
    fftw_free(blank_in);
    fftw_free(blank_out);
}

// Foward FFT along second axis (columns) for 3D array (M rows, N cols, P pgs). 
// Assumes that arrays are in column major order
void fty(const fftw_complex *s, fftw_complex *fs, int M, int N, int P) {

    // Create blank arrays to use for the FFTW plan
    fftw_complex *blank_in = fftw_malloc(M * sizeof(*blank_in));
    fftw_complex *blank_out = fftw_malloc(M * sizeof(*blank_out));
    if ( blank_in == NULL || blank_out == NULL ) {
        perror("Error allocating memory for blank arrays in ftx\n");
        fflush(stdout);
        exit(EXIT_FAILURE);
    }
    
    // Create the "plan" for the FFT, that makes an optimized FFT code for current machine.
    fftw_plan plan = fftw_plan_dft_1d(N, blank_in, blank_out, FFTW_FORWARD, FFTW_ESTIMATE);

    #pragma omp parallel
    {

        // Allocate buffers
        fftw_complex *vec = fftw_malloc(N * sizeof(*vec));
        fftw_complex *vec_out = fftw_malloc(N * sizeof(*vec_out));
        if ( vec != NULL && vec_out != NULL ) {
            
            // Loop through each row and page to ensure all columns are processed
            #pragma omp for collapse(2)
            for ( int k  = 0 ; k < P ; k++ ) {
                for ( int i = 0 ; i < M ; i++ ) {

                    // Extract a col vector from parameter to FFT input array
                    for ( int j = 0 ; j < N ; j++ )
                        vec[j] = s[i + (j * M) + (k * M * N)];

                    // ifftshift    
                    ifftshift_1d(vec, N);

                    // FFT
                    fftw_execute_dft(plan, vec, vec_out);

                    // fftshift
                    fftshift_1d(vec_out, N);
                    
                    // Store result
                    for ( int j = 0 ; j < N ; j++ )
                        fs[i + (j * M) + (k * M * N)] = vec_out[j];
                }
            }
        }
        // Free memory no longer needed 
        fftw_free(vec);
        fftw_free(vec_out);
    }
    // Free memory no longer needed 
    fftw_destroy_plan(plan);
    fftw_free(blank_in);
    fftw_free(blank_out);
}

// Forward FFT along axis 3 (pages) for 3D array (M rows, N cols, P pgs). 
// Assumes that arrays are in column major order
void ftz(const fftw_complex *s, fftw_complex *fs, int M, int N, int P) {

    // Create blank arrays to use for the FFTW plan
    fftw_complex *blank_in = fftw_malloc(M * sizeof(*blank_in));
    fftw_complex *blank_out = fftw_malloc(M * sizeof(*blank_out));
    if ( blank_in == NULL || blank_out == NULL ) {
        perror("Error allocating memory for blank arrays in ftx\n");
        fflush(stdout);
        exit(EXIT_FAILURE);
    }

    // Create the "plan" for the FFT, that makes an optimized FFT code for current machine.
    fftw_plan plan = fftw_plan_dft_1d(P, blank_in, blank_out, FFTW_FORWARD, FFTW_ESTIMATE);

    #pragma omp parallel
    {

        // Allocate buffers
        fftw_complex *page = fftw_malloc(P * sizeof(*page));
        fftw_complex *page_out = fftw_malloc(P * sizeof(*page_out));
        if ( page != NULL && page_out != NULL ) {

            // Loop through each row and column to ensure all pages are processed
            #pragma omp for collapse(2)
            for ( int j = 0 ; j < N ; j++ ) {
                for ( int i = 0 ; i < M ; i++ ) {
                    
                    // Extracta a page from parameter to FFT input array
                    for ( int k = 0 ; k < P ; k++ )
                        page[k] = s[i + (j * M) + (k * M * N)];
                    
                    // FFT
                    fftw_execute_dft(plan, page, page_out);

                    // Store result
                    for ( int k = 0 ; k < P ; k++ )
                        fs[i + (j * M) + (k * M * N)] = page_out[k];
                }
            }
        }
        // Free memory no longer needed
        fftw_free(page);
        fftw_free(page_out);
    }
    // Free memory no longer needed
    fftw_destroy_plan(plan);
    fftw_free(blank_in);
    fftw_free(blank_out);
}

// Backward FFT along first axis (rows) for 2D array (M rows, N cols). Can also take in 3D arrays where N = # columns * # pages.
// Assumes that arrays are in column major order
void iftx(const fftw_complex *fs, fftw_complex *s, int M, int N) {

    // Create blank arrays to use for the FFTW plan
    fftw_complex *blank_in = fftw_malloc(M * sizeof(*blank_in));
    fftw_complex *blank_out = fftw_malloc(M * sizeof(*blank_out));
    if ( blank_in == NULL || blank_out == NULL ) {
        perror("Error allocating memory for blank arrays in ftx\n");
        fflush(stdout);
        exit(EXIT_FAILURE);
    }

    // Create the "plan" for the FFT, that makes an optimized FFT code for current machine.
    fftw_plan plan = fftw_plan_dft_1d(M, blank_in, blank_out, FFTW_BACKWARD, FFTW_ESTIMATE);

    #pragma omp parallel
    {

        // Allocate buffers
        fftw_complex *arr = fftw_malloc(M * sizeof(*arr));
        fftw_complex *arr_out = fftw_malloc(M * sizeof(*arr_out));
        if ( arr != NULL && arr_out != NULL ) {

            // Loop through each column (and technically also pages) to ensure all rows are processed
            #pragma omp for
            for ( int j = 0 ; j < N ; j++ ) {
                // Copy column data to row buffer
                memcpy(arr, &fs[j * M], M * sizeof(*fs));

                // ifftshift
                ifftshift_1d(arr, M);

                // IFFT
                fftw_execute_dft(plan, arr, arr_out);

                // fftshift
                fftshift_1d(arr_out, M);

                // Store result ("normalize" by dividing by M)
                for ( int i = 0 ; i < M ; i++ )
                    s[i + (j * M)] = arr_out[i] / M;
            }
        }
        // Free memory no longer needed
        fftw_free(arr);
        fftw_free(arr_out);
    }
    // Free memory no longer needed
    fftw_destroy_plan(plan);
    fftw_free(blank_in);
    fftw_free(blank_out);    
}

// Backward FFT along second axis (columns) for 3D array (M rows, N cols, P pages).
// Assumes that arrays are in column major order
void ifty(const fftw_complex *fs, fftw_complex *s, int M, int N, int P) {

    // Create blank arrays to use for the FFTW plan
    fftw_complex *blank_in = fftw_malloc(M * sizeof(*blank_in));
    fftw_complex *blank_out = fftw_malloc(M * sizeof(*blank_out));
    if ( blank_in == NULL || blank_out == NULL ) {
        perror("Error allocating memory for blank arrays in ftx\n");
        fflush(stdout);
        exit(EXIT_FAILURE);
    }

    // Create the "plan" for the FFT, that makes an optimized FFT code for current machine.
    fftw_plan plan = fftw_plan_dft_1d(N, blank_in, blank_out, FFTW_BACKWARD, FFTW_ESTIMATE);

    #pragma omp parallel
    {

        // Allocate buffers
        fftw_complex *vec = fftw_malloc(N * sizeof(*vec));
        fftw_complex *vec_out = fftw_malloc(N * sizeof(*vec_out));
        if ( vec != NULL && vec_out != NULL ) {
            
            // Loop through each row and page to ensure all columns are processed
            #pragma omp for collapse(2)
            for ( int k = 0 ; k < P ; k++ ) {
                for ( int i = 0 ; i < M ; i++ ) {

                    // Extract vector along axis 2 (columns)
                    for ( int j = 0 ; j < N ; j++ )
                        vec[j] = fs[i + (j * M) + (k * M * N)];

                    // ifftshift
                    ifftshift_1d(vec, N);

                    // IFFT
                    fftw_execute_dft(plan, vec, vec_out);

                    // fftshift
                    fftshift_1d(vec_out, N);

                    // Store result ("normalize" by dividing by N)
                    for ( int j = 0 ; j < N ; j++ )
                        s[i + (j * M) + (k * M * N)] = vec_out[j] / N;
                }
            }
        }
        // Free memory no longer needed
        fftw_free(vec);
        fftw_free(vec_out);
    }
    // Free memory no longer needed
    fftw_destroy_plan(plan);
    fftw_free(blank_in);
    fftw_free(blank_out);
}

// Helper: inverse shift (ifftshift) with 1D array
static void ifftshift_1d(fftw_complex *data, int M) {

    int half = M / 2;
    
    if ( M % 2 == 0 ) { // In-place swap for efficiency for even M

        for ( int i = 0; i < half; i++ ) {
            
            fftw_complex tmp;
            tmp = data[i];
            data[i] = data[i + half];
            data[i + half] = tmp;
        }
    } else {  // For odd M, use temp buffer strategy

        fftw_complex *tmp = malloc(M * sizeof(*tmp));

        if (tmp != NULL) {

            // Copy second half (excluding the middle element) of data to the start of tmp
            for ( int i = 0 ; i < M - half ; i++ )
                tmp[i] = data[i + half];

            // Copy first half (including the middle element) of data to the end of tmp
            for ( int i = 0 ; i < half ; i++ )
                tmp[M - half + i] = data[i];

            // Copy back
            memcpy(data, tmp, M * sizeof(*tmp));
        }
        free(tmp);
    }
}

/* In a 3D array, a fftshift causes rows to rotate by floored(#rows/2) and pages to rotate floored(#pages/2) */

// Helper: forward shift (fftshift) with 1D array 
static void fftshift_1d(fftw_complex *data, int M) {

    int half = M / 2;
    
    if ( M % 2 == 0 ) { // In-place swap for efficiency for even M
        for ( int i = 0; i < half; i++ ) {

            fftw_complex tmp;
            tmp = data[i];
            data[i] = data[i + half];
            data[i + half] = tmp;
        }
    } else { // For odd M, use temp buffer strategy

        
        fftw_complex *tmp = malloc(M * sizeof(*tmp));

        if (tmp != NULL) {

            // Copy second half (excluding the middle element) of data to the start of tmp
            memcpy(tmp, &data[M - half], half * sizeof(*data));

            // Copy first half (including the middle element) of data to the end of tmp
            memcpy(&tmp[half], data, (M - half) * sizeof(*data));

            // Copy back
            memcpy(data, tmp, M * sizeof(*tmp)); 
        }
        free(tmp);
    }
}

#endif
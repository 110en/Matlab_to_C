/*

    :Created: July 2025
    :Author: Ali Nilforoushan 
        - With The Help Of: Ali Ghazizadeh (University of Michigan, Ann Arbor, ECE Department)
        - Original Matlab Code: Aditya Muppala (University of Michigan, Ann Arbor, ECE Department)
        - Paper in which this is based on: https://ieeexplore.ieee.org/abstract/document/10549956
        - Credit Where Credit is Due: FFTW, GSL, Matio
        
*/


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <complex.h>
#include <math.h>
#include <omp.h>
#include <gsl/gsl_interp.h>
#include "linspace.h"
#include "coloncolon.h"
#include "arr_size.h"
#include "pole_array.h"
#include "matio.h"
#include "fftw3.h"
#include "fft.h"


static int nextpow2(int x);
static double tand(double degrees);

#define PI          (3.141592653589793)

/***************************************/
/* Radar System Operational Parameters */
/***************************************/
                                            
#define FREQ_BW     (24e9)                                                   // Bandwidth 
#define FC          (15.9875*16e9)                                           // Carrier Frequency
#define LIGHT       (3e8)                                                    // Speed of Light 
#define WC          (2 * PI * FC)
#define WBW         (2 * PI * FREQ_BW)
#define LAMBDA      (LIGHT / FC) 
#define LAMBDA_MIN  (LIGHT / (FC + (FREQ_BW / 2)))                           // Wavelength at Highest Frequency
#define LAMBDA_MAX  (LIGHT / (FC - (FREQ_BW / 2)))                           // Wavelength at Lowest Frequency
#define FREQ_MIN    (2 * PI / LAMBDA_MAX)                                    // Wavenumber at Lowest Frequency
#define FREQ_MAX    (2 * PI / LAMBDA_MIN)                                    // Wavenumber at Highest Frequency

/**************************************/
/* Fast-Time Domain Params and Arrays */
/**************************************/

#define FS          (500e6)                                                  // Sampling Rate...
#define T           (100e-6)                                                 // ...in Microseconds
#define N           (50000)                                                  // Number of Samples
#define GAMMA       (WBW / T)                                                // Chirp Rate
#define DT          (1 / FS)                                                 // Time Domain Sampling
#define T_END       ((N - 1) * DT)                                           // End Time of Sampling

/**************************************/
/* Slow-Time Domain Params and Arrays */
/**************************************/

#define STEP        (0.5e-3)                                                 // Step Size
#define FINISH_X    (49.5e-3)                                                // End of X Axis
#define FINISH_Y    (24.5e-3)                                                // End of Y Axis


/**************************************/
/* Main Function aka Start of Program */
/**************************************/

int main(void) {

    /**************************************/
    /*     Fast Time Domain Variables     */
    /**************************************/

    double *t_arr = NULL;                                                    // Time array for data acquisition

    /**************************************/
    /*    File Manipulation Variables     */
    /**************************************/

    mat_t *img_fp = NULL;                                                    // Pointer to the .mat file
    matvar_t *X_POS = NULL;                                                  // Variable to hold X_POS variable, content of variable doesn't matter only dimensions
    matvar_t *Y_POS = NULL;                                                  // Variable to hold Y_POS variable, content of variable doesn't matter only dimensions
    matvar_t *myStruct = NULL;                                               // Variable to hold myStruct variable
    size_t file_rows;                                                        // Number of rows for the final data tensor, also the number of data arrays in each .mat file's myStruct
    size_t file_cols;                                                        // Number of columns for the final data tensor, also the number of .mat files
    size_t struct_fld_len;                                                   // Number of pages for the final data tensor, also the length of an array field in the data of a myStruct
    double *final_data = NULL;                                               // Tensor that holds all data from .mat files. Every column holds data for one .mat file, 
    //                                                                          in which (i,j,:) is one array from a file's myStruct and (:,j,:) is one file
    matvar_t *struct_fld= NULL;                                              // Variable to hold a field from myStruct, which is an array of arrays
    double *struct_fld_data = NULL;                                          // Pointer to an array in the struct_fld
    char f_naming[300] = {'\0'};                                             // The naming convention for the .mat files  
    char f_path[32768] = {'\0'};                                             // Path to files with data. 32768 chars is max for file path. "C:\Users\zimmy\Code\Mat_2_C\Input\" Is for me
    char f_name[33068] = {'\0'};                                             // Path and name of the .mat file to read. 33068 chars is max for file path + file name + safety space
    char f_bkg_naming[300] = {'\0'};                                         // The name of the background files 
    char f_bkg_name[33068] = {'\0'};                                         // Path and name of the background file to read. 33068 chars is max for file path + file name + safety space
    FILE *fp_bkg = NULL;                                                     // Pointer to the file with background data   
    double *bkg_data = NULL;                                                 // Array to hold data read from the background file

}
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


    printf("Program started\n"); 
    fflush(stdout);

    /**************************************/
    /*       Give Variables Values        */
    /**************************************/

    t_arr = linspace(0, T_END, N);
    if ( t_arr == NULL ) {                                                   // Check if memory allocation failed
        perror("Error allocating memory for t_arr\n");
        fflush(stdout);
        exit(EXIT_FAILURE);
    }

    /**************************************/
    /*     Acquiring Data from Files      */
    /**************************************/

    /*
    .mat file architecture:

    FILE        VARIABLE                                        FIELD                                       DATA

                |-- X_POS (1xfile_rows double array)
                |-- Y_POS (1xfile_cols double array)
    .mat file --|
                |-- myStruct (1xfile_rows struct array) --------|-- data (1xfile_rows double array) ----|-- data(0) (1xstruct_fld_len double array) ; upto data(file_rows-1)
    */

    // Obtain number of rows, columns, and pages for the final data tensor

    // Prompt user for file path and naming convention
    printf("Enter Path to Files with Data (e.g. C:\\Users\\zimmy\\Code\\Mat_2_C\\Input\\): ");
    scanf("%32767s", f_path);

    printf("Enter Naming Convention for Files with Data (e.g. Img_Row, where all files named Img_Row1.mat, Img_Row2.mat etc.): ");
    scanf("%299s", f_naming);

    sprintf(f_name, "%s%s%d.mat", f_path, f_naming, 1);                      // Which of the files is read here doesn't matter.

    // Open .mat file
    img_fp = Mat_Open(f_name, MAT_ACC_RDONLY);
    if ( img_fp == NULL ) {
        printf("Error opening file '%s'\n", f_name);
        fflush(stdout);
        exit(EXIT_FAILURE);
    }

    // Read variables in .mat file
    X_POS = Mat_VarRead(img_fp, "X_POS");
    Y_POS = Mat_VarRead(img_fp, "Y_POS");
    myStruct = Mat_VarRead(img_fp, "myStruct");
    if ( X_POS == NULL || Y_POS == NULL || myStruct == NULL ) {              // Check if reading variables failed
        printf("Error reading variable(s) from file.\n");
        fflush(stdout);
        exit(EXIT_FAILURE);
    }

    // Extract info from variables
    file_rows = X_POS -> dims[1];
    file_cols = Y_POS -> dims[1];
    struct_fld = Mat_VarGetStructField(myStruct, "data", MAT_BY_NAME, 1);    // Which of the files is read here doesn't matter.
    if ( struct_fld == NULL ) {
        printf("Error: 'data' field missing in myStruct element 1\n"); 
        fflush(stdout);
        exit(EXIT_FAILURE);
    }
    struct_fld_len = struct_fld -> dims[1];

    // Close .mat file and variables that are no longer used
    Mat_VarFree(X_POS);
    Mat_VarFree(Y_POS);
    Mat_VarFree(myStruct);
    Mat_Close(img_fp);

    // Prompt user for background file name
    printf("Enter Name of File with the Background Data (Assuming it is in the same folder as the rest of the data files): ");
    scanf("%299s", f_bkg_naming);

    sprintf(f_bkg_name, "%s%s", f_path, f_bkg_naming);

    fp_bkg = fopen(f_bkg_name, "r");
    if ( fp_bkg == NULL ) {                                                  // Check if file opening failed
        perror("Error opening background file\n");
        fflush(stdout);
        exit(EXIT_FAILURE);
    }

    // Allocate buffer to hold data from file_bkg
    bkg_data = malloc(struct_fld_len * sizeof(*bkg_data));
    if ( bkg_data == NULL ) {                                                // Check if memory allocation failed
        perror("Error allocating memory for data_bkg\n");
        fflush(stdout);
        exit(EXIT_FAILURE);
    }    

    // Read data from background file
    for ( size_t i = 0 ; i < struct_fld_len ; i++ ) {

        if ( fscanf(fp_bkg, " %lf,", &bkg_data[i]) != 1 ) {
            fprintf(stderr, "Error reading value %zu from background file\n", i);
            exit(EXIT_FAILURE);
        }
    }

    // Close file no longer needed
    fclose(fp_bkg);


}
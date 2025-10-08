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
    /*     Data Manipulation Variables    */
    /**************************************/

    double complex *time2_data = NULL;
    double ele1_math;                                                        // Number used for element-wise math operations on t_arr
    double ele2_math;                                                        // Number used for element-wise math operations on t_arr
    double complex expo;                                                     // Exponential of the element-wise math operations on t_arr
    
    /**************************************/
    /*           FFT Variables            */
    /**************************************/

    fftw_complex *fft_arr = NULL;                                            // Data array that has gone through fft 
    fftw_complex *fft2_arr = NULL;                                           // Data array that has gone through fft twice
    double Qt_arr[N] = {0};
    double dkx;
    double dky;
    float ks = PI / STEP;
    double *kx = NULL;                                                       // Length will be equivalent to file_rows
    double *ky = NULL;                                                       // Length will be equivalent to file_cols
    double *trunc_kx = NULL;                                                 // Column vector containing all the values of kx which are within range, terminated by NAN. 
    double *trunc_ky = NULL;                                                 // Row vector containing all the values of ky which are within range, terminated by NAN.
    int *kx_index = NULL;                                                    // Has all the indices of kx that were copied over to trunc_kx.
    int *ky_index = NULL;                                                    // Has all the indices of ky that were copied over to trunc_ky.
    int num_trunc_kx;                                                        // Length of trunc_kx
    int num_trunc_ky;                                                        // Length of trunc_ky
    int center_kx;                                                           // Row index of the center of kz_3D
    int center_ky;                                                           // Column index of the center of kz_3D
    float *kz_3D = NULL;
    double *kz_1D = NULL;                                                    // A 1x1xN array that holds a slice of kz_3D at the very center index (center_kx, center_ky, :)
    double *part_kz_3D = NULL;                                               // A 1x1xN array that holds a slice of kz_3D for a (num_trunc_kx/2,num_trunc_ky/2,:), terminated by NAN

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

    /**************************************/
    /*      Element-Wise Arithmetics      */
    /**************************************/

    time2_data = malloc(file_rows * file_cols * struct_fld_len * sizeof(*time2_data));

    // Check if memory allocation failed
    if ( time2_data == NULL ) { 
        perror("Error allocating memory for time2_data\n");
        fflush(stdout);
        exit(EXIT_FAILURE);
    }

    ele1_math = (2.0 * GAMMA * 32.85e-2) / LIGHT;
    ele2_math = (WC * 2.0 * 32.85e-2) / LIGHT;
    // Get data from all .mat files and copy them to final_data

    // Opens a file, reads the myStruct variable, then stores all its arrays from the data field into final_data in indicies (i,j,:)
    #pragma omp parallel for
    for ( size_t j = 0 ; j < file_cols ; j++ ) {                              // Loop through all .mat files

        char thread_f_name[33068] = {'\0'};                                   // Thread safe file name variable
        mat_t *thread_img_fp = NULL;                                          // Thread safe pointer to the .mat file
        matvar_t *thread_myStruct = NULL;                                     // Thread safe variable to hold myStruct variable

        sprintf(thread_f_name, "%s%s%d.mat", f_path, f_naming, j + 1);        // +1 is due to 0-based indexing in C, but files start at 1

        thread_img_fp = Mat_Open(thread_f_name, MAT_ACC_RDONLY);
        if ( thread_img_fp == NULL ) {
            printf("Error opening file '%s'\n", thread_f_name);
            fflush(stdout);
            exit(EXIT_FAILURE);
        }

        thread_myStruct = Mat_VarRead(thread_img_fp, "myStruct");
        if ( thread_myStruct == NULL ) {
            printf("Error reading myStruct from file #%d.\n", j + 1);
            fflush(stdout);
            exit(EXIT_FAILURE);
        }

        for ( size_t i = 0 ; i < file_rows ; i++ ) {                         // Loop through all rows aka arrays in data field of the myStruct

            matvar_t *thread_struct_fld = NULL;                              // Variable to hold a field from myStruct for this thread
            double *struct_fld_data = NULL;                                  // Pointer to an array in the struct

            thread_struct_fld = Mat_VarGetStructField(thread_myStruct, "data", MAT_BY_NAME, i);

            if ( thread_struct_fld == NULL ) {
                printf("Error: 'data' field missing in myStruct element %zu\n", i + 1);
                fflush(stdout);
                exit(EXIT_FAILURE);
            }

            struct_fld_data = thread_struct_fld -> data;
            if ( struct_fld_data == NULL ) {
                printf("Error: 'data' field in myStruct element %zu is NULL\n", i + 1);
                fflush(stdout);
                exit(EXIT_FAILURE);
            }

            // Value at (i,j,k) in tensor = tensor[i + (j*rows) + (k*rows*cols)]
            // Multiply each element in t_arr by ele1_math, then add ele2_math to each element.
            // Next, make each element complex, then take the exponential of each complex element. 
            // Then, subtract each corresponding element in bkg_data from each element in struct_fld_data
            // Finally, multiply exponentiated number with "clean" data and store in time2_data.
            // Visually: final_data is a file_rows x file_cols matrix, which holds arrays (1xstruct_fld_len). 
            // We are subtracting bkg_data(1xstruct_fld_len) from each array in final_data before multiplying with the exponentiated number.
            for ( size_t k = 0 ; k < struct_fld_len ; k++ ) {
                expo = cexp(((t_arr[k] * ele1_math) + ele2_math) * I);
                time2_data[i + (j * file_rows) + (k * file_rows * file_cols)] = (struct_fld_data[k] - bkg_data[k]) * expo;
            }    
            // Don't need to free struct_fld, as it is freed when myStruct is destroyed.
        }

        // Close .mat file and variables that are no longer used
        Mat_VarFree(thread_myStruct);
        Mat_Close(thread_img_fp);
    }

    free(bkg_data);                                                          // Free memory no longer needed

    printf("reached post file reading\n");
    fflush(stdout);

    /**************************************/
    /*       FFT and Reconstruction       */
    /**************************************/

    fft_arr = fftw_malloc( file_rows * file_cols * struct_fld_len * sizeof(*fft_arr));
    if ( fft_arr == NULL ) {                                                 // Check if memory allocation failed
        perror("Error allocating memory for fft_arr\n");
        fflush(stdout);
        exit(EXIT_FAILURE);
    }
    fft_arr[0] = nan("");                                                    // Set first element to NAN to check if FFT failed later'

    ftx(time2_data, fft_arr, file_rows, file_cols * struct_fld_len);         // Perform FFT along the first axis (rows) of time2_data
    if ( isnan((double) fft_arr[0]) ) {                                      // Check if FFT failed
        printf("FFTx failed, first element is NAN\n");
        fflush(stdout);
        exit(EXIT_FAILURE);
    }
    free(time2_data);                                                        // Free memory no longer needed

    fft2_arr = fftw_malloc( file_rows * file_cols * struct_fld_len * sizeof(*fft2_arr));
    if ( fft2_arr == NULL ) {                                                // Check if memory allocation failed
        perror("Error allocating memory for fft2_arr\n");
        fflush(stdout);
        exit(EXIT_FAILURE);
    }
    fft2_arr[0] = nan("");                                                   // Set first element to NAN to check if FFT failed later

    fty(fft_arr, fft2_arr, file_rows, file_cols, struct_fld_len);            // Perform FFT along the second axis (columns) of fft_arr
    if ( isnan((double) fft2_arr[0]) ) {                                     // Check if FFT failed
        printf("FFTy failed, first element is NAN\n");
        fflush(stdout);
        exit(EXIT_FAILURE);
    }
    fftw_free(fft_arr);                                                      // Free memory no longer needed

    // Do some math to t_arr and get values for Qt_arr
    for ( int a = 0 ; a < N ; a++ )
        Qt_arr[a] = (t_arr[a] * ((-2.0 * GAMMA) / LIGHT)) + ((WC * 2.0) / LIGHT);  

    dkx = PI / FINISH_X;
    dky = PI / FINISH_Y;

    kx = coloncolon(((file_rows - 1)) / -2.0f, 1.0f, ((file_rows - 1)) / 2.0f);
    ky = coloncolon(((file_cols - 1)) / -2.0f, 1.0f, ((file_cols - 1)) / 2.0f);
    if ( kx == NULL || ky == NULL ) {                                        // Check if memory allocation failed
        perror("Error allocating memory for kx or ky\n");
        fflush(stdout);
        exit(EXIT_FAILURE);
    }

    // Will possibly be larger than needed
    trunc_kx = calloc(file_rows + 1, sizeof(*trunc_kx));                        
    trunc_ky = calloc(file_cols + 1, sizeof(*trunc_ky));
    kx_index = malloc(file_rows * sizeof(*kx_index)); 
    ky_index = malloc(file_cols * sizeof(*ky_index));
    if ( trunc_kx == NULL || trunc_ky == NULL || kx_index == NULL || ky_index == NULL ) { // Check if memory allocation failed
        perror("Error allocating memory for truncated and index arrays of kx/ky\n");
        fflush(stdout);
        exit(EXIT_FAILURE);
    }
    
    int count = 0;
    // Multiply values of kx by dkx
    // Copy over values of kx to trunc_kx if thier absolute value is within the range of ks * tan(30 degrees).
    // Also copy over the index of kx to kx_index
    for ( size_t a = 0 ; a < file_rows ; a++ ) {

        if ( fabs(kx[a] * dkx) <= ks * tand(30) ) {                          // 60 is beamwidth, 30 is half of it. 

            trunc_kx[count] = kx[a] * dkx;
            kx_index[count++] = a;
        }            
    }
    trunc_kx[count] = nan("");                                               // terminate trunc_kx with NAN
    num_trunc_kx = count - 1;                                                // -1 to exclude NAN 

    count = 0;                                                               // Reset count for trunc_ky
    // Multiply values of ky by dky
    // Copy over values of ky to trunc_ky if thier absolute value is within the range of ks * tan(30 degrees).
    // Also copy over the index of ky to ky_index
    for ( size_t a = 0 ; a < file_cols ; a++ ) {

        if ( fabs(ky[a] * dky) <= ks * tand(30) ) {                          // 60 is beamwidth, 30 is half of it.

            trunc_ky[count] = ky[a] * dky;
            ky_index[count++] = a;
        }            
    }
    trunc_ky[count] = nan("");                                               // terminate trunc_ky with NAN
    num_trunc_ky = count - 1;                                                // -1 to exclude NAN

    // Free memory no longer needed
    free(t_arr);
    free(kx);
    free(ky);

    kz_3D = malloc(num_trunc_kx * num_trunc_ky * N * sizeof(*kz_3D));
    part_kz_3D = malloc((N + 1) * sizeof(*part_kz_3D)); 
    kz_1D = malloc(N * sizeof(*kz_1D)); 
    if ( kz_3D == NULL || part_kz_3D == NULL || kz_1D == NULL ) {            // Check if memory allocation failed
        perror("Error allocating memory for kz3D array(s)\n");
        fflush(stdout);
        exit(EXIT_FAILURE);
    }

    // Calculate center indices of kz_3D
    center_kx = num_trunc_kx / 2;
    center_ky = num_trunc_ky / 2;

    // Fill in values for kz_3D, kz_1D, and part_kz_3D
    for ( int k = 0 ; k < N ; k++ ) {
        for ( int j = 0 ; j < num_trunc_ky ; j++ ) {
            for ( int i = 0 ; i < num_trunc_kx ; i++ )
                kz_3D[i + (j * num_trunc_kx) + (k * num_trunc_kx * num_trunc_ky)] = sqrt((Qt_arr[k]*Qt_arr[k]) - (trunc_kx[i] * trunc_kx[i]) - (trunc_ky[j] * trunc_ky[j]));
        }

        kz_1D[k] = (double) kz_3D[center_kx + (center_ky * num_trunc_kx) + (k * num_trunc_kx * num_trunc_ky)];
        part_kz_3D[k] = kz_3D[(num_trunc_kx/2) + ((num_trunc_ky/2) * num_trunc_kx) + (k * num_trunc_kx * num_trunc_ky)];            
    }
    part_kz_3D[N] = nan("");                                                 // terminate part_kz_3D with NAN as the last element

}
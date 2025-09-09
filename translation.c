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

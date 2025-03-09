#ifndef DSP_H
#define DSP_H

#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>
#include <complex.h>
#include <math.h>


// ============================================================
//                 Digital Signal Processing Header
// ============================================================

/*
 * File: dsp.h
 * Description: This header file contains the declarations and 
 *              function prototypes for the Digital Signal Processing (DSP) 
 *              functions used in the project.
 * Author: Thane Gallo
 * Date: 3/8/25
 */


// ============================================================
//                      Function Declarations
// ============================================================

void debug_print_fft(complex float *data, uint8_t size);
complex float twiddle_factor(int k, int N);
complex float* fft(complex float *data, uint8_t size);


// ============================================================
//                       Constants & Macros
// ============================================================

#define DSP_MAX_BUFFER_SIZE 32

#define SAMPLE_RATE 100
#define BLOCK_LENGTH 0x00

#endif // DSP_H

#ifndef DSP_H
#define DSP_H

#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
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

//fft functions

complex float twiddle_factor(int k, int N, int8_t sign);
complex float* fft(complex float *data, uint8_t size);
complex float* inverse_fft(complex float *data, uint8_t size);

//filters 
void low_pass_filter(complex float *data, float freq, float rolloff, uint8_t size);
void high_pass_filter(complex float *data, float freq, float rolloff, uint8_t size);
void band_pass_filter(complex float *data, float upper_cutoff, float lower_cutoff, float hpf_rolloff, float lpf_rolloff, uint8_t size);
void band_stop_filter(complex float *data, float upper_cutoff, float lower_cutoff, float hpf_rolloff, float lpf_rolloff, uint8_t size);


//other pedals to be implimented
/*
delay 
gain 




*/

//misc utilities


// ============================================================
//                       Constants & Macros
// ============================================================

#define DSP_MAX_BUFFER_SIZE 32

#define SAMPLE_RATE 32
#define BLOCK_LENGTH 0x00

#endif // DSP_H

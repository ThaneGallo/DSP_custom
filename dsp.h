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


//z-transform 
// complex float* z_transform( float *data, uint8_t size);

//fft functions
complex float twiddle_factor(int k, int N, int8_t sign);
complex float* turn_to_imag(float *data, uint8_t size);
complex float* fft(complex float *data, uint8_t size);
complex float* inverse_fft(complex float *data, uint8_t size);

//filters in freq domain
void low_pass_filter_freq(complex float *data, float freq, float rolloff, uint8_t size);
void high_pass_filter_freq(complex float *data, float freq, float rolloff, uint8_t size);
void band_pass_filter_freq(complex float *data, float upper_cutoff, float lower_cutoff, float hpf_rolloff, float lpf_rolloff, uint8_t size);
void band_stop_filter_freq(complex float *data, float upper_cutoff, float lower_cutoff, float hpf_rolloff, float lpf_rolloff, uint8_t size);
// void resonance_filter_freq(complex float *data, uint8_t size);

//filter in time domain 
void low_pass_filter_time(float *data, float freq, uint8_t order, uint8_t size);
void high_pass_filter_time(float *data, float freq, uint8_t order, uint8_t size);






//other pedals to be implimented
/*
delay 
gain 
wah 
frog 
clip / distortion
equilizer
*/


//misc utilities
float find_peak_to_peak(complex float *data, uint8_t size);
float find_fundamental_freq(complex float *data, uint8_t size);

void clip_signal(complex float *data, float upper_clip_limit, float lower_clip_limit, uint8_t size);




// ============================================================
//                       Constants & Macros
// ============================================================

#define DSP_MAX_BUFFER_SIZE 32
#define SAMPLE_RATE 32

#define IDEAL_FILTER -1



#endif // DSP_H

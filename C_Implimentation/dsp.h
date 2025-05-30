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


 enum filter_type
 {
     LOW_PASS_FILTER,
     HIGH_PASS_FILTER,
     RESONANCE_FILTER
 
 };

 
struct first_order_coeff
{
    float alpha_0;
    float alpha_1;
};

struct second_order_coeff
{
    float alpha_0;
    float alpha_1;

    float beta_0;
    float beta_1;
    float beta_2;
};


 union filter_coefficients
 {
     struct first_order_coeff first_order;
     struct second_order_coeff second_order;
 };

 typedef struct filter {
    enum filter_type type;
    uint8_t order;
    float cutoff;
    float Q;
    union filter_coefficients coeff;

}filter; 






// ============================================================
//                      Function Declarations
// ============================================================

// z-transform
//  complex float* z_transform( float *data, uint8_t size);

// fft functions
complex float twiddle_factor(int k, int N, int8_t sign);
complex float *turn_to_imag(float *data, uint8_t size);
complex float *fft(complex float *data, uint8_t size);
complex float *inverse_fft(complex float *data, uint8_t size);

void create_filter(filter *filter, enum filter_type type, uint8_t order, float Q, float freq);

// filters in freq domain
void low_pass_filter_freq(complex float *data, filter *filter, uint8_t size);
void high_pass_filter_freq(complex float *data, filter *filter, uint8_t size);
void band_pass_filter_freq(complex float *data, filter *filter, uint8_t size);
void band_stop_filter_freq(complex float *data, filter *filter, uint8_t size);
// void resonance_filter_freq(complex float *data, uint8_t size);

// filter in time domain
void low_pass_filter_time(float *data,  filter *filter, uint8_t size);
void high_pass_filter_time(float *data,  filter *filter, uint8_t size);
void band_pass_filter_time(float *data, filter *low_pass_filter,  filter *high_pass_filter, uint8_t size);
void band_stop_filter_time(float *data, filter *low_pass_filter,  filter *high_pass_filter, uint8_t size);

// other pedals to be implimented
/*
delay
gain
wah
frog
clip / distortion
equilizer
*/

// misc utilities
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

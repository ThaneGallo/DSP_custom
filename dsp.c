#include "dsp.h"




/**
 * @brief uses low pass filter to adjust signal (mimics analog filter)
 * @param data pointer to transformed signal freq response (post fft)
 * @param cutoff desired cutoff frequency
 * @param rolloff desired rate of rolloff (-1 for ideal)
 * @param size size of signal array 
 */
void low_pass_filter(complex float *data, float cutoff, float rolloff, uint8_t size)
{
    uint8_t i = 0;
   
    //in case of ideal filter
    if (rolloff == -1)
    {
        for (i = 0; i < size; i++)
        {
            float freq_bin = ((float)i * SAMPLE_RATE) / size;

            if (freq_bin >= cutoff)
            {
                data[i] = 0;
            }
        }
    }

 
    for (i = 0; i < size; i++)
    {
        float freq_bin = ((float)i * SAMPLE_RATE) / size;

        if (freq_bin >= cutoff)
        {
           float amplitude_ratio = pow(10, (-rolloff * log10f(freq_bin / cutoff)) / 20.0f);
            data[i] = amplitude_ratio * data[i];
        }
    }
}

void high_pass_filter(complex float *data, float cutoff, float rolloff, uint8_t size){
    uint8_t i = 0;
    float magnitude_db;
   
    //in case of ideal filter
    if (rolloff == -1)
    {
        for (i = 0; i < size; i++)
        {
            float freq_bin = ((float)i * SAMPLE_RATE) / size;

            if (freq_bin <= cutoff)
            {
                data[i] = 0;
            }
        }
    }

 
    for (i = 0; i < size; i++)
    {
        float freq_bin = ((float)i * SAMPLE_RATE) / size;

        if (freq_bin <= cutoff)
        {
           float amplitude_ratio = pow(10, (-rolloff * log10f(cutoff / freq_bin)) / 20.0f);
            data[i] = amplitude_ratio * data[i];
        }
    }


}

/**
 * @brief calculates twiddle factor
 * @param k sample number
 * @param N total number of samples
 * @param sign 1 or -1 depending on fft or ifft
 * @return twiddle factor
 */
complex float twiddle_factor(int k, int N, int8_t sign)
{
    float angle = (sign) * 2.0 * M_PI * k / N;
    return cosf(angle) + I * sinf(angle);
}

/**
 * @brief cooley-tukey FFT divide and conquer
 * @param data digitized audio data
 * @param size size of data array (2^n)
 * @return pointer to transformed array
 */
complex float *fft(complex float *data, uint8_t size)
{

    // base case
    if (size == 1)
    {
        return data;
    }

    else
    {

        complex float *even_data = malloc(sizeof(complex float) * size / 2);
        complex float *odd_data = malloc(sizeof(complex float) * size / 2);

        uint8_t i;

        // moves data into even and odd parts
        for (i = 0; i < size / 2; i++)
        {
            even_data[i] = data[2 * i];
            odd_data[i] = data[2 * i + 1];
        }

        even_data = fft(even_data, size / 2);
        odd_data = fft(odd_data, size / 2);

        // apply twiddle factor and recombine
        for (i = 0; i < size / 2; i++)
        {
            data[i] = even_data[i] + twiddle_factor(i, size, 1) * odd_data[i];
            data[i + size / 2] = even_data[i] - twiddle_factor(i, size, 1) * odd_data[i];
        }

        free(even_data);
        free(odd_data);

        return data;
    }
}

/**
 * @brief inverts fft
 * @param data digitized audio data
 * @param size size of data array (2^n)
 * @return pointer to reverted array
 */
complex float *inverse_fft(complex float *data, uint8_t size)
{

    // base case
    if (size == 1)
    {
        return data;
    }

    else
    {

        complex float *even_data = malloc(sizeof(complex float) * size / 2);
        complex float *odd_data = malloc(sizeof(complex float) * size / 2);

        uint8_t i;

        // moves data into even and odd parts
        for (i = 0; i < size / 2; i++)
        {
            even_data[i] = data[2 * i];
            odd_data[i] = data[2 * i + 1];
        }

        even_data = inverse_fft(even_data, size / 2);
        odd_data = inverse_fft(odd_data, size / 2);

        // apply twiddle factor and recombine
        for (i = 0; i < size / 2; i++)
        {
            data[i] = even_data[i] + twiddle_factor(i, size, -1) * odd_data[i];
            data[i + size / 2] = even_data[i] - twiddle_factor(i, size, -1) * odd_data[i];
        }

        // normallizes
        for (i = 0; i < size; i++)
        {
            data[i] /= (size / 2);
        }

        free(even_data);
        free(odd_data);

        return data;
    }

    return;
}
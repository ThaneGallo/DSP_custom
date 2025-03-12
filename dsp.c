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

/**
 * @brief uses high pass filter to adjust signal (mimics analog filter)
 * @param data pointer to transformed signal freq response (post fft)
 * @param cutoff desired cutoff frequency
 * @param rolloff desired rate of rolloff (-1 for ideal)
 * @param size size of signal array 
 */
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
 * @brief mimics analog band pass filter
 * @param data pointer to transformed signal freq response (post fft)
 * @param upper_cutoff desired cutoff for lpf
 * @param lower_cutoff desired cutoff for hpf
 * @param hpf_rolloff desired rate of upper rolloff (-1 for ideal)
 * @param lpf_rolloff desired rate of lower rolloff (-1 for ideal)
 * @param size size of signal array 
 */
void band_pass_filter(complex float *data, float upper_cutoff, float lower_cutoff, float hpf_rolloff, float lpf_rolloff, uint8_t size){
        //does not need any special ordering like bsp bc the filters do not interact with one another
        high_pass_filter(data, lower_cutoff, hpf_rolloff, size);
        low_pass_filter(data, upper_cutoff, lpf_rolloff, size);

}

/**
 * @brief mimics analog band stop filter
 * @param data pointer to transformed signal freq response (post fft)
 * @param upper_cutoff desired cutoff for hpf
 * @param lower_cutoff desired cutoff for lpf
 * @param hpf_rolloff desired rate of lower rolloff (-1 for ideal)
 * @param lpf_rolloff desired rate of higher rolloff (-1 for ideal)
 * @param size size of signal array 
 */
void band_stop_filter(complex float *data, float upper_cutoff, float lower_cutoff, float hpf_rolloff, float lpf_rolloff, uint8_t size){
        uint8_t i;

        // split then sum
        complex float *low_pass_data = malloc(sizeof(complex float) * size);
        complex float *high_pass_data = malloc(sizeof(complex float) * size);

        if(!low_pass_data || !high_pass_data){
            printf("Low or high pass data is NULL");
            return;
        }

        //copies data to new arrays
        memcpy(low_pass_data, data, sizeof(complex float) * size);
        memcpy(high_pass_data, data, sizeof(complex float) * size);


        low_pass_filter(low_pass_data, lower_cutoff, lpf_rolloff, size);
        high_pass_filter(high_pass_data, upper_cutoff, hpf_rolloff, size);

        //thanks to handy dandy superposition i can sum to get filter results
        for(i = 0; i < size; i++){
            data[i] = low_pass_data[i] + high_pass_data[i];
        }

        free(low_pass_data);
        free(high_pass_data);

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
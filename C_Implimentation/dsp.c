#include "dsp.h"

/**
 * @brief returns approximate peak to peak voltage for time domain signal
 * @param data pointer to transformed signal freq response (post fft)
 * @param size size of signal array
 */
float find_peak_to_peak(complex float *data, uint8_t size)
{
    uint8_t i;

    // temp storage for min and max
    float min = 0;
    float max = 0;

    // use magnitude? or real value
    for (i = 0; i < size; i++)
    {
        data[i];
    }

    return max - min;
}

/**
 * @brief clips time domain signal
 * @param data pointer to transformed signal freq response (post fft)
 * @param upper_clip_limit where clipping occurs on + side
 * @param lower_clip_limit where clipping occurs on - side
 * @param size size of signal array
 */
void clip_signal(complex float *data, float upper_clip_limit, float lower_clip_limit, uint8_t size)
{
    uint8_t i;

    // use magnitude? or real value
    for (i = 0; i < size; i++)
    {

        if (cabs(data[i]) > upper_clip_limit)
        {
            data[i] = upper_clip_limit;
        }

        if (cabs(data[i]) < lower_clip_limit)
        {
            data[i] = lower_clip_limit;
        }
    }
};

/**
 * @brief finds f0 frequency for finding harmonics later
 * @param data pointer to transformed signal freq response (post fft)
 * @param size size of signal array
 */
float find_fundamental_freq(complex float *data, uint8_t size)
{

    float fund_freq = 0;
    uint8_t i;

    for (i = 0; i < size; i++)
    {
        float freq_bin = ((float)i * SAMPLE_RATE) / size;

        if (cabsf(data[i]) > fund_freq)
        {
            fund_freq = freq_bin;
        }
    }

    return fund_freq;
}

/**
 * @brief turns float into complex  float for fft
 * @param data pointer to real signal
 * @param size size of signal array
 */
complex float *turn_to_imag(float *data, uint8_t size)
{

    int i = 0;
    complex float *return_data = malloc(sizeof(complex float) * size);

    for (i = 0; i < size; i++)
    {

        return_data[i] = data[i] + 0.0 * I;
    }

    free(data);

    return return_data;
}

/**
 * @brief computes filter coefficients for easier application later
 * @param filter pointer to pilter union
 * @param type enum for type
 * @param order order number for attenuation rate
 * @param Q q-factor (if applicable (order > 1))
 * @param freq cutoff frequency point
 */
void create_filter(filter *filter, enum filter_type type, uint8_t order, float Q, float freq)
{
    filter->type = type;
    filter->cutoff = freq;
    filter->Q = Q;
    
    
    switch (order)
    {
    case -1:
        break;


    case 1:
        
        float tau = 1.0f / (2.0f * M_PI * freq);

        switch (type)
        {
        case LOW_PASS_FILTER:
            filter->coeff.first_order.alpha_0 = SAMPLE_RATE / (tau + SAMPLE_RATE);
            filter->coeff.first_order.alpha_1 = 1 - filter->coeff.first_order.alpha_0;

            break;
        case HIGH_PASS_FILTER:
            filter->coeff.first_order.alpha_0 = tau / (tau + SAMPLE_RATE);
            filter->coeff.first_order.alpha_1 = 1 - filter->coeff.first_order.alpha_0;

            break;

        default:
            printf("Unexpected Filter type: %d\n", filter->type);
            break;
        }

        break;

    case 2:

        float ang_freq = (2 * M_PI * freq) / SAMPLE_RATE;
        float alpha = sinf(ang_freq) / (2 * Q);

        filter->coeff.second_order.alpha_0 = (1 - cosf(ang_freq)) / 2;
        filter->coeff.second_order.alpha_1 = 1 - cosf(ang_freq);

        filter->coeff.second_order.beta_0 = 1 + alpha;
        filter->coeff.second_order.beta_1 = -2 * cosf(ang_freq);
        filter->coeff.second_order.beta_2 = 1 - alpha;

        break;

    default:
        printf("Unexpected Filter Order: %d\n", filter->order);
        break;
    }
}

/**
 * @brief uses low pass filter to adjust signal (mimics analog filter) in time domain
 * @param data pointer to transformed signal freq response
 * @param cutoff desired cutoff frequency
 * @param order decides the order of the filter (higher num = sharper response)
 * @param size size of signal array
 */
void low_pass_filter_time(float *data, filter *filter, uint8_t size)
{
    uint8_t i;
    float previous_value;

    // initializes first ouput sample
    previous_value = data[0];

    for (i = 1; i < size; i++)
    {

        data[i] = filter->coeff.first_order.alpha_0 * data[i] + filter->coeff.first_order.alpha_1 * previous_value;
        previous_value = data[i];
    }
}

/**
 * @brief uses low pass filter to adjust signal (mimics analog filter) in time domain
 * @param data pointer to transformed signal freq response
 * @param cutoff desired cutoff frequency
 * @param order decides the order of the filter (higher num = sharper response)
 * @param size size of signal array
 */
void high_pass_filter_time(float *data, filter *filter, uint8_t size)
{
    uint8_t i;
    float previous_value;
    float temp_value;

    // initializes first ouput sample
    previous_value = data[0];

    for (i = 1; i < size; i++)
    {

        // storage of prev input since we write in place
        temp_value = data[i];

        // Y[n] = A (y[n-1] + x[n] - x[n-1]) is general form
        data[i] = filter->coeff.first_order.alpha_0 * data[i] + filter->coeff.first_order.alpha_1 * previous_value;

        previous_value = temp_value;
    }
}

/**
 * @brief uses filters to adjust signal (mimics analog filter) in time domain
 * @param data pointer to transformed signal freq response
 * @param low_pass_filter low pass filter coefficients and information
 * @param high_pass_filter high pass filter coefficients and information
 * @param size size of signal array
 */
void band_pass_filter_time(float *data, filter *low_pass_filter, filter *high_pass_filter, uint8_t size)
{
    low_pass_filter_time(data, low_pass_filter, size);
    high_pass_filter_time(data, high_pass_filter, size);
}

/**
 * @brief uses low pass filter to adjust signal (mimics analog filter) in time domain
 * @param data pointer to transformed signal freq response
 * @param low_pass_filter low pass filter coefficients and information
 * @param high_pass_filter high pass filter coefficients and information
 * @param size size of signal array
 */
void band_stop_filter_time(float *data, filter *low_pass_filter, filter *high_pass_filter, uint8_t size)
{

    uint8_t i;

    // split then sum
    float *low_pass_data = malloc(sizeof(float) * size);
    float *high_pass_data = malloc(sizeof(float) * size);

    if (!low_pass_data || !high_pass_data)
    {
        printf("Low or high pass data is NULL");
        return;
    }

    // copies data to new arrays
    memcpy(low_pass_data, data, sizeof(float) * size);
    memcpy(high_pass_data, data, sizeof(float) * size);

    low_pass_filter_time(low_pass_data, low_pass_filter, size);
    high_pass_filter_time(high_pass_data, high_pass_filter, size);

    // thanks to handy dandy superposition i can sum to get filter results
    for (i = 0; i < size; i++)
    {
        data[i] = (low_pass_data[i] + high_pass_data[i]);
    }

    free(low_pass_data);
    free(high_pass_data);
}

// /**
//  * @brief uses low pass filter to adjust signal (mimics analog filter) in frequency domain
//  * @param data pointer to transformed signal freq response (post fft)
//  * @param cutoff desired cutoff frequency
//  * @param rolloff desired rate of rolloff (-1 for ideal)
//  * @param size size of signal array
//  */
// void low_pass_filter_freq(complex float *data, filter *filter, uint8_t size)
// {
//     uint8_t i = 0;


    
    
//     // in case of ideal filter
//     if (filter->order == -1)
//     {
//         for (i = 0; i < size; i++)
//         {
//             float freq_bin = ((float)i * SAMPLE_RATE) / size;

//             if (freq_bin >= cutoff)
//             {
//                 data[i] = 0;
//             }
//         }
//     }

//     for (i = 0; i < size; i++)
//     {
//         float freq_bin = ((float)i * SAMPLE_RATE) / size;

//         if (freq_bin >= cutoff)
//         {
//             float amplitude_ratio = pow(10, (-rolloff * log10f(freq_bin / cutoff)) / 20.0f);
//             data[i] = amplitude_ratio * data[i];
//         }
//     }
// }

// /**
//  * @brief uses high pass filter to adjust signal (mimics analog filter)in frequency domain
//  * @param data pointer to transformed signal freq response (post fft)
//  * @param cutoff desired cutoff frequency
//  * @param rolloff desired rate of rolloff (-1 for ideal)
//  * @param size size of signal array
//  */
// void high_pass_filter_freq(complex float *data, filter *filter, uint8_t size)
// {
//     uint8_t i = 0;
//     float magnitude_db;

//     float cutoff = (1 - filter->coeff.first_order.alpha_0) / (filter->coeff.first_order.alpha_0 * 2.0f * M_PI * SAMPLE_RATE);

//     // in case of ideal filter
//     if (filter->order == -1)
//     {
//         for (i = 0; i < size; i++)
//         {
//             float freq_bin = ((float)i * SAMPLE_RATE) / size;

//             if (freq_bin <= cutoff)
//             {
//                 data[i] = 0;
//             }
//         }
//     }

//     for (i = 0; i < size; i++)
//     {
//         float freq_bin = ((float)i * SAMPLE_RATE) / size;

//         if (freq_bin <= cutoff)
//         {
//             float amplitude_ratio = pow(10, (-rolloff * log10f(cutoff / freq_bin)) / 20.0f);
//             data[i] = amplitude_ratio * data[i];
//         }
//     }
// }

// /**
//  * @brief mimics analog band pass filter in frequency domain
//  * @param data pointer to transformed signal freq response (post fft)
//  * @param upper_cutoff desired cutoff for lpf
//  * @param lower_cutoff desired cutoff for hpf
//  * @param hpf_rolloff desired rate of upper rolloff (-1 for ideal)
//  * @param lpf_rolloff desired rate of lower rolloff (-1 for ideal)
//  * @param size size of signal array
//  */
// void band_pass_filter_freq(complex float *data, filter *filter, uint8_t size)
// {
//     // does not need any special ordering like bsp bc the filters do not interact with one another
//     high_pass_filter_freq(data, filter, size);
//     low_pass_filter_freq(data, filter, size);
// }

// /**
//  * @brief mimics analog band stop filter in frequency domain
//  * @param data pointer to transformed signal freq response (post fft)
//  * @param upper_cutoff desired cutoff for hpf
//  * @param lower_cutoff desired cutoff for lpf
//  * @param hpf_rolloff desired rate of lower rolloff (-1 for ideal)
//  * @param lpf_rolloff desired rate of higher rolloff (-1 for ideal)
//  * @param size size of signal array
//  */
// void band_stop_filter_freq(complex float *data, filter *filter, uint8_t size)
// {
//     uint8_t i;

//     // split then sum
//     complex float *low_pass_data = malloc(sizeof(complex float) * size);
//     complex float *high_pass_data = malloc(sizeof(complex float) * size);

//     if (!low_pass_data || !high_pass_data)
//     {
//         printf("Low or high pass data is NULL");
//         return;
//     }

//     // copies data to new arrays
//     memcpy(low_pass_data, data, sizeof(complex float) * size);
//     memcpy(high_pass_data, data, sizeof(complex float) * size);

//     low_pass_filter_freq(low_pass_data, lower_cutoff, lpf_rolloff, size);
//     high_pass_filter_freq(high_pass_data, upper_cutoff, hpf_rolloff, size);

//     // thanks to handy dandy superposition i can sum to get filter results
//     for (i = 0; i < size; i++)
//     {
//         data[i] = (low_pass_data[i] + high_pass_data[i]);
//     }

//     free(low_pass_data);
//     free(high_pass_data);
// }

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
}
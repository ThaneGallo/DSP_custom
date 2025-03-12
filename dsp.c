#include "dsp.h"




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
        for (i = 0; i < size/2; i++)
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
            data[i + size/2] = even_data[i] - twiddle_factor(i, size, 1) * odd_data[i];
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
complex float *inverse_fft(complex float *data, uint8_t size){
    
    
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
           for (i = 0; i < size/2; i++)
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
               data[i + size/2] = even_data[i] - twiddle_factor(i, size, -1) * odd_data[i];
           }

        //normallizes 
           for (i = 0; i < size; i++)
           {
              data[i] /= (size/2);
           }
   
           free(even_data);
           free(odd_data);
   
           return data;
       }
    
    
    return;
}
#include <assert.h>
#include "dsp.h"

/**
 * @brief prints fft of data before and after for debugging
 * @param data digitized audio data
 * @param size size of data array (2^n)
 */
void debug_print_fft(complex float *data, float size, float sample_rate)
{
    uint8_t i;

    printf("buffer after fft:");
    for (i = 0; i < size; i++)
    {
        float freq_bin = ((float)i * sample_rate) / size;
        printf("freq %.4f HZ  magnitude = %.4f\n", freq_bin, cabsf(data[i]));
    }
}

complex float *sine_generator(float freq, float size, float sample_rate)
{

    complex float *data = malloc(sizeof(complex float) * size);
    uint8_t i = 0;

    for (i = 0; i < size; i++)
    {

        float t = (float)i / sample_rate;

        data[i] = sinf(2 * M_PI * freq * t) + 0.0f * I;
    }

    return data;
}

complex float *sine_sum(complex float *sine1, complex float *sine2, int size1, int size2)
{
    int i;

    if (size1 > size2)
    {
        for (i = 0; i < size2; i++)
        {
            sine1[i] += sine2[i];
        }

        return sine1;
    }
    else
    {

        for (i = 0; i < size1; i++)
        {
            sine2[i] += sine1[i];
        }

        return sine2;
    }
}

void main()
{
    uint8_t i = 0;
    uint8_t size = 8;
    uint8_t sine_size = 64;

    // initial dataset (periodic sine)
    complex float data_all_real[8] = {1.0, 0.0, -1.0, 0.0, 1.0, 0.0, -1.0, 0.0};
    // for (i = 0; i < size; i++)
    // {
    //     printf("%f\n", crealf(data_all_real[i]));
    // }

    // only dc
    complex float data_dc_only[8] = {1, 1, 1, 1, 1, 1, 1, 1};
    // for testing dc offset in bin 0 [freq 0Hz]
    complex float data_dc_offset[8] = {2.0, 1.0, 0.0, 1.0, 2.0, 1.0, 0.0, 1.0};

    // generates sine wave

    complex float *data_sine;
    complex float *data_sine_result = malloc(sizeof(complex float) * sine_size);

    // for (i = 0; i < sine_size; i++)
    // {
    //     printf("sine sample %f\n", crealf(data_sine[i]));
    // }

    // creates a wave with magnitue 32 at frequwencies 0 --> 31.5
    for (i = 0; i < 32; i++)
    {
        data_sine = sine_generator(i, sine_size, 64);
        // data_sine_result = sine_generator(, sine_size, SAMPLE_RATE);

        data_sine_result = sine_sum(data_sine_result, data_sine, sine_size, sine_size);
    }

    // perform ffts
    fft(data_dc_offset, size);
    fft(data_all_real, size);
    fft(data_dc_only, size);
    fft(data_sine_result, sine_size);

    // low_pass_filter(data_sine_result, 20, 20, sine_size);
    // high_pass_filter(data_sine_result, 20, 10, sine_size);

    debug_print_fft(data_sine_result, sine_size, SAMPLE_RATE);

    // FFT ASSERTS

    // all real signals are periodic around half of size (nyquist freq)
    for (i = 0; i > size / 2; i++)
    {
        assert(cabs(data_all_real[i + 1]) == cabs(data_all_real[i + (size / 2)]));
    }

    assert(cabs(data_dc_only[0]) == (1.0 * size));
    // shoul be all zero is the signal is not periodic
    for (i = 1; i < size; i++)
    {
        assert(cabs(data_dc_only[i]) == 0.0);
    }

    // verifies dc offset bin 0 [freq 0Hz]
    assert(cabs(data_dc_offset[0]) == (1.0 * size));

    // inverse_fft(data_dc_only, size);
    // inverse_fft(data_dc_offset, size);
    // inverse_fft(data_all_real, size);

    // for (i = 0; i < size; i++)
    // {
    //     printf("%f\n", crealf(data_dc_offset[i]));
    // }

    // debug_print_fft(data_sine, 64, SAMPLE_RATE);
}
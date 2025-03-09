#include "dsp.h"

void main()
{
    complex float data[8] = {1.0, 0.0, -1.0, 0.0, 1.0, 0.0, -1.0, 0.0};

    debug_print_fft(data, 8);
}
#include <math.h>
#include <stdio.h>

#include "BnFFT.h"

int main(int argc, char** argv){

    {
        //Create instance of fft object.
        //Sample length is 2^11 = 2048.
        bn::ComplexFFT2n<11> fft;
        
        /*
        Or you can use ComplexFFT2n3 class for length of (2^n)*3 fft.
        bn::ComplexFFT2n3<10> fft;
        In this case, sample length is 2^10 * 3= 3072.
        */

        //Input/output sample data of fft.
        bn::Complex data[2048];

        //Initiate with input data.
        for (int i= 0; i<2048; i++){
            data[i] = bn::Complex(sin(i*0.2), cos(i*0.4));
        }

        //Execute fft.
        //Output data are orverwritten on input data array.
        fft.execute(data);

        //To execute inverse FFT, set seccond template parameter to true. 
        bn::ComplexFFT2n<11, true> ifft;
        //Data will be original value with scale (2*pi*n).
        ifft.execute(data);


        ////////////
        //Or, you can use prain double array as I/O data with evil C cast
        //in normal c++ compiler
        {
            bn::Complex complex_array[2];
            double* double_pointer = (double*)complex_array;
            bool condition0 = &double_pointer[0] == &complex_array[0].re;
            bool condition1 = &double_pointer[1] == &complex_array[0].im;
            bool condition2 = &double_pointer[2] == &complex_array[1].re;
            bool condition3 = &double_pointer[3] == &complex_array[1].im;

            if (condition0 && condition1 && condition2 && condition3){
                printf("Yes, you can cast *bn::Complex into *double\n");
            }
        }

        double double_data[2048*2];
        for (int i= 0; i<2048; i++){
            double_data[i] = sin(i*0.2);
            double_data[i] = cos(i*0.4);
        }
        fft.execute((bn::Complex*)data);

    }
    return 0;
}

//main.c


#include 		"fft.h"
#include 		"fft-org.h"
#include        <stdio.h>
#include        <math.h>


#define M       3

//number of points
#define N       (1<<M)

fixed real[N], imag[N];
fixed real_org[N], imag_org[N];

int main()
{
    int i;

    for(i=0; i<N; i++)
    {
        real[i] = 1000*cos(i*2*3.1415926535/N);
        real_org[i] = 1000*cos(i*2*3.1415926535/N);
        imag[i] = 0;
        imag_org[i] = 0;
    }

    printf("\nInput Data\n");
    for (i=0; i<N; i++)
    {

        printf("%d: %d, %d\n", i, real[i], imag[i]);
    }


    //FFT
    fix_fft(real, imag, M, 0);
    fix_fft(real_org, imag_org, M, 0);

    printf("\nFFT\n");
    for (i=0; i<N; i++)
    {
        printf("%d: %d, %d", i, real[i], imag[i]);
		if (real[i] != real_org[i] || imag[i] != imag_org[i]) {
			printf(" expected (%d, %d)", real_org[i], imag_org[i]);
		}
        printf("\n");
    }


    //IFFT
    fix_fft(real, imag, M, 1);
    fix_fft(real_org, imag_org, M, 1);

    printf("\nIFFT\n");
    for (i=0; i<N; i++)
    {
        printf("%d: %d, %d", i, real[i], imag[i]);
		if (real[i] != real_org[i] || imag[i] != imag_org[i]) {
			printf(" expected (%d, %d)", real_org[i], imag_org[i]);
		}
        printf("\n");
    }


    return 0;
}


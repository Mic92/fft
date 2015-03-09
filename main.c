//main.c


#include        "fft.h"
#include        "fft-org.h"
#include        <stdio.h>
#include        <math.h>


#define M       8

//number of points
#define N       (1<<M)

#define M2    100
#define N2    (3<<M)

fixed real[N], imag[N];
fixed real_org[N], imag_org[N];

void show_result(fixed* real, fixed* real_org, fixed* imag, fixed* imag_org, int n)
{
    int i;
    for (i=0; i<n; i++)
    {
        printf("%d: %d, %d", i, real[i], imag[i]);
        if (real[i] != real_org[i] || imag[i] != imag_org[i]) {
            printf(" expected (%d, %d)", real_org[i], imag_org[i]);
        }
        printf("\n");
    }
}

int main()
{
    int i;

    for(i=0; i<N; i++)
    {
        real[i] = 1000*cos(i*2*3.1415926535/N);
        real_org[i] = real[i];
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
    fix_fft_org(real_org, imag_org, M, 0);

    printf("\nFFT\n");
    show_result(real, real_org, imag, imag_org, N);

    //IFFT
    fix_fft(real, imag, M, 1);
    fix_fft_org(real_org, imag_org, M, 1);

    printf("\nIFFT\n");
    show_result(real, real_org, imag, imag_org, N);

    for(i=0; i<N; i++)
    {
        real[i] = 1000*(cos(i*2*3.1415926535/N) + sin(i*2*3.1415926535/N));
        real_org[i] = real[i];
        imag[i] = 0;
        imag_org[i] = 0;
    }

    /*
    //FFT
    fix_fft(real, imag, M, 0);
    fix_fft_org(real_org, imag_org, M, 0);

    printf("\nFFT2\n");
    show_result(real, real_org, imag, imag_org, N);

    //IFFT
    fix_fft(real, imag, M, 1);
    fix_fft_org(real_org, imag_org, M, 1);

    printf("\nIFFT2\n");
    show_result(real, real_org, imag, imag_org, N);
    */

    return 0;
}


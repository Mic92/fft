//main.c

#include        "fft.h"
#include        "fft-org.h"
#include        <stdio.h>
#include        <math.h>

//defined by Compile Target
//#define M       3

//number of points
#define N       (1<<M)

//#define DEBUG 1

fixed real[N], imag[N];
fixed real_org[N], imag_org[N];

void show_result(fixed* real, fixed* real_org, fixed* imag, fixed* imag_org, int n)
{
    int i;
    for (i=0; i<n; i++)
    {
#ifdef DEBUG
        printf("%d: %d, %d", i, real[i], imag[i]);
        if (real[i] != real_org[i] || imag[i] != imag_org[i]) {
            printf(" expected (%d, %d)", real_org[i], imag_org[i]);
        }
        printf("\n");
#else
        if (real[i] != real_org[i] || imag[i] != imag_org[i])
            printf("got (%d, %d), expected (%d, %d)\n", real[i], imag[i], real_org[i], imag_org[i]);
#endif
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

    int scale, scale_org;
    
    //FFT
    scale = fix_fft(real, imag, M, 0);
    scale_org = fix_fft_org(real_org, imag_org, M, 0);

    printf("\nFFT\n");
    show_result(real, real_org, imag, imag_org, N);
    if (scale != scale_org) {
    	printf("got scale: %d, expected: %d\n", scale, scale_org);
    }

    //IFFT
    scale = fix_fft(real, imag, M, 1);
    scale_org = fix_fft_org(real_org, imag_org, M, 1);

    printf("\nIFFT\n");
    show_result(real, real_org, imag, imag_org, N);
    if (scale != scale_org) {
    	printf("got scale: %d, expected: %d\n", scale, scale_org);
    }

    printf("done\n");

    return 0;
}


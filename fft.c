/*  fft.c - Fixed-point Fast Fourier Transform  */
/*
    fix_fft()       perform FFT or inverse FFT
    fix_mpy()       perform fixed-point multiplication.
    Sinewave[1024]  sinewave normalized to 32767 (= 1.0).

    All data are fixed-point short integers, in which
    -32768 to +32768 represent -1.0 to +1.0. Integer arithmetic
    is used for speed, instead of the more natural floating-point.

    For the forward FFT (time -> freq), fixed scaling is
    performed to prevent arithmetic overflow, and to map a 0dB
    sine/cosine wave (i.e. amplitude = 32767) to two -6dB freq
    coefficients; the one in the lower half is reported as 0dB.

    For the inverse FFT (freq -> time), fixed scaling cannot be
    done, as two 0dB coefficients would sum to a peak amplitude of
    64K, overflowing the 32k range of the fixed-point integers.
    Thus, the fix_fft() routine performs variable scaling, and
    returns a value which is the number of bits LEFT by which
    the output must be shifted to get the actual amplitude
    (i.e. if fix_fft() returns 3, each value of fr[] and fi[]
    must be multiplied by 8 (2^3) for proper scaling.
    Clearly, this cannot be done within the fixed-point short
    integers. In practice, if the result is to be used as a
    filter, the scale_shift can usually be ignored, as the
    result will be approximately correctly normalized as is.


    Source Code taken by http://www.jjj.de/crs4: integer_fft.c
    Last Modified by Sebastian Haas at Oct. 2014.
*/

#include "fft.h"

#include <xtensa/tie/fft_inst.h>

/*
 *	fix_fft() - perform fast Fourier transform.
 *
 *  if n>0 FFT is done, if n<0 inverse FFT is done
 *	fr[n],fi[n] are real,imaginary arrays, INPUT AND RESULT.
 *	size of data = 2^m
 *  set inverse to 0=dft, 1=idft
 */
int fix_fft(fixed fr[], fixed fi[], int m, int inverse)
{
    int mr,nn,i,j,l,k,istep, n, scale, shift;

    fixed qr,qi;		//even input
    fixed tr,ti;		//odd input
    fixed wr,wi;		//twiddle factor

    //number of input data
    n = 1<<m;

    if(n > N_WAVE) return -1;

    mr = 0;
    nn = n - 1;
    scale = 0;
    
    int mm = m;
    xtbool t;

    /* decimation in time - re-order data */
    for(m=1; m<=nn; ++m) {
        mr = FFT_bit_reverse(m, mm);

        if(mr <= m) continue;
        tr = fr[m];
        fr[m] = fr[mr];
        fr[mr] = tr;

        ti = fi[m];
        fi[m] = fi[mr];
        fi[mr] = ti;
    }


    l = 1;
    k = LOG2_N_WAVE-1;
    while(l < n)
    {
        if(inverse)
        {
        	/* variable scaling, depending upon data */
        	shift = 0;
        	for(i=0; i<n; ++i)
        	{
        	    j = fr[i];
        	    if(j < 0) j = -j;

        	    m = fi[i];
        	    if(m < 0) m = -m;

        	    if(j > 16383 || m > 16383)
        	    {
        	        shift = 1;
        	        break;
        	    }
        	}
        	if(shift) ++scale;        	
        }
        else
        {
            /* fixed scaling, for proper normalization -
               there will be log2(n) passes, so this
               results in an overall factor of 1/n,
               distributed to maximize arithmetic accuracy. */
            shift = 1;
        }

        /* it may not be obvious, but the shift will be performed
           on each data point exactly once, during this pass. */
        istep = l << 1;		//step width of current butterfly
        for(m=0; m<l; ++m)
        {
            j = m << k;
            
            FFT_twiddle(wr, wi, j, shift, inverse);
            
            FFT_reg reg;
            fixed *reg_s = ((fixed*) &reg);
            
            WUR_FFT_loop(m);
            while(1)
            {
            	i = RUR_FFT_loop();
            	
            	FFT_loop_check(n, istep, t, j);
            	if(!t)
            		break;
                
                reg_s[3] = fr[i];
                reg_s[2] = fr[j];
                reg_s[1] = fi[i];
                reg_s[0] = fi[j];
                
                FFT_calc(reg, wr, wi, shift);
                
                fr[i] = reg_s[3];
                fr[j] = reg_s[2];
                fi[i] = reg_s[1];
                fi[j] = reg_s[0];
            }
        }
        --k;
        l = istep;
    }

    return scale;
}

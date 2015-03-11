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

#include <stdio.h>

#define FFT_COMBINED_STORE(_fr, _fi, _i, _simd_r, _simd_i) \
    asm ("{"                                  "\n" \
         "    fft_simd_store %1, %0, %2"      "\n" \
         "    nop"                            "\n" \
         "    fft_simd_store %3, %0, %4"      "\n" \
         "}" \
         :: "r" (_i), "r" (_fr), "r" (_simd_r), "r" (_fi), "r" (_simd_i));<e
		        		
#define FUSED_STORE(_p1, _v1, _p2, _v2) \
	asm ("{"                                  "\n" \
         "    s16i %1, %0, 0"      "\n" \
         "    nop"                            "\n" \
         "    s16i %3, %2, 0"      "\n" \
         "}" \
         :: "r" (_p1), "r" (_v1), "r" (_p2), "r" (_v2));

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
    int mr,nn,i,j,l,k,istep, n, scale;
    xtbool shift;

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

    /* decimation in time - re-order data */
    for(m=0;;) {
        FFT_BIT_REVERSE(m, mr, mm);

        if(m >= nn) break;
        if(mr <= m) continue;

        tr = fr[m];
        ti = fi[m];
        //fixed foo = fi[mr];
        //fixed bar = fr[mr];
        
        fi[m] = fi[mr];
        fr[m] = fr[mr];     
        //FUSED_STORE(fi + m, foo, fr + m, bar);
        
        FUSED_STORE(fi + mr, ti, fr + mr, tr);
    }

    l = 1;
    k = LOG2_N_WAVE-1;
    while(l < n)
    {
        if(inverse)
        {
        	/* variable scaling, depending upon data */
        	shift = 0;
        	for(i=0; i<n/8; i+=8)
        	{
        	    if(FFT_SHIFT_CHECK(fr, i) | FFT_SHIFT_CHECK(fi, i))
        	    {
        	        shift = 1;
        	        ++scale;
        	        break;
        	    }
        	}
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

        FFT_REG reg;
        FFT_REG_SIMD simd_r, simd_i, simd_r2, simd_i2;
        fixed *reg_s = ((fixed*) &reg);
        
        switch(l)
        {
	        case 1:
	        	for(i=0; i<n; i+=8)
	        	{
	        		simd_r = FFT_SIMD_LOAD(fr, i);
	        		simd_i = FFT_SIMD_LOAD(fi, i);
	    			FFT_SIMD_FIRST(simd_r, simd_i, shift);
	    			FFT_COMBINED_STORE(fr, fi, i, simd_r, simd_i);
	        	}
	        	break;
	        
	        case 2:
	        	for(i=0; i<n; i+=8)
	        	{
	        		simd_r = FFT_SIMD_LOAD(fr, i);
	        		simd_i = FFT_SIMD_LOAD(fi, i);
	    			FFT_SIMD_SECOND(simd_r, simd_i, shift, inverse);
	    			FFT_COMBINED_STORE(fr, fi, i, simd_r, simd_i);
	        	}
	        	break;
	        
	        case 4:
	        	WUR_FFT_SIMD_K(7);
	        	for(i=0; i<n; i+=8)
	        	{
	        		simd_r = FFT_SIMD_LOAD(fr, i);
	        		simd_i = FFT_SIMD_LOAD(fi, i);
	    			FFT_SIMD_THIRD(simd_r, simd_i, i, shift, inverse);
	    			FFT_COMBINED_STORE(fr, fi, i, simd_r, simd_i);
	        	}
	        	break;
	        
	        default:
	        	WUR_FFT_SIMD_K(k);
        		xtbool b0 = 0, b1 = 1;
		        for(m=0; m<n; m+=istep)
		        {
		        	for(i=m; i<m+l; i+=8)
		            {
		                j = i + l;
		                
		        		FFT_SIMD_LOAD_SHUFFLE(fr, i, simd_r, simd_r2, 0);
		        		FFT_SIMD_LOAD_SHUFFLE(fi, i, simd_i, simd_i2, 0);
		        		FFT_SIMD_LOAD_SHUFFLE(fr, j, simd_r, simd_r2, 1);
		        		FFT_SIMD_LOAD_SHUFFLE(fi, j, simd_i, simd_i2, 1);
		        		FFT_SIMD_THIRD(simd_r,  simd_i,  i,   shift, inverse);
		        		FFT_SIMD_THIRD(simd_r2, simd_i2, i+4, shift, inverse);
		        		
		        	    asm ("{"                                  "\n"
		        	         "    fft_simd_store_shuffle %1, %0, %2, %5, %7"  "\n"
		        	         "    nop"                            "\n"
		        	         "    fft_simd_store_shuffle %3, %0, %4, %6, %7"      "\n"
		        	         "}"                                  "\n"
		        	    :: "r" (i), "r" (fr), "r" (simd_r), "r" (fi), "r" (simd_i), "r" (simd_r2), "r" (simd_i2), "r" (b0));
		        	    asm ("{"                                  "\n"
		        	         "    fft_simd_store_shuffle %1, %0, %2, %5, %7"  "\n"
		        	         "    nop"                            "\n"
		        	         "    fft_simd_store_shuffle %3, %0, %4, %6, %7"      "\n"
		        	         "}"
		        	    :: "r" (j), "r" (fr), "r" (simd_r), "r" (fi), "r" (simd_i), "r" (simd_r2), "r" (simd_i2), "r" (b1));

		        	    /*
		        		FFT_SIMD_STORE_SHUFFLE(fr, i, simd_r, simd_r2, 0);
		        		FFT_SIMD_STORE_SHUFFLE(fi, i, simd_i, simd_i2, 0);
		        		FFT_SIMD_STORE_SHUFFLE(fr, j, simd_r, simd_r2, 1);
		        		FFT_SIMD_STORE_SHUFFLE(fi, j, simd_i, simd_i2, 1);
		        		*/
		        	}
		        }
	        	/*
		        for(m=0; m<n; m+=istep)
		        {
		        	for(i=m; i<m+l; ++i)
		            {
		                j = i + l;
		                
		                reg_s[3] = fr[i];
		                reg_s[2] = fr[j];
		                reg_s[1] = fi[i];
		                reg_s[0] = fi[j];
		                
		                FFT_CALC(reg, i << k, (xtbool) shift, inverse);
		                
		                fr[i] = reg_s[3];
		                fr[j] = reg_s[2];
		                fi[i] = reg_s[1];
		                fi[j] = reg_s[0];
		            }
		        }
		        */
		        break;
        }
        --k;
        l = istep;
    }

    return scale;
}

#include "dif.h"
#include "dit-flix.h"
#include "dif-flix.h"

#include <xtensa/tie/dit.h>
#include <xtensa/tie/dif.h>


int fix_dif_fft(fixed fr[], fixed fi[], int size, int _inverse)
{
    int i, j, l, k, istep, n, m;
    xtbool shift, inverse = _inverse;
    register FFT_REG_SIMD simd_r, simd_i, simd_r2, simd_i2;
    
    fixed tr,ti,qr,qi,wr,wi;

    //number of input data
    n = 1 << size;
    if(n > N_WAVE) return -1;
    
    int mr = 0, nn = n - 1;
    
    int scale = 0;

    scale = 0;

    l = n>>1;
    k = LOG2_N_WAVE-size;
    while(l > 0)
    {
        if(inverse)
        {
            /* variable scaling, depending upon data */
            shift = 0;
        	for (i = 0; i < (n / 8); i += 8)
        	{
        	    if (FFT_SHIFT_CHECK(fr, i) | FFT_SHIFT_CHECK(fi, i))	
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
        
        switch(l)
        {
	        case 1:
	        	for (i = 0; i < n; i += 8)
	        	{
	        		simd_r = FFT_SIMD_LOAD(fr, i);
	        		simd_i = FFT_SIMD_LOAD(fi, i);
	    			DIF_FFT_SIMD_FIRST(simd_r, simd_i, shift);
	    			FFT_FLIX_SIMD_STORE(fr, fi, i, simd_r, simd_i);
	        	}
	        	break;
	        
	        case 2:
	        	for (i = 0; i < n; i += 8)
	        	{
	        		simd_r = FFT_SIMD_LOAD(fr, i);
	        		simd_i = FFT_SIMD_LOAD(fi, i);
	        		DIF_FFT_SIMD_SECOND(simd_r, simd_i, shift, inverse);
	    			FFT_FLIX_SIMD_STORE(fr, fi, i, simd_r, simd_i);
	        	}
	        	break;
	        
	        case 4:
	        	WUR_FFT_SIMD_K(7);
	        	for (i = 0; i < n; i += 8)
	        	{
	        		simd_r = FFT_SIMD_LOAD(fr, i);
	        		simd_i = FFT_SIMD_LOAD(fi, i);
	        		DIF_FFT_SIMD_THIRD(i, simd_r, simd_i, shift, inverse);
	    			FFT_FLIX_SIMD_STORE(fr, fi, i, simd_r, simd_i);
	        	}
	        	break;
	        
	        default:
	        	WUR_FFT_SIMD_K(k);
		        for (m = 0; m < n; m += istep)
		        {
		        	j = m + l;
		        	for (i = m; i < m + l;)
		            {
		        		FFT_FLIX_SIMD_LOAD_SHUFFLE(i, fr, simd_r, simd_r2, fi, simd_i, simd_i2, 0);
		        		FFT_FLIX_SIMD_LOAD_SHUFFLE(j, fr, simd_r, simd_r2, fi, simd_i, simd_i2, 1);

		        		DIF_FFT_FLIX_SIMD_THIRD(i, simd_r, simd_r2, simd_i, simd_i2, shift, inverse);
		        		
		        		// inlines also i++ and j = i + l
		        		FFT_FLIX_SIMD_STORE_SHUFFLE(i, j, fr, simd_r, simd_r2, fi, simd_i, simd_i2);
		        	}
		        }
		        break;
        }

        ++k;
        l >>= 1;
    }
    
    /* decimation in frequency - re-order data */
    for (m = 1; m < n; m++) {
        int mr = FFT_BIT_REVERSE(m, size);

        if(mr <= m) continue;
        int ti = fi[m];
        int tr = fr[m];
        
        FLIX_S16I(&fr[m], fr[mr], &fi[m], fi[mr]);
        FLIX_S16I(&fi[mr], ti, &fr[mr], tr);   
     }

    return scale;
}

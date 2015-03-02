
#ifndef FFT_H
#define FFT_H

/* FIX_MPY() - fixed-point multiplication macro.
   This macro is a statement, not an expression (uses asm).
   BEWARE: make sure _DX is not clobbered by evaluating (A) or DEST.
   args are all of type fixed.
   Scaling ensures that 32767*32767 = 32767. */
#define FIX_MPY(DEST,A,B)       DEST = ((long)(A) * (long)(B))>>15

#define N_WAVE          1024    /* dimension of Sinewave[] */
#define LOG2_N_WAVE     10      /* log2(N_WAVE) */

#ifndef fixed
#define fixed short
#endif

extern fixed Sinewave_org[N_WAVE];

//function prototypes
fixed fix_mpy_org(fixed a, fixed b);
int fix_fft_org(fixed *fr, fixed *fi, int m, int inverse);


#endif	//FFT_H

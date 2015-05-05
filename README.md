# FFT
FFT Implementation for Tensilica DSP Processors

Written for the [Lab Hardware/Software-Codesign](https://mns.ifn.et.tu-dresden.de/Teaching/Courses/Pages/P-HWSW-Codesign.aspx) in 2015 at TU Dresden.

Cycles needed for N-Point FFT with TIE instructions, SIMD, pipelining and FLIX
(-O2, compiled with feedback optimization):
- N=8: 327
- N=256: 7 612
- N=1024: 32 050

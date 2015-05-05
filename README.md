# FFT
FFT Implementation for Tensilica DSP Processors

Written for the [Lab Hardware/Software-Codesign](https://mns.ifn.et.tu-dresden.de/Teaching/Courses/Pages/P-HWSW-Codesign.aspx) in 2015 in the TU Dresden.

Cycles needed with TIE instructions, SIMD, pipelining and FLIX
(-O2, compiled with feedback optimization):
- M=3 (N=8): 327
- M=8 (N=256): 7 612
- M=10 (N=1024): 32 050

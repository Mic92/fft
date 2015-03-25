#ifdef FLIX
#define DIF_FFT_FLIX_SIMD_THIRD(i, simd_r, simd_r2, simd_i, simd_i2, shift, inverse) \
	do { \
		/* also i += 4 in the next cycle */ \
		asm ("{                                                                     \n" \
		     "    dif_fft_simd_third %0, %1, %2, %3, %4                                 \n" \
		     "    addi %0, %0, 4                                                    \n" \
		     "    nop                                                               \n" \
		     "}" \
		     :: "r" (i), "r" (simd_r), "r" (simd_i), "r" (shift), "r" (inverse)); \
		/* also i -= 4 in the next cycle */ \
		asm ("{                                                                     \n" \
		     "    dif_fft_simd_third %0, %1, %2, %3, %4                                 \n" \
		     "    addi %0, %0, -4                                                   \n" \
		     "    nop                                                               \n" \
		     "}" \
		     :: "r" (i), "r" (simd_r2), "r" (simd_i2), "r" (shift), "r" (inverse)); \
    } while(0)
#else
#define DIF_FFT_FLIX_SIMD_THIRD(i, simd_r, simd_r2, simd_i, simd_i2, shift, inverse) \
	DIF_FFT_SIMD_THIRD(i, simd_r, simd_i, shift, inverse); \
	DIF_FFT_SIMD_THIRD(i + 4, simd_r2, simd_i2, shift, inverse);
#endif

#ifndef FFT_ASM_H
#define FFT_ASM_H

#define FLIX

#ifdef FLIX
#define FFT_FLIX_SIMD_STORE(_fr, _fi, _i, _simd_r, _simd_i) \
        asm ("{                                         \n" \
             "    fft_simd_store %1, %0, %2             \n" \
             "    nop                                   \n" \
             "    fft_simd_store %3, %0, %4             \n" \
             "}" \
             :: "r" (_i), "r" (_fr), "r" (_simd_r), "r" (_fi), "r" (_simd_i));

#define FLIX_S16I(_p1, _v1, _p2, _v2) \
        asm ("{                         \n" \
             "    s16i %1, %0, 0        \n" \
             "    nop                   \n" \
             "    s16i %3, %2, 0        \n" \
             "}" \
             :: "r" (_p1), "r" (_v1), "r" (_p2), "r" (_v2) \
             : "memory");

#define FFT_FLIX_SIMD_THIRD(i, simd_r, simd_r2, simd_i, simd_i2, shift, inverse) \
        do { \
                /* also i += 4 in the next cycle */ \
                asm ("{                                                      \n" \
                                "    fft_simd_third %0, %1, %2, %3, %4       \n" \
                                "    addi %0, %0, 4                          \n" \
                                "    nop                                     \n" \
                                "}" \
                                :: "r" (i), "r" (simd_r), "r" (simd_i), "r" (shift), "r" (inverse)); \
                /* also i -= 4 in the next cycle */ \
                asm ("{                                                                          \n" \
                                "    fft_simd_third %0, %1, %2, %3, %4                           \n" \
                                "    addi %0, %0, -4                                             \n" \
                                "    nop                                                         \n" \
                                "}" \
                                :: "r" (i), "r" (simd_r2), "r" (simd_i2), "r" (shift), "r" (inverse)); \
        } while(0)
#define FFT_FLIX_SIMD_STORE_SHUFFLE(i, j, fr, simd_r, simd_r2, fi, simd_i, simd_i2) \
        do { \
                /* also i++ */ \
                asm ("{                                                         \n" \
                     "    fft_simd_store_shuffle %1, %0, %2, %5, 0              \n" \
                     "    addi %0, %0, 8                                        \n" \
                     "    fft_simd_store_shuffle %3, %0, %4, %6, 0              \n" \
                     "}" \
                     /*  0         1         2             3         4*/ \
                     : "=r" (i) : "r" (fr), "r" (simd_r), "r" (fi), "r" (simd_i), \
                     "r" (simd_r2), "r" (simd_i2)); \
                   /* 5              6*/ \
                /* also j = i + l; */ \
                asm ("{                                                        \n" \
                     "    fft_simd_store_shuffle %3, %0, %4, %7, 1             \n" \
                     "    add.n %0, %2, %9                                     \n" \
                     "    fft_simd_store_shuffle %5, %0, %6, %8, 1             \n" \
                     "}" \
                     /*  0            1         2              3              4*/ \
                     : "=r" (j)    : "r" (j),  "r" (i),       "r" (fr),      "r" (simd_r), \
                     "r" (fi), "r" (simd_i), "r" (simd_r2), "r" (simd_i2), "r" (l)); \
                         /* 5             6         7              8              9*/ \
        } while(0)
#else

#define FFT_FLIX_SIMD_STORE(_fr, _fi, _i, _simd_r, _simd_i) \
        FFT_SIMD_STORE(fr, i, simd_r); \
FFT_SIMD_STORE(fi, i, simd_i); \

#define FFT_FLIX_SIMD_THIRD(i, simd_r, simd_r2, simd_i, simd_i2, shift, inverse) \
        FFT_SIMD_THIRD(i, simd_r, simd_i, shift, inverse); \
FFT_SIMD_THIRD(i + 4, simd_r2, simd_i2, shift, inverse);
#define FLIX_S16I(_p1, _v1, _p2, _v2) \
        do { \
                *(_p1) = (_v1); \
                *(_p2) = (_v2); \
        } while(0)
#define FFT_FLIX_SIMD_STORE_SHUFFLE(i, j, fr, simd_r, simd_r2, fi, simd_i, simd_i2) \
        do { \
                FFT_SIMD_STORE_SHUFFLE(fr, i, simd_r, simd_r2, 0); \
                FFT_SIMD_STORE_SHUFFLE(fi, i, simd_i, simd_i2, 0); \
                FFT_SIMD_STORE_SHUFFLE(fr, j, simd_r, simd_r2, 1); \
                FFT_SIMD_STORE_SHUFFLE(fi, j, simd_i, simd_i2, 1); \
                i++; \
                j = i + l; \
        } while(0)

#endif // FLIX

#define FFT_FLIX_SIMD_LOAD_SHUFFLE(i, fr, simd_r, simd_r2, fi, simd_i, simd_i2, high) \
        do { \
                FFT_SIMD_LOAD_SHUFFLE(fr, i, simd_r, simd_r2, high); \
                FFT_SIMD_LOAD_SHUFFLE(fi, i, simd_i, simd_i2, high); \
        } while(0)

#endif // FFT_ASM_H

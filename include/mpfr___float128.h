// convert __float128 <--> mpfr_t
//http://permalink.gmane.org/gmane.comp.lib.mpfr.general/769
//

#include <stdint.h>
#include <mpfr.h>

// Frankly, if you have __float128, you have 64-bit integers, right?
//XXX
#if GMP_LIMB_BITS != 64
# error "cannot compile"
#endif

//XXX
#define __LITTLE_ENDIAN__ 1	//Caution! MacOSX PPC is Big endian

typedef union {
    __float128 x;

    struct {
#if __BIG_ENDIAN__
	unsigned negative:1;
	unsigned exponent:15;
	uint64_t mant_high:48;
	uint64_t mant_low:64;
#endif
#if __LITTLE_ENDIAN__
	uint64_t mant_low:64;
	uint64_t mant_high:48;
	unsigned exponent:15;
	unsigned negative:1;
#endif
    } ieee;

    struct {
	uint8_t b[16];
    } bytes;

    struct {
#if __BIG_ENDIAN__
	uint64_t high;
	uint64_t low;
#endif
#if __LITTLE_ENDIAN__
	uint64_t low;
	uint64_t high;
#endif
    } words64;

    struct {
#if __BIG_ENDIAN__
	unsigned negative:1;
	unsigned exponent:15;
	unsigned quiet_nan:1;
	uint64_t mant_high:47;
	uint64_t mant_low:64;
#endif
#if __LITTLE_ENDIAN__
	uint64_t mant_low:64;
	uint64_t mant_high:47;
	unsigned quiet_nan:1;
	unsigned exponent:15;
	unsigned negative:1;
#endif
    } nan;

} ieee754_float128;

#define IEEE754_FLOAT128_BIAS 0x3fff

inline static void mpfr_set_float128(mpfr_t * a, const __float128 val)
{
    ieee754_float128 f;
    f.x = val;

    int exp = f.ieee.exponent;
    int subnormal = 0;
    int numoflimbs;

    uint64_t h = f.ieee.mant_high;
    uint64_t l = f.ieee.mant_low;
    (*a)-> _mpfr_sign = f.ieee.negative ? -1 : +1;
    // Special case for zero exponent
    if (exp == 0) {
	if (h == 0 && l == 0) {
	    (*a)->_mpfr_exp = __MPFR_EXP_ZERO;
	    return;
	} else
	    subnormal = 1;
    } else if (exp == 0x7fff) {
	if (h == 0 && l == 0) {
	    (*a)->_mpfr_exp = __MPFR_EXP_INF;
	    return;
	} else {
	    (*a)->_mpfr_exp = __MPFR_EXP_NAN;
	    return;
	}
    } else
	(*a)->_mpfr_exp = exp - (IEEE754_FLOAT128_BIAS - 1);
    numoflimbs = ((*a)->_mpfr_prec - 1 )/ GMP_LIMB_BITS + 1;
    switch ( numoflimbs - 1 ) {
        case 0:
            (*a)->_mpfr_d[0] =  1UL << 63 | h << 15 | l >> 49;
        break;
        case 1:
            if (subnormal)
	        (*a)->_mpfr_d[1] = h << 15 | l >> 49;
            else
                (*a)->_mpfr_d[1] = 1UL << 63 | h << 15 | l >> 49;
            (*a)->_mpfr_d[0] = l << 15;
        break;
        default:
           for (int i = 0; i < numoflimbs - 2 ; i++ ) (*a)->_mpfr_d[i] = 0L;
           (*a)->_mpfr_d[numoflimbs-2] = l << 15;
           (*a)->_mpfr_d[numoflimbs-1] = 1UL << 63 | h << 15 | l >> 49;
        break;
   }
}

inline static __float128 mpfr_get_float128(const mpfr_t a)
{
    ieee754_float128 f;
    int numoflimbs;

    if (mpfr_nan_p(a)) {
	f.ieee.exponent = 0x7fff;
	f.ieee.mant_high = f.ieee.mant_low = 1;
	return f.x;
    }

    if (mpfr_inf_p(a)) {
	f.ieee.exponent = 0x7fff;
	f.ieee.mant_high = f.ieee.mant_low = 0;
	f.ieee.negative = (MPFR_SIGN(a) == -1) ? 1 : 0;
	return f.x;
    }

    if (mpfr_zero_p(a)) {
	f.ieee.exponent = 0;
	f.ieee.mant_high = f.ieee.mant_low = 0;
	f.ieee.negative = (MPFR_SIGN(a) == -1) ? 1 : 0;
	return f.x;
    }
    // BIG TODO: deal with denormals

    f.ieee.negative = (MPFR_SIGN(a) == -1) ? 1 : 0;
    f.ieee.exponent = a->_mpfr_exp + (IEEE754_FLOAT128_BIAS - 1);
    numoflimbs = (a->_mpfr_prec - 1 )/ GMP_LIMB_BITS + 1;
    uint64_t x = a->_mpfr_d[numoflimbs-1];
    uint64_t y = a->_mpfr_d[numoflimbs-2];
    f.ieee.mant_high = x >> 15;	//f.ieee.mant_high = (x << 1) >> 16;
    f.ieee.mant_low = (x << 49) | (y >> 15);

    return f.x;
}

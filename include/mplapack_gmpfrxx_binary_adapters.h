#ifndef MPLAPACK_GMPFRXX_BINARY_ADAPTERS_H
#define MPLAPACK_GMPFRXX_BINARY_ADAPTERS_H

#ifndef __STDC_WANT_IEC_60559_TYPES_EXT__
#define __STDC_WANT_IEC_60559_TYPES_EXT__ 1
#endif
#ifndef __STDC_WANT_IEC_60559_FUNCS_EXT__
#define __STDC_WANT_IEC_60559_FUNCS_EXT__ 1
#endif

#include <mplapack_config.h>

#include <mplapack_gmpfrxx_mkII_config.h>

#include <gmpfrxx_mkII/adapters/binary80_complex.hpp>
#include <gmpfrxx_mkII/adapters/binary128_complex.hpp>

#include <cstdio>
#include <cstdlib>

#if MPLAPACK_BINARY128_IO == MPLAPACK_BINARY128_IO_QUADMATH_SNPRINTF
#include <quadmath.h>
#endif

namespace mplapack {
namespace gmpfrxx_adapter {

inline mpfr_rnd_t default_mpfr_rounding()
{
    return mpfrxx::mpfr_class::default_rounding();
}

inline mpfrxx::mpfr_class make_default_mpfr()
{
    return mpfrxx::mpfr_class::with_precision(
        mpfrxx::mpfr_class::default_precision());
}

#if MPLAPACK_BINARY80_MODE != MPLAPACK_BINARY80_MODE_DISABLED
using binary80_native_type = mplapack_binary80_t;
using binary80_source_type =
    gmpfrxx_mkII::adapters::binary80_source<binary80_native_type>;
using binary80_complex_source_type =
    gmpfrxx_mkII::adapters::binary80_complex_source<binary80_native_type>;

static_assert(gmpfrxx_mkII::detail::is_binary80_native_v<binary80_native_type>,
              "MPLAPACK binary80 type has no gmpfrxx_mkII import path");

inline binary80_source_type make_binary80_source(binary80_native_type value)
{
    return gmpfrxx_mkII::adapters::make_binary80_source(value);
}

inline binary80_complex_source_type
make_binary80_complex_source(binary80_native_type real,
                             binary80_native_type imag)
{
    return gmpfrxx_mkII::adapters::make_binary80_complex_source(real, imag);
}

inline binary80_native_type export_binary80(const mpfrxx::mpfr_class &value)
{
    const long double native = mpfr_get_ld(
        value.mpfr_data(), default_mpfr_rounding());
    return static_cast<binary80_native_type>(native);
}
#endif

#if MPLAPACK_BINARY128_MODE != MPLAPACK_BINARY128_MODE_DISABLED
using binary128_native_type = mplapack_binary128_t;

#if MPLAPACK_BINARY128_MODE == MPLAPACK_BINARY128_MODE_LDBL || \
    GMPFRXX_MKII_ADAPTERS_HAVE_MPFR_FLOAT128
#define MPLAPACK_GMPFRXX_BINARY128_DIRECT 1
#else
#define MPLAPACK_GMPFRXX_BINARY128_DIRECT 0
#endif

#if MPLAPACK_GMPFRXX_BINARY128_DIRECT
using binary128_source_type =
    gmpfrxx_mkII::adapters::binary128_source<binary128_native_type>;
using binary128_complex_source_type =
    gmpfrxx_mkII::adapters::binary128_complex_source<binary128_native_type>;

static_assert(gmpfrxx_mkII::detail::is_binary128_native_v<binary128_native_type>,
              "MPLAPACK binary128 type has no gmpfrxx_mkII import path");

inline binary128_source_type make_binary128_source(binary128_native_type value)
{
    return gmpfrxx_mkII::adapters::make_binary128_source(value);
}

inline binary128_complex_source_type
make_binary128_complex_source(binary128_native_type real,
                              binary128_native_type imag)
{
    return gmpfrxx_mkII::adapters::make_binary128_complex_source(real, imag);
}
#endif

inline void binary128_to_decimal(char *buf, size_t buflen,
                                 binary128_native_type value)
{
#if MPLAPACK_BINARY128_IO == MPLAPACK_BINARY128_IO_SNPRINTF_LDBL
    std::snprintf(buf, buflen, "%.40Le", static_cast<long double>(value));
#elif MPLAPACK_BINARY128_IO == MPLAPACK_BINARY128_IO_STRFROMF128
    strfromf128(buf, buflen, "%.40e", static_cast<_Float128>(value));
#elif MPLAPACK_BINARY128_IO == MPLAPACK_BINARY128_IO_QUADMATH_SNPRINTF
    quadmath_snprintf(buf, buflen, "%.40Qe", static_cast<__float128>(value));
#else
#error "Unsupported MPLAPACK_BINARY128_IO value"
#endif
}

inline binary128_native_type binary128_from_decimal(const char *buf)
{
#if MPLAPACK_BINARY128_IO == MPLAPACK_BINARY128_IO_SNPRINTF_LDBL
    return static_cast<binary128_native_type>(strtold(buf, nullptr));
#elif MPLAPACK_BINARY128_IO == MPLAPACK_BINARY128_IO_STRFROMF128
    return static_cast<binary128_native_type>(strtof128(buf, nullptr));
#elif MPLAPACK_BINARY128_IO == MPLAPACK_BINARY128_IO_QUADMATH_SNPRINTF
    return static_cast<binary128_native_type>(strtoflt128(buf, nullptr));
#else
#error "Unsupported MPLAPACK_BINARY128_IO value"
#endif
}

inline mpfrxx::mpfr_class binary128_to_mpfr(binary128_native_type value)
{
#if MPLAPACK_GMPFRXX_BINARY128_DIRECT
    return mpfrxx::mpfr_class(make_binary128_source(value));
#else
    char buf[128];
    binary128_to_decimal(buf, sizeof(buf), value);
    mpfrxx::mpfr_class result = make_default_mpfr();
    mpfr_set_str(result.mpfr_data(), buf, 10, default_mpfr_rounding());
    return result;
#endif
}

inline binary128_native_type export_binary128(const mpfrxx::mpfr_class &value)
{
#if MPLAPACK_BINARY128_MODE == MPLAPACK_BINARY128_MODE_LDBL
    return mpfr_get_ld(value.mpfr_data(),
                       default_mpfr_rounding());
#elif GMPFRXX_MKII_ADAPTERS_HAVE_MPFR_FLOAT128
    return static_cast<binary128_native_type>(mpfr_get_float128(
        value.mpfr_data(), default_mpfr_rounding()));
#else
    char buf[128];
    mpfr_snprintf(buf, sizeof(buf), "%.40Re", value.mpfr_data());
    return binary128_from_decimal(buf);
#endif
}
#endif

} // namespace gmpfrxx_adapter
} // namespace mplapack

#if MPLAPACK_BINARY80_MODE != MPLAPACK_BINARY80_MODE_DISABLED
inline mplapack_binary80_t cast2binary80_t(const mpfrxx::mpfr_class &value)
{
    return mplapack::gmpfrxx_adapter::export_binary80(value);
}
#endif

#if MPLAPACK_BINARY128_MODE != MPLAPACK_BINARY128_MODE_DISABLED
inline mplapack_binary128_t cast2binary128_t(const mpfrxx::mpfr_class &value)
{
    return mplapack::gmpfrxx_adapter::export_binary128(value);
}
#endif

#endif // MPLAPACK_GMPFRXX_BINARY_ADAPTERS_H

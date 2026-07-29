#ifndef MPLAPACK_GMPFRXX_BINARY_ADAPTERS_H
#define MPLAPACK_GMPFRXX_BINARY_ADAPTERS_H

#include <mplapack_config.h>

#include <mplapack_gmpfrxx_mkII_config.h>

#include <gmpfrxx_mkII/adapters/binary80_complex.hpp>
#include <gmpfrxx_mkII/adapters/binary128_complex.hpp>

namespace mplapack {
namespace gmpfrxx_adapter {

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
        value.mpfr_data(), mpfrxx::mpfr_class::default_rounding());
    return static_cast<binary80_native_type>(native);
}
#endif

#if MPLAPACK_BINARY128_MODE != MPLAPACK_BINARY128_MODE_DISABLED
using binary128_native_type = mplapack_binary128_t;
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

inline binary128_native_type export_binary128(const mpfrxx::mpfr_class &value)
{
#if MPLAPACK_BINARY128_MODE == MPLAPACK_BINARY128_MODE_LDBL
    return mpfr_get_ld(value.mpfr_data(),
                       mpfrxx::mpfr_class::default_rounding());
#else
    return static_cast<binary128_native_type>(mpfr_get_float128(
        value.mpfr_data(), mpfrxx::mpfr_class::default_rounding()));
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

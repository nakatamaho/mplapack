#ifndef FEM_READ_HPP
#define FEM_READ_HPP

// Ensure MPLAPACK utils expose sprintnum()/sprintnum_short() and MPLAPACK_BUFLEN.
// In mplapack_utils_*.h these are guarded by MPLAPACK_INTERNAL.
#ifndef MPLAPACK_INTERNAL
#define MPLAPACK_INTERNAL 1
#endif

// MPLAPACK backend utilities (printnum/sprintnum, precision, buffer length, etc.)
#if defined(MPLAPACK_BUILD_WITH_GMP)
#include "mplapack_utils_gmp.h"
#elif defined(MPLAPACK_BUILD_WITH_MPFR)
#include "mplapack_utils_mpfr.h"
#elif defined(MPLAPACK_BUILD_WITH_BINARY128)
#include "mplapack_utils_binary128.h"
#elif defined(MPLAPACK_BUILD_WITH_BINARY80)
#include "mplapack_utils_binary80.h"
#elif defined(MPLAPACK_BUILD_WITH_DD)
#include "mplapack_utils_dd.h"
#elif defined(MPLAPACK_BUILD_WITH_QD)
#include "mplapack_utils_qd.h"
#elif defined(MPLAPACK_BUILD_WITH_DOUBLE)
#include "mplapack_utils_double.h"
#else
#error "No MPLAPACK backend macro is defined (MPLAPACK_BUILD_WITH_*)."
#endif

#include <fem/common.hpp>
#include <fem/format.hpp>
#include <fem/star.hpp>
#include <fem/str_arr_ref.hpp>
#include <fem/utils/misc.hpp>
#include <fem/utils/string_to_double_fmt.hpp>
#include <string>
#include <cstdint>
#include <cstdlib>
#include <type_traits>
#include <mplapack_gmpfrxx_mkII_config.h>
#if defined(MPLAPACK_BUILD_WITH_GMP)
#if __has_include(<gmpxx_mkII.h>)
#include <gmpxx_mkII.h>
using namespace gmpxx;
#endif
#endif
#if defined(MPLAPACK_BUILD_WITH_MPFR)
#if __has_include(<mpfrxx_mkII.h>) && __has_include(<mpcxx_mkII.h>)
#include <mpfrxx_mkII.h>
#include <mpcxx_mkII.h>
using namespace mpfrxx;
#endif
#endif
#if defined(MPLAPACK_BUILD_WITH_QD) || defined(MPLAPACK_BUILD_WITH_DD)
// QD headers define and use qd::nint (and other short identifiers) inside the headers.
// MPLAPACK (or other code) may define macros like `nint`, which would macro-expand
// `qd::nint` into `qd::__dd_nint` and break the QD headers.
// Temporarily disable such macros while including QD headers.
#if defined(nint)
#pragma push_macro("nint")
#undef nint
#define FEM_RESTORE_nint 1
#endif
#if defined(min)
#pragma push_macro("min")
#undef min
#define FEM_RESTORE_min 1
#endif
#if defined(max)
#pragma push_macro("max")
#undef max
#define FEM_RESTORE_max 1
#endif
#if defined(abs)
#pragma push_macro("abs")
#undef abs
#define FEM_RESTORE_abs 1
#endif
#if defined(sign)
#pragma push_macro("sign")
#undef sign
#define FEM_RESTORE_sign 1
#endif
#if __has_include(<qd/dd_real.h>)
#include <qd/dd_real.h>
#endif
#if __has_include(<qd/qd_real.h>)
#include <qd/qd_real.h>
#endif
#if __has_include(<dd_complex.h>)
#include <dd_complex.h>
#endif
#if __has_include(<qd_complex.h>)
#include <qd_complex.h>
#endif
#if defined(FEM_RESTORE_sign)
#pragma pop_macro("sign")
#undef FEM_RESTORE_sign
#endif
#if defined(FEM_RESTORE_abs)
#pragma pop_macro("abs")
#undef FEM_RESTORE_abs
#endif
#if defined(FEM_RESTORE_max)
#pragma pop_macro("max")
#undef FEM_RESTORE_max
#endif
#if defined(FEM_RESTORE_min)
#pragma pop_macro("min")
#undef FEM_RESTORE_min
#endif
#if defined(FEM_RESTORE_nint)
#pragma pop_macro("nint")
#undef FEM_RESTORE_nint
#endif

#endif
#define IOSTAT_OK 0
#define IOSTAT_ERROR 1
#define IOSTAT_END -1
namespace fem {
class read_loop // TODO copy-constructor potential performance problem
{
  private:
    utils::slick_ptr<utils::simple_istream> inp;
    bool first_inp_get;
    format::token_loop fmt_loop;
    bool blanks_zero;
    int exp_scale;
    io_modes io_mode;
    int *iostat_ptr;
    static inline void normalize_fortran_exponent(std::string &s) {
        for (char &ch: s) {
            if (ch == 'D' || ch == 'd') {
                ch = 'E';
            }
        }
    }
    template <class T> static inline void assign_from_token_string(T &val, std::string s) {
        normalize_fortran_exponent(s);
        if constexpr (std::is_constructible<T, const char *>::value) {
            val = T(s.c_str());
        } else if constexpr (std::is_assignable<T &, double>::value) {
            long double ld = std::strtold(s.c_str(), nullptr);
            val = static_cast<double>(ld);
        } else if constexpr (std::is_constructible<T, double>::value) {
            long double ld = std::strtold(s.c_str(), nullptr);
            val = T(static_cast<double>(ld));
        } else {
            // As a last resort, try long double via double conversion.
            long double ld = std::strtold(s.c_str(), nullptr);
            val = T(static_cast<double>(ld));
        }
    }

  public:
    read_loop(common &cmn, int const &unit, unformatted_type const &) : inp(cmn.io.simple_istream(unit)), first_inp_get(true), blanks_zero(false), exp_scale(0), io_mode(io_unformatted), iostat_ptr(0) {}
    read_loop(common &cmn, int const &unit, star_type const &) : inp(cmn.io.simple_istream(unit)), first_inp_get(true), blanks_zero(false), exp_scale(0), io_mode(io_list_directed), iostat_ptr(0) {}
    read_loop(common &cmn, int const &unit, str_cref fmt) : inp(cmn.io.simple_istream(unit)), first_inp_get(true), fmt_loop(fmt), blanks_zero(false), exp_scale(0), io_mode(io_formatted), iostat_ptr(0) {}
    read_loop(str_cref const &internal_file, star_type const &) : inp(utils::slick_ptr<utils::simple_istream>(new utils::simple_istream_from_char_ptr_and_size(internal_file.elems(), internal_file.len()))), first_inp_get(true), blanks_zero(false), exp_scale(0), io_mode(io_list_directed), iostat_ptr(0) {}
    read_loop(str_cref const &internal_file, str_cref fmt) : inp(utils::slick_ptr<utils::simple_istream>(new utils::simple_istream_from_char_ptr_and_size(internal_file.elems(), internal_file.len()))), first_inp_get(true), fmt_loop(fmt), blanks_zero(false), exp_scale(0), io_mode(io_formatted), iostat_ptr(0) {}
    read_loop &rec(int const &) {
        inp.reset();
        throw TBXX_NOT_IMPLEMENTED();
    }
    read_loop &iostat(int &iostat) {
        this->iostat_ptr = &iostat;
        iostat = IOSTAT_OK;
        return *this;
    }
    std::string const &next_edit_descriptor() {
        while (true) {
            utils::token const *t = fmt_loop.next_executable_token();
            std::string const &tv = t->value;
            if (t->type == "string") {
                inp.reset();
                throw TBXX_NOT_IMPLEMENTED();
            } else if (t->type == "op") {
                if (tv[0] == ':') {
                    // ignored
                } else if (tv[0] == '/') {
                    skip_to_end_of_line();
                    int c = inp_get();
                    if (utils::is_stream_end(c)) {
                        inp.reset();
                        if (this->iostat_ptr != 0)
                            *iostat_ptr = IOSTAT_END;
                        throw read_end("End of input while reading string");
                    }
                    inp->backup();
                } else if (tv[0] == '$') {
                    inp.reset();
                    throw TBXX_NOT_IMPLEMENTED();
                } else {
                    inp.reset();
                    throw TBXX_UNREACHABLE_ERROR();
                }
            } else if (t->type == "format") {
                if (utils::ends_with_char(tv, 'x')) {
                    if (tv.size() == 1) {
                        process_fmt_x(1);
                    } else {
                        process_fmt_x(utils::signed_integer_value(tv.data(), 0, tv.size() - 1));
                    }
                } else if (std::strchr("adefgilz", tv[0]) != 0) {
                    return tv;
                } else if (utils::ends_with_char(tv, 'p')) {
                    if (tv.size() == 1) {
                        exp_scale = 1;
                    } else {
                        exp_scale = utils::signed_integer_value(tv.data(), 0, tv.size() - 1);
                    }
                } else if (tv[0] == 't') {
                    inp.reset();
                    throw TBXX_NOT_IMPLEMENTED();
                } else if (tv[0] == 's') {
                    inp.reset();
                    throw TBXX_NOT_IMPLEMENTED();
                } else if (tv[0] == 'b') {
                    blanks_zero = (tv[1] == 'z');
                } else {
                    inp.reset();
                    throw TBXX_UNREACHABLE_ERROR();
                }
            } else {
                inp.reset();
                throw TBXX_UNREACHABLE_ERROR();
            }
        }
    }
    int inp_get() {
        int result = 0;
        result = inp->get();
        if (utils::is_stream_err(result)) {
            inp.reset();
            if (this->iostat_ptr != 0)
                *iostat_ptr = IOSTAT_ERROR;
            throw io_err("Error during read");
        }
        if (first_inp_get || io_mode == io_unformatted) {
            first_inp_get = false;
            if (utils::is_stream_end(result)) {
                inp.reset();
                if (this->iostat_ptr != 0)
                    *iostat_ptr = IOSTAT_END;
                throw read_end("End of input during read");
            }
        }
        if (io_mode == io_formatted && result == '\r') {
            int next_char = inp->get();
            if (next_char == '\n') {
                result = '\n';
            } else {
                inp->backup();
            }
        }
        return result;
    }
    void process_fmt_x(unsigned n) {
        for (unsigned i = 0; i < n; i++) {
            int c = inp_get();
            if (c == utils::stream_end) {
                return;
            }
            if (utils::is_end_of_line(c)) {
                inp->backup();
                return;
            }
        }
    }
    read_loop &operator,(char &val) {
        inp.reset();
        throw TBXX_NOT_IMPLEMENTED();
        return *this;
    }
    read_loop &operator,(bool &val) {
        if (io_mode == io_unformatted) {
            from_stream_unformatted(reinterpret_cast<char *>(&val), sizeof(bool));
            return *this;
        }

        auto parse_logical = [&](std::string s) -> bool {
            // Trim spaces
            auto l = s.find_first_not_of(' ');
            if (l == std::string::npos) {
                throw io_err("Empty token while reading logical value");
            }
            auto r = s.find_last_not_of(' ');
            s = s.substr(l, r - l + 1);

            // Remove surrounding dots, e.g. ".TRUE." -> "TRUE"
            if (!s.empty() && s.front() == '.') {
                s.erase(s.begin());
            }
            if (!s.empty() && s.back() == '.') {
                s.pop_back();
            }

            // Uppercase
            for (char &c: s) {
                c = static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
            }

            if (s == "T" || s == "TRUE")
                return true;
            if (s == "F" || s == "FALSE")
                return false;

            // Also accept "1"/"0" defensively
            if (s == "1")
                return true;
            if (s == "0")
                return false;

            throw io_err("Invalid token while reading logical value: " + s);
        };

        if (io_mode == io_list_directed) {
            std::string s = read_star_token_string();
            val = parse_logical(s);
            return *this;
        }

        // io_formatted: Lw (logical) or fall back to list-directed tokenization
        std::string const &ed = next_edit_descriptor();
        int n = static_cast<int>(ed.size());
        if (n >= 2 && ed[0] == 'l') {
            int iw = utils::unsigned_integer_scan(ed.data(), 1, ed.size());
            int w = utils::unsigned_integer_value(ed.data(), 1, iw);
            std::string s = read_fmt_token_string(w);
            val = parse_logical(s);
        } else {
            std::string s = read_star_token_string();
            val = parse_logical(s);
        }
        return *this;
    }
    read_loop &operator,(integer_star_1 &val) {
        if (io_mode == io_unformatted) {
            inp.reset();
            throw TBXX_NOT_IMPLEMENTED();
        } else if (io_mode == io_list_directed) {
            inp.reset();
            throw TBXX_NOT_IMPLEMENTED();
        } else {
            inp.reset();
            throw TBXX_NOT_IMPLEMENTED();
        }
        return *this;
    }
    read_loop &operator,(integer_star_2 &val) {
        if (io_mode == io_unformatted) {
            from_stream_unformatted(reinterpret_cast<char *>(&val), sizeof(integer_star_2));
        } else if (io_mode == io_list_directed) {
            inp.reset();
            throw TBXX_NOT_IMPLEMENTED();
        } else {
            inp.reset();
            throw TBXX_NOT_IMPLEMENTED();
        }
        return *this;
    }
    read_loop &operator,(integer_star_4 &val) {
        if (io_mode == io_unformatted) {
            from_stream_unformatted(reinterpret_cast<char *>(&val), sizeof(integer_star_4));
        } else if (io_mode == io_list_directed) {
            val = static_cast<int>(read_star_long());
        } else {
            std::string const &ed = next_edit_descriptor();
            int n = ed.size();
            if (ed[0] == 'i' && n > 1) {
                n = utils::unsigned_integer_value(ed.data(), 1, n);
                val = static_cast<int>(read_fmt_long(n));
            } else {
                val = static_cast<int>(read_star_long());
            }
        }
        return *this;
    }
    // ---------------------------------------------------------------------
    // Explicit overloads for 'long' and 'unsigned long'.
    //
    // In C++, 'long' is always a distinct type from 'int' and 'long long',
    // even when sizeof(long)==sizeof(int) (ILP32, LLP64 / Windows).
    // Without these overloads the compiler cannot match 'long &' to the
    // integer_star_4 (int) or integer_star_8 (long long) overloads above,
    // causing compile errors or silent wrong conversions.
    //
    // The enable_if guard prevents a duplicate-overload error on LP64 Linux
    // where integer_star_8 is already typedef'd to 'long'.
    // ---------------------------------------------------------------------
    template <typename T, std::enable_if_t<std::is_same_v<T, long> && !std::is_same_v<long, integer_star_4> && !std::is_same_v<long, integer_star_8>, int> = 0> read_loop &operator,(T &val) {
        if constexpr (sizeof(long) == sizeof(integer_star_4)) {
            integer_star_4 tmp{};
            (*this), tmp;
            val = static_cast<long>(tmp);
        } else {
            integer_star_8 tmp{};
            (*this), tmp;
            val = static_cast<long>(tmp);
        }
        return *this;
    }
    template <typename T, std::enable_if_t<std::is_same_v<T, long long> && !std::is_same_v<long long, integer_star_4> && !std::is_same_v<long long, integer_star_8>, int> = 0> read_loop &operator,(T &val) {
        integer_star_8 tmp{};
        (*this), tmp;
        val = static_cast<long long>(tmp);
        return *this;
    }
    template <typename T, std::enable_if_t<std::is_same_v<T, unsigned long> && !std::is_same_v<unsigned long, std::uint32_t> && !std::is_same_v<unsigned long, std::uint64_t>, int> = 0> read_loop &operator,(T &val) {
        if constexpr (sizeof(unsigned long) == sizeof(integer_star_4)) {
            integer_star_4 tmp{};
            (*this), tmp;
            val = static_cast<unsigned long>(tmp);
        } else {
            integer_star_8 tmp{};
            (*this), tmp;
            val = static_cast<unsigned long>(tmp);
        }
        return *this;
    }
    read_loop &operator,(integer_star_8 &val) {
        if (io_mode == io_unformatted) {
            from_stream_unformatted(reinterpret_cast<char *>(&val), sizeof(integer_star_8));
        } else if (io_mode == io_list_directed) {
            // List-directed integer read (Fortran "*")
            val = static_cast<integer_star_8>(read_star_long());
        } else {
            // Formatted integer read (Iw) or fallback to list-directed rules
            std::string const &ed = next_edit_descriptor();
            int n = static_cast<int>(ed.size());
            if (ed[0] == 'i' && n > 1) {
                n = utils::unsigned_integer_value(ed.data(), 1, n);
                val = static_cast<integer_star_8>(read_fmt_long(n));
            } else {
                val = static_cast<integer_star_8>(read_star_long());
            }
        }
        return *this;
    }
    read_loop &operator,(float &val) {
        if (io_mode == io_unformatted) {
            from_stream_unformatted(reinterpret_cast<char *>(&val), sizeof(float));
        } else {
            val = static_cast<float>((io_mode == io_formatted ? read_fmt_double() : read_star_double()));
        }
        return *this;
    }
    read_loop &operator,(double &val) {
        if (io_mode == io_unformatted) {
            from_stream_unformatted(reinterpret_cast<char *>(&val), sizeof(double));
        } else {
            val = (io_mode == io_formatted ? read_fmt_double() : read_star_double());
        }
        return *this;
    }
    // --- Selected multiprecision/back-end numeric types: read as string token, normalize, then assign ---
    //
    // Rationale: LAPACK test inputs are within binary64, but MPLAPACK back-ends may use
    // mpf_class/mpfr_class/binary128/binary80/QD/DD. Reading as a raw token and then
    // constructing/assigning avoids heavy formatted parsing (g2.16 etc.) and handles
    // Fortran 'D' exponents by normalization to 'E'.
    //
#if defined(MPLAPACK_BUILD_WITH_GMP)
    read_loop &operator,(mpf_class &val) {
        std::string s = read_numeric_as_string();
        normalize_fortran_exponent(s);
        val = mpf_class(s.c_str());
        return *this;
    }
#endif
#if defined(MPLAPACK_BUILD_WITH_MPFR)
    read_loop &operator,(mpfrxx::mpfr_class &val) {
        std::string s = read_numeric_as_string();
        normalize_fortran_exponent(s);
        val = mpfrxx::mpfr_class(s.c_str());
        return *this;
    }
#endif
#if defined(MPLAPACK_BUILD_WITH_BINARY128)
    read_loop &operator,(mplapack_binary128_t &val) {
        std::string s = read_numeric_as_string();
        normalize_fortran_exponent(s);
        long double ld = std::strtold(s.c_str(), nullptr);
        val = static_cast<mplapack_binary128_t>(ld);
        return *this;
    }
#endif
#if defined(MPLAPACK_BUILD_WITH_BINARY80)
    read_loop &operator,(mplapack_binary80_t &val) {
        std::string s = read_numeric_as_string();
        normalize_fortran_exponent(s);
        long double ld = std::strtold(s.c_str(), nullptr);
        val = static_cast<mplapack_binary80_t>(ld);
        return *this;
    }
#endif
#if defined(MPLAPACK_BUILD_WITH_DD)
    // QD library dd_real (double-double)
    read_loop &operator,(dd_real &val) {
        std::string s = read_numeric_as_string();
        assign_from_token_string(val, s);
        return *this;
    }
#endif
#if defined(MPLAPACK_BUILD_WITH_QD)
    // QD library qd_real (quad-double)
    read_loop &operator,(qd_real &val) {
        std::string s = read_numeric_as_string();
        assign_from_token_string(val, s);
        return *this;
    }
#endif
// double is already handled by the existing arithmetic path.
#if defined(MPLAPACK_BUILD_WITH_DOUBLE)
// no-op
#endif

    // --- Helper function to read a complex literal "(real, imag)" as a string ---
    // Returns the full token including parentheses, with D->E normalization applied.
    std::string read_complex_token_string() {
        // Skip leading whitespace and separators (list-directed: newlines are separators)
        int c;
        while (true) {
            c = inp_get();
            if (utils::is_stream_end(c)) {
                if (this->iostat_ptr != 0)
                    *iostat_ptr = IOSTAT_END;
                throw read_end("End of input while reading complex value");
            }
            if (utils::is_end_of_line(c) || c == ',') {
                continue; // skip separators
            }
            if (!utils::is_whitespace(c)) {
                break; // found start of token
            }
        }

        // Expect '('
        if (c != '(') {
            throw io_err("Expected '(' at start of complex value, got: " + utils::format_char_for_display(c));
        }

        // Read until ')' to get the full complex literal "(real,imag)"
        std::string token;
        token.push_back('(');
        while (true) {
            c = inp_get();
            if (utils::is_stream_end(c)) {
                throw io_err("Unexpected end of input while reading complex value");
            }
            token.push_back(static_cast<char>(c));
            if (c == ')') {
                break;
            }
        }

        return normalize_fortran_numeric_string(token);
    }

    // --- Helper to parse "(real,imag)" string into two component strings ---
    static void parse_complex_components(const std::string &normalized, std::string &real_str, std::string &imag_str) {
        // normalized is "(real,imag)" with D->E conversion already done
        std::string inside = normalized.substr(1, normalized.size() - 2);
        size_t comma = inside.find(',');
        if (comma == std::string::npos) {
            throw io_err("Invalid complex format, missing comma: " + normalized);
        }
        real_str = inside.substr(0, comma);
        imag_str = inside.substr(comma + 1);

        // Trim whitespace
        auto trim = [](std::string &s) {
            size_t b = s.find_first_not_of(' ');
            if (b == std::string::npos) {
                s.clear();
                return;
            }
            size_t e = s.find_last_not_of(' ');
            s = s.substr(b, e - b + 1);
        };
        trim(real_str);
        trim(imag_str);
    }

    read_loop &operator,(std::complex<float> &val) {
        if (io_mode == io_unformatted) {
            float re, im;
            from_stream_unformatted(reinterpret_cast<char *>(&re), sizeof(float));
            from_stream_unformatted(reinterpret_cast<char *>(&im), sizeof(float));
            val = std::complex<float>(re, im);
        } else {
            std::string token = read_complex_token_string();
            std::string real_str, imag_str;
            parse_complex_components(token, real_str, imag_str);
            float re = static_cast<float>(std::strtod(real_str.c_str(), nullptr));
            float im = static_cast<float>(std::strtod(imag_str.c_str(), nullptr));
            val = std::complex<float>(re, im);
        }
        return *this;
    }
    read_loop &operator,(std::complex<double> &val) {
        if (io_mode == io_unformatted) {
            double re, im;
            from_stream_unformatted(reinterpret_cast<char *>(&re), sizeof(double));
            from_stream_unformatted(reinterpret_cast<char *>(&im), sizeof(double));
            val = std::complex<double>(re, im);
        } else {
            std::string token = read_complex_token_string();
            std::string real_str, imag_str;
            parse_complex_components(token, real_str, imag_str);
            double re = std::strtod(real_str.c_str(), nullptr);
            double im = std::strtod(imag_str.c_str(), nullptr);
            val = std::complex<double>(re, im);
        }
        return *this;
    }
#if defined(MPLAPACK_BUILD_WITH_GMP)
    read_loop &operator,(mpfc_class &val) {
        if (io_mode == io_unformatted) {
            throw TBXX_NOT_IMPLEMENTED();
        }
        std::string token = read_complex_token_string();
        std::string real_str, imag_str;
        parse_complex_components(token, real_str, imag_str);
        mpf_class re(real_str.c_str());
        mpf_class im(imag_str.c_str());
        val = mpfc_class(re, im);
        return *this;
    }
#endif
#if defined(MPLAPACK_BUILD_WITH_MPFR)
    read_loop &operator,(mpc_class &val) {
        if (io_mode == io_unformatted) {
            throw TBXX_NOT_IMPLEMENTED();
        }
        std::string token = read_complex_token_string();
        std::string real_str, imag_str;
        parse_complex_components(token, real_str, imag_str);
        mpfrxx::mpfr_class re(real_str.c_str());
        mpfrxx::mpfr_class im(imag_str.c_str());
        val = mpc_class(re, im);
        return *this;
    }
#endif
#if defined(MPLAPACK_BUILD_WITH_DD)
    read_loop &operator,(dd_complex &val) {
        if (io_mode == io_unformatted) {
            throw TBXX_NOT_IMPLEMENTED();
        }
        std::string token = read_complex_token_string();
        std::string real_str, imag_str;
        parse_complex_components(token, real_str, imag_str);
        dd_real re(real_str.c_str());
        dd_real im(imag_str.c_str());
        val = dd_complex(re, im);
        return *this;
    }
#endif
#if defined(MPLAPACK_BUILD_WITH_QD)
    read_loop &operator,(qd_complex &val) {
        if (io_mode == io_unformatted) {
            throw TBXX_NOT_IMPLEMENTED();
        }
        std::string token = read_complex_token_string();
        std::string real_str, imag_str;
        parse_complex_components(token, real_str, imag_str);
        qd_real re(real_str.c_str());
        qd_real im(imag_str.c_str());
        val = qd_complex(re, im);
        return *this;
    }
#endif
#if defined(MPLAPACK_BUILD_WITH_BINARY128)
    read_loop &operator,(std::complex<mplapack_binary128_t> &val) {
        if (io_mode == io_unformatted) {
            throw TBXX_NOT_IMPLEMENTED();
        }
        std::string token = read_complex_token_string();
        std::string real_str, imag_str;
        parse_complex_components(token, real_str, imag_str);
        long double ld_re = std::strtold(real_str.c_str(), nullptr);
        long double ld_im = std::strtold(imag_str.c_str(), nullptr);
        val = std::complex<mplapack_binary128_t>(static_cast<mplapack_binary128_t>(ld_re), static_cast<mplapack_binary128_t>(ld_im));
        return *this;
    }
#endif
#if defined(MPLAPACK_BUILD_WITH_BINARY80)
    read_loop &operator,(std::complex<mplapack_binary80_t> &val) {
        if (io_mode == io_unformatted) {
            throw TBXX_NOT_IMPLEMENTED();
        }
        std::string token = read_complex_token_string();
        std::string real_str, imag_str;
        parse_complex_components(token, real_str, imag_str);
        long double ld_re = std::strtold(real_str.c_str(), nullptr);
        long double ld_im = std::strtold(imag_str.c_str(), nullptr);
        val = std::complex<mplapack_binary80_t>(static_cast<mplapack_binary80_t>(ld_re), static_cast<mplapack_binary80_t>(ld_im));
        return *this;
    }
#endif
    read_loop &operator,(str_ref const &val) {
        if (io_mode == io_unformatted) {
            from_stream_unformatted(val.elems(), val.len());
        } else {
            int vl = val.len();
            int n = vl;
            if (io_mode == io_formatted) {
                std::string const &ed = next_edit_descriptor();
                if (ed[0] == 'a' && ed.size() > 1) {
                    n = utils::unsigned_integer_value(ed.data(), 1, ed.size());
                }
            }
            for (int i = 0; i < n - vl; i++) {
                int c = inp_get();
                if (utils::is_stream_end(c)) {
                    if (i == 0) {
                        inp.reset();
                        if (this->iostat_ptr != 0)
                            *iostat_ptr = IOSTAT_END;
                        throw read_end("End of input while reading string");
                    }
                    break;
                }
                if (utils::is_end_of_line(c)) {
                    inp->backup();
                }
            }
            for (int i = 0; i < vl; i++) {
                int c = inp_get();
                if (utils::is_stream_end(c)) {
                    if (i == 0) {
                        inp.reset();
                        if (this->iostat_ptr != 0)
                            *iostat_ptr = IOSTAT_END;
                        throw read_end("End of input while reading string");
                    }
                    val[i] = ' ';
                } else if (utils::is_end_of_line(c)) {
                    inp->backup();
                    val[i] = ' ';
                } else {
                    val[i] = c;
                }
            }
        }
        return *this;
    }
    template <typename T, size_t Ndims> read_loop &operator,(arr_ref<T, Ndims> const &val) {
        T *v = val.begin();
        size_t n = val.size_1d();
        for (size_t i = 0; i < n; i++) {
            (*this), v[i];
        }
        return *this;
    }
    template <size_t Ndims> read_loop &operator,(str_arr_ref<Ndims> const &val) {
        size_t n = val.size_1d();
        int l = val.len();
        char *val_begin = val.begin();
        for (size_t i = 0; i < n; i++) {
            (*this), str_ref(&val_begin[i * l], l);
        }
        return *this;
    }
    void skip_to_end_of_line() {
        while (true) {
            int c = inp_get();
            if (utils::is_stream_end(c) || utils::is_end_of_line(c)) {
                break;
            }
        }
    }
    ~read_loop() NOEXCEPT_FALSE {
        if (inp.get() == 0)
            return;
        if (io_mode == io_unformatted) {
            skip_to_end_of_unformatted_record();
        } else {
            skip_to_end_of_line();
        }
    }
    long read_fmt_long(unsigned n) {
        bool had_non_blank = false;
        bool negative = false;
        long result = 0;
        for (unsigned i = 0; i < n; i++) {
            int c = inp_get();
            if (utils::is_stream_end(c)) {
                break;
            }
            if (c == ',') {
                break;
            }
            if (utils::is_end_of_line(c)) {
                inp->backup();
                break;
            }
            if (c != ' ') {
                if (!had_non_blank) {
                    had_non_blank = true;
                    if (c == '-') {
                        negative = true;
                        continue;
                    }
                    if (c == '+') {
                        continue;
                    }
                }
                if (!utils::is_digit(c)) {
                    inp.reset();
                    if (this->iostat_ptr != 0)
                        *iostat_ptr = IOSTAT_ERROR;
                    throw io_err("Invalid character while reading integer value.");
                }
                result *= 10;
                result += utils::digit_as_int(c);
            }
        }
        if (negative)
            result *= -1;
        return result;
    }
    long read_star_long() {
        while (true) { // loop scanning for first token character
            int c = inp_get();
            if (utils::is_stream_end(c)) {
                inp.reset();
                if (this->iostat_ptr != 0)
                    *iostat_ptr = IOSTAT_END;
                throw read_end("End of input while reading integer value");
            }
            // List-directed ("*") input: end-of-record is a separator, not data.
            if (utils::is_end_of_line(c) || c == ',') {
                continue;
            }
            if (!utils::is_whitespace(c)) {
                bool negative = (c == '-');
                if (negative || c == '+') {
                    c = inp_get();
                    if (utils::is_stream_end(c)) {
                        throw read_end("End of input while reading integer value");
                    }
                }
                long result = 0;
                while (true) { // loop collecting digits
                    if (!utils::is_digit(c)) {
                        throw io_err("Invalid character while reading integer value.");
                    }
                    result *= 10;
                    result += utils::digit_as_int(c);
                    c = inp_get();
                    if (utils::is_stream_end(c) || utils::is_whitespace(c) || c == ',' || utils::is_end_of_line(c)) {
                        if (negative)
                            result *= -1;
                        if (utils::is_end_of_line(c))
                            inp->backup();
                        return result;
                    }
                }
            }
        }
    }
    void throw_if_conv_error_message(utils::string_to_double const &conv) {
        if (conv.error_message) {
            inp.reset();
            if (conv.stream_end) {
                if (this->iostat_ptr != 0)
                    *iostat_ptr = IOSTAT_END;
                throw read_end(*conv.error_message);
            }
            if (this->iostat_ptr != 0)
                *iostat_ptr = IOSTAT_ERROR;
            throw io_err(*conv.error_message);
        }
    }
    double read_fmt_double() {
        std::string const &ed = next_edit_descriptor();
        int n = ed.size();
        if (n < 2 || std::strchr("defg", ed[0]) == 0) {
            return read_star_double();
        }
        int iw = utils::unsigned_integer_scan(ed.data(), 1, ed.size());
        int w = utils::unsigned_integer_value(ed.data(), 1, iw);
        int d = 0;
        if (iw + 1 != ed.size()) {
            d = utils::unsigned_integer_value(ed.data(), iw + 1, ed.size());
        }
        utils::string_to_double_fmt conv(*inp, w, d, blanks_zero, exp_scale);
        throw_if_conv_error_message(conv);
        first_inp_get = false;
        return conv.result;
    }
    double read_star_double() {
        utils::string_to_double conv(*inp);
        throw_if_conv_error_message(conv);
        first_inp_get = false;
        int c = inp_get();
        if (utils::is_stream_end(c) || utils::is_whitespace(c) || c == ',') {
            if (utils::is_end_of_line(c))
                inp->backup();
            return conv.result;
        }
        inp.reset();
        if (this->iostat_ptr != 0)
            *iostat_ptr = IOSTAT_ERROR;
        throw io_err("Invalid character while reading floating-point value: " + utils::format_char_for_display(c));
    }
    // --- Helpers for reading numeric tokens as raw strings ---
    // These are used for selected multiprecision types (mpf_class, mpfrxx::mpfr_class, mplapack_binary128_t, ...).
    // For LAPACK test inputs, values are at most binary64, so string->ctor assignment is sufficient.
    static std::string normalize_fortran_numeric_string(std::string s) {
        // Convert Fortran 'D' exponent to 'E' (e.g., 1.0D+00 -> 1.0E+00).
        // Also supports complex literals like "(-1.0D+00, 9.0D+00)" by
        // normalizing each component.
        auto trim_inplace = [](std::string &x) {
            size_t b = x.find_first_not_of(' ');
            if (b == std::string::npos) {
                x.clear();
                return;
            }
            size_t e = x.find_last_not_of(' ');
            x = x.substr(b, e - b + 1);
        };
        trim_inplace(s);
        if (s.size() >= 2 && s.front() == '(' && s.back() == ')') {
            std::string inside = s.substr(1, s.size() - 2);
            size_t comma = inside.find(',');
            if (comma != std::string::npos) {
                std::string a = inside.substr(0, comma);
                std::string b = inside.substr(comma + 1);
                trim_inplace(a);
                trim_inplace(b);
                for (char &c: a) {
                    if (c == 'D' || c == 'd')
                        c = 'E';
                }
                for (char &c: b) {
                    if (c == 'D' || c == 'd')
                        c = 'E';
                }
                return "(" + a + "," + b + ")";
            }
        }
        for (char &c: s) {
            if (c == 'D' || c == 'd')
                c = 'E';
        }
        return s;
    }
    std::string read_star_token_string() {
        // Skip leading whitespace and commas
        int c = 0;
        while (true) {
            c = inp_get();
            if (utils::is_stream_end(c)) {
                inp.reset();
                if (this->iostat_ptr != 0) {
                    *iostat_ptr = IOSTAT_END;
                }
                throw read_end("End of input while reading token");
            }
            if (utils::is_end_of_line(c)) {
                // List-directed ("*") input: end-of-record is a separator, not an error.
                // If we have not started a token yet, skip the record boundary and keep scanning.
                continue;
            }
            if (!utils::is_whitespace(c) && c != ',') {
                break;
            }
        }
        std::string s;
        s.push_back(static_cast<char>(c));
        while (true) {
            c = inp_get();
            if (utils::is_stream_end(c) || utils::is_whitespace(c) || c == ',') {
                if (utils::is_end_of_line(c)) {
                    inp->backup();
                }
                return s;
            }
            if (utils::is_end_of_line(c)) {
                inp->backup();
                return s;
            }
            s.push_back(static_cast<char>(c));
        }
    }
    std::string read_fmt_token_string(int w) {
        std::string s;
        s.reserve(static_cast<size_t>(w));
        for (int i = 0; i < w; i++) {
            int c = inp_get();
            if (utils::is_stream_end(c)) {
                break;
            }
            if (utils::is_end_of_line(c)) {
                inp->backup();
                break;
            }
            s.push_back(static_cast<char>(c));
        }
        // Trim spaces on both ends
        size_t b = s.find_first_not_of(' ');
        if (b == std::string::npos) {
            return std::string();
        }
        size_t e = s.find_last_not_of(' ');
        return s.substr(b, e - b + 1);
    }
    std::string read_numeric_as_string() {
        if (io_mode == io_unformatted) {
            inp.reset();
            throw TBXX_NOT_IMPLEMENTED();
        }
        if (io_mode == io_formatted) {
            std::string const &ed = next_edit_descriptor();
            int n = static_cast<int>(ed.size());
            if (n >= 2 && std::strchr("defg", ed[0]) != 0) {
                int iw = utils::unsigned_integer_scan(ed.data(), 1, ed.size());
                int w = utils::unsigned_integer_value(ed.data(), 1, iw);
                return normalize_fortran_numeric_string(read_fmt_token_string(w));
            }
            // Fallback to list-directed tokenization.
        }
        return normalize_fortran_numeric_string(read_star_token_string());
    }
    void from_stream_unformatted(char *target, unsigned target_size) {
        for (unsigned i = 0; i < target_size; i++) {
            int ic = inp_get();
            char c = static_cast<char>(ic);
            if (c == end_of_unformatted_record) {
                if (inp_get() != ic) {
                    inp.reset();
                    if (this->iostat_ptr != 0)
                        *iostat_ptr = IOSTAT_END;
                    throw read_end("End of record during unformatted read");
                }
            }
            target[i] = c;
        }
    }
    void skip_to_end_of_unformatted_record() {
        while (true) {
            char c = static_cast<char>(inp_get());
            if (c == end_of_unformatted_record) {
                if (inp_get() == 0) {
                    break;
                }
            }
        }
    }
};
struct read_from_string : public read_loop {
    read_from_string(std::string const &internal_file, str_cref fmt) : read_loop(str_cref(internal_file.data(), static_cast<int>(internal_file.size())), fmt) {}
};
struct common_read {
    common &cmn;
    common_read(common &cmn_) : cmn(cmn_) {}
    read_loop operator()(int unit, unformatted_type const &) {
        read_loop result(cmn, unit, unformatted);
        return result;
    }
    read_loop operator()(int unit, star_type const &) {
        read_loop result(cmn, unit, star);
        return result;
    }
    read_loop operator()(int const &unit, str_cref fmt) {
        read_loop result(cmn, unit, fmt);
        return result;
    }
    read_loop operator()(str_cref const &internal_file, star_type const &) {
        read_loop result(internal_file, star);
        return result;
    }
    read_loop operator()(str_cref const &internal_file, str_cref fmt) {
        read_loop result(internal_file, fmt);
        return result;
    }
};
} // namespace fem
#endif // GUARD

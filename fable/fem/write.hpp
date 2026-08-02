#ifndef FEM_WRITE_HPP
#define FEM_WRITE_HPP

#include <limits>
#include <cmath>
#include <cerrno>

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

#if defined(MPLAPACK_BUILD_WITH_QD) || defined(MPLAPACK_BUILD_WITH_DD)
// QD headers define and use qd::nint (and other short identifiers).
// Temporarily disable macros that would interfere.
#if defined(nint)
#pragma push_macro("nint")
#undef nint
#define FEM_WRITE_RESTORE_nint 1
#endif
#if defined(min)
#pragma push_macro("min")
#undef min
#define FEM_WRITE_RESTORE_min 1
#endif
#if defined(max)
#pragma push_macro("max")
#undef max
#define FEM_WRITE_RESTORE_max 1
#endif
#if defined(abs)
#pragma push_macro("abs")
#undef abs
#define FEM_WRITE_RESTORE_abs 1
#endif
#if defined(sign)
#pragma push_macro("sign")
#undef sign
#define FEM_WRITE_RESTORE_sign 1
#endif
#if __has_include(<qd/dd_real.h>)
#include <qd/dd_real.h>
#endif
#if __has_include(<qd/qd_real.h>)
#include <qd/qd_real.h>
#endif
#if __has_include(<qd/dd_complex.h>)
#include <qd/dd_complex.h>
#endif
#if __has_include(<qd/qd_complex.h>)
#include <qd/qd_complex.h>
#endif
#if defined(FEM_WRITE_RESTORE_sign)
#pragma pop_macro("sign")
#undef FEM_WRITE_RESTORE_sign
#endif
#if defined(FEM_WRITE_RESTORE_abs)
#pragma pop_macro("abs")
#undef FEM_WRITE_RESTORE_abs
#endif
#if defined(FEM_WRITE_RESTORE_max)
#pragma pop_macro("max")
#undef FEM_WRITE_RESTORE_max
#endif
#if defined(FEM_WRITE_RESTORE_min)
#pragma pop_macro("min")
#undef FEM_WRITE_RESTORE_min
#endif
#if defined(FEM_WRITE_RESTORE_nint)
#pragma pop_macro("nint")
#undef FEM_WRITE_RESTORE_nint
#endif
#endif

// Fallback buffer length (should be provided by mplapack_utils_*.h).
#ifndef MPLAPACK_BUFLEN
#define MPLAPACK_BUFLEN 1024
#endif

#include <noexcept_false.hpp>
#include <fem/common.hpp>
#include <fem/format.hpp>
#include <fem/star.hpp>
#include <fem/str_arr_ref.hpp>
#include <fem/utils/double_to_string.hpp>
#include <fem/utils/misc.hpp>
#include <type_traits>
#include <utility>
#include <fem/utils/real_as_string.hpp>
#if defined(_MSC_VER)
#define FEM_WRITE_CRLF true
#else
#define FEM_WRITE_CRLF false
#endif
namespace fem {
// Helper for printing user-defined numeric types (e.g. multiprecision reals).
// If sprintnum_short(char*, T const&) is available, we convert the value to a
// short string and emit it as text, ignoring the numeric edit descriptor.
namespace detail {
    template <typename T> struct make_void {
        typedef void type;
    };
    template <typename T, typename = void> struct has_sprintnum_short : std::false_type {};
    template <typename T> struct has_sprintnum_short<T, typename make_void<decltype(sprintnum_short(std::declval<char *>(), std::declval<T const &>()))>::type> : std::true_type {};
} // namespace detail
struct write_loop_base {
    bool write_crlf;
    unsigned pos;
    bool prev_was_string;
    int exp_scale;
    unsigned number_of_x_held;
    bool suppress_new_line_at_end;
    write_loop_base(bool write_crlf_ = false) : write_crlf(write_crlf_), pos(0), prev_was_string(false), exp_scale(0), number_of_x_held(0), suppress_new_line_at_end(false) {}
};
class write_loop : write_loop_base
// TODO copy-constructor potential performance problem
{
  private:
    utils::slick_ptr<utils::simple_ostream> out;
    int internal_file_len;
    io_modes io_mode;
    format::token_loop fmt_loop;
    bool terminated_by_colon;

  public:
    write_loop(common &cmn, int const &unit, unformatted_type const &) : write_loop_base(FEM_WRITE_CRLF && !is_std_io_unit(unit)), out(cmn.io.simple_ostream(unit)), internal_file_len(-1), io_mode(io_unformatted), terminated_by_colon(false) {}
    write_loop(common &cmn, int const &unit, star_type const &) : write_loop_base(FEM_WRITE_CRLF && !is_std_io_unit(unit)), out(cmn.io.simple_ostream(unit)), internal_file_len(-1), io_mode(io_list_directed), terminated_by_colon(false) {}
    write_loop(common &cmn, int const &unit, str_cref fmt) : write_loop_base(FEM_WRITE_CRLF && !is_std_io_unit(unit)), out(cmn.io.simple_ostream(unit)), internal_file_len(-1), io_mode(io_formatted), fmt_loop(fmt), terminated_by_colon(false) {}
    write_loop(str_ref const &internal_file, star_type const &) : out(utils::slick_ptr<utils::simple_ostream>(new utils::simple_ostream_to_char_ptr_and_size(internal_file.elems(), internal_file.len()))), internal_file_len(internal_file.len()), io_mode(io_list_directed), terminated_by_colon(false) {}
    write_loop(str_ref const &internal_file, str_cref fmt) : out(utils::slick_ptr<utils::simple_ostream>(new utils::simple_ostream_to_char_ptr_and_size(internal_file.elems(), internal_file.len()))), internal_file_len(internal_file.len()), io_mode(io_formatted), fmt_loop(fmt), terminated_by_colon(false) {}
    template <size_t Ndims> write_loop(str_arr_ref<Ndims> const &internal_file, star_type const &) : out(utils::slick_ptr<utils::simple_ostream>(new utils::simple_ostream_to_char_ptr_and_size(internal_file.begin(), internal_file.len()))), internal_file_len(internal_file.len()), io_mode(io_list_directed), terminated_by_colon(false) {}
    template <size_t Ndims> write_loop(str_arr_ref<Ndims> const &internal_file, str_cref fmt) : out(utils::slick_ptr<utils::simple_ostream>(new utils::simple_ostream_to_char_ptr_and_size(internal_file.begin(), internal_file.len()))), internal_file_len(internal_file.len()), io_mode(io_formatted), fmt_loop(fmt), terminated_by_colon(false) {}
    std::string const &next_edit_descriptor(bool final = false) {
        while (true) {
            if (terminated_by_colon) {
                static const std::string empty("");
                return empty;
            }
            utils::token const *t = fmt_loop.next_executable_token(final);
            if (t == 0) {
                static const std::string empty("");
                return empty;
            }
            std::string const &tv = t->value;
            if (t->type == "string") {
                to_stream_fmt(tv.data(), tv.size());
            } else if (t->type == "op") {
                if (tv[0] == ':') {
                    if (final)
                        terminated_by_colon = true;
                } else if (tv[0] == '/') {
                    if (write_crlf) {
                        to_stream_fmt("\r\n", 2);
                    } else {
                        to_stream_fmt("\n", 1);
                    }
                } else if (tv[0] == '$') {
                    suppress_new_line_at_end = true;
                } else {
                    out.reset();
                    throw TBXX_UNREACHABLE_ERROR();
                }
            } else if (t->type == "format") {
                if (utils::ends_with_char(tv, 'x')) {
                    unsigned n = tv.size();
                    if (n != 1)
                        n = utils::unsigned_integer_value(tv.data(), n - 1);
                    number_of_x_held += n;
                } else if (std::strchr("adefgilz", tv[0]) != 0) {
                    return tv;
                } else if (utils::ends_with_char(tv, 'p')) {
                    if (tv.size() == 1) {
                        exp_scale = 1;
                    } else {
                        exp_scale = utils::signed_integer_value(tv.data(), 0, tv.size() - 1);
                    }
                } else if (tv[0] == 't') {
                    out.reset();
                    throw TBXX_NOT_IMPLEMENTED();
                } else if (tv[0] == 's') {
                    out.reset();
                    throw TBXX_NOT_IMPLEMENTED();
                } else if (tv[0] == 'b') {
                    out.reset();
                    throw TBXX_NOT_IMPLEMENTED();
                } else {
                    out.reset();
                    throw TBXX_UNREACHABLE_ERROR();
                }
            } else {
                out.reset();
                throw TBXX_UNREACHABLE_ERROR();
            }
        }
    }
    // C-style array: expand all elements into the write loop (Fortran array I/O semantics)
    template <typename T, std::size_t N, typename = typename std::enable_if<!std::is_same<typename std::remove_cv<T>::type, char>::value>::type> write_loop &operator,(T const (&arr)[N]) {
        for (std::size_t i = 0; i < N; i++) {
            (*this), arr[i];
        }
        return *this;
    }
    write_loop &operator,(char const &val) {
        if (io_mode == io_unformatted) {
            to_stream_unformatted(&val, 1);
        } else if (io_mode == io_list_directed) {
            to_stream(&val, 1, /*space*/ !prev_was_string);
            prev_was_string = true;
        } else {
            out.reset();
            throw TBXX_NOT_IMPLEMENTED();
        }
        return *this;
    }
    write_loop &operator,(char const *val) {
        if (io_mode == io_unformatted) {
            to_stream_unformatted(val, std::strlen(val));
        } else if (io_mode == io_list_directed) {
            to_stream_star(val, std::strlen(val), /*space*/ !prev_was_string);
            prev_was_string = true;
        } else {
            std::string const &ed = next_edit_descriptor();
            int l = std::strlen(val);
            if (ed[0] == 'a') {
                int n = ed.size();
                if (n > 1) {
                    n = utils::unsigned_integer_value(ed.data() + 1, n - 1);
                    int b = n - l;
                    if (b < 0) {
                        l = n;
                    } else {
                        for (int i = 0; i < b; i++)
                            to_stream_fmt(" ", 1);
                    }
                }
                to_stream_fmt(val, l);
            } else {
                to_stream_fmt(val, l);
            }
        }
        return *this;
    }
    write_loop &operator,(bool const &val) {
        if (io_mode == io_unformatted) {
            to_stream_unformatted(reinterpret_cast<char const *>(&val), sizeof(bool));
        } else if (io_mode == io_list_directed) {
            to_stream_star((val ? "T" : "F"), 1);
            prev_was_string = false;
        } else {
            std::string const &ed = next_edit_descriptor();
            int n = ed.size();
            if (ed[0] == 'l' && n > 1) {
                n = utils::unsigned_integer_value(ed.data() + 1, n - 1);
            } else {
                n = 1;
            }
            for (int i = 1; i < n; i++)
                to_stream_fmt(" ", 1);
            to_stream_fmt((val ? "T" : "F"), 1);
        }
        return *this;
    }
    write_loop &operator,(integer_star_1 const &val) {
        if (io_mode == io_unformatted) {
            to_stream_unformatted(reinterpret_cast<char const *>(&val), sizeof(integer_star_1));
        } else if (io_mode == io_list_directed) {
            char buf[64];
            int n = std::snprintf(buf, sizeof(buf), "%4d", static_cast<int>(val));
            to_stream_star(buf, n);
            prev_was_string = false;
        } else {
            out.reset();
            throw TBXX_NOT_IMPLEMENTED();
        }
        return *this;
    }
    write_loop &operator,(integer_star_2 const &val) {
        if (io_mode == io_unformatted) {
            to_stream_unformatted(reinterpret_cast<char const *>(&val), sizeof(integer_star_2));
        } else if (io_mode == io_list_directed) {
            char buf[64];
            int n = std::snprintf(buf, sizeof(buf), "%6d", static_cast<int>(val));
            to_stream_star(buf, n);
            prev_was_string = false;
        } else {
            out.reset();
            throw TBXX_NOT_IMPLEMENTED();
        }
        return *this;
    }
    write_loop &operator,(integer_star_4 const &val) {
        if (io_mode == io_unformatted) {
            to_stream_unformatted(reinterpret_cast<char const *>(&val), sizeof(integer_star_4));
        } else if (io_mode == io_list_directed) {
            char buf[64];
            int n = std::snprintf(buf, sizeof(buf), "%11d", val);
            to_stream_star(buf, n);
            prev_was_string = false;
        } else {
            std::string const &ed = next_edit_descriptor();
            if (ed[0] == 'i') {
                int n = ed.size();
                TBXX_ASSERT(n + 2 < 64);
                char fmt[64];
                fmt[0] = '%';
                std::strncpy(fmt + 1, ed.data() + 1, n - 1);
                fmt[n] = 'd';
                fmt[n + 1] = '\0';
                char buf[64];
                n = std::snprintf(buf, sizeof(buf), fmt, val);
                to_stream_fmt(buf, n);
            } else {
                char buf[64];
                int n = std::snprintf(buf, sizeof(buf), " %d", val);
                to_stream_fmt(buf, n);
            }
        }
        return *this;
    }
    write_loop &operator,(integer_star_8 const &val) {
        if (io_mode == io_unformatted) {
            to_stream_unformatted(reinterpret_cast<char const *>(&val), sizeof(integer_star_8));
        } else if (io_mode == io_list_directed) {
            // TODO faster implementation
            std::ostringstream o;
            o.width(21);
            o << val;
            std::string s = o.str();
            to_stream_star(s.data(), s.size());
            prev_was_string = false;
        } else {
            std::string const &ed = next_edit_descriptor();
            if (ed[0] == 'i') {
                // Use long long formatting for INTEGER*8.
                std::string fmt = "%" + ed.substr(1) + "lld";
                char buf[64];
                int n = std::snprintf(buf, sizeof(buf), fmt.c_str(), static_cast<long long>(val));
                to_stream_fmt(buf, n);
            } else {
                char buf[64];
                int n = std::snprintf(buf, sizeof(buf), " %lld", static_cast<long long>(val));
                to_stream_fmt(buf, n);
            }
        }
        return *this;
    }
    // Explicit overload for plain 'long int' to resolve ambiguity on platforms
    // where long int is distinct from both int32_t and int64_t.
    //
    // ILP32  (32-bit Linux/etc.) : int=32b, long=32b, ll=64b
    //                              int32_t=int, int64_t=long long  => long is distinct
    // LLP64  (Windows 64-bit)   : int=32b, long=32b, ll=64b
    //                              int32_t=int, int64_t=long long  => long is distinct
    // LP64   (64-bit Linux/etc.) : int=32b, long=64b, ll=64b
    // int32_t=int, int64_t=long       => long == integer_star_8 (already covered)
    //
    // SFINAE: only instantiate when 'long' is not already covered by the
    // integer_star_4 or integer_star_8 overloads above.
    template <typename U = long int,
              typename std::enable_if<
                  !std::is_same<U, integer_star_4>::value &&
                  !std::is_same<U, integer_star_8>::value,
              int>::type = 0>
    write_loop &operator,(long int const &val) {
        // Delegate to the appropriately-sized fixed-width overload.
        if (sizeof(long int) == sizeof(integer_star_4)) {
            return (*this), static_cast<integer_star_4>(val);
        } else {
            return (*this), static_cast<integer_star_8>(val);
        }
    }
    // Ditto for 'unsigned long int'.
    template <typename U = unsigned long int,
              typename std::enable_if<
                  !std::is_same<U, uint32_t>::value &&
                  !std::is_same<U, uint64_t>::value,
              int>::type = 0>
    write_loop &operator,(unsigned long int const &val) {
        if (sizeof(unsigned long int) == sizeof(integer_star_4)) {
            return (*this), static_cast<integer_star_4>(static_cast<long int>(val));
        } else {
            return (*this), static_cast<integer_star_8>(static_cast<long long>(val));
        }
    }
    // Helper function to format a string according to an edit descriptor.
    // Must be defined before dd_real/qd_real overloads that use it.
    void to_stream_fmt_double_given_string(std::string const &val_str, std::string const &ed) {
        // Parse edit descriptor (e.g., "d20.10", "e25.15", "f15.8")
        int n = static_cast<int>(ed.size());
        if (n < 2) {
            // Fallback: just output the string
            for (char c: val_str) {
                out->put(c);
            }
            pos += val_str.size();
            return;
        }
        int iw = utils::unsigned_integer_scan(ed.data(), 1, ed.size());
        int w = utils::unsigned_integer_value(ed.data(), 1, iw);

        // Right-justify in field of width w
        int padding = w - static_cast<int>(val_str.size());
        for (int i = 0; i < padding; i++) {
            out->put(' ');
            pos++;
        }
        for (char c: val_str) {
            out->put(c);
            pos++;
        }
    }

#if defined(MPLAPACK_BUILD_WITH_DD) || defined(MPLAPACK_BUILD_WITH_QD)
    //
    // Explicit overload for dd_real to prevent infinite recursion.
    // Without this, dd_real goes through generic template -> char buffer ->
    // dd_real(const char*) -> infinite loop -> stack overflow.
    //
    write_loop &operator,(dd_real const &val) {
        if (io_mode == io_list_directed) {
            // List-directed: convert to string and write
            char buf[64];
            val.write(buf, sizeof(buf), 32);
            // Trim and write
            std::string s(buf);
            size_t b = s.find_first_not_of(' ');
            size_t e = s.find_last_not_of(' ');
            if (b != std::string::npos && e != std::string::npos) {
                s = s.substr(b, e - b + 1);
            }
            // Add space separator for list-directed
            if (pos != 0) {
                out->put(' ');
                pos++;
            }
            for (char c: s) {
                out->put(c);
                pos++;
            }
        } else if (io_mode == io_formatted) {
            // Formatted output using edit descriptor
            char buf[64];
            val.write(buf, sizeof(buf), 32);
            std::string s(buf);
            size_t b = s.find_first_not_of(' ');
            size_t e = s.find_last_not_of(' ');
            if (b != std::string::npos && e != std::string::npos) {
                s = s.substr(b, e - b + 1);
            }
            // Get edit descriptor and format accordingly
            std::string const &ed = next_edit_descriptor();
            to_stream_fmt_double_given_string(s, ed);
        } else {
            // Unformatted: write raw bytes
            to_stream_unformatted(reinterpret_cast<char const *>(&val), sizeof(val));
        }
        return *this;
    }

    //
    // Explicit overload for dd_complex
    //
    write_loop &operator,(dd_complex const &val) { return (*this), val.real(), val.imag(); }
#endif

#if defined(MPLAPACK_BUILD_WITH_QD)
    //
    // Explicit overload for qd_real
    //
    write_loop &operator,(qd_real const &val) {
        if (io_mode == io_list_directed) {
            char buf[128];
            val.write(buf, sizeof(buf), 64);
            std::string s(buf);
            size_t b = s.find_first_not_of(' ');
            size_t e = s.find_last_not_of(' ');
            if (b != std::string::npos && e != std::string::npos) {
                s = s.substr(b, e - b + 1);
            }
            if (pos != 0) {
                out->put(' ');
                pos++;
            }
            for (char c: s) {
                out->put(c);
                pos++;
            }
        } else if (io_mode == io_formatted) {
            char buf[128];
            val.write(buf, sizeof(buf), 64);
            std::string s(buf);
            size_t b = s.find_first_not_of(' ');
            size_t e = s.find_last_not_of(' ');
            if (b != std::string::npos && e != std::string::npos) {
                s = s.substr(b, e - b + 1);
            }
            std::string const &ed = next_edit_descriptor();
            to_stream_fmt_double_given_string(s, ed);
        } else {
            to_stream_unformatted(reinterpret_cast<char const *>(&val), sizeof(val));
        }
        return *this;
    }

    //
    // Explicit overload for qd_complex
    //
    write_loop &operator,(qd_complex const &val) { return (*this), val.real(), val.imag(); }
#endif

  protected: // implementation detail
    void to_stream_fmt_double(double const &val) {
        std::string const &ed = next_edit_descriptor();
        if (ed[0] == 'f') {
            int n = ed.size();
            TBXX_ASSERT(n + 2 < 64);
            char fmt[64];
            fmt[0] = '%';
            std::strncpy(fmt + 1, ed.data() + 1, n - 1);
            fmt[n] = 'f';
            fmt[n + 1] = '\0';
            char buf[64];
            n = std::snprintf(buf, sizeof(buf), fmt, val);
            to_stream_fmt(buf, n);
        } else if ((ed[0] == 'd' || ed[0] == 'e') && ed.size() > 1) {
            int es = ed.size();
            int nw = utils::unsigned_integer_scan(ed.data(), 1, es);
            TBXX_ASSERT(nw > 0);
            int w = utils::unsigned_integer_value(ed.data(), 1, nw);
            int d = 0;
            if (nw != es) {
                TBXX_ASSERT(ed[nw] == '.');
                TBXX_ASSERT(nw + 1 != es);
                d = utils::unsigned_integer_value(ed.data(), nw + 1, es);
            }
            utils::double_to_string_scientific_notation conv(val, w, d, exp_scale, utils::to_upper(ed[0]));
            to_stream_fmt(conv.buffer, w);
        } else {
            char buf[64];
            int n = std::snprintf(buf, sizeof(buf), " %.6g", val);
            to_stream_fmt(buf, n);
        }
    }

  public:
    write_loop &operator,(float const &val) {
        if (io_mode == io_unformatted) {
            to_stream_unformatted(reinterpret_cast<char const *>(&val), sizeof(float));
        } else if (io_mode == io_list_directed) {
            utils::float_as_string_list_directed conv(val);
            to_stream(conv.begin, conv.n);
            prev_was_string = false;
        } else {
            to_stream_fmt_double(static_cast<double>(val));
        }
        return *this;
    }
    write_loop &operator,(double const &val) {
        if (io_mode == io_unformatted) {
            to_stream_unformatted(reinterpret_cast<char const *>(&val), sizeof(double));
        } else if (io_mode == io_list_directed) {
            utils::double_as_string_list_directed conv(val);
            to_stream(conv.begin, conv.n);
            prev_was_string = false;
        } else {
            to_stream_fmt_double(val);
        }
        return *this;
    }
    // Convert long double to a reasonable Fortran-ish textual form.
    // This is not a perfect Fortran formatting emulation, but it avoids crashes
    // and preserves enough precision for test logs.
    static inline void long_double_to_chars(char *buf, std::size_t buf_sz, long double v) {
        // Handle NaN/Inf explicitly for stable output across libcs.
        if (std::isnan(v)) {
            std::snprintf(buf, buf_sz, "NaN");
            return;
        }
        if (std::isinf(v)) {
            if (v > 0)
                std::snprintf(buf, buf_sz, "Inf");
            else
                std::snprintf(buf, buf_sz, "-Inf");
            return;
        }

        // Use max_digits10 to round-trip when possible.
        // %.*Lg uses either fixed or scientific depending on magnitude (like G format).
        int prec = std::numeric_limits<long double>::max_digits10;
        if (prec < 1)
            prec = 1;
        // Leave room; if snprintf fails, fall back to something minimal.
        int n = std::snprintf(buf, buf_sz, "%.*Lg", prec, v);
        if (n <= 0 || static_cast<std::size_t>(n) >= buf_sz) {
            // Fallback: shorter precision
            std::snprintf(buf, buf_sz, "%.18Lg", v);
        }
    }
    write_loop &operator,(long double const &val) {
        if (io_mode == io_unformatted) {
            to_stream_unformatted(reinterpret_cast<char const *>(&val), actual_sizeof_long_double);
        } else if (io_mode == io_list_directed) {
            // List-directed: emit a space separator if needed, then value.
            // This follows the existing behavior for doubles (prev_was_string).
            char tmp[256];
            long_double_to_chars(tmp, sizeof(tmp), val);
            to_stream(tmp, std::strlen(tmp));
            prev_was_string = false;
        } else {
            // Formatted: for binary80 long double builds that use snprintf("%Lg") I/O,
            // fall back to sprintnum_short() so edit-descriptor width handling is at least consistent.
#if defined(MPLAPACK_BINARY80_IO) && (MPLAPACK_BINARY80_IO == MPLAPACK_BINARY80_IO_SNPRINTF_LDBL)
            char buf[MPLAPACK_BUFLEN];
            buf[0] = '\0';
            // Avoid sprintnum_short(...) overload ambiguity (qd/gmp/etc).
            // Generate a "short" long double representation directly.
            // NOTE: This is intentionally minimal; strict Fortran formatting is not implemented here.
            std::snprintf(buf, sizeof(buf), "%+24.20Le", val);
            buf[sizeof(buf) - 1] = '\0';
            std::string s(buf);
            size_t b = s.find_first_not_of(' ');
            size_t e = s.find_last_not_of(' ');
            if (b != std::string::npos && e != std::string::npos) {
                s = s.substr(b, e - b + 1);
            }
            std::string const &ed = next_edit_descriptor();
            to_stream_fmt_double_given_string(s, ed);
#elif defined(MPLAPACK_BINARY128_IO) && (MPLAPACK_BINARY128_IO == MPLAPACK_BINARY128_IO_SNPRINTF_LDBL)
            char buf[MPLAPACK_BUFLEN];
            buf[0] = '\0';
            // Avoid sprintnum_short(...) overload ambiguity (qd/gmp/etc).
            // Generate a "short" long double (binary128 on arm64) representation directly.
            // binary128: 112-bit mantissa => ~33.97 decimal digits; use 35 for safe round-trip.
            // Exponent range: +-16383 (5 digits) => field width 46 is sufficient.
            // NOTE: This is intentionally minimal; strict Fortran formatting is not implemented here.
            std::snprintf(buf, sizeof(buf), "%+46.35Le", val);
            buf[sizeof(buf) - 1] = '\0';
            std::string s(buf);
            size_t b = s.find_first_not_of(' ');
            size_t e = s.find_last_not_of(' ');
            if (b != std::string::npos && e != std::string::npos) {
                s = s.substr(b, e - b + 1);
            }
            std::string const &ed = next_edit_descriptor();
            to_stream_fmt_double_given_string(s, ed);
#else
            // Generic fallback: emit the textual value (does not honor edit-descriptor width/scale).
            char tmp[256];
            long_double_to_chars(tmp, sizeof(tmp), val);
            to_stream(tmp, std::strlen(tmp));
#endif
            prev_was_string = false;
        }
        return *this;
    }
    // Generic output path for MPLAPACK MP types:
    // If ::sprintnum(char*, T) exists, stringify and forward to existing string output.
    //
    // Fallback for MPLAPACK numeric types providing sprintnum_short().
    // Converts value to string and outputs with edit descriptor formatting.
    template <typename T, typename = typename std::enable_if<!std::is_array<T>::value, decltype(sprintnum_short(static_cast<char *>(nullptr), std::declval<T const &>()))>::type> write_loop &operator,(T const &val) {
        if (io_mode == io_unformatted) {
            out.reset();
            throw TBXX_NOT_IMPLEMENTED();
        }
        char buf[MPLAPACK_BUFLEN];
        buf[0] = '\0';
        sprintnum_short(buf, val);
        buf[MPLAPACK_BUFLEN - 1] = '\0';
        if (io_mode == io_list_directed) {
            to_stream(buf, std::strlen(buf));
            prev_was_string = false;
        } else {
            // Use to_stream_fmt_double_given_string() to properly handle edit descriptors,
            // consistent with dd_real/qd_real explicit overloads.
            std::string s(buf);
            size_t b = s.find_first_not_of(' ');
            size_t e = s.find_last_not_of(' ');
            if (b != std::string::npos && e != std::string::npos) {
                s = s.substr(b, e - b + 1);
            }
            std::string const &ed = next_edit_descriptor();
            to_stream_fmt_double_given_string(s, ed);
        }
        return *this;
    }
    write_loop &operator,(std::complex<float> const &val) {
        if (io_mode == io_unformatted) {
            float re = val.real();
            float im = val.imag();
            to_stream_unformatted(reinterpret_cast<char const *>(&re), sizeof(float));
            to_stream_unformatted(reinterpret_cast<char const *>(&im), sizeof(float));
        } else if (io_mode == io_list_directed) {
            utils::float_as_string_list_directed conv_re(val.real());
            utils::float_as_string_list_directed conv_im(val.imag());
            to_stream_star_complex(conv_re.begin, conv_re.n, conv_im.begin, conv_im.n);
        } else {
            out.reset();
            throw TBXX_NOT_IMPLEMENTED();
        }
        return *this;
    }
    write_loop &operator,(std::complex<double> const &val) {
        if (io_mode == io_unformatted) {
            double re = val.real();
            double im = val.imag();
            to_stream_unformatted(reinterpret_cast<char const *>(&re), sizeof(double));
            to_stream_unformatted(reinterpret_cast<char const *>(&im), sizeof(double));
        } else if (io_mode == io_list_directed) {
            utils::double_as_string_list_directed conv_re(val.real());
            utils::double_as_string_list_directed conv_im(val.imag());
            to_stream_star_complex(conv_re.begin, conv_re.n, conv_im.begin, conv_im.n);
        } else {
            out.reset();
            throw TBXX_NOT_IMPLEMENTED();
        }
        return *this;
    }
    write_loop &operator,(std::complex<long double> const &val) {
        if (io_mode == io_unformatted) {
            long double re = val.real();
            long double im = val.imag();
            to_stream_unformatted(reinterpret_cast<char const *>(&re), actual_sizeof_long_double);
            to_stream_unformatted(reinterpret_cast<char const *>(&im), actual_sizeof_long_double);
        } else if (io_mode == io_list_directed) {
            // List-directed complex: reuse the same "* complex" emitter as double.
            char re_buf[256];
            char im_buf[256];
            long_double_to_chars(re_buf, sizeof(re_buf), val.real());
            long_double_to_chars(im_buf, sizeof(im_buf), val.imag());
            to_stream_star_complex(re_buf, std::strlen(re_buf), im_buf, std::strlen(im_buf));
        } else {
            // Formatted complex output is not fully implemented.
            // Fallback to a reasonable textual representation to avoid abort.
            char re_buf[256];
            char im_buf[256];
            long_double_to_chars(re_buf, sizeof(re_buf), val.real());
            long_double_to_chars(im_buf, sizeof(im_buf), val.imag());
            to_stream_star_complex(re_buf, std::strlen(re_buf), im_buf, std::strlen(im_buf));
        }
        return *this;
    }
    write_loop &operator,(str_cref const &val) {
        if (io_mode == io_unformatted) {
            to_stream_unformatted(val.elems(), val.len());
        } else if (io_mode == io_list_directed) {
            to_stream(val.elems(), val.len(), /*space*/ !prev_was_string);
            prev_was_string = true;
        } else {
            std::string const &ed = next_edit_descriptor();
            int n = ed.size();
            if (ed[0] == 'a' && n > 1) {
                n = utils::unsigned_integer_value(ed.data() + 1, n - 1);
                to_stream(val.elems(), std::min(val.len(), n));
                for (int i = val.len(); i < n; i++) {
                    to_stream(" ", 1);
                }
            } else {
                to_stream(val.elems(), val.len());
            }
        }
        return *this;
    }
    write_loop &operator,(str_addends const &val) {
        int ll = val.lhs.len();
        int rl = val.rhs.len();
        utils::simple_buffer<char> buffer(ll + rl);
        char *b = buffer.space;
        std::memcpy(b, val.lhs.elems(), ll);
        std::memcpy(b + ll, val.rhs.elems(), rl);
        (*this), str_cref(b, ll + rl);
        return *this;
    }
    template <typename T, size_t Ndims> write_loop &operator,(arr_cref<T, Ndims> const &val) {
        size_t n = val.size_1d();
        T const *val_begin = val.begin();
        for (size_t i = 0; i < n; i++) {
            (*this), val_begin[i];
        }
        return *this;
    }
    template <size_t Ndims> write_loop &operator,(str_arr_cref<Ndims> const &val) {
        size_t n = val.size_1d();
        int l = val.len();
        char const *val_begin = val.begin();
        for (size_t i = 0; i < n; i++) {
            (*this), str_cref(&val_begin[i * l], l);
        }
        return *this;
    }
    ~write_loop() NOEXCEPT_FALSE {
        if (out.get() == 0)
            return;
        if (internal_file_len < 0) {
            if (io_mode == io_unformatted) {
                out->put(end_of_unformatted_record);
                out->put('\0');
            } else {
                if (io_mode == io_list_directed) {
                    if (pos == 0)
                        out->put(' ');
                } else {
                    next_edit_descriptor(/*final*/ true);
                }
                if (!suppress_new_line_at_end) {
                    if (write_crlf) {
                        out->put("\r\n", 2);
                    } else {
                        out->put('\n');
                    }
                }
            }
            out->flush();
        } else {
            if (io_mode == io_unformatted) {
                out.reset();
                throw TBXX_NOT_IMPLEMENTED();
            } else {
                if (io_mode == io_list_directed) {
                    if (pos == 0) {
                        out->put(' ');
                        pos++;
                    }
                } else {
                    next_edit_descriptor(/*final*/ true);
                }
                while (pos < internal_file_len) {
                    out->put(' ');
                    pos++;
                }
            }
        }
    }

  private:
    void to_stream(char const *buf, unsigned n, bool space = true) {
        switch (io_mode) {
        case io_unformatted:
            to_stream_unformatted(buf, n);
            break;
        case io_list_directed:
            to_stream_star(buf, n, space);
            break;
        default:
            to_stream_fmt(buf, n);
        }
    }
    void to_stream_unformatted(char const *buf, unsigned n) {
        for (unsigned i = 0; i < n; i++) {
            char c = buf[i];
            out->put(c);
            if (c == end_of_unformatted_record) {
                out->put(c);
            }
        }
    }
    void to_stream_fmt(char const *buf, unsigned n) {
        while (number_of_x_held != 0) {
            out->put(" ", 1);
            number_of_x_held--;
        }
        out->put(buf, n);
    }
    void to_stream_star(char const *buf, unsigned n, bool space = true) {
        if (pos == 0) {
            out->put(' ');
            pos = 1;
        } else if (pos + (space ? 1 : 0) + n > 80) {
            if (write_crlf) {
                out->put("\r\n ", 3);
            } else {
                out->put("\n ", 2);
            }
            pos = 1;
        } else if (space) {
            out->put(' ');
            pos++;
        }
        out->put(buf, n);
        pos += n;
    }
    void to_stream_star_strip_leading_and_trailing_blank_padding(char const *buf, unsigned n, bool space = true) {
        size_t_2 indices = utils::find_leading_and_trailing_blank_padding(buf, n);
        to_stream_star(buf + indices.elems[0], indices.elems[1] - indices.elems[0],
                       /*space*/ false);
    }
    void to_stream_star_complex(char const *buf_re, unsigned n_re, char const *buf_im, unsigned n_im) {
        to_stream_star("(", 1);
        {
            to_stream_star_strip_leading_and_trailing_blank_padding(buf_re, n_re, /*space*/ false);
        }
        to_stream_star(",", 1, /*space*/ false);
        {
            to_stream_star_strip_leading_and_trailing_blank_padding(buf_im, n_im, /*space*/ false);
        }
        to_stream_star(")", 1, /*space*/ false);
        prev_was_string = false;
    }
};
struct common_write {
    common &cmn;
    common_write(common &cmn_) : cmn(cmn_) {}
    write_loop operator()(int unit, unformatted_type const &) {
        write_loop result(cmn, unit, unformatted);
        return result;
    }
    write_loop operator()(int unit, star_type const &) {
        write_loop result(cmn, unit, star);
        return result;
    }
    write_loop operator()(int const &unit, str_cref fmt) {
        write_loop result(cmn, unit, fmt);
        return result;
    }
    write_loop operator()(str_ref const &internal_file, star_type const &) {
        write_loop result(internal_file, star);
        return result;
    }
    write_loop operator()(str_ref const &internal_file, str_cref fmt) {
        write_loop result(internal_file, fmt);
        return result;
    }
    template <size_t Ndims> write_loop operator()(str_arr_ref<Ndims> const &internal_file, star_type const &) {
        write_loop result(internal_file, star);
        return result;
    }
    template <size_t Ndims> write_loop operator()(str_arr_ref<Ndims> const &internal_file, str_cref fmt) {
        write_loop result(internal_file, fmt);
        return result;
    }
};
} // namespace fem
#endif // GUARD

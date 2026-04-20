#ifndef FEM_UTILS_INT_TYPES_HPP
#define FEM_UTILS_INT_TYPES_HPP
namespace fem {
namespace utils {
    // int type names as used in boost/cstdint.hpp
    typedef signed char int8_t;
    typedef signed short int16_t;
    typedef signed int int32_t;
#if defined(_WIN32) || defined(__MINGW32__) || defined(__MINGW64__) || defined(__i386__) || defined(__ppc__) || defined(_MSC_VER)
    // Windows uses LLP64: long remains 32-bit even on 64-bit MinGW/MSVC.
    typedef signed long long int64_t;
#else
    typedef signed long int64_t;
#endif
} // namespace utils
} // namespace fem
#endif // GUARD

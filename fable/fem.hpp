#ifndef FEM_HPP
#define FEM_HPP
#include <fem/do.hpp>
#include <fem/data_of_type.hpp>
#include <fem/error_utils.hpp>
#include <fem/format.hpp>
#include <fem/intrinsics.hpp>
#include <fem/main.hpp>
#include <fem/major_types.hpp>
#include <chrono>
// LAPACK-like timer: seconds since an arbitrary epoch (first call).
// Use steady_clock for monotonic behavior.
inline double dsecnd() noexcept {
    using clock = std::chrono::steady_clock;
    static const auto t0 = clock::now();
    const auto t1 = clock::now();
    return std::chrono::duration<double>(t1 - t0).count();
}
#endif // GUARD

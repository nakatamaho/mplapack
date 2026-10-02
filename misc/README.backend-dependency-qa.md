# Backend dependency smoke tests (2.4)

The 2.4 public GMP API uses GMP's C++ wrapper. Its link interface therefore
includes gmpxx and gmp. The old MPFR API uses mpfrc++ and MPC and also exposes
GMP C++ types. Its interface includes mpc, mpfr, gmpxx, and gmp. DD and QD both
use the QD package. These dependencies are expressed by library LIBADD and
the public pkg-config compile/link interface, without requiring dependency
pkg-config modules. Static QD consumers also receive libm. Optimized static
consumers receive the detected OpenMP runtime through Libs.private.

Configure with OPENMP_LIBS to override runtime discovery for unusual compiler
configurations. --disable-openmp leaves the OpenMP link interface empty.
Compiler identification uses predefined macros, so wrappers and absolute
compiler paths do not change runtime selection.

Use an isolated installation prefix for these checks. Do not preload backend
libraries. A command-local loader path is only for testing the temporary
installation, not a production workaround.

For each of gmp, mpfr, qd, and dd, and each built reference/optimized variant:

```sh
PKG_CONFIG_PATH="$prefix/lib/pkgconfig" \
  sh misc/check_backend_pkgconfig.sh shared mplapack_gmp \
  misc/backend_gmp_pkgconfig_consumer.cpp
PKG_CONFIG_PATH="$prefix/lib/pkgconfig" \
  sh misc/check_backend_pkgconfig.sh static mplapack_gmp \
  misc/backend_gmp_pkgconfig_consumer.cpp
python3 misc/check_backend_shared_load.py "$library"
sh misc/check_backend_elf_dependencies.sh "$library" gmp
```

Substitute the backend and actual library filename. For MPFR check the mpc,
mpfr, and gmp dependency stems; for DD/QD check qd. Do not require unused
libraries to survive the linker's as-needed processing. The ELF check is
Linux-specific and exits 77 when its tools are unavailable. The static
pkg-config probe uses GNU linker archive selection and should be run on
Linux; other platforms need their equivalent static consumer checks.

Run the consumers against shared and static Autotools/CMake installations.
For CMake additionally test an installed package consumer using exported
mplapack targets. Inspect support-library relocations and libtool dependency
metadata as well. These smoke tests do not replace release numerical QA.

EIG and LIN support libraries additionally link their matching MATGEN,
override, and primary backend libraries. FABLE common-block globals such as
fs and iparms are intentionally defined by the test driver. Their unresolved
references in an isolated support-library check must not be mistaken for
missing GMP, MPC, MPFR, or QD dependencies. Check these libraries with the
driver's common-block definitions present, or classify those references
separately. Public backend libraries must have no unresolved relocations.
Use `sh misc/check_backend_support_libraries.sh <top-builddir>` for this
classification; unexpected unresolved component/backend symbols still fail.

After building the primary, MATGEN, LIN, and EIG libraries, run
`sh misc/check_backend_support_install.sh <top-builddir>` to exercise their
first installation into a fresh temporary DESTDIR. Use a configured prefix
without previously installed support libraries, since libtool may otherwise
find an old override there and mask an installation-order failure. The test
installs the prerequisites and then the support libraries using the generated
Automake rules, without running numerical tests. Each override must appear
before its dependent library in `lib_LTLIBRARIES` so install-time relinking
can find it. The temporary installation is removed when the test exits.

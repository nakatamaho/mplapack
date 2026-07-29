/*
 * MPLAPACK gmpfrxx_mkII integration contract.
 *
 * MPLAPACK fixes the symbolic external-provider mode before any gmpfrxx_mkII
 * header is included so every translation unit uses the provider implemented by
 * mplapackinit.cpp.
 */

#ifndef MPLAPACK_GMPFRXX_MKII_CONFIG_H
#define MPLAPACK_GMPFRXX_MKII_CONFIG_H

#if defined(GMPFRXX_MKII_DETAIL_CONFIG_HPP)
#if GMPXX_MKII_DEFAULT_CONTEXT_MODE != GMPXX_MKII_DEFAULT_CONTEXT_EXTERNAL_PROVIDER
#error "MPLAPACK requires the gmpfrxx_mkII external default-context provider"
#endif
#else
#ifdef GMPXX_MKII_DEFAULT_CONTEXT_MODE
#undef GMPXX_MKII_DEFAULT_CONTEXT_MODE
#endif
#define GMPXX_MKII_DEFAULT_CONTEXT_MODE GMPXX_MKII_DEFAULT_CONTEXT_EXTERNAL_PROVIDER
#endif

#endif

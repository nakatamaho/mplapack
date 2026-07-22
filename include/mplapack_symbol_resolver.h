/*
 * Copyright (c) 2008-2026
 *      Nakata, Maho
 *      All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE AUTHOR AND CONTRIBUTORS ``AS IS'' AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED.  IN NO EVENT SHALL THE AUTHOR OR CONTRIBUTORS BE LIABLE
 * FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
 * DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS
 * OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
 * HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
 * LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
 * OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF
 * SUCH DAMAGE.
 *
 */

// Dynamic symbol resolver for MPLAPACK benchmark
// Eliminates the need for compiler/version-specific mangled symbol headers.
//
// Requires: C++17 or later (uses inline variables/functions for ODR safety
//           across multiple translation units).
//
// Supported platforms:
//   - Linux:  ELF dynamic symbol table parsing via DT_HASH / DT_GNU_HASH
//   - macOS:  Mach-O symbol table parsing via LC_SYMTAB + dyld APIs
//   - MinGW:  PE export table parsing via IMAGE_EXPORT_DIRECTORY
//             using dlfcn-win32 for dlopen/dlsym/dlclose
//
// Note: This resolver does not use an external `nm` fallback.
//
// Usage:
//   void *sym = mplapack_resolver::resolve_symbol(handle, "Rgemm");         // verbose (default)
//   void *sym = mplapack_resolver::resolve_symbol(handle, "Rgemm", false);  // silent
//
// Verbose output (stderr):
//   resolve_symbol: "Rgemm" -> "_Z5RgemmPKcS0_lllgPglS1_lgS1_l"
//                   in /usr/local/lib/libmplapack_dd.so [via elf, scanned 342 syms]
//   resolve_symbol: "Raxpy" NOT FOUND in /usr/local/lib/libmplapack_dd.so [scanned 342 syms via elf]
//
// Notes:
//   - Matches by demangled pattern "FuncName(" with boundary detection,
//     so overloaded functions resolve to whichever is found first
//     (non-deterministic).
//   - Namespace-qualified names (e.g. mplapack::Rgemm) ARE matched
//     via boundary-aware substring search.
//   - MSVC-demangled names with return type / calling convention prefixes
//     (e.g. "void __cdecl Rgemm(") are also handled correctly.
//   - Results are cached per (handle, func_name) pair for the process
//     lifetime (single cache shared across all translation units).
//     Call resolve_cache_clear() to invalidate.

#ifndef MPLAPACK_SYMBOL_RESOLVER_H
#define MPLAPACK_SYMBOL_RESOLVER_H

// ---- Build-time guard ----
// This resolver relies on the GCC/Clang C++ ABI demangler (cxxabi.h).
// It is intended for MinGW / GCC / Clang toolchains, not MSVC.
#if defined(_WIN32) && defined(_MSC_VER) && !defined(__MINGW32__) && !defined(__MINGW64__)
#error "mplapack_symbol_resolver: MSVC is not supported (requires <cxxabi.h> / Itanium demangling). Use MinGW/Clang, or add an MSVC-only demangler path."
#endif

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cstdint>
#include <string>
#include <unordered_map>
#include <utility>

// ---- Platform headers ----
#if defined(_WIN32)
// dlfcn-win32 provides POSIX-compatible dlopen/dlsym/dlclose/dlerror.
// On MinGW, dlfcn-win32's dlopen() returns the HMODULE directly cast
// to void*, so (HMODULE)handle is a valid conversion.
// We do NOT statically link dbghelp; it is loaded on demand at runtime
// to support MSVC-mangled symbols when present.
#include <windows.h>
#include <dlfcn.h>

// Some MinGW headers may not define these search flags.
#ifndef LOAD_LIBRARY_SEARCH_SYSTEM32
#define LOAD_LIBRARY_SEARCH_SYSTEM32 0x00000800
#endif
#elif defined(__APPLE__)
#include <dlfcn.h>
#if defined(__cplusplus) && defined(__GNUC__) && !defined(__clang__) && !defined(_Static_assert)
#define MPLAPACK_RESTORE_STATIC_ASSERT_FOR_MACOS26_SDK 1
#define _Static_assert static_assert
#endif
#include <mach-o/dyld.h>
#include <mach-o/loader.h>
#include <mach-o/nlist.h>
#if defined(MPLAPACK_RESTORE_STATIC_ASSERT_FOR_MACOS26_SDK)
#undef _Static_assert
#undef MPLAPACK_RESTORE_STATIC_ASSERT_FOR_MACOS26_SDK
#endif
#else
#include <dlfcn.h>
#include <elf.h>
#include <link.h>
#endif

#include <cxxabi.h>

namespace mplapack_resolver {

// ============================================================
//  Result cache (process-wide via C++17 inline)
// ============================================================
// Key:   (handle address as uintptr_t, func_name)
// Value: resolved pointer (nullptr means "confirmed not found")
//
// All functions and variables in this header are 'inline' so that
// multiple TUs including this header share a single definition
// (and thus a single cache instance) per C++17 ODR guarantees.

inline const void *kNotCached = reinterpret_cast<const void *>(~uintptr_t(0));

struct CacheKeyHash {
    size_t operator()(const std::pair<uintptr_t, std::string> &k) const {
        size_t h1 = std::hash<uintptr_t>()(k.first);
        size_t h2 = std::hash<std::string>()(k.second);
        return h1 ^ (h2 * 0x9e3779b97f4a7c15ULL + 0x517cc1b727220a95ULL);
    }
};

inline std::unordered_map<std::pair<uintptr_t, std::string>, void *, CacheKeyHash> &resolve_cache_map() {
    static std::unordered_map<std::pair<uintptr_t, std::string>, void *, CacheKeyHash> cache;
    return cache;
}

inline void *cache_lookup(void *handle, const char *func_name) {
    auto &cache = resolve_cache_map();
    auto key = std::make_pair(reinterpret_cast<uintptr_t>(handle), std::string(func_name));
    auto it = cache.find(key);
    if (it != cache.end())
        return it->second;
    return const_cast<void *>(kNotCached);
}

inline void cache_store(void *handle, const char *func_name, void *result) {
    auto &cache = resolve_cache_map();
    auto key = std::make_pair(reinterpret_cast<uintptr_t>(handle), std::string(func_name));
    cache[key] = result;
}

// Clear the entire resolution cache (e.g. after dlclose + re-dlopen).
inline void resolve_cache_clear() { resolve_cache_map().clear(); }

// ============================================================
//  Windows: dlfcn-win32 handle -> HMODULE conversion
// ============================================================
#if defined(_WIN32)

inline HMODULE handle_to_hmodule(void *handle) {
    if (!handle)
        return nullptr;

    // Fast path: dlfcn-win32 returns HMODULE directly.
    // Validate by probing the memory region with VirtualQuery.
    // This is safe on MinGW (no SEH __try/__except needed).
    HMODULE hmod = reinterpret_cast<HMODULE>(handle);
    MEMORY_BASIC_INFORMATION mbi;
    if (VirtualQuery(hmod, &mbi, sizeof(mbi)) != 0) {
        // The allocation base of a loaded module equals the HMODULE.
        if (mbi.AllocationBase == hmod) {
            const IMAGE_DOS_HEADER *dos = reinterpret_cast<const IMAGE_DOS_HEADER *>(hmod);
            if (dos->e_magic == IMAGE_DOS_SIGNATURE)
                return hmod;
        }
    }

    // Slow path: resolve any symbol from the library, then use
    // GetModuleHandleExA(FROM_ADDRESS) to recover the HMODULE.
    static const char *probe_syms[] = {"_init", "_fini", "DllMain", nullptr};
    for (const char **p = probe_syms; *p; ++p) {
        void *addr = dlsym(handle, *p);
        if (addr) {
            HMODULE recovered = nullptr;
            if (GetModuleHandleExA(GET_MODULE_HANDLE_EX_FLAG_FROM_ADDRESS | GET_MODULE_HANDLE_EX_FLAG_UNCHANGED_REFCOUNT, reinterpret_cast<LPCSTR>(addr), &recovered) && recovered) {
                return recovered;
            }
        }
    }

    // Last resort: trust the dlfcn-win32 contract
    return hmod;
}

// ============================================================
//  Windows: dbghelp UnDecorateSymbolName (dynamically loaded)
// ============================================================

typedef DWORD(WINAPI *UnDecorateSymbolName_t)(const char *, char *, DWORD, DWORD);

inline UnDecorateSymbolName_t get_undecorate_fn() {
    static UnDecorateSymbolName_t fn = nullptr;
    static bool tried = false;
    if (!tried) {
        tried = true;
        // Prefer loading from %SystemRoot%\System32 to reduce DLL preloading risk.
        HMODULE hDbg = LoadLibraryExA("dbghelp.dll", nullptr, LOAD_LIBRARY_SEARCH_SYSTEM32);
        if (!hDbg)
            hDbg = LoadLibraryA("dbghelp.dll");
        if (hDbg)
            fn = reinterpret_cast<UnDecorateSymbolName_t>(GetProcAddress(hDbg, "UnDecorateSymbolName"));
        // Intentionally leak hDbg; it stays loaded for process lifetime.
    }
    return fn;
}

#ifndef UNDNAME_COMPLETE
#define UNDNAME_COMPLETE 0x0000
#endif
#ifndef UNDNAME_NO_ACCESS_SPECIFIERS
#define UNDNAME_NO_ACCESS_SPECIFIERS 0x0080
#endif

#endif // _WIN32

// ============================================================
//  Helper: get library file path from a dlopen handle
// ============================================================
inline const char *get_library_path(void *handle) {
    // Note: this buffer is shared (not thread-safe). Acceptable for
    // benchmark use; if thread safety is needed, use thread_local.
    static char path_buf[4096];
    path_buf[0] = '\0';

#if defined(__linux__)
    struct link_map *lm = nullptr;
    if (dlinfo(handle, RTLD_DI_LINKMAP, &lm) == 0 && lm && lm->l_name[0]) {
        strncpy(path_buf, lm->l_name, sizeof(path_buf) - 1);
        path_buf[sizeof(path_buf) - 1] = '\0';
    }
#elif defined(__APPLE__)
    uint32_t image_count = _dyld_image_count();
    for (uint32_t i = 0; i < image_count; ++i) {
        const char *name = _dyld_get_image_name(i);
        if (!name)
            continue;
        void *test = dlopen(name, RTLD_LAZY | RTLD_NOLOAD);
        if (test == handle) {
            strncpy(path_buf, name, sizeof(path_buf) - 1);
            path_buf[sizeof(path_buf) - 1] = '\0';
            dlclose(test);
            break;
        }
        if (test)
            dlclose(test);
    }
#elif defined(_WIN32)
    HMODULE hmod = handle_to_hmodule(handle);
    if (hmod)
        GetModuleFileNameA(hmod, path_buf, sizeof(path_buf));
#endif

    if (!path_buf[0])
        return "(unknown library)";
    return path_buf;
}

// ============================================================
//  Helper: boundary-aware substring match for "FuncName("
// ============================================================
//
// Checks whether the demangled string contains "FuncName(" at a
// valid word boundary.  This correctly handles:
//
//   "Rgemm("                     -- plain global function
//   "mplapack::Rgemm("           -- namespace-qualified
//   "void __cdecl Rgemm("        -- MSVC with return type prefix
//   "dd_real Rgemm("             -- return type before name
//
// A "boundary" is: start of string, or a preceding character that
// is not alphanumeric and not '_' (so "fooRgemm(" does NOT match).
//
inline bool contains_func_call(const char *demangled, const char *target_with_paren, size_t target_len) {
    const char *p = demangled;
    while ((p = strstr(p, target_with_paren)) != nullptr) {
        if (p == demangled) {
            return true;
        }
        char prev = *(p - 1);
        // Valid boundary: anything that is NOT part of a C++ identifier
        if (prev != '_' && !(prev >= 'a' && prev <= 'z') && !(prev >= 'A' && prev <= 'Z') && !(prev >= '0' && prev <= '9')) {
            return true;
        }
        p += target_len;
    }
    return false;
}

// ============================================================
//  Helper: demangle and check if name matches "FuncName("
//  If matched, copies the original mangled name to out_mangled.
// ============================================================
inline bool demangled_matches(const char *mangled, const char *target, size_t target_len, char *out_mangled, size_t out_mangled_size) {
    // ---- Itanium ABI demangling (GCC / Clang / MinGW) ----
    int status = -1;
    char *demangled = abi::__cxa_demangle(mangled, nullptr, nullptr, &status);
    if (status == 0 && demangled) {
        bool match = contains_func_call(demangled, target, target_len);
        if (match && out_mangled && out_mangled_size > 0) {
            strncpy(out_mangled, mangled, out_mangled_size - 1);
            out_mangled[out_mangled_size - 1] = '\0';
        }
        free(demangled);
        return match;
    }
    free(demangled);

    // Some MinGW toolchains / PE exports prefix Itanium-mangled names
    // with an extra underscore: "__Z..." instead of "_Z...".
    // Strip the leading '_' and retry demangling.
    const char *itanium = nullptr;
    if (mangled[0] == '_' && mangled[1] == '_' && mangled[2] == 'Z')
        itanium = mangled + 1; // "__Z..." -> "_Z..."
    else if (mangled[0] == '_' && mangled[1] == 'Z')
        itanium = mangled; // already "_Z..." (no extra prefix)

    if (itanium && itanium != mangled) {
        status = -1;
        demangled = abi::__cxa_demangle(itanium, nullptr, nullptr, &status);
        if (status == 0 && demangled) {
            bool match = contains_func_call(demangled, target, target_len);
            if (match && out_mangled && out_mangled_size > 0) {
                // Store the ORIGINAL export name (with extra '_') so that
                // dlsym / GetProcAddress uses the correct lookup key.
                strncpy(out_mangled, mangled, out_mangled_size - 1);
                out_mangled[out_mangled_size - 1] = '\0';
            }
            free(demangled);
            return match;
        }
        free(demangled);
    }

#if defined(_WIN32)
    // ---- MSVC demangling via dbghelp (for DLLs built with MSVC) ----
    // MSVC-demangled names typically look like:
    //   "void __cdecl Rgemm(char const *,...)"
    //   "class dd_real __cdecl Rgemm(char const *,...)"
    // so we use boundary-aware substring matching.
    UnDecorateSymbolName_t undecorate = get_undecorate_fn();
    if (undecorate) {
        char buf[1024];
        DWORD result = undecorate(mangled, buf, sizeof(buf), UNDNAME_COMPLETE | UNDNAME_NO_ACCESS_SPECIFIERS);
        if (result > 0) {
            bool match = contains_func_call(buf, target, target_len);
            if (match && out_mangled && out_mangled_size > 0) {
                strncpy(out_mangled, mangled, out_mangled_size - 1);
                out_mangled[out_mangled_size - 1] = '\0';
            }
            return match;
        }
    }
#endif

    return false;
}

// ============================================================
//  Linux: ELF dynamic symbol table
// ============================================================
#if defined(__linux__)

inline void *resolve_via_elf(void *handle, const char *func_name, char *out_mangled, size_t out_mangled_size, size_t *out_nsyms) {
    *out_nsyms = 0;

    struct link_map *lm = nullptr;
    if (dlinfo(handle, RTLD_DI_LINKMAP, &lm) != 0 || !lm)
        return nullptr;

    ElfW(Dyn) *dyn = lm->l_ld;
    if (!dyn)
        return nullptr;

    ElfW(Sym) *symtab = nullptr;
    const char *strtab = nullptr;
    ElfW(Word) *hashtab = nullptr;
    uint32_t *gnu_hashtab = nullptr;

    for (ElfW(Dyn) *d = dyn; d->d_tag != DT_NULL; ++d) {
        switch (d->d_tag) {
        case DT_SYMTAB:
            symtab = reinterpret_cast<ElfW(Sym) *>(d->d_un.d_ptr);
            break;
        case DT_STRTAB:
            strtab = reinterpret_cast<const char *>(d->d_un.d_ptr);
            break;
        case DT_HASH:
            hashtab = reinterpret_cast<ElfW(Word) *>(d->d_un.d_ptr);
            break;
        case DT_GNU_HASH:
            gnu_hashtab = reinterpret_cast<uint32_t *>(d->d_un.d_ptr);
            break;
        }
    }

    if (!symtab || !strtab)
        return nullptr;

    // Determine symbol count.
    // Prefer DT_HASH (nchain is exact) over DT_GNU_HASH (estimated).
    size_t nsyms = 0;

    if (hashtab) {
        // DT_HASH: [nbucket, nchain, ...]; nchain == total symbol count
        nsyms = hashtab[1];
    } else if (gnu_hashtab) {
        // DT_GNU_HASH layout:
        //   uint32_t nbuckets, symoffset, bloom_size, bloom_shift
        //   ElfW(Addr) bloom[bloom_size]
        //   uint32_t buckets[nbuckets]
        //   uint32_t chains[...]
        uint32_t nbuckets = gnu_hashtab[0];
        uint32_t symoffset = gnu_hashtab[1];
        uint32_t bloom_size = gnu_hashtab[2];

        const uint32_t *buckets = &gnu_hashtab[4 + bloom_size * (sizeof(ElfW(Addr)) / 4)];
        const uint32_t *chains = &buckets[nbuckets];

        // Find the maximum symbol index across all buckets
        uint32_t last_sym = 0;
        for (uint32_t i = 0; i < nbuckets; ++i) {
            if (buckets[i] > last_sym)
                last_sym = buckets[i];
        }

        if (last_sym >= symoffset) {
            // Walk chain from last_sym until terminator bit (bit 0) is set
            uint32_t idx = last_sym - symoffset;
            while (!(chains[idx] & 1))
                ++idx;
            nsyms = last_sym + (idx - (last_sym - symoffset)) + 1;
        }
    }

    *out_nsyms = nsyms;
    if (nsyms == 0)
        return nullptr;

    char target[256];
    snprintf(target, sizeof(target), "%s(", func_name);
    size_t target_len = strlen(target);

    for (size_t i = 0; i < nsyms; ++i) {
        if (symtab[i].st_shndx == SHN_UNDEF)
            continue;
        unsigned char st_type = ELF64_ST_TYPE(symtab[i].st_info);
        if (st_type != STT_FUNC && st_type != STT_GNU_IFUNC)
            continue;

        const char *mangled = strtab + symtab[i].st_name;
        if (!mangled[0])
            continue;

        if (demangled_matches(mangled, target, target_len, out_mangled, out_mangled_size))
            return dlsym(handle, mangled);
    }
    return nullptr;
}

#endif // __linux__

// ============================================================
//  macOS: Mach-O symbol table via dyld APIs
// ============================================================
#if defined(__APPLE__)

inline void *resolve_via_macho(void *handle, const char *func_name, char *out_mangled, size_t out_mangled_size, size_t *out_nsyms) {
    *out_nsyms = 0;

    char target[256];
    snprintf(target, sizeof(target), "%s(", func_name);
    size_t target_len = strlen(target);

    uint32_t image_count = _dyld_image_count();

    for (uint32_t img = 0; img < image_count; ++img) {
        // Only scan the image that corresponds to this handle.
        const char *img_name = _dyld_get_image_name(img);
        if (!img_name)
            continue;
        void *img_handle = dlopen(img_name, RTLD_LAZY | RTLD_NOLOAD);
        bool is_target = (img_handle == handle);
        if (img_handle)
            dlclose(img_handle);
        if (!is_target)
            continue;

        const struct mach_header_64 *hdr = reinterpret_cast<const struct mach_header_64 *>(_dyld_get_image_header(img));
        if (!hdr || hdr->magic != MH_MAGIC_64)
            continue;

        intptr_t slide = _dyld_get_image_vmaddr_slide(img);

        const uint8_t *cmd_ptr = reinterpret_cast<const uint8_t *>(hdr + 1);
        const struct symtab_command *symtab_cmd = nullptr;
        const struct segment_command_64 *linkedit_seg = nullptr;
        const struct segment_command_64 *text_seg = nullptr;

        for (uint32_t c = 0; c < hdr->ncmds; ++c) {
            const struct load_command *lc = reinterpret_cast<const struct load_command *>(cmd_ptr);

            if (lc->cmd == LC_SYMTAB) {
                symtab_cmd = reinterpret_cast<const struct symtab_command *>(lc);
            } else if (lc->cmd == LC_SEGMENT_64) {
                const struct segment_command_64 *seg = reinterpret_cast<const struct segment_command_64 *>(lc);
                if (strcmp(seg->segname, SEG_LINKEDIT) == 0)
                    linkedit_seg = seg;
                else if (strcmp(seg->segname, SEG_TEXT) == 0)
                    text_seg = seg;
            }
            cmd_ptr += lc->cmdsize;
        }

        if (!symtab_cmd || !linkedit_seg || !text_seg)
            continue;

        uintptr_t linkedit_base = slide + linkedit_seg->vmaddr - linkedit_seg->fileoff;

        const struct nlist_64 *symtab = reinterpret_cast<const struct nlist_64 *>(linkedit_base + symtab_cmd->symoff);
        const char *strtab = reinterpret_cast<const char *>(linkedit_base + symtab_cmd->stroff);

        *out_nsyms = symtab_cmd->nsyms;

        for (uint32_t s = 0; s < symtab_cmd->nsyms; ++s) {
            const struct nlist_64 *nl = &symtab[s];

            if (nl->n_type & N_STAB)
                continue;
            if ((nl->n_type & N_TYPE) == N_UNDF)
                continue;

            uint32_t strx = nl->n_un.n_strx;
            const char *name = strtab + strx;
            if (!name[0])
                continue;

            // Mach-O external symbols have a leading underscore.
            // Strip it for demangling, but keep the original for dlsym fallback.
            const char *mangled = name;
            if (name[0] == '_')
                mangled = name + 1;

            if (demangled_matches(mangled, target, target_len, out_mangled, out_mangled_size)) {
                // Try the Itanium mangled name (without leading '_')
                void *sym = dlsym(handle, mangled);
                if (sym)
                    return sym;
                // Try the raw Mach-O symbol name (with leading '_')
                sym = dlsym(handle, name);
                if (sym)
                    return sym;
            }
        }
        // We found and scanned the target image; stop iterating.
        break;
    }
    return nullptr;
}

#endif // __APPLE__

// ============================================================
//  MinGW / Windows: PE export table parsing (dlfcn-win32)
// ============================================================
#if defined(_WIN32)

inline void *resolve_via_pe(void *handle, const char *func_name, char *out_mangled, size_t out_mangled_size, size_t *out_nsyms) {
    *out_nsyms = 0;

    HMODULE hmod = handle_to_hmodule(handle);
    if (!hmod)
        return nullptr;

    const IMAGE_DOS_HEADER *dos_hdr = reinterpret_cast<const IMAGE_DOS_HEADER *>(hmod);
    if (dos_hdr->e_magic != IMAGE_DOS_SIGNATURE)
        return nullptr;

    const IMAGE_NT_HEADERS *nt_hdr = reinterpret_cast<const IMAGE_NT_HEADERS *>(reinterpret_cast<const uint8_t *>(hmod) + dos_hdr->e_lfanew);
    if (nt_hdr->Signature != IMAGE_NT_SIGNATURE)
        return nullptr;

    DWORD export_rva, export_size;
#if defined(_WIN64)
    const IMAGE_OPTIONAL_HEADER64 *opt = reinterpret_cast<const IMAGE_OPTIONAL_HEADER64 *>(&nt_hdr->OptionalHeader);
#else
    const IMAGE_OPTIONAL_HEADER32 *opt = reinterpret_cast<const IMAGE_OPTIONAL_HEADER32 *>(&nt_hdr->OptionalHeader);
#endif
    export_rva = opt->DataDirectory[IMAGE_DIRECTORY_ENTRY_EXPORT].VirtualAddress;
    export_size = opt->DataDirectory[IMAGE_DIRECTORY_ENTRY_EXPORT].Size;

    if (export_rva == 0 || export_size == 0)
        return nullptr;

    const uint8_t *base = reinterpret_cast<const uint8_t *>(hmod);
    const IMAGE_EXPORT_DIRECTORY *exports = reinterpret_cast<const IMAGE_EXPORT_DIRECTORY *>(base + export_rva);

    const DWORD *names = reinterpret_cast<const DWORD *>(base + exports->AddressOfNames);
    const WORD *ordinals = reinterpret_cast<const WORD *>(base + exports->AddressOfNameOrdinals);
    const DWORD *funcs = reinterpret_cast<const DWORD *>(base + exports->AddressOfFunctions);

    *out_nsyms = exports->NumberOfNames;

    char target[256];
    snprintf(target, sizeof(target), "%s(", func_name);
    size_t target_len = strlen(target);

    for (DWORD i = 0; i < exports->NumberOfNames; ++i) {
        const char *mangled = reinterpret_cast<const char *>(base + names[i]);
        if (!mangled || !mangled[0])
            continue;

        if (demangled_matches(mangled, target, target_len, out_mangled, out_mangled_size)) {
            WORD ord = ordinals[i];

            // Bounds check: ordinal must index into the functions array.
            if (ord >= exports->NumberOfFunctions)
                continue;

            DWORD func_rva = funcs[ord];

            // RVA 0 means "not exported" (placeholder entry).
            if (func_rva == 0)
                continue;

            // Prefer dlsym (= GetProcAddress) — it handles forwarders,
            // API sets, and other special export forms correctly.
            void *sym = dlsym(handle, mangled);
            if (sym)
                return sym;

            // Forwarder exports have RVA within the export directory.
            // If dlsym already failed on a forwarder, skip it.
            if (func_rva >= export_rva && func_rva < export_rva + export_size)
                continue;

            // Last resort: compute address from base + RVA directly.
            return reinterpret_cast<void *>(const_cast<uint8_t *>(base) + func_rva);
        }
    }
    return nullptr;
}

#endif // _WIN32

// ============================================================
//  Main entry point
// ============================================================

// Resolve a C++ function by its unmangled base name (e.g. "Rgemm").
// Searches the dynamic/exported symbol table of the library opened
// via dlopen(). Returns the function pointer, or nullptr on failure.
//
// If verbose == true (default), prints diagnostic info to stderr:
//   resolve_symbol: "Rgemm" -> "_Z5RgemmPKcS0_lllgPglS1_lgS1_l"
//                   in /usr/local/lib/libmplapack_dd.so [via elf, scanned 342 syms]
//   resolve_symbol: "Raxpy" NOT FOUND in /usr/local/lib/libmplapack_dd.so [scanned 342 syms via elf]
//
// Resolution order:
//   1. Cache lookup
//   2. Native table parsing (ELF / Mach-O / PE)
inline void *resolve_symbol(void *handle, const char *func_name, bool verbose = true) {
    if (!handle || !func_name)
        return nullptr;

    // ---- Cache lookup ----
    void *cached = cache_lookup(handle, func_name);
    if (cached != kNotCached) {
        if (verbose && cached) {
            const char *lib = get_library_path(handle);
            fprintf(stderr, "resolve_symbol: \"%s\" (cached) in %s\n", func_name, lib);
        }
        return cached;
    }

    char mangled_buf[1024];
    mangled_buf[0] = '\0';

    void *sym = nullptr;
    const char *primary_method = "none";
    size_t nsyms = 0;

#if defined(__linux__)
    sym = resolve_via_elf(handle, func_name, mangled_buf, sizeof(mangled_buf), &nsyms);
    primary_method = "elf";
#elif defined(__APPLE__)
    sym = resolve_via_macho(handle, func_name, mangled_buf, sizeof(mangled_buf), &nsyms);
    primary_method = "macho";
#elif defined(_WIN32)
    sym = resolve_via_pe(handle, func_name, mangled_buf, sizeof(mangled_buf), &nsyms);
    primary_method = "pe";
#endif

    // Method string for logging.
    char method_str[32];
    snprintf(method_str, sizeof(method_str), "%s", primary_method);

    // ---- Store in cache ----
    cache_store(handle, func_name, sym);

    if (verbose) {
        const char *lib = get_library_path(handle);
        if (sym) {
            fprintf(stderr, "resolve_symbol: \"%s\" -> \"%s\" in %s [via %s, scanned %zu syms]\n", func_name, mangled_buf[0] ? mangled_buf : "(unknown mangled name)", lib, method_str, nsyms);
        } else {
            fprintf(stderr, "resolve_symbol: \"%s\" NOT FOUND in %s [scanned %zu syms via %s]\n", func_name, lib, nsyms, method_str);
        }
    }

    return sym;
}

} // namespace mplapack_resolver

#endif // MPLAPACK_SYMBOL_RESOLVER_H

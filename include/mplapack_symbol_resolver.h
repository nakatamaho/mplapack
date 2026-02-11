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
// Supported platforms:
//   - Linux:  ELF dynamic symbol table parsing via DT_HASH / DT_GNU_HASH
//   - macOS:  Mach-O symbol table parsing via LC_SYMTAB + dyld APIs
//   - MinGW:  PE export table parsing via IMAGE_EXPORT_DIRECTORY
//
// All platforms share a common nm-based fallback.
//
// Usage:
//   void *sym = mplapack_resolver::resolve_symbol(handle, "Rgemm");         // verbose (default)
//   void *sym = mplapack_resolver::resolve_symbol(handle, "Rgemm", false);  // silent
//
// Verbose output (stderr):
//   resolve_symbol: "Rgemm" -> "_Z5RgemmPKcS0_lllgPglS1_lgS1_l"
//                   in /usr/local/lib/libmpblas_dd.so [via elf]
//   resolve_symbol: "Raxpy" NOT FOUND in /usr/local/lib/libmpblas_dd.so

#ifndef MPLAPACK_SYMBOL_RESOLVER_H
#define MPLAPACK_SYMBOL_RESOLVER_H

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cstdint>

// ---- Platform headers ----
#if defined(_WIN32)
#include <windows.h>
#include <dbghelp.h>
#include <dlfcn.h>
#elif defined(__APPLE__)
#include <dlfcn.h>
#include <mach-o/dyld.h>
#include <mach-o/loader.h>
#include <mach-o/nlist.h>
#else
#include <dlfcn.h>
#include <elf.h>
#include <link.h>
#endif

#include <cxxabi.h>

namespace mplapack_resolver {

// ============================================================
//  Helper: get library file path from a dlopen handle
// ============================================================
static const char *get_library_path(void *handle) {
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
    HMODULE hmod = reinterpret_cast<HMODULE>(handle);
    GetModuleFileNameA(hmod, path_buf, sizeof(path_buf));
#endif

    if (!path_buf[0])
        return "(unknown library)";
    return path_buf;
}

// ============================================================
//  Helper: demangle and check if name matches "FuncName("
//  If matched, copies the original mangled name to out_mangled.
// ============================================================
static bool demangled_matches(const char *mangled, const char *target, size_t target_len, char *out_mangled, size_t out_mangled_size) {
    int status = -1;
    char *demangled = abi::__cxa_demangle(mangled, nullptr, nullptr, &status);
    if (status != 0 || !demangled) {
#if defined(_WIN32)
        // Try MSVC-style demangling (for MSVC-compiled DLLs loaded from MinGW)
        char buf[1024];
        if (UnDecorateSymbolName(mangled, buf, sizeof(buf), UNDNAME_COMPLETE | UNDNAME_NO_ACCESS_SPECIFIERS)) {
            bool match = (strncmp(buf, target, target_len) == 0);
            if (match && out_mangled && out_mangled_size > 0) {
                strncpy(out_mangled, mangled, out_mangled_size - 1);
                out_mangled[out_mangled_size - 1] = '\0';
            }
            free(demangled);
            return match;
        }
#endif
        free(demangled);
        return false;
    }
    bool match = (strncmp(demangled, target, target_len) == 0);
    if (match && out_mangled && out_mangled_size > 0) {
        strncpy(out_mangled, mangled, out_mangled_size - 1);
        out_mangled[out_mangled_size - 1] = '\0';
    }
    free(demangled);
    return match;
}

// ============================================================
//  Linux: ELF dynamic symbol table
// ============================================================
#if defined(__linux__)

static void *resolve_via_elf(void *handle, const char *func_name, char *out_mangled, size_t out_mangled_size) {
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

    // Determine symbol count
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

static void *resolve_via_macho(void *handle, const char *func_name, char *out_mangled, size_t out_mangled_size) {
    char target[256];
    snprintf(target, sizeof(target), "%s(", func_name);
    size_t target_len = strlen(target);

    uint32_t image_count = _dyld_image_count();

    for (uint32_t img = 0; img < image_count; ++img) {
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

            // Mach-O symbols have a leading underscore; skip it for demangling
            const char *mangled = name;
            if (name[0] == '_')
                mangled = name + 1;

            if (demangled_matches(mangled, target, target_len, out_mangled, out_mangled_size)) {
                void *sym = dlsym(handle, mangled);
                if (sym)
                    return sym;
                sym = dlsym(handle, name + 1);
                if (sym)
                    return sym;
            }
        }
    }
    return nullptr;
}

#endif // __APPLE__

// ============================================================
//  MinGW / Windows: PE export table parsing
// ============================================================
#if defined(_WIN32)

static void *resolve_via_pe(void *handle, const char *func_name, char *out_mangled, size_t out_mangled_size) {
    HMODULE hmod = reinterpret_cast<HMODULE>(handle);

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

    char target[256];
    snprintf(target, sizeof(target), "%s(", func_name);
    size_t target_len = strlen(target);

    for (DWORD i = 0; i < exports->NumberOfNames; ++i) {
        const char *mangled = reinterpret_cast<const char *>(base + names[i]);
        if (!mangled || !mangled[0])
            continue;

        if (demangled_matches(mangled, target, target_len, out_mangled, out_mangled_size)) {
            WORD ord = ordinals[i];
            DWORD func_rva = funcs[ord];

            // Check for forwarder (RVA falls within export directory range)
            if (func_rva >= export_rva && func_rva < export_rva + export_size) {
                void *sym = dlsym(handle, mangled);
                if (sym)
                    return sym;
                continue;
            }
            return reinterpret_cast<void *>(const_cast<uint8_t *>(base) + func_rva);
        }
    }
    return nullptr;
}

#endif // _WIN32

// ============================================================
//  Fallback: nm command (all platforms)
// ============================================================
static void *resolve_via_nm(void *handle, const char *func_name, char *out_mangled, size_t out_mangled_size) {
    const char *lib_path = nullptr;
    const char *nm_flags = "";

#if defined(__linux__)
    struct link_map *lm = nullptr;
    if (dlinfo(handle, RTLD_DI_LINKMAP, &lm) != 0 || !lm || !lm->l_name[0])
        return nullptr;
    lib_path = lm->l_name;
    nm_flags = "-D";
#elif defined(__APPLE__)
    static char macho_path[4096];
    macho_path[0] = '\0';
    uint32_t image_count = _dyld_image_count();
    for (uint32_t i = 0; i < image_count; ++i) {
        const char *name = _dyld_get_image_name(i);
        if (!name)
            continue;
        void *test_handle = dlopen(name, RTLD_LAZY | RTLD_NOLOAD);
        if (test_handle == handle) {
            strncpy(macho_path, name, sizeof(macho_path) - 1);
            macho_path[sizeof(macho_path) - 1] = '\0';
            dlclose(test_handle);
            break;
        }
        if (test_handle)
            dlclose(test_handle);
    }
    if (!macho_path[0])
        return nullptr;
    lib_path = macho_path;
    nm_flags = "-gU";
#elif defined(_WIN32)
    static char win_path[MAX_PATH];
    HMODULE hmod = reinterpret_cast<HMODULE>(handle);
    if (GetModuleFileNameA(hmod, win_path, sizeof(win_path)) == 0)
        return nullptr;
    lib_path = win_path;
    nm_flags = "";
#else
    return nullptr;
#endif

    if (!lib_path)
        return nullptr;

    char cmd[4096];
    snprintf(cmd, sizeof(cmd), "nm %s '%s' 2>/dev/null", nm_flags, lib_path);

    FILE *pipe = popen(cmd, "r");
    if (!pipe)
        return nullptr;

    char target[256];
    snprintf(target, sizeof(target), "%s(", func_name);
    size_t target_len = strlen(target);

    char line[1024];
    void *result = nullptr;
    while (fgets(line, sizeof(line), pipe)) {
        char *p = line;
        while (*p == ' ' || (*p >= '0' && *p <= '9') || (*p >= 'a' && *p <= 'f') || (*p >= 'A' && *p <= 'F'))
            p++;
        if (*p == ' ')
            p++;
        char sym_type = *p;
        if (sym_type == 'U' || sym_type == 'u' || sym_type == 'w')
            continue;
        p++;
        while (*p == ' ')
            p++;
        char *name_pos = p;

        size_t len = strlen(name_pos);
        while (len > 0 && (name_pos[len - 1] == '\n' || name_pos[len - 1] == '\r' || name_pos[len - 1] == ' '))
            name_pos[--len] = '\0';
        if (!name_pos[0])
            continue;

        const char *mangled = name_pos;
#if defined(__APPLE__)
        if (name_pos[0] == '_')
            mangled = name_pos + 1;
#endif

        if (demangled_matches(mangled, target, target_len, out_mangled, out_mangled_size)) {
            result = dlsym(handle, mangled);
#if defined(__APPLE__)
            if (!result && name_pos[0] == '_')
                result = dlsym(handle, name_pos + 1);
#endif
            if (result)
                break;
        }
    }
    pclose(pipe);
    return result;
}

// ============================================================
//  Main entry point
// ============================================================

// Resolve a C++ function by its unmangled base name (e.g. "Rgemm").
// Searches the dynamic/exported symbol table of the library opened
// via dlopen(). Returns the function pointer, or nullptr on failure.
//
// If verbose == true (default), prints diagnostic info to stderr:
//   resolve_symbol: "Rgemm" -> "_Z5RgemmPKcS0_lllgPglS1_lgS1_l"
//                   in /usr/local/lib/libmpblas_dd.so [via elf]
//   resolve_symbol: "Raxpy" NOT FOUND in /usr/local/lib/libmpblas_dd.so
//
// Resolution order:
//   1. Native table parsing (ELF / Mach-O / PE)
//   2. nm command fallback
static void *resolve_symbol(void *handle, const char *func_name, bool verbose = true) {
    if (!handle || !func_name)
        return nullptr;

    char mangled_buf[1024];
    mangled_buf[0] = '\0';

    void *sym = nullptr;
    const char *method = "none";

#if defined(__linux__)
    sym = resolve_via_elf(handle, func_name, mangled_buf, sizeof(mangled_buf));
    if (sym)
        method = "elf";
#elif defined(__APPLE__)
    sym = resolve_via_macho(handle, func_name, mangled_buf, sizeof(mangled_buf));
    if (sym)
        method = "macho";
#elif defined(_WIN32)
    sym = resolve_via_pe(handle, func_name, mangled_buf, sizeof(mangled_buf));
    if (sym)
        method = "pe";
#endif

    if (!sym) {
        mangled_buf[0] = '\0';
        sym = resolve_via_nm(handle, func_name, mangled_buf, sizeof(mangled_buf));
        if (sym)
            method = "nm";
    }

    if (verbose) {
        const char *lib = get_library_path(handle);
        if (sym) {
            fprintf(stderr, "resolve_symbol: \"%s\" -> \"%s\" in %s [via %s]\n", func_name, mangled_buf[0] ? mangled_buf : "(unknown mangled name)", lib, method);
        } else {
            fprintf(stderr, "resolve_symbol: \"%s\" NOT FOUND in %s\n", func_name, lib);
        }
    }

    return sym;
}

} // namespace mplapack_resolver

#endif // MPLAPACK_SYMBOL_RESOLVER_H

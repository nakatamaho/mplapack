/* Load exactly one library without linking or preloading its dependencies. */
#include <dlfcn.h>
#include <stdio.h>

int main(int argc, char **argv) {
    if (argc != 2) {
        fprintf(stderr, "usage: %s <shared library>\n", argv[0]);
        return 2;
    }
    void *handle = dlopen(argv[1], RTLD_NOW | RTLD_LOCAL);
    if (!handle) {
        fprintf(stderr, "FAIL: unable to load %s: %s\n", argv[1], dlerror());
        return 1;
    }
    printf("PASS: clean-process load of %s\n", argv[1]);
    dlclose(handle);
    return 0;
}

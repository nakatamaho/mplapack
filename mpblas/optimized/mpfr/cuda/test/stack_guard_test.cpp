// No device: deterministically exercise otherwise hard-to-inject API errors.
#include "../mpfr_cuda_stack.h"
#include <cstdio>
#include <initializer_list>
using namespace mplapack_mpfr_cuda;
static size_t limit=1024;
static int sets=0,clears=0;
static bool fail_set=false,fail_get=false,pending=false;
cudaError_t cudaDeviceGetLimit(size_t *p,cudaLimit) {
    if(fail_get) return cudaErrorInvalidValue;
    *p=limit;return cudaSuccess;
}
cudaError_t cudaDeviceSetLimit(cudaLimit,size_t n) {
    ++sets;if(fail_set) {pending=true;return cudaErrorInvalidValue;}
    limit=n;return cudaSuccess;
}
cudaError_t cudaGetLastError() {++clears;pending=false;return cudaSuccess;}
static void check(bool ok,const char *what) {
    if(!ok) {std::fprintf(stderr,"FAIL stack guard: %s\n",what);std::exit(1);}
}
static bool early_return() {mpfr_stack_guard g(8192);return false;}
int main() {
    unsetenv("MPLAPACK_MPFR_CUDA_STACK_RESTORE");
    {mpfr_stack_guard g(4096);check(!g.error && g.changed,"raise");}
    check(limit==4096 && sets==1,"default retains raised limit");
    for(const char *v:{"0","other","11","1 "}) {
        setenv("MPLAPACK_MPFR_CUDA_STACK_RESTORE",v,1);limit=1024;
        {mpfr_stack_guard g(4096);}
        check(limit==4096,"only exact 1 restores");
    }
    setenv("MPLAPACK_MPFR_CUDA_STACK_RESTORE","1",1);limit=1024;
    {mpfr_stack_guard g(4096);check(!g.restore() && !g.restore(),"explicit restore idempotent");}
    check(limit==1024,"restore original");
    limit=7777;int before=sets;
    {mpfr_stack_guard g(4096);check(!g.changed,"no raise");}
    check(limit==7777 && sets==before,"no change means no set or restore");
    limit=1024;check(!early_return() && limit==1024,"early return destructor restores");
    {mpfr_stack_guard outer(4096);{mpfr_stack_guard inner(8192);}check(limit==4096,"nested original");}
    check(limit==1024,"outer original");
    {mpfr_stack_guard g(4096);unsetenv("MPLAPACK_MPFR_CUDA_STACK_RESTORE");}
    check(limit==1024,"policy sampled at call entry");
    {mpfr_stack_guard g(4096);setenv("MPLAPACK_MPFR_CUDA_STACK_RESTORE","1",1);}
    check(limit==4096,"default call not retroactively restored");limit=1024;
    {mpfr_stack_guard g(4096);fail_set=true;check(!g.restore(),"restore failure cannot change computation status");}
    check(clears==1 && !pending,"failed restore clears CUDA error once");fail_set=false;
    fail_get=true;before=sets;{mpfr_stack_guard g(8192);check(g.error && !g.changed,"get failure");}
    check(sets==before,"get failure does not set/restore");fail_get=false;
    limit=1024;fail_set=true;{mpfr_stack_guard g(8192);check(g.error && !g.changed,"raise failure");}
    fail_set=false;cudaGetLastError();
    std::puts("stack guard host: all policies, nesting, early return, restore error cleanup PASS");
}

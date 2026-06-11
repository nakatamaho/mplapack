#include <assert.h>
#include <malloc.h>
#include <mpblas__Float128.h>
#include "OclHelper.hpp"

static __CPP::__CPP_Helper *p;
static cl_mem b_x, b_y, b_z;
static cl_command_queue q;
static cl_kernel ker, ker_nn, ker_tn, ker_nt, ker_tt;
static size_t buffer_size_a_byte_global;
static size_t buffer_size_b_byte_global;
static size_t buffer_size_c_byte_global;
//static _Float128 *AB;
static size_t bs_opencl = 1;
static int __firstcall__ = 1;

void init_gpu(int bs = 4)
{
  char options[128];
  bs_opencl = bs;
  sprintf(options, "-D TS=%i", bs);
  
  p = new __CPP::__CPP_Helper();
  
  cl_context ctx = p->getctx();
  q = p->getq();
  p->build(options);

  ker_nn = p->getkernel("gemm_NN_X");
  ker_tn = p->getkernel("gemm_TN_X");
  ker_nt = p->getkernel("gemm_NT_X");
  ker_tt = p->getkernel("gemm_TT_X");

  buffer_size_a_byte_global = 8*1024*1024*sizeof(_Float128);
  buffer_size_b_byte_global = 8*1024*1024*sizeof(_Float128);
  buffer_size_c_byte_global = 8*1024*1024*sizeof(_Float128);
 
  cl_int result;
  b_x  = clCreateBuffer(ctx, CL_MEM_READ_ONLY,  buffer_size_a_byte_global, NULL, &result);
  assert(result == CL_SUCCESS);
  b_y  = clCreateBuffer(ctx, CL_MEM_READ_ONLY,  buffer_size_b_byte_global, NULL, &result);
  assert(result == CL_SUCCESS);
  b_z  = clCreateBuffer(ctx, CL_MEM_READ_WRITE, buffer_size_c_byte_global, NULL, &result);
  assert(result == CL_SUCCESS);

  //  AB = (_Float128 *)memalign(64, buffer_size_c_byte_global);
  //  assert(AB != NULL);
  __firstcall__ = 0;

  std::cout << "block size for OpenCL kernels = " << bs_opencl << "\n";
}

void Rgemm_internal(bool btransa, bool btransb,
		    int m, int n, int k,
		    _Float128 alpha,
		    _Float128 *A, int lda,
		    _Float128 *B, int ldb,
		    _Float128 beta,
		    _Float128 *C, int ldc)
{
  cl_event event;

  if (__firstcall__ == 1) {
    //    std::cout << ">>>First call (gpu_init) \n";
    char* bs_str = getenv("RGEMM_OPENCL_BLOCKSIZE");
    if (bs_str == NULL) {
      init_gpu(8);
    } else {
      init_gpu(atoi(bs_str));
    } 
  }else{
  }

  //  std::cerr << "called Rgemm for OpenCL\n";
  
  // Matrix A : m * k  
  // Matrix B : k * n
  // Matrix C : m * n

  // Column major layout
  //
  // if A is transposed     : shape k-by-m, buffer lda*m
  // if A is not transposed : shape m-by-k, buffer lda*k
  //
  // if B is transposed     : shape n-by-k, buffer ldb*k
  // if B is not transposed : shape k-by-n, buffer ldb*n
  //
  // C : shape m-by-n, buffer ldc*n
  //
  //
  //  printf("RRR %i %i : %i %i %i\n", btransa, btransb, m, n, k);

  size_t buffer_size_a_byte;
  size_t buffer_size_b_byte;
  size_t buffer_size_c_byte;
  const size_t w_byte = sizeof(_Float128);
  
  if (btransa) buffer_size_a_byte = lda*m*w_byte;
  else         buffer_size_a_byte = lda*k*w_byte;

  if (btransb) buffer_size_b_byte = ldb*k*w_byte;
  else         buffer_size_b_byte = ldb*n*w_byte;
  
  buffer_size_c_byte = ldc*n*w_byte;

  cl_int result;
  cl_context ctx = p->getctx();

  //  printf("%i %i\n", buffer_size_a_byte , buffer_size_a_byte_global);

  if (buffer_size_a_byte > buffer_size_a_byte_global) {
    clReleaseMemObject(b_x);
    buffer_size_a_byte_global = buffer_size_a_byte;
    b_x  = clCreateBuffer(ctx, CL_MEM_READ_ONLY, buffer_size_a_byte_global, NULL, &result);
    assert(result == CL_SUCCESS);
  }
  if (buffer_size_b_byte > buffer_size_b_byte_global) {
    clReleaseMemObject(b_y);
    buffer_size_b_byte_global = buffer_size_b_byte;
    b_y  = clCreateBuffer(ctx, CL_MEM_READ_ONLY,  buffer_size_b_byte_global, NULL, &result);
    assert(result == CL_SUCCESS);
  }
  if (buffer_size_c_byte > buffer_size_c_byte_global) {
    clReleaseMemObject(b_z);
    buffer_size_c_byte_global = buffer_size_c_byte;
    b_z  = clCreateBuffer(ctx, CL_MEM_READ_WRITE,  buffer_size_c_byte_global, NULL, &result);
    assert(result == CL_SUCCESS);
    //    free(AB);
    //    AB = (_Float128 *)memalign(64, buffer_size_c_byte_global);
    //    assert(AB != NULL);
  }

  if (btransa == true  && btransb == true ) ker = ker_tt;
  if (btransa == false && btransb == true ) ker = ker_nt;
  if (btransa == true  && btransb == false) ker = ker_tn;
  if (btransa == false && btransb == false) ker = ker_nn;

  result = clSetKernelArg(ker, 0, sizeof(cl_int), &m);        assert(result == CL_SUCCESS);
  result = clSetKernelArg(ker, 1, sizeof(cl_int), &n);        assert(result == CL_SUCCESS);
  result = clSetKernelArg(ker, 2, sizeof(cl_int), &k);        assert(result == CL_SUCCESS);
  result = clSetKernelArg(ker, 3, sizeof(_Float128), &alpha); assert(result == CL_SUCCESS);
  result = clSetKernelArg(ker, 4, sizeof(cl_mem), &b_x);      assert(result == CL_SUCCESS);
  result = clSetKernelArg(ker, 5, sizeof(cl_int), &lda);      assert(result == CL_SUCCESS);
  result = clSetKernelArg(ker, 6, sizeof(cl_mem), &b_y);      assert(result == CL_SUCCESS);
  result = clSetKernelArg(ker, 7, sizeof(cl_int), &ldb);      assert(result == CL_SUCCESS);
  result = clSetKernelArg(ker, 8, sizeof(_Float128), &beta);  assert(result == CL_SUCCESS);
  result = clSetKernelArg(ker, 9, sizeof(cl_mem), &b_z);      assert(result == CL_SUCCESS);
  result = clSetKernelArg(ker, 10, sizeof(cl_int), &ldc);     assert(result == CL_SUCCESS);

#define roundup(n,nt) ((n) % (nt) == 0 ? (n) : ((n)/(nt) + 1)*(nt))
  cl_event e;
  const size_t lt[2] = {bs_opencl, bs_opencl};
  const size_t ni = roundup(m, lt[0]);
  const size_t nj = roundup(n, lt[1]);
  const size_t nt[2] = {ni, nj};

  clEnqueueWriteBuffer(q, b_x, CL_TRUE, 0, buffer_size_a_byte, A,  0, NULL, &event);
  clEnqueueWriteBuffer(q, b_y, CL_TRUE, 0, buffer_size_b_byte, B,  0, NULL, &event);
  clEnqueueWriteBuffer(q, b_z, CL_TRUE, 0, buffer_size_c_byte, C,  0, NULL, &event);
  
  result = clEnqueueNDRangeKernel(q, ker, 2, NULL, nt, lt, 0, NULL, &e);
  assert(result == CL_SUCCESS);
  clFinish(q);
  
  clEnqueueReadBuffer(q, b_z, CL_TRUE, 0, buffer_size_c_byte, C, 0, NULL, &event);  
  clFinish(q);
}

void Rgemm_ref(const char *transa, const char *transb,
	       mplapackint const m, mplapackint const n, mplapackint const k,
	       _Float128 const alpha,
	       _Float128 * A, mplapackint const lda,
	       _Float128 * B, mplapackint const ldb,
	       _Float128 const beta,
	       _Float128 * C, mplapackint const ldc);

void Rgemm(const char *transa, const char *transb,
	   mplapackint const m, mplapackint const n, mplapackint const k,
	   _Float128 const alpha,
	   _Float128 * A, mplapackint const lda,
	   _Float128 * B, mplapackint const ldb,
	   _Float128 const beta,
	   _Float128 * C, mplapackint const ldc)
{
  bool b_transa, b_transb;
  const _Float128 Zero = 0.0, One = 1.0;
  
  //Quick return if possible.
  if ((m == 0) || (n == 0) || (((alpha == Zero) || (k == 0)) && (beta == One)))
    return;

  //And when alpha == 0.0
  if (alpha == Zero) {
    if (beta == Zero) {
      for (int j = 0; j < n; j++) {
	for (int i = 0; i < m; i++) {
	  C[i + j * ldc] = Zero;
	}
      }
    } else {
      for (int j = 0; j < n; j++) {
	for (int i = 0; i < m; i++) {
	  C[i + j * ldc] = beta * C[i + j * ldc];
	}
      }
    }
    return;
  }

  if ( (m <= bs_opencl) || (n <= bs_opencl) ) {
    //    std::cout << "call Rgemm for small matrices\n";
    Rgemm_ref(transa, transb, m, n, k, alpha, A, lda, B, ldb, beta, C, ldc);
    return;
  }

  b_transa = false;
  b_transb = false;

  if (transa[0] == 'T' || transa[0] == 't') b_transa = true;
  if (transb[0] == 'T' || transb[0] == 't') b_transb = true;

  Rgemm_internal(b_transa, b_transb, m, n, k, alpha, A, lda, B, ldb, beta, C, ldc);
          
  return;
}

#ifdef PEZYCL
#include <pzcl/pzcl_ocl_wrapper.h>
#define cl_mem_flags pzcl_mem_flags
#include "base64.h"
#else
#ifdef __XILLINX__
#include "xcl2.hpp"
#else
#define __CL_ENABLE_EXCEPTIONS
#define CL_USE_DEPRECATED_OPENCL_1_1_APIS
#ifdef __APPLE__
#include <OpenCL/opencl.h>
#else
#include <CL/cl.h>
#endif
#endif
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <vector>
#include <iostream>
#include <assert.h>

namespace __CPP {
#ifdef PEZYCL
#include "cpp_kernel_pz.h"
#else

#ifdef __XILLINX__
#else
#ifdef __AOCL__
#else
#include "cpp_kernel.h"
#endif
#endif
#endif
  class __CPP_Helper0 {
  public:
    __CPP_Helper0() { clinfo(); }; 
    ~__CPP_Helper0();

    void clinfo();
    //    void build();
    //    void build(std::string);
    void build(const char []);
    void setup(unsigned int ip, unsigned int id);      
    void setup();
    cl_kernel getkernel(const char *kernel_name);

  public:
    cl_context ctx;
    cl_command_queue q;
    std::vector<std::pair<int,int> > list;

  private:
    cl_program program;
    cl_platform_id platform_id[16];
    cl_kernel ker;
    cl_uint npl;
    cl_device_id device_id[16];
    cl_device_id dev;
    char pname[128];
    char dname[128];
    char pver[128];
    void dumperror();
  };

  __CPP_Helper0::~__CPP_Helper0 () {
    clReleaseKernel(ker);
    clReleaseProgram(program);
    clReleaseCommandQueue(q);
    clReleaseContext(ctx);
  }
  
  cl_kernel __CPP_Helper0::getkernel(const char *kernel_name)
  {
    cl_int status = CL_SUCCESS;
    cl_kernel res = clCreateKernel(program, kernel_name, &status);
    ker = res;

    if (status != CL_SUCCESS) {  
      fprintf(stderr, "Create kernel %s failed %d\n", kernel_name, status);
      exit(-1);
    }

    // https://arrayfire.com/blog/generating-ptx-files-from-opencl-code/
    {
      // Query binary (PTX file) size
      size_t bin_sz;
      cl_int err = clGetProgramInfo(program, CL_PROGRAM_BINARY_SIZES, sizeof(size_t), &bin_sz, NULL);

      // Read binary (PTX file) to memory buffer
      unsigned char *bin = (unsigned char *)malloc(bin_sz);
      err = clGetProgramInfo(program, CL_PROGRAM_BINARIES, sizeof(unsigned char *), &bin, NULL);

      // Save PTX
      char buf[1024];
      sprintf(buf, "%s_ocl.ptx", kernel_name);  
      FILE *fp = fopen(buf, "wb");
      fwrite(bin, sizeof(char), bin_sz, fp);
      fclose(fp);
      free(bin);
    }

    return res;
  }

  void __CPP_Helper0::setup()
  {
    char* p_type_str = getenv("RGEMM_OPENCL_PLATFORM");
    char* d_type_str = getenv("RGEMM_OPENCL_DEVICE");

    if (p_type_str == NULL || d_type_str == NULL) {
      setup(0, 0);
    } else {
      setup(atoi(p_type_str), atoi(d_type_str));
    }
  }

  void __CPP_Helper0::setup(unsigned int ip, unsigned int id)
  {
    cl_int result = CL_SUCCESS;
    std::cerr << "Selected: "; 
    if (npl <= ip) {
      fprintf(stderr, "FATAL: the specifed platform does not exist");
      exit(-1);      
    }
    cl_uint ndev = 0;
    if ((result = clGetDeviceIDs(platform_id[ip], CL_DEVICE_TYPE_ALL, 16, device_id, &ndev)) != CL_SUCCESS) {
      fprintf(stderr, "clGetDeviceIDs() failed : %d\n", result);
      exit(-1);
    }

    if (ndev <= id) {
      fprintf(stderr, "FATAL: the specifed device does not exist");
      exit(-1);      
    }

    clGetPlatformInfo(platform_id[ip],  CL_PLATFORM_NAME, sizeof(pname), pname, NULL);
    clGetPlatformInfo(platform_id[ip],  CL_PLATFORM_VERSION, sizeof(pver), pver, NULL);
    clGetDeviceInfo(device_id[id], CL_DEVICE_NAME, sizeof(dname), dname, NULL);

    std::cerr << pname << " " << pver << "::"	<< dname << "\n";

    // Create Context
    if( (ctx = clCreateContext(NULL, 1, &device_id[id], NULL, NULL, &result)) == NULL ) {
      fprintf(stderr, "Create context failed %d : dev %i\n", result, id);
      exit(-1);
    }

    // use the first device only
    cl_int status = CL_SUCCESS;
    q = clCreateCommandQueue(ctx, device_id[id], CL_QUEUE_PROFILING_ENABLE, &status);
    if (status != CL_SUCCESS) {  
      fprintf(stderr, "Create commandq failed %d\n", status);
      exit(-1);
    }
    dev = device_id[id];
  }

#ifdef __AOCL__
  // from AOCLUtils
  unsigned char *loadBinaryFile(const char *file_name, size_t *size) {
    fprintf(stderr, "loading %s\n", file_name);

    // Open the File
    FILE* fp;
    fp = fopen(file_name, "rb");
    if(fp == 0) {
      *size = 0;
      return NULL;
    }

    // Get the size of the file
    fseek(fp, 0, SEEK_END);
    *size = ftell(fp);
    //std::cout << *size << "\n";

    // Allocate space for the binary
    unsigned char *binary = new unsigned char[*size];

    // Go back to the file start
    rewind(fp);

    // Read the file into the binary
    if(fread((void*)binary, *size, 1, fp) == 0) {
      delete[] binary;
      fclose(fp);
      return NULL;
    }
    
    return binary;
  }
#endif
  void __CPP_Helper0::build(const char options[] = NULL)
  {
    cl_int status = CL_SUCCESS;
#ifdef PEZYCL
    {
      // Create program object
      std::vector<cl_device_id> device_id_lists;
      device_id_lists.push_back( dev );

      int blen;
      const unsigned char *b = unbase64(kernel_pz_str, sizeof(kernel_pz_str), &blen);
      blen--;
      std::cerr << "GOOSE: unbase len " << blen << " byte \n";
      {
	const unsigned char* listBin[1];
	listBin[0] = (unsigned char*)b;
	cl_int binary_status = CL_SUCCESS;
	size_t length = blen;
	cl_int result;

	program = clCreateProgramWithBinary( ctx, (cl_uint)device_id_lists.size(), &device_id_lists[0],
					  &length, listBin, &binary_status, &result);
	if(program == NULL)
	  {
	    printf("clCreateProgramWithBinary failed, %d\n", result);
	    exit(-1);
	  }
      }
    }
#else
#ifdef __XILLINX__
    std::string device_name(dname);

    std::string binaryFile = xcl::find_binary_file(device_name, kernel_name);
    cl::Program::Binaries bins = xcl::import_binary_file(binaryFile);
    //    devices.resize(1);
  
    //    std::vector<cl::Device> devices;
    //    devices.push_back(dev);
    //    cl::Program program(ctx, devices, bins);  

    cl_int binary_status = CL_SUCCESS;

    //    const unsigned char* listBin[1];
    //    const unsigned char *b = (const unsigned char *)
    //    listBin[0] = (unsigned char**)&bins[0].first;
    size_t length = bins[0].second;
    cl_int result;

    std::cout << length << "\n";

    std::vector<cl_device_id> device_id_lists;
    device_id_lists.push_back( dev );
    program = clCreateProgramWithBinary(ctx, (cl_uint)device_id_lists.size(), &device_id_lists[0],
					&length, (const unsigned char**)&bins[0].first, 
					&binary_status, &result);
    
    if(program == NULL) {
      printf("clCreateProgramWithBinary failed, %d\n", result);
      exit(-1);
    }
#else
#ifdef __AOCL__
    size_t binary_size;
    unsigned char *binary = loadBinaryFile(kernel_name, &binary_size);
    program = clCreateProgramWithBinary(ctx, 1, device_id, &binary_size, 
					(const unsigned char **)&binary, NULL, &status);
    if(program == NULL) {
      printf("clCreateProgramWithBinary failed, %d\n", status);
      exit(-1);
    }
#else
    {
      char *prog;
      prog = (char *)malloc(sizeof(__CPP::cpp_kernel_str));
      size_t ss[1];
      ss[0] = sizeof(__CPP::cpp_kernel_str);
      strcpy(prog, __CPP::cpp_kernel_str);
      program = clCreateProgramWithSource(ctx, 1, (const char **)&prog, ss, &status);
      if (status != CL_SUCCESS) {  
	fprintf(stderr, "cl create program failed %d\n", status);
	exit(-1);
      }
      free(prog);
    }
#endif
#endif
#endif

    status = clBuildProgram(program, 1, &dev, options, NULL, NULL);
    if(status != CL_SUCCESS) {
      fprintf(stderr, "build failed\n");
      dumperror();
      exit(-1);
    }
  }

  void __CPP_Helper0::dumperror()
  {
#ifdef PEZYCL
#else
    cl_int logStatus;
    char * buildLog = NULL;
    size_t buildLogSize = 0;

    logStatus = clGetProgramBuildInfo (program,
                                       dev,
                                       CL_PROGRAM_BUILD_LOG,
                                       buildLogSize,
                                       buildLog,
                                       &buildLogSize);

    buildLog = (char*)malloc(buildLogSize);
    memset(buildLog, 0, buildLogSize);

    logStatus = clGetProgramBuildInfo (program,
                                       dev,
                                       CL_PROGRAM_BUILD_LOG,
                                       buildLogSize,
                                       buildLog,
                                       NULL);

    fprintf(stderr, "%s\n", buildLog);
    std::cout << logStatus << "\n";

    free(buildLog);
#endif
  }

  void __CPP_Helper0::clinfo() 
  {
    cl_int result = CL_SUCCESS;

    if ((result = clGetPlatformIDs(16, platform_id, &npl)) != CL_SUCCESS) {
      fprintf(stderr, "clGetPlatformIDs() failed : %dn", result);
      exit(-1);
    }

    for(unsigned int i = 0; i < npl; i++) {
      clGetPlatformInfo(platform_id[i],  CL_PLATFORM_NAME, sizeof(pname), pname, NULL);
      clGetPlatformInfo(platform_id[i],  CL_PLATFORM_VERSION, sizeof(pver), pver, NULL);
      fprintf(stderr, "platform %d %s %s\n", i, pname, pver);

      // Get Device ID
      cl_uint ndev = 0;
      if ((result = clGetDeviceIDs(platform_id[i], CL_DEVICE_TYPE_ALL, 16, device_id, &ndev)) != CL_SUCCESS) {
	fprintf(stderr, "clGetDeviceIDs() failed : %d\n", result);
	exit(-1);
      }

      for(unsigned int j = 0; j < ndev; j++) {
	clGetDeviceInfo(device_id[j], CL_DEVICE_NAME, sizeof(dname), dname, NULL);
	fprintf(stderr, "\tdevice %d %s\n", j, dname);
	list.push_back(std::make_pair(i, j));
      }
    }
  }

  static bool cl_first       = true;
  static bool cl_first_build = true;
  static __CPP_Helper0 *acc0;

  class __CPP_Helper {
  public:
    __CPP_Helper(int ip = 0, int id = 0) {
      if (cl_first) {
	acc0 = new __CPP_Helper0;
	acc0->setup();
	cl_first = false;
      }
    }
    ~__CPP_Helper() {
      delete acc0;
    };

    void build() { 
      if (cl_first_build) {
	acc0->build();
	cl_first_build = false;
      }
    }

    void build(const char options[]) { 
      if (cl_first_build) {
	acc0->build(options);
	cl_first_build = false;
      }
    }

    cl_context getctx() {
      return acc0->ctx;
    }

    cl_command_queue getq() {
      return acc0->q;
    }

    cl_kernel getkernel(const char *kernel_name) {
      return acc0->getkernel(kernel_name);
    }
  };
}

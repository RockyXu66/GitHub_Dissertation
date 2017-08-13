/***** GCL Generated File *********************/
/* Automatically generated file, do not edit! */
/**********************************************/

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <dispatch/dispatch.h>
#include <OpenCL/opencl.h>
#include <OpenCL/gcl_priv.h>
#include "mykernel.cl.h"

static void initBlocks(void);

// Initialize static data structures
static block_kernel_pair pair_map[2] = {
      { NULL, NULL },
      { NULL, NULL }
};

static block_kernel_map bmap = { 0, 2, initBlocks, pair_map };

// Block function
void (^addVector_kernel)(const cl_ndrange *ndrange, cl_float* output, cl_float* score, cl_float* eigenvector) =
^(const cl_ndrange *ndrange, cl_float* output, cl_float* score, cl_float* eigenvector) {
  int err = 0;
  cl_kernel k = bmap.map[0].kernel;
  if (!k) {
    initBlocks();
    k = bmap.map[0].kernel;
  }
  if (!k)
    gcl_log_fatal("kernel addVector does not exist for device");
  kargs_struct kargs;
  gclCreateArgsAPPLE(k, &kargs);
  err |= gclSetKernelArgMemAPPLE(k, 0, output, &kargs);
  err |= gclSetKernelArgMemAPPLE(k, 1, score, &kargs);
  err |= gclSetKernelArgMemAPPLE(k, 2, eigenvector, &kargs);
  gcl_log_cl_fatal(err, "setting argument for addVector failed");
  err = gclExecKernelAPPLE(k, ndrange, &kargs);
  gcl_log_cl_fatal(err, "Executing addVector failed");
  gclDeleteArgsAPPLE(k, &kargs);
};

void (^addMean_kernel)(const cl_ndrange *ndrange, cl_float* output, cl_float* mean) =
^(const cl_ndrange *ndrange, cl_float* output, cl_float* mean) {
  int err = 0;
  cl_kernel k = bmap.map[1].kernel;
  if (!k) {
    initBlocks();
    k = bmap.map[1].kernel;
  }
  if (!k)
    gcl_log_fatal("kernel addMean does not exist for device");
  kargs_struct kargs;
  gclCreateArgsAPPLE(k, &kargs);
  err |= gclSetKernelArgMemAPPLE(k, 0, output, &kargs);
  err |= gclSetKernelArgMemAPPLE(k, 1, mean, &kargs);
  gcl_log_cl_fatal(err, "setting argument for addMean failed");
  err = gclExecKernelAPPLE(k, ndrange, &kargs);
  gcl_log_cl_fatal(err, "Executing addMean failed");
  gclDeleteArgsAPPLE(k, &kargs);
};

// Initialization functions
static void initBlocks(void) {
  const char* build_opts = "";
  static dispatch_once_t once;
  dispatch_once(&once,
    ^{ int err = gclBuildProgramBinaryAPPLE("OpenCL/mykernel.cl", "", &bmap, build_opts);
       if (!err) {
          assert(bmap.map[0].block_ptr == addVector_kernel && "mismatch block");
          bmap.map[0].kernel = clCreateKernel(bmap.program, "addVector", &err);
          assert(bmap.map[1].block_ptr == addMean_kernel && "mismatch block");
          bmap.map[1].kernel = clCreateKernel(bmap.program, "addMean", &err);
       }
     });
}

__attribute__((constructor))
static void RegisterMap(void) {
  gclRegisterBlockKernelMap(&bmap);
  bmap.map[0].block_ptr = addVector_kernel;
  bmap.map[1].block_ptr = addMean_kernel;
}


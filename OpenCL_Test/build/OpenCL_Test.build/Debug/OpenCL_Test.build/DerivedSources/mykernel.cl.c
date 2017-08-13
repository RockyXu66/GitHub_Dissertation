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
void (^square_kernel)(const cl_ndrange *ndrange, cl_float* input1, cl_float* input2, cl_float* output) =
^(const cl_ndrange *ndrange, cl_float* input1, cl_float* input2, cl_float* output) {
  int err = 0;
  cl_kernel k = bmap.map[0].kernel;
  if (!k) {
    initBlocks();
    k = bmap.map[0].kernel;
  }
  if (!k)
    gcl_log_fatal("kernel square does not exist for device");
  kargs_struct kargs;
  gclCreateArgsAPPLE(k, &kargs);
  err |= gclSetKernelArgMemAPPLE(k, 0, input1, &kargs);
  err |= gclSetKernelArgMemAPPLE(k, 1, input2, &kargs);
  err |= gclSetKernelArgMemAPPLE(k, 2, output, &kargs);
  gcl_log_cl_fatal(err, "setting argument for square failed");
  err = gclExecKernelAPPLE(k, ndrange, &kargs);
  gcl_log_cl_fatal(err, "Executing square failed");
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
          assert(bmap.map[0].block_ptr == square_kernel && "mismatch block");
          bmap.map[0].kernel = clCreateKernel(bmap.program, "square", &err);
          assert(bmap.map[1].block_ptr == addMean_kernel && "mismatch block");
          bmap.map[1].kernel = clCreateKernel(bmap.program, "addMean", &err);
       }
     });
}

__attribute__((constructor))
static void RegisterMap(void) {
  gclRegisterBlockKernelMap(&bmap);
  bmap.map[0].block_ptr = square_kernel;
  bmap.map[1].block_ptr = addMean_kernel;
}


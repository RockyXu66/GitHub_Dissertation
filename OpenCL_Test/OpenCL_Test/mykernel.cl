// Simple OpenCL kernel that squares an input array.
// This code is stored in a file called mykernel.cl.
// You can name your kernel file as you would name any other
// file.  Use .cl as the file extension for all kernel
// source files.

// Kernel block.                                      //   1
kernel void square(                                   //   2
                   global float* input1,               //   3
                   global float* input2,
                   global float* output)
{
    size_t i = get_global_id(0);
    output[i] = output[i] + input1[i] * input2[0];
}


kernel void addMean(global float* output, global float* mean){
    size_t i = get_global_id(0);
    output[i] = output[i] + mean[i];
}

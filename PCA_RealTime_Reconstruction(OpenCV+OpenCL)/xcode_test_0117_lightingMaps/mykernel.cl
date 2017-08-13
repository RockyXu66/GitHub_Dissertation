// Kernel block.


kernel void addVector(global float* output, global float* score, global float* eigenvector){
    size_t i = get_global_id(0);
    output[i] = 1.0;//output[i] + eigenvector[i] * score[0];
}




kernel void addMean(global float* output, global float* mean){
    size_t i = get_global_id(0);
    output[i] = output[i] + mean[i];
}

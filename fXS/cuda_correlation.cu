#include <sastbx/fXS/cuda_correlation.cuh>

namespace sastbx {
namespace fXS {

  const int blocks_per_grid = 128;
  const int threads_per_block = 1024;
  const int chunk_size = blocks_per_grid*threads_per_block;

  /* ==========================================================================
     Basic cumulative addition of complex arrays using streams
  */
  __global__ void add_images_kernel_streams
  (const cuDoubleComplex* image_chunk, const int begin, const int image_size,
   cuDoubleComplex* summed_image) {
    int chunk_index = blockDim.x * blockIdx.x + threadIdx.x;
    int image_index = begin + chunk_index;
    if (image_index < image_size) {
      summed_image[image_index] = cuCadd(summed_image[image_index],
                                         image_chunk[chunk_index]);
    }
  }

  scitbx::af::shared<std::complex<double> > cuda_add_images_streams
    (const scitbx::af::const_ref<std::complex<double> >& images,
     const int& image_size, const int& n_images, const int& gpu_id) {

    // start GPU
    cudaSafeCall( cudaSetDevice(gpu_id) );

    // initialize timer
    // cudaEvent_t start, stop;
    // cudaSafeCall( cudaEventCreate(&start) );
    // cudaSafeCall( cudaEventCreate(&stop) );
    // cudaSafeCall( cudaEventRecord(start) );

    // allocate and initialize arrays
    cuDoubleComplex * h_images, * h_summed_image;
    cuDoubleComplex * d_images0, * d_images1, * d_summed_image;
    cudaSafeCall( cudaHostAlloc( (void**)&h_images,
                                 images.size() * sizeof(cuDoubleComplex),
                                 cudaHostAllocDefault ) );
    cudaSafeCall( cudaHostAlloc( (void**)&h_summed_image,
                                 image_size * sizeof(cuDoubleComplex),
                                 cudaHostAllocDefault ) );
    cudaSafeCall( cudaMalloc( (void**)&d_images0,
                              chunk_size * sizeof(cuDoubleComplex) ) );
    cudaSafeCall( cudaMalloc( (void**)&d_images1,
                              chunk_size * sizeof(cuDoubleComplex) ) );
    cudaSafeCall( cudaMalloc( (void**)&d_summed_image,
                              image_size * sizeof(cuDoubleComplex) ) );
    for (int i=0; i<images.size(); i++) {
      h_images[i] = make_cuDoubleComplex(images[i].real(),images[i].imag());
    }
    for (int i=0; i<image_size; i++) {
      h_summed_image[i] = make_cuDoubleComplex(0.0,0.0);
    }
    cudaSafeCall( cudaMemcpy( d_summed_image, h_summed_image,
                              image_size * sizeof(cuDoubleComplex),
                              cudaMemcpyHostToDevice ) );

    // initialize streams
    cudaStream_t s0, s1;
    cudaSafeCall( cudaStreamCreate( &s0 ) );
    cudaSafeCall( cudaStreamCreate( &s1 ) );

    int begin0, begin1, image_offset;
    for (int i=0; i<n_images; i++) {
      image_offset = i*image_size;
      for (int j=0; j<int(ceil(image_size/chunk_size)); j+=2) {
        // copy images one chunk at a time
        begin0 = j*chunk_size;
        begin1 = begin0 + chunk_size;
        cudaSafeCall( cudaMemcpyAsync
                      ( d_images0,
                        h_images + image_offset + begin0,
                        chunk_size * sizeof(cuDoubleComplex),
                        cudaMemcpyHostToDevice, s0 ) );
        cudaSafeCall( cudaMemcpyAsync
                      ( d_images1,
                        h_images + image_offset + begin1,
                        chunk_size * sizeof(cuDoubleComplex),
                        cudaMemcpyHostToDevice, s1 ) );
        // add to sum
        add_images_kernel_streams<<<blocks_per_grid,threads_per_block,0,s0>>>
          (d_images0,begin0,image_size,d_summed_image);
        add_images_kernel_streams<<<blocks_per_grid,threads_per_block,0,s1>>>
          (d_images1,begin1,image_size,d_summed_image);
      }
    }
    // copy result from GPU to host
    cudaSafeCall( cudaStreamSynchronize( s0 ) );
    cudaSafeCall( cudaStreamSynchronize( s1 ) );
    cudaSafeCall( cudaMemcpy( h_summed_image, d_summed_image,
                              image_size * sizeof(cuDoubleComplex),
                              cudaMemcpyDeviceToHost ) );
    scitbx::af::shared<std::complex<double> > summed_image(image_size);
    for (int i=0; i<image_size; i++) {
      summed_image[i] = std::complex<double>(cuCreal(h_summed_image[i]),
                                             cuCimag(h_summed_image[i]));
    }

    // clean up
    cudaSafeCall( cudaFreeHost( h_images ) );
    cudaSafeCall( cudaFreeHost( h_summed_image ) );
    cudaSafeCall( cudaFree( d_images0 ) );
    cudaSafeCall( cudaFree( d_images1 ) );
    cudaSafeCall( cudaFree( d_summed_image ) );
    cudaSafeCall( cudaStreamDestroy( s0 ) );
    cudaSafeCall( cudaStreamDestroy( s1 ) );

    // end timer
    // cudaSafeCall( cudaEventRecord(stop) );
    // cudaSafeCall( cudaEventSynchronize(stop) );
    // float elapsedTime;
    // cudaSafeCall( cudaEventElapsedTime(&elapsedTime,start,stop) );
    // std::cout << elapsedTime << "\n";
    // cudaSafeCall( cudaEventDestroy(start) );
    // cudaSafeCall( cudaEventDestroy(stop) );

    return summed_image;
  }

  /* ==========================================================================
     Fast implementation without streams (probably due to excessive transfers
     to and from global memory, and no copying of original data)
  */
  __global__ void add_images_kernel
  (const cuDoubleComplex* images, const int image_size, const int n_images,
   cuDoubleComplex* summed_image) {

    int i = blockDim.x * blockIdx.x + threadIdx.x;

    // each thead sums over all images for one element
    if (i < image_size) {
      cuDoubleComplex sum = make_cuDoubleComplex(0.0,0.0);
      for (int c_i=0; c_i<n_images; c_i++) {
        sum = cuCadd(sum,images[c_i*image_size + i]);
      }

      // transfer sum to global
      summed_image[i] = sum;
    }
  }

  scitbx::af::shared<std::complex<double> > cuda_add_images
  (const scitbx::af::const_ref<std::complex<double> >& images,
   const int& image_size, const int& n_images, const int& gpu_id) {

    // start GPU
    cudaSafeCall( cudaSetDevice(gpu_id) );

    // initialize timer
    // cudaEvent_t start, stop;
    // cudaSafeCall( cudaEventCreate(&start) );
    // cudaSafeCall( cudaEventCreate(&stop) );
    // cudaSafeCall( cudaEventRecord(start) );

    // allocate and initialize arrays
    cuDoubleComplex * h_images, * d_images;
    cuDoubleComplex * h_summed_image, * d_summed_image;
    h_images = (cuDoubleComplex*)&images[0];
    h_summed_image = new cuDoubleComplex[image_size];
    cudaSafeCall( cudaMalloc( (void**)&d_images,
                              images.size() * sizeof(cuDoubleComplex) ) );
    cudaSafeCall( cudaMalloc( (void**)&d_summed_image,
                              image_size * sizeof(cuDoubleComplex) ) );
    cudaSafeCall( cudaMemcpy( d_images, h_images,
                              images.size() * sizeof(cuDoubleComplex),
                              cudaMemcpyHostToDevice ) );

    // run kernel
    int bpg = (image_size + threads_per_block - 1)/threads_per_block;
    add_images_kernel<<<bpg,threads_per_block>>>
      (d_images,image_size,n_images,d_summed_image);

    // copy result from GPU
    cudaSafeCall( cudaMemcpy( h_summed_image, d_summed_image,
                              image_size * sizeof(cuDoubleComplex),
                              cudaMemcpyDeviceToHost ) );
    scitbx::af::shared<std::complex<double> > summed_image
      ((std::complex<double>*)&h_summed_image[0],
       (std::complex<double>*)&h_summed_image[0] + image_size);

    // clean up
    cudaSafeCall( cudaFree( d_images ) );
    cudaSafeCall( cudaFree( d_summed_image ) );

    // end timer
    // cudaSafeCall( cudaEventRecord(stop) );
    // cudaSafeCall( cudaEventSynchronize(stop) );
    // float elapsedTime;
    // cudaSafeCall( cudaEventElapsedTime(&elapsedTime,start,stop) );
    // std::cout << elapsedTime << "\n";
    // cudaSafeCall( cudaEventDestroy(start) );
    // cudaSafeCall( cudaEventDestroy(stop) );

    return summed_image;
  }

}
}

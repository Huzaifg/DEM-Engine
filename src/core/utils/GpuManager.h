#ifndef SGPS_GPU_MANAGER
#define SGPS_GPU_MANAGER

#include <cuda_runtime_api.h>
#include <vector>
#include <core/ApiVersion.h>

class GpuManager {
  public:
    GpuManager(unsigned int total_streams = 1);
    ~GpuManager();

    struct StreamInfo {
      public:
        int device;
        cudaStream_t stream;
        
        bool _impl_active; // Reserved for the implementation
    };

    // Returns the LEAST number of streams available on any device
    unsigned int getStreamsPerDevice();
    // Returns the HIGHEST number of streams per device
    unsigned int getMaxStreamsPerDevice();

    unsigned int getNumDevices();

    // DO NOT USE UNLESS YOU INTEND TO MANUALLY HANDLE YOUR STREAMS
    const std::vector<StreamInfo>& getStreamsFromDevice(int index);
    
    // DO NOT USE UNLESS YOU INTEND TO MANUALLY HANDLE YOUR STREAMS
    const std::vector<StreamInfo>& getStreamsFromDevice(const StreamInfo&);

    // Get a stream which hasn't been used yet and mark it as used
    const StreamInfo& getAvailableStream();
    const StreamInfo& getAvailableStreamFromDevice(int index);

    // Mark a stream as unused
    void setStreamAvailable(const StreamInfo&);

  private:
    std::vector<std::vector<StreamInfo>> streams;
};

#endif

//  Copyright (c) 2021, SBEL GPU Development Team
//  Copyright (c) 2021, University of Wisconsin - Madison
//  All rights reserved.

#pragma once

#include <granular/DataStructs.h>
#include <granular/GranularStructs.h>
#include <granular/GranularDefines.h>
#include <core/utils/GpuManager.h>

namespace sgps {

// This file should not be visible to gcc so it's difficult to make functions here templated. We probably have to bear
// with writing each version of the same functions individually.
void cubPrefixScan_binSphere(binsSphereTouches_t* d_in,
                             binSphereTouchPairs_t* d_out,
                             size_t n,
                             cudaStream_t& this_stream,
                             DEMSolverStateDataKT& scratchPad);

void cubPrefixScan_contacts(spheresBinTouches_t* d_in,
                            contactPairs_t* d_out,
                            size_t n,
                            cudaStream_t& this_stream,
                            DEMSolverStateDataKT& scratchPad);

// template <typename T1, typename T2>
// void cubPrefixScan(T1* d_in,
//                    T2* d_out,
//                    size_t n,
//                    cudaStream_t& this_stream,
//                    DEMSolverStateData& scratchPad);

void cubSortByKeys(binID_t* d_keys_in,
                   binID_t* d_keys_out,
                   bodyID_t* d_vals_in,
                   bodyID_t* d_vals_out,
                   size_t n,
                   cudaStream_t& this_stream,
                   DEMSolverStateDataKT& scratchPad);

void cubUnique(binID_t* d_in,
               binID_t* d_out,
               size_t* d_num_out,
               size_t n,
               cudaStream_t& this_stream,
               DEMSolverStateDataKT& scratchPad);

void cubRunLengthEncode(binID_t* d_in,
                        binID_t* d_unique_out,
                        spheresBinTouches_t* d_counts_out,
                        size_t* d_num_out,
                        size_t n,
                        cudaStream_t& this_stream,
                        DEMSolverStateDataKT& scratchPad);

void cubCollectForces(clumpBodyInertiaOffset_t* inertiaPropOffsets,
                      bodyID_t* idA,
                      bodyID_t* idB,
                      float3* contactForces,
                      float3* contactPointA,
                      float3* contactPointB,
                      float* clump_h2aX,
                      float* clump_h2aY,
                      float* clump_h2aZ,
                      float* clump_h2AlphaX,
                      float* clump_h2AlphaY,
                      float* clump_h2AlphaZ,
                      bodyID_t* ownerClumpBody,
                      float* massClumpBody,
                      float* mmiXX,
                      float* mmiYY,
                      float* mmiZZ,
                      double h,
                      size_t nContactPairs,
                      size_t nClumps,
                      double l,
                      bool contactPairArr_isFresh,
                      cudaStream_t& this_stream,
                      DEMSolverStateDataDT& scratchPad,
                      clumpBodyInertiaOffset_t nDistinctClumpBodyTopologies);

}  // namespace sgps
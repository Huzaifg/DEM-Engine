//  Copyright (c) 2021, SBEL GPU Development Team
//  Copyright (c) 2021, University of Wisconsin - Madison
//
//	SPDX-License-Identifier: BSD-3-Clause

#include <cub/cub.cuh>
// #include <thrust/sort.h>
#include <core/utils/JitHelper.h>
#include <nvmath/helper_math.cuh>

#include <algorithms/DEMCubBasedSubroutines.h>
#include <DEM/HostSideHelpers.hpp>

#include <algorithms/DEMCubWrappers.cu>

#include <core/utils/GpuError.h>

namespace sgps {

void collectContactForces(std::shared_ptr<jitify::Program>& collect_force_kernels,
                          DEMDataDT* granData,
                          const size_t nContactPairs,
                          const size_t nClumps,
                          bool contactPairArr_isFresh,
                          cudaStream_t& this_stream,
                          DEMSolverStateData& scratchPad,
                          SolverTimers& timers) {
    // Preparation: allocate enough temp array memory and chop it to pieces, for the usage of cub operations. Note that
    // if contactPairArr_isFresh is false, then this allocation should not alter the size and content of the temp array
    // space, so the information in it can be used in the next iteration.
    size_t cachedArraySizeOwner = (size_t)2 * nContactPairs * sizeof(bodyID_t);
    // Use temp vector 0 to store the flattened owner IDs. So the rest of the temp arrays start from number 1.
    bodyID_t* idAOwner = (bodyID_t*)scratchPad.allocateTempVector(0, cachedArraySizeOwner);
    bodyID_t* idBOwner = (bodyID_t*)(idAOwner + nContactPairs);
    // size_t cachedArraySizeMass = (size_t)2 * nContactPairs * sizeof(float);
    // size_t cachedArraySizeMOI = (size_t)2 * nContactPairs * sizeof(float3);
    // float* massAOwner = (float*)scratchPad.allocateCachedMass(cachedArraySizeMass);
    // float* massBOwner = (float*)(massAOwner + nContactPairs);
    // float3* moiAOwner = (float3*)scratchPad.allocateCachedMOI(cachedArraySizeMOI);
    // float3* moiBOwner = (float3*)(moiAOwner + nContactPairs);

    // size_t blocks_needed_for_twice_contacts = (2 * nContactPairs + SGPS_DEM_NUM_BODIES_PER_BLOCK - 1) /
    // SGPS_DEM_NUM_BODIES_PER_BLOCK;
    size_t blocks_needed_for_contacts =
        (nContactPairs + SGPS_DEM_NUM_BODIES_PER_BLOCK - 1) / SGPS_DEM_NUM_BODIES_PER_BLOCK;
    if (contactPairArr_isFresh) {
        // First step, prepare the owner ID array (nContactPairs * bodyID_t) for usage in final reduction by key (do it
        // for both A and B)
        // Note for A, it is always a sphere or a triangle
        collect_force_kernels->kernel("cashInOwnerIndexA")
            .instantiate()
            .configure(dim3(blocks_needed_for_contacts), dim3(SGPS_DEM_NUM_BODIES_PER_BLOCK), 0, this_stream)
            .launch(idAOwner, granData->idGeometryA, granData->ownerClumpBody, granData->contactType, nContactPairs);
        GPU_CALL(cudaStreamSynchronize(this_stream));

        // But for B, it can be sphere, triangle or some analytical geometries
        collect_force_kernels->kernel("cashInOwnerIndexB")
            .instantiate()
            .configure(dim3(blocks_needed_for_contacts), dim3(SGPS_DEM_NUM_BODIES_PER_BLOCK), 0, this_stream)
            .launch(idBOwner, granData->idGeometryB, granData->ownerClumpBody, granData->contactType, nContactPairs);
        GPU_CALL(cudaStreamSynchronize(this_stream));
        // displayArray<bodyID_t>(idAOwner, nContactPairs);
        // displayArray<bodyID_t>(idBOwner, nContactPairs);
    }

    // ==============================================
    // 2nd, combine mass and force to get (contact pair-wise) acceleration, which will be reduced...
    // Note here allocated is temp vector, since unlike cached vectors, they cannot be reused in the next iteration
    size_t tempArraySizeAcc = (size_t)2 * nContactPairs * sizeof(float3);
    size_t tempArraySizeAcc_sorted = (size_t)2 * nContactPairs * sizeof(float3);
    size_t tempArraySizeOwnerAcc = (size_t)nClumps * sizeof(float3);
    size_t tempArraySizeOwner = (size_t)nClumps * sizeof(bodyID_t);
    float3* acc_A = (float3*)scratchPad.allocateTempVector(1, tempArraySizeAcc);
    float3* acc_B = (float3*)(acc_A + nContactPairs);
    float3* acc_A_sorted = (float3*)scratchPad.allocateTempVector(2, tempArraySizeAcc_sorted);
    // float3* acc_B_sorted = (float3*)(acc_A_sorted  + nContactPairs);
    bodyID_t* idAOwner_sorted = (bodyID_t*)scratchPad.allocateTempVector(3, cachedArraySizeOwner);
    // bodyID_t* idBOwner_sorted = (bodyID_t*)(idAOwner_sorted + nContactPairs);
    float3* accOwner = (float3*)scratchPad.allocateTempVector(
        4, tempArraySizeOwnerAcc);  // can store both linear and angular acceleration
    bodyID_t* uniqueOwner = (bodyID_t*)scratchPad.allocateTempVector(5, tempArraySizeOwner);
    // Collect accelerations for body A (modifier used to be h * h / l when we stored acc as h^2*acc)
    // NOTE!! If you pass floating point number to kernels, the number needs to be something like 1.f, not 1.0.
    // Somtimes 1.0 got converted to 0.f with the kernel call.
    collect_force_kernels->kernel("forceToAcc")
        .instantiate()
        .configure(dim3(blocks_needed_for_contacts), dim3(SGPS_DEM_NUM_BODIES_PER_BLOCK), 0, this_stream)
        .launch(acc_A, granData->contactForces, idAOwner, 1.f, nContactPairs, granData);
    GPU_CALL(cudaStreamSynchronize(this_stream));
    // and don't forget body B
    collect_force_kernels->kernel("forceToAcc")
        .instantiate()
        .configure(dim3(blocks_needed_for_contacts), dim3(SGPS_DEM_NUM_BODIES_PER_BLOCK), 0, this_stream)
        .launch(acc_B, granData->contactForces, idBOwner, -1.f, nContactPairs, granData);
    GPU_CALL(cudaStreamSynchronize(this_stream));
    // displayFloat3(acc_A, 2 * nContactPairs);
    // displayFloat3(granData->contactForces, nContactPairs);

    // Reducing the acceleration (2 * nContactPairs for both body A and B)
    // Note: to do this, idAOwner needs to be sorted along with acc_A. So we sort first.
    cubDEMSortByKeys<bodyID_t, float3, DEMSolverStateData>(idAOwner, idAOwner_sorted, acc_A, acc_A_sorted,
                                                           nContactPairs * 2, this_stream, scratchPad);
    // Then we reduce by key
    // This variable stores the cub output of how many cub runs it executed for collecting forces
    size_t* pForceCollectionRuns = scratchPad.pTempSizeVar1;
    CubFloat3Add float3_add_op;
    cubDEMReduceByKeys<bodyID_t, float3, CubFloat3Add, DEMSolverStateData>(
        idAOwner_sorted, uniqueOwner, acc_A_sorted, accOwner, pForceCollectionRuns, float3_add_op, nContactPairs * 2,
        this_stream, scratchPad);
    // Then we stash acceleration
    size_t blocks_needed_for_stashing =
        (*pForceCollectionRuns + SGPS_DEM_NUM_BODIES_PER_BLOCK - 1) / SGPS_DEM_NUM_BODIES_PER_BLOCK;
    collect_force_kernels->kernel("stashElem")
        .instantiate()
        .configure(dim3(blocks_needed_for_stashing), dim3(SGPS_DEM_NUM_BODIES_PER_BLOCK), 0, this_stream)
        .launch(granData->aX, granData->aY, granData->aZ, uniqueOwner, accOwner, *pForceCollectionRuns);
    GPU_CALL(cudaStreamSynchronize(this_stream));
    // displayArray<float>(granData->aX, nClumps);
    // displayArray<float>(granData->aY, nClumps);
    // displayArray<float>(granData->aZ, nClumps);

    // =====================================================
    // Then take care of angular accelerations
    float3* alpha_A = (float3*)(acc_A);  // Memory spaces for accelerations can be reused
    float3* alpha_B = (float3*)(acc_B);
    float3* alpha_A_sorted = (float3*)(acc_A_sorted);
    // float3* alpha_B_sorted = (float3*)(acc_B_sorted);
    // collect angular accelerations for body A (modifier used to be h * h when we stored acc as h^2*acc)
    collect_force_kernels->kernel("forceToAngAcc")
        .instantiate()
        .configure(dim3(blocks_needed_for_contacts), dim3(SGPS_DEM_NUM_BODIES_PER_BLOCK), 0, this_stream)
        .launch(alpha_A, granData->contactPointGeometryA, granData->oriQw, granData->oriQx, granData->oriQy,
                granData->oriQz, granData->contactForces, granData->contactTorque_convToForce, idAOwner, 1.f,
                nContactPairs, granData);
    GPU_CALL(cudaStreamSynchronize(this_stream));
    // and don't forget body B
    collect_force_kernels->kernel("forceToAngAcc")
        .instantiate()
        .configure(dim3(blocks_needed_for_contacts), dim3(SGPS_DEM_NUM_BODIES_PER_BLOCK), 0, this_stream)
        .launch(alpha_B, granData->contactPointGeometryB, granData->oriQw, granData->oriQx, granData->oriQy,
                granData->oriQz, granData->contactForces, granData->contactTorque_convToForce, idBOwner, -1.f,
                nContactPairs, granData);
    GPU_CALL(cudaStreamSynchronize(this_stream));
    // Reducing the angular acceleration (2 * nContactPairs for both body A and B)
    // Note: to do this, idAOwner needs to be sorted along with alpha_A. So we sort first.
    cubDEMSortByKeys<bodyID_t, float3, DEMSolverStateData>(idAOwner, idAOwner_sorted, alpha_A, alpha_A_sorted,
                                                           nContactPairs * 2, this_stream, scratchPad);
    // Then we reduce
    cubDEMReduceByKeys<bodyID_t, float3, CubFloat3Add, DEMSolverStateData>(
        idAOwner_sorted, uniqueOwner, alpha_A_sorted, accOwner, pForceCollectionRuns, float3_add_op, nContactPairs * 2,
        this_stream, scratchPad);
    // Then we stash angular acceleration
    blocks_needed_for_stashing =
        (*pForceCollectionRuns + SGPS_DEM_NUM_BODIES_PER_BLOCK - 1) / SGPS_DEM_NUM_BODIES_PER_BLOCK;
    collect_force_kernels->kernel("stashElem")
        .instantiate()
        .configure(dim3(blocks_needed_for_stashing), dim3(SGPS_DEM_NUM_BODIES_PER_BLOCK), 0, this_stream)
        .launch(granData->alphaX, granData->alphaY, granData->alphaZ, uniqueOwner, accOwner, *pForceCollectionRuns);
    GPU_CALL(cudaStreamSynchronize(this_stream));
}

}  // namespace sgps

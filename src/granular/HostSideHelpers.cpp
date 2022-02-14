//  Copyright (c) 2021, SBEL GPU Development Team
//  Copyright (c) 2021, University of Wisconsin - Madison
//  All rights reserved.

#ifndef SGPS_DEM_HOST_HELPERS
#define SGPS_DEM_HOST_HELPERS

#pragma once
#include <iostream>
#include <list>
#include <cmath>
#include <vector>
#include <algorithm>
#include <helper_math.cuh>

namespace sgps {

template <typename T1>
inline void displayArray(T1* arr, size_t n) {
    for (size_t i = 0; i < n; i++) {
        std::cout << arr[i] << " ";
    }
    std::cout << std::endl;
}

inline void displayFloat3(float3* arr, size_t n) {
    for (size_t i = 0; i < n; i++) {
        std::cout << "(" << arr[i].x << ", " << arr[i].y << ", " << arr[i].z << "), ";
    }
    std::cout << std::endl;
}

template <typename T1>
inline void hostPrefixScan(T1* arr, size_t n) {
    T1 buffer_previous = arr[0];
    arr[0] = 0;
    for (size_t i = 1; i < n; i++) {
        T1 item_to_add = buffer_previous;
        buffer_previous = arr[i];
        arr[i] = arr[i - 1] + item_to_add;
    }
}

template <typename T1>
inline void elemSwap(T1* x, T1* y) {
    T1 tmp = *x;
    *x = *y;
    *y = tmp;
}

template <typename T1, typename T2>
inline void hostSortByKey(T1* keys, T2* vals, size_t n) {
    // Just bubble sort it
    bool swapped;
    for (size_t i = 0; i < n - 1; i++) {
        swapped = false;
        for (size_t j = 0; j < n - i - 1; j++) {
            if (keys[j] > keys[j + 1]) {
                elemSwap<T1>(&keys[j], &keys[j + 1]);
                elemSwap<T2>(&vals[j], &vals[j + 1]);
                swapped = true;
            }
        }
        if (!swapped) {
            break;
        }
    }
}

template <typename T1>
inline void hostScanForJumpsNum(T1* arr, size_t n, unsigned int minSegLen, size_t& total_found) {
    size_t i = 0;
    total_found = 0;
    while (i < n - 1) {
        size_t thisIndx = i;
        T1 thisItem = arr[i];
        do {
            i++;
        } while (arr[i] == thisItem && i < n - 1);
        if (i - thisIndx >= minSegLen || (i == n - 1 && i - thisIndx + 1 >= minSegLen && arr[i] == thisItem)) {
            total_found++;
        }
    }
}

// Tell each active bin where to find its touching spheres
template <typename T1, typename T2, typename T3>
inline void hostScanForJumps(T1* arr, T1* arr_elem, T2* jump_loc, T3* jump_len, size_t n, unsigned int minSegLen) {
    size_t total_found = 0;
    T2 i = 0;
    unsigned int thisSegLen;
    while (i < n - 1) {
        thisSegLen = 0;
        T2 thisIndx = i;
        T1 thisItem = arr[i];
        do {
            i++;
            thisSegLen++;
        } while (arr[i] == thisItem && i < n - 1);
        if (i - thisIndx >= minSegLen || (i == n - 1 && i - thisIndx + 1 >= minSegLen && arr[i] == thisItem)) {
            jump_loc[total_found] = thisIndx;
            if (i == n - 1 && i - thisIndx + 1 >= minSegLen && arr[i] == thisItem)
                thisSegLen++;
            jump_len[total_found] = thisSegLen;
            arr_elem[total_found] = arr[thisIndx];
            total_found++;
        }
    }
}

// Note we assume the ``contact force'' here is actually acceleration written in terms of multiples of l
// And we collect h2a, not force
inline void hostCollectForces(bodyID_t* idA,
                              bodyID_t* idB,
                              float3* contactForces,
                              float* clump_h2aX,
                              float* clump_h2aY,
                              float* clump_h2aZ,
                              bodyID_t* ownerClumpBody,
                              double h,
                              size_t n) {
    for (size_t i = 0; i < n; i++) {
        bodyID_t bodyA = idA[i];
        bodyID_t bodyB = idB[i];
        bodyID_t AOwner = ownerClumpBody[bodyA];
        clump_h2aX[AOwner] += (double)contactForces[i].x * h * h;
        clump_h2aY[AOwner] += (double)contactForces[i].y * h * h;
        clump_h2aZ[AOwner] += (double)contactForces[i].z * h * h;

        bodyID_t BOwner = ownerClumpBody[bodyB];
        clump_h2aX[BOwner] += -(double)contactForces[i].x * h * h;
        clump_h2aY[BOwner] += -(double)contactForces[i].y * h * h;
        clump_h2aZ[BOwner] += -(double)contactForces[i].z * h * h;
    }
}

// Note we assume the ``contact force'' here is actually acceleration written in terms of multiples of l
// And we collect h2Alpha, not torque
inline void hostCollectTorques(bodyID_t* idA,
                               bodyID_t* idB,
                               float3* contactForces,
                               float3* contactLocA,
                               float3* contactLocB,
                               float* clump_h2AlphaX,
                               float* clump_h2AlphaY,
                               float* clump_h2AlphaZ,
                               bodyID_t* ownerClumpBody,
                               double h,
                               size_t n,
                               double l) {
    float mass = 1.;
    float MOI = 1.;
    for (size_t i = 0; i < n; i++) {
        // First, recover true force
        float3 F = contactForces[i] * l * mass;
        // Then, compute alpha as torque/moi
        float3 CPA = contactLocA[i];
        float3 CPB = contactLocB[i];
        float3 alphaA = cross(CPA, F) / 1.;
        float3 alphaB = cross(CPB, -F) / 1.;

        bodyID_t bodyA = idA[i];
        bodyID_t bodyB = idB[i];
        bodyID_t AOwner = ownerClumpBody[bodyA];
        clump_h2AlphaX[AOwner] += (double)alphaA.x * h * h;
        clump_h2AlphaY[AOwner] += (double)alphaA.y * h * h;
        clump_h2AlphaZ[AOwner] += (double)alphaA.z * h * h;

        bodyID_t BOwner = ownerClumpBody[bodyB];
        clump_h2AlphaX[BOwner] += (double)alphaB.x * h * h;
        clump_h2AlphaY[BOwner] += (double)alphaB.y * h * h;
        clump_h2AlphaZ[BOwner] += (double)alphaB.z * h * h;
    }
}

/// A light-weight grid sampler that can be used to generate the initial stage of the granular system
inline std::vector<float3> DEMBoxGridSampler(float3 BoxCenter, float3 HalfDims, float GridSize) {
    std::vector<float3> points;
    for (float z = BoxCenter.z - HalfDims.z; z <= BoxCenter.z + HalfDims.z; z += GridSize) {
        for (float y = BoxCenter.y - HalfDims.y; y <= BoxCenter.y + HalfDims.y; y += GridSize) {
            for (float x = BoxCenter.x - HalfDims.x; x <= BoxCenter.x + HalfDims.x; x += GridSize) {
                float3 xyz;
                xyz.x = x;
                xyz.y = y;
                xyz.z = z;
                points.push_back(xyz);
            }
        }
    }
    return points;
}

/// Host version of applying a quaternion to a vector
template <typename T1, typename T2>
inline void hostApplyOriQ2Vector3(T1& X, T1& Y, T1& Z, const T2& Q0, const T2& Q1, const T2& Q2, const T2& Q3) {
    T1 oldX = X;
    T1 oldY = Y;
    T1 oldZ = Z;
    X = ((T2)2.0 * (Q0 * Q0 + Q1 * Q1) - (T2)1.0) * oldX + ((T2)2.0 * (Q1 * Q2 - Q0 * Q3)) * oldY +
        ((T2)2.0 * (Q1 * Q3 + Q0 * Q2)) * oldZ;
    Y = ((T2)2.0 * (Q1 * Q2 + Q0 * Q3)) * oldX + ((T2)2.0 * (Q0 * Q0 + Q2 * Q2) - (T2)1.0) * oldY +
        ((T2)2.0 * (Q2 * Q3 - Q0 * Q1)) * oldZ;
    Z = ((T2)2.0 * (Q1 * Q3 - Q0 * Q2)) * oldX + ((T2)2.0 * (Q2 * Q3 + Q0 * Q1)) * oldY +
        ((T2)2.0 * (Q0 * Q0 + Q3 * Q3) - (T2)1.0) * oldZ;
}

}  // namespace sgps

#endif

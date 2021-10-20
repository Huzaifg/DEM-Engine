//  Copyright (c) 2021, SBEL GPU Development Team
//  Copyright (c) 2021, University of Wisconsin - Madison
//  All rights reserved.

#pragma once
#include <climits>
#include <stdint.h>

namespace sgps {
#define BITS_PER_BYTE 8
#define N_MANUFACTURED_ITEMS 4
#define N_INPUT_ITEMS 4
#define VOXEL_RES_POWER2 16
#define WAIT_GRANULARITY_MS 1

typedef uint64_t voxelID_default_t;
typedef unsigned int bodyID_default_t;
typedef unsigned int materialsOffset_default_t;
typedef unsigned int clumpBodyInertiaOffset_default_t;
typedef unsigned int stateVectors_default_t;
typedef unsigned int clumpComponentOffset_t;

// typedef unsigned int distinctSphereRelativePositions_default_t;
// typedef unsigned int distinctSphereRadiiOffset_default_t;

// somehow add array materialsArray and radiiArray??

}  // namespace sgps

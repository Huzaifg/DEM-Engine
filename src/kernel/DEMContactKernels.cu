// DEM contact detection-related custom kernels
#include <DEM/DEMDefines.h>
#include <kernel/DEMHelperKernels.cu>

#include <cub/util_ptx.cuh>

// If clump templates are jitified, they will be below
_clumpTemplateDefs_;
// Family mask, _nFamilyMaskEntries_ elements are in this array
// __constant__ __device__ bool familyMasks[] = {_familyMasks_};

__global__ void getNumberOfContactsEachBin(sgps::DEMSimParams* simParams,
                                           sgps::DEMDataKT* granData,
                                           sgps::bodyID_t* sphereIDsEachBinTouches_sorted,
                                           sgps::binID_t* activeBinIDs,
                                           sgps::spheresBinTouches_t* numSpheresBinTouches,
                                           sgps::binSphereTouchPairs_t* sphereIDsLookUpTable,
                                           sgps::spheresBinTouches_t* numContactsInEachBin,
                                           size_t nActiveBins) {
    // Only active bins got execute this...
    sgps::binID_t myActiveID = blockIdx.x * blockDim.x + threadIdx.x;
    // I need to store all the sphereIDs that I am supposed to look into
    // A100 has about 164K shMem... these arrays really need to be small, or we can only fit a small number of bins in
    // one block
    sgps::bodyID_t ownerIDs[SGPS_DEM_MAX_SPHERES_PER_BIN];
    float radii[SGPS_DEM_MAX_SPHERES_PER_BIN];
    double bodyX[SGPS_DEM_MAX_SPHERES_PER_BIN];
    double bodyY[SGPS_DEM_MAX_SPHERES_PER_BIN];
    double bodyZ[SGPS_DEM_MAX_SPHERES_PER_BIN];
    sgps::family_t ownerFamily[SGPS_DEM_MAX_SPHERES_PER_BIN];
    if (myActiveID < nActiveBins) {
        // I got a true bin ID
        sgps::binID_t binID = activeBinIDs[myActiveID];

        sgps::spheresBinTouches_t contact_count = 0;
        // Grab the bodies that I care, put into local memory
        sgps::spheresBinTouches_t nBodiesMeHandle = numSpheresBinTouches[myActiveID];
        if (nBodiesMeHandle > SGPS_DEM_MAX_SPHERES_PER_BIN) {
            SGPS_DEM_ABORT_KERNEL("Bin %u contains %u sphere components, exceeding maximum allowance (%u)\n",
                                  myActiveID, nBodiesMeHandle, SGPS_DEM_MAX_SPHERES_PER_BIN);
        }

        sgps::binSphereTouchPairs_t myBodiesTableEntry = sphereIDsLookUpTable[myActiveID];
        // printf("nBodies: %u\n", nBodiesMeHandle);
        for (sgps::spheresBinTouches_t i = 0; i < nBodiesMeHandle; i++) {
            sgps::bodyID_t sphereID = sphereIDsEachBinTouches_sorted[myBodiesTableEntry + i];
            ownerIDs[i] = granData->ownerClumpBody[sphereID];
            ownerFamily[i] = granData->familyID[ownerIDs[i]];
            double ownerX, ownerY, ownerZ;
            float myRelPosX, myRelPosY, myRelPosZ, myRadius;

            // Get my component offset info from either jitified arrays or global memory
            // Outputs myRelPosXYZ, myRadius (in CD kernels, radius needs to be expanded)
            // Use an input named exactly `sphereID' which is the id of this sphere component
            {
                _componentAcqStrat_;
                myRadius += simParams->beta;
            }

            voxelID2Position<double, sgps::voxelID_t, sgps::subVoxelPos_t>(
                ownerX, ownerY, ownerZ, granData->voxelID[ownerIDs[i]], granData->locX[ownerIDs[i]],
                granData->locY[ownerIDs[i]], granData->locZ[ownerIDs[i]], _nvXp2_, _nvYp2_, _voxelSize_, _l_);
            float myOriQ0 = granData->oriQ0[ownerIDs[i]];
            float myOriQ1 = granData->oriQ1[ownerIDs[i]];
            float myOriQ2 = granData->oriQ2[ownerIDs[i]];
            float myOriQ3 = granData->oriQ3[ownerIDs[i]];
            applyOriQToVector3<float, sgps::oriQ_t>(myRelPosX, myRelPosY, myRelPosZ, myOriQ0, myOriQ1, myOriQ2,
                                                    myOriQ3);
            bodyX[i] = ownerX + (double)myRelPosX;
            bodyY[i] = ownerY + (double)myRelPosY;
            bodyZ[i] = ownerZ + (double)myRelPosZ;
            radii[i] = myRadius;
        }

        for (sgps::spheresBinTouches_t bodyA = 0; bodyA < nBodiesMeHandle; bodyA++) {
            for (sgps::spheresBinTouches_t bodyB = bodyA + 1; bodyB < nBodiesMeHandle; bodyB++) {
                // For 2 bodies to be considered in contact, the contact point must be in this bin (to avoid
                // double-counting), and they do not belong to the same clump
                if (ownerIDs[bodyA] == ownerIDs[bodyB])
                    continue;

                // Grab family number from memory (not jitified: b/c family number can change frequently in a sim)
                unsigned int bodyAFamily = ownerFamily[bodyA];
                unsigned int bodyBFamily = ownerFamily[bodyB];
                unsigned int maskMatID = locateMaskPair<unsigned int>(bodyAFamily, bodyBFamily);
                // If marked no contact, skip ths iteration
                if (granData->familyMasks[maskMatID] != sgps::DEM_DONT_PREVENT_CONTACT) {
                    continue;
                }

                double contactPntX;
                double contactPntY;
                double contactPntZ;
                bool in_contact;
                in_contact = checkSpheresOverlap<double>(bodyX[bodyA], bodyY[bodyA], bodyZ[bodyA], radii[bodyA],
                                                         bodyX[bodyB], bodyY[bodyB], bodyZ[bodyB], radii[bodyB],
                                                         contactPntX, contactPntY, contactPntZ);
                sgps::binID_t contactPntBin = getPointBinID<sgps::binID_t>(
                    contactPntX, contactPntY, contactPntZ, simParams->binSize, simParams->nbX, simParams->nbY);

                /*
                printf("contactPntBin: %u, %u, %u\n", (unsigned int)(contactPntX/_binSize_),
                                                        (unsigned int)(contactPntY/_binSize_),
                                                        (unsigned int)(contactPntZ/_binSize_));
                unsigned int ZZ = binID/(_nbX_*_nbY_);
                unsigned int YY = binID%(_nbX_*_nbY_)/_nbX_;
                unsigned int XX = binID%(_nbX_*_nbY_)%_nbX_;
                printf("binID: %u, %u, %u\n", XX,YY,ZZ);
                printf("bodyA: %f, %f, %f\n", bodyX[bodyA], bodyY[bodyA], bodyZ[bodyA]);
                printf("bodyB: %f, %f, %f\n", bodyX[bodyB], bodyY[bodyB], bodyZ[bodyB]);
                printf("contactPnt: %f, %f, %f\n", contactPntX, contactPntY, contactPntZ);
                printf("contactPntBin: %u\n", contactPntBin);
                */

                if (in_contact && (contactPntBin == binID)) {
                    contact_count++;
                }
            }
        }
        numContactsInEachBin[myActiveID] = contact_count;
    }
}

__global__ void populateContactPairsEachBin(sgps::DEMSimParams* simParams,
                                            sgps::DEMDataKT* granData,
                                            sgps::bodyID_t* sphereIDsEachBinTouches_sorted,
                                            sgps::binID_t* activeBinIDs,
                                            sgps::spheresBinTouches_t* numSpheresBinTouches,
                                            sgps::binSphereTouchPairs_t* sphereIDsLookUpTable,
                                            sgps::contactPairs_t* contactReportOffsets,
                                            sgps::bodyID_t* idSphA,
                                            sgps::bodyID_t* idSphB,
                                            size_t nActiveBins) {
    // Only active bins got to execute this...
    sgps::binID_t myActiveID = blockIdx.x * blockDim.x + threadIdx.x;
    // I need to store all the sphereIDs that I am supposed to look into
    // A100 has about 164K shMem... these arrays really need to be small, or we can only fit a small number of bins in
    // one block
    sgps::bodyID_t ownerIDs[SGPS_DEM_MAX_SPHERES_PER_BIN];
    sgps::bodyID_t bodyIDs[SGPS_DEM_MAX_SPHERES_PER_BIN];
    float radii[SGPS_DEM_MAX_SPHERES_PER_BIN];
    double bodyX[SGPS_DEM_MAX_SPHERES_PER_BIN];
    double bodyY[SGPS_DEM_MAX_SPHERES_PER_BIN];
    double bodyZ[SGPS_DEM_MAX_SPHERES_PER_BIN];
    sgps::family_t ownerFamily[SGPS_DEM_MAX_SPHERES_PER_BIN];
    if (myActiveID < nActiveBins) {
        // But I got a true bin ID
        sgps::binID_t binID = activeBinIDs[myActiveID];

        // Grab the bodies that I care, put into local memory
        sgps::spheresBinTouches_t nBodiesMeHandle = numSpheresBinTouches[myActiveID];
        sgps::binSphereTouchPairs_t myBodiesTableEntry = sphereIDsLookUpTable[myActiveID];
        for (sgps::spheresBinTouches_t i = 0; i < nBodiesMeHandle; i++) {
            sgps::bodyID_t sphereID = sphereIDsEachBinTouches_sorted[myBodiesTableEntry + i];
            ownerIDs[i] = granData->ownerClumpBody[sphereID];
            ownerFamily[i] = granData->familyID[ownerIDs[i]];
            bodyIDs[i] = sphereID;
            double ownerX, ownerY, ownerZ;
            float myRelPosX, myRelPosY, myRelPosZ, myRadius;

            // Get my component offset info from either jitified arrays or global memory
            // Outputs myRelPosXYZ, myRadius (in CD kernels, radius needs to be expanded)
            // Use an input named exactly `sphereID' which is the id of this sphere component
            {
                _componentAcqStrat_;
                myRadius += simParams->beta;
            }

            voxelID2Position<double, sgps::voxelID_t, sgps::subVoxelPos_t>(
                ownerX, ownerY, ownerZ, granData->voxelID[ownerIDs[i]], granData->locX[ownerIDs[i]],
                granData->locY[ownerIDs[i]], granData->locZ[ownerIDs[i]], _nvXp2_, _nvYp2_, _voxelSize_, _l_);
            float myOriQ0 = granData->oriQ0[ownerIDs[i]];
            float myOriQ1 = granData->oriQ1[ownerIDs[i]];
            float myOriQ2 = granData->oriQ2[ownerIDs[i]];
            float myOriQ3 = granData->oriQ3[ownerIDs[i]];
            applyOriQToVector3<float, sgps::oriQ_t>(myRelPosX, myRelPosY, myRelPosZ, myOriQ0, myOriQ1, myOriQ2,
                                                    myOriQ3);
            bodyX[i] = ownerX + (double)myRelPosX;
            bodyY[i] = ownerY + (double)myRelPosY;
            bodyZ[i] = ownerZ + (double)myRelPosZ;
            radii[i] = myRadius;
        }

        // Get my offset for writing back to the global arrays that contain contact pair info
        sgps::contactPairs_t myReportOffset = contactReportOffsets[myActiveID];

        for (sgps::spheresBinTouches_t bodyA = 0; bodyA < nBodiesMeHandle; bodyA++) {
            for (sgps::spheresBinTouches_t bodyB = bodyA + 1; bodyB < nBodiesMeHandle; bodyB++) {
                // For 2 bodies to be considered in contact, the contact point must be in this bin (to avoid
                // double-counting), and they do not belong to the same clump
                if (ownerIDs[bodyA] == ownerIDs[bodyB])
                    continue;

                // Grab family number from memory (not jitified: b/c family number can change frequently in a sim)
                unsigned int bodyAFamily = ownerFamily[bodyA];
                unsigned int bodyBFamily = ownerFamily[bodyB];
                unsigned int maskMatID = locateMaskPair<unsigned int>(bodyAFamily, bodyBFamily);
                // If marked no contact, skip ths iteration
                if (granData->familyMasks[maskMatID] != sgps::DEM_DONT_PREVENT_CONTACT) {
                    continue;
                }

                double contactPntX;
                double contactPntY;
                double contactPntZ;
                bool in_contact;
                in_contact = checkSpheresOverlap<double>(bodyX[bodyA], bodyY[bodyA], bodyZ[bodyA], radii[bodyA],
                                                         bodyX[bodyB], bodyY[bodyB], bodyZ[bodyB], radii[bodyB],
                                                         contactPntX, contactPntY, contactPntZ);
                sgps::binID_t contactPntBin = getPointBinID<sgps::binID_t>(
                    contactPntX, contactPntY, contactPntZ, simParams->binSize, simParams->nbX, simParams->nbY);

                if (in_contact && (contactPntBin == binID)) {
                    idSphA[myReportOffset] = bodyIDs[bodyA];
                    idSphB[myReportOffset] = bodyIDs[bodyB];
                    myReportOffset++;
                }
            }
        }
    }
}

//  Copyright (c) 2021, SBEL GPU Development Team
//  Copyright (c) 2021, University of Wisconsin - Madison
//
//	SPDX-License-Identifier: BSD-3-Clause

// =============================================================================
// A demo that is basically the hello-world script for DEME.
// =============================================================================

#include <core/ApiVersion.h>
#include <core/utils/ThreadManager.h>
#include <DEM/API.h>
#include <DEM/HostSideHelpers.hpp>
#include <DEM/utils/Samplers.hpp>

#include <filesystem>
#include <cstdio>
#include <time.h>
#include <filesystem>

using namespace deme;
using namespace std::filesystem;

int main() {
    DEMSolver DEMSim;
    DEMSim.SetVerbosity(STEP_DEBUG);
    DEMSim.SetOutputFormat(OUTPUT_FORMAT::CSV);
    DEMSim.SetContactOutputContent({"OWNER", "FORCE", "POINT", "COMPONENT", "NORMAL", "TORQUE"});
    DEMSim.EnsureKernelErrMsgLineNum();
    DEMSim.SetNoForceRecord();

    // srand(time(NULL));
    srand(4150);

    // Special material: has a cohesion param
    auto mat_type_1 =
        DEMSim.LoadMaterial({{"E", 1e9}, {"nu", 0.3}, {"CoR", 0.8}, {"mu", 0.3}, {"Crr", 0.01}, {"Cohesion", 50}});
    auto mat_type_2 =
        DEMSim.LoadMaterial({{"E", 2e9}, {"nu", 0.4}, {"CoR", 0.6}, {"mu", 0.3}, {"Crr", 0.01}, {"Cohesion", 50}});
    std::shared_ptr<DEMMaterial> mat_type_3 = DEMSim.Duplicate(mat_type_2);
    // If you don't have this line, then CoR between thw 2 materials will take average when they are in contact
    DEMSim.SetMaterialPropertyPair("CoR", mat_type_1, mat_type_2, 0.6);
    // Even though set elsewhere to be 50, this pairwise value takes precedence when two different materials are in
    // contact.
    DEMSim.SetMaterialPropertyPair("Cohesion", mat_type_1, mat_type_2, 100.);
    DEMSim.SetMaterialPropertyPair("Cohesion", mat_type_1, mat_type_3, 100.);
    DEMSim.SetMaterialPropertyPair("Cohesion", mat_type_2, mat_type_3, 100.);

    auto sph_type_1 = DEMSim.LoadSphereType(11728., 1., mat_type_1);
    // Test clump template duplication...
    auto sph_type_2 = DEMSim.Duplicate(sph_type_1);

    std::vector<float3> input_xyz1, input_xyz2;
    std::vector<float3> input_vel1, input_vel2;
    std::vector<std::shared_ptr<DEMClumpTemplate>> input_clump_type1(1, sph_type_1);
    std::vector<std::shared_ptr<DEMClumpTemplate>> input_clump_type2(1, sph_type_2);

    // Inputs are just 2 spheres
    float sphPos = 1.2f;
    input_xyz1.push_back(make_float3(-sphPos, 0, 0));
    input_xyz2.push_back(make_float3(sphPos, 0, 0));
    input_vel1.push_back(make_float3(1.f, 0, 0));
    input_vel2.push_back(make_float3(-1.f, 0, 0));

    auto particles1 = DEMSim.AddClumps(input_clump_type1, input_xyz1);
    particles1->SetVel(input_vel1);
    particles1->SetFamily(0);
    particles1->AddOwnerWildcard("mu_custom", 0.5);
    // This one is never used in force model, yet it should not create an error
    particles1->AddOwnerWildcard("some_property", 1.0);
    auto tracker1 = DEMSim.Track(particles1);

    // DEMSim.DisableContactBetweenFamilies(0, 1);

    // Add bottom plane mesh
    auto bot_plane = DEMSim.AddWavefrontMeshObject((GET_DATA_PATH() / "mesh/plane_20by20.obj").string(), mat_type_2);
    bot_plane->SetInitPos(make_float3(0, 0, -1.25));
    bot_plane->SetMass(10000.);
    // Just testing adding another mesh...
    auto another_plane =
        DEMSim.AddWavefrontMeshObject((GET_DATA_PATH() / "mesh/plane_20by20.obj").string(), mat_type_2);
    another_plane->SetInitPos(make_float3(0, 0, -1.5));
    another_plane->SetFamily(100);
    DEMSim.DisableFamilyOutput(100);
    DEMSim.SetFamilyFixed(100);

    // Create a inspector to find out stuff
    auto max_z_finder = DEMSim.CreateInspector("clump_max_z");
    float max_z;
    auto max_v_finder = DEMSim.CreateInspector("clump_max_absv");
    float max_v;
    auto KE_finder = DEMSim.CreateInspector("clump_kinetic_energy");
    float KE;

    // A custom force model can be read in through a file and used by the simulation. Magic, right?
    auto my_force_model = DEMSim.ReadContactForceModel("ForceModelWithCohesion.cu");
    // This custom force model still uses contact history arrays, so let's define it
    my_force_model->SetPerContactWildcards({"delta_time", "delta_tan_x", "delta_tan_y", "delta_tan_z"});
    my_force_model->SetMustPairwiseMatProp({"CoR", "mu", "Crr", "Cohesion"});

    DEMSim.SetInitTimeStep(2e-5);
    DEMSim.SetGravitationalAcceleration(make_float3(0, 0, -9.8));
    DEMSim.SetCDUpdateFreq(10);
    DEMSim.SetMaxVelocity(6.);
    DEMSim.SetExpandSafetyType("auto");
    DEMSim.SetExpandSafetyMultiplier(1.2);
    DEMSim.SetIntegrator("centered_difference");

    DEMSim.Initialize();

    DEMSim.UpdateSimParams();  // Not needed; just testing if this function works...

    // You can add more clumps to simulation after initialization, like this...
    // DEMSim.ClearCache();  // Clearing cache is no longer needed
    auto particles2 = DEMSim.AddClumps(input_clump_type2, input_xyz2);
    particles2->SetVel(input_vel2);
    particles2->SetFamily(1);
    auto tracker2 = DEMSim.Track(particles2);
    DEMSim.UpdateClumps();

    // Ready simulation
    path out_dir = current_path();
    out_dir += "/DemoOutput_SingleSphereCollide";
    create_directory(out_dir);
    // bool changed_family = false;

    for (int i = 0; i < 100; i++) {
        std::cout << "Frame: " << i << std::endl;

        // if ((!changed_family) && i >= 10) {
        //     // DEMSim.ChangeFamily(1, 0);
        //     DEMSim.EnableContactBetweenFamilies(0, 1);
        //     changed_family = true;
        // }

        char filename[200];
        sprintf(filename, "%s/DEMdemo_output_%04d.csv", out_dir.c_str(), i);
        DEMSim.WriteSphereFile(std::string(filename));

        char cnt_filename[200];
        sprintf(cnt_filename, "%s/Contact_pairs_%04d.csv", out_dir.c_str(), i);
        // DEMSim.WriteContactFile(std::string(cnt_filename));

        char meshfilename[200];
        sprintf(meshfilename, "%s/DEMdemo_mesh_%04d.vtk", out_dir.c_str(), i);
        DEMSim.WriteMeshFile(std::string(meshfilename));

        // Testing persistent contact functionality...
        DEMSim.MarkPersistentContact();

        DEMSim.DoDynamicsThenSync(1e-2);
        max_z = max_z_finder->GetValue();
        max_v = max_v_finder->GetValue();
        KE = KE_finder->GetValue();

        // Test if family changer works
        float3 pos1 = tracker1->Pos();
        float3 pos2 = tracker2->Pos();
        DEMSim.ChangeClumpFamily(i % 10, std::pair<float, float>(pos1.x - 0.1, pos1.x + 0.1),
                                 std::pair<float, float>(pos1.y - 0.1, pos1.y + 0.1),
                                 std::pair<float, float>(pos1.z - 0.1, pos1.z + 0.1));
        tracker2->SetFamily(i % 10 + 1);
        unsigned int fam1 = tracker1->GetFamily(0);
        unsigned int fam2 = tracker2->GetFamily();

        std::cout << "Max Z coord is " << max_z << std::endl;
        std::cout << "Max velocity of any point is " << max_v << std::endl;
        std::cout << "Total kinetic energy is " << KE << std::endl;
        std::cout << "Particle 1 X coord is " << pos1.x << std::endl;
        std::cout << "Particle 2 X coord is " << pos2.x << std::endl;
        std::cout << "Particle 1 family is " << fam1 << std::endl;
        std::cout << "Particle 2 family is " << fam2 << std::endl;
        std::cout << "Average contacts each sphere has: " << DEMSim.GetAvgSphContacts() << std::endl;

        // Test changing material type on-the-fly...
        DEMSim.SetFamilyClumpMaterial(1, mat_type_3);
    }

    DEMSim.ShowThreadCollaborationStats();
    DEMSim.ShowTimingStats();
    std::cout << "DEMdemo_SingleSphereCollide exiting..." << std::endl;
    return 0;
}

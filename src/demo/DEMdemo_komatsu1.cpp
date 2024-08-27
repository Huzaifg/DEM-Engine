//  Copyright (c) 2021, SBEL GPU Development Team
//  Copyright (c) 2021, University of Wisconsin - Madison
//
//	SPDX-License-Identifier: BSD-3-Clause

// =============================================================================
// A bowl plowing in a pile of granular material.
// =============================================================================

#include <core/ApiVersion.h>
#include <core/utils/ThreadManager.h>
#include <DEM/API.h>
#include <DEM/HostSideHelpers.hpp>
#include <DEM/utils/Samplers.hpp>

#include <cstdio>
#include <chrono>
#include <filesystem>

using namespace deme;
using namespace std::filesystem;

const double math_PI = 3.1415927;

void AdvanceSimulation(DEMSolver& DEMSim,
                       double time,
                       float step_size,
                       unsigned int out_steps,
                       const path& out_dir,
                       unsigned int& curr_step,
                       unsigned int& currframe) {
    for (double t = 0; t < time; t += step_size, curr_step++) {
        if (curr_step % out_steps == 0) {
            char filename[200], meshfile[200];
            sprintf(filename, "%s/DEMdemo_output_%04d.csv", out_dir.c_str(), currframe);
            sprintf(meshfile, "%s/DEMdemo_excavator_%04d.vtk", out_dir.c_str(), currframe);
            DEMSim.WriteSphereFile(std::string(filename));
            DEMSim.WriteMeshFile(std::string(meshfile));
            std::cout << "Frame: " << currframe << std::endl;
            currframe++;
            DEMSim.ShowThreadCollaborationStats();
        }

        DEMSim.DoStepDynamics();
    }
}

int main() {
    DEMSolver DEMSim;
    DEMSim.UseFrictionalHertzianModel();
    DEMSim.SetVerbosity(INFO);
    DEMSim.SetOutputFormat(OUTPUT_FORMAT::CSV);

    // Scale-defining numbers of this simulation.
    float world_halfsize = 0.5;
    float bowl_bottom = -world_halfsize;

    auto mat_type_walls = DEMSim.LoadMaterial({{"E", 1e8}, {"nu", 0.3}, {"CoR", 0.3}, {"mu", 0.5}, {"Crr", 0.3}});
    auto mat_type_particles = DEMSim.LoadMaterial({{"E", 1e8}, {"nu", 0.3}, {"CoR", 0.05}, {"mu", 0.6}, {"Crr", 0.05}});
    // auto mat_type_particles = DEMSim.LoadMaterial({{"E", 1e7}, {"nu", 0.3}, {"CoR", 0.3}, {"mu", 0.5}});
    // If you don't have this line, then CoR between wall material and granular material will be 0.5 (average of the
    // two).
    DEMSim.SetMaterialPropertyPair("CoR", mat_type_walls, mat_type_particles, 0.3);
    DEMSim.SetMaterialPropertyPair("Crr", mat_type_walls, mat_type_particles, 0.3);
    DEMSim.SetMaterialPropertyPair("mu", mat_type_walls, mat_type_particles, 0.5);

    // Define the terrain particle templates
    // Calculate its mass and MOI
    float mass = 2.6e3 * 4. / 3. * math_PI * 2 * 1 * 1;
    float3 MOI = make_float3(1. / 5. * mass * (1 * 1 + 2 * 2), 1. / 5. * mass * (1 * 1 + 2 * 2),
                             1. / 5. * mass * (1 * 1 + 1 * 1));
    // We can scale this general template to make it smaller, like a DEM particle that you would actually use
    float scaling = 0.003;  // Pretty much 6 mm particle
    std::shared_ptr<DEMClumpTemplate> my_template =
        DEMSim.LoadClumpType(mass, MOI, GetDEMEDataFile("clumps/ellipsoid_2_1_1.csv"), mat_type_particles);
    my_template->Scale(scaling);

    auto mat_type_bucket = DEMSim.LoadMaterial({{"E", 1e8}, {"nu", 0.3}, {"CoR", 0.3}, {"mu", 0.5}, {"Crr", 0.0}});
    DEMSim.SetMaterialPropertyPair("CoR", mat_type_particles, mat_type_bucket, 0.0);
    DEMSim.SetMaterialPropertyPair("Crr", mat_type_particles, mat_type_bucket, 0.0);
    // Add the excavator mesh
    auto excavator = DEMSim.AddWavefrontMeshObject(GetDEMEDataFile("mesh/komatsu_bucket_2.obj"), mat_type_walls);
    excavator->Move(make_float3(0, 0, 0), make_float4(0, 0, 0, 1));
    excavator->Scale(1. / 15);
    float bucket_height_approx = 0.105;  // This is half height
    float3 init_pos = make_float3(-0.1, 0, -world_halfsize + bucket_height_approx + 2 * scaling);
    float4 init_Q = make_float4(0.7071, 0, 0, 0.7071);  // 90 deg about x
    excavator->SetInitPos(init_pos);
    excavator->SetInitQuat(init_Q);
    excavator->SetFamily(10);

    // auto excavator = DEMSim.AddWavefrontMeshObject(GetDEMEDataFile("mesh/excavator.obj"), mat_type_bucket);
    // excavator->Move(make_float3(0, 0, 0), make_float4(0, 0, 0, 1));
    // excavator->Scale(1. / 200);
    // float bucket_height_approx = 0.105;  // This is half height
    // float3 init_pos = make_float3(-0.1, 0, -0.2);
    // float4 init_Q = make_float4(0, 0, 0.7071, 0.7071);  // 90 deg about x
    // excavator->SetInitPos(init_pos);
    // excavator->SetInitQuat(init_Q);
    // excavator->SetFamily(10);

    // DEMSim.DisableContactBetweenFamilies(0, 10);
    // DEMSim.EnableContactBetweenFamilies(0, 10);

    DEMSim.SetFamilyFixed(10);

    // Excavator

    // Generate initial clumps for piling
    float spacing = 2. * scaling;
    float fill_halfwidth = world_halfsize - 4. * scaling;
    float fill_height = 0.173 * 1.5 - 3 * scaling;
    float height_needed = 0.173;
    float fill_bottom = bowl_bottom + 3. * scaling;
    float max_length_at_layer;
    float max_width_at_layer;
    float tan_alpha = tan(30. * 3.14 / 180.);
    PDSampler sampler(spacing);
    // Use a PDSampler-based clump generation process. For PD sampler it is better to do it layer by layer.
    std::vector<float3> input_pile_xyz;
    float layer_z = 0;
    float theoritical_layer_z = 0;
    while (theoritical_layer_z < height_needed) {
        max_width_at_layer = fill_halfwidth;
        max_length_at_layer = (height_needed - theoritical_layer_z) / (tan_alpha);
        std::cout << max_length_at_layer << std::endl;
        float3 sample_center = make_float3(0.5, 0, fill_bottom + layer_z);
        auto layer_xyz = sampler.SampleBox(sample_center, make_float3(max_length_at_layer, max_width_at_layer, 0));
        // for (const auto& point : layer_xyz) {
        //     std::cout << "Layer Point: (" << point.x << ", " << point.y << ", " << point.z << ")" << std::endl;
        // }
        for (auto& xyz : layer_xyz) {
            if (xyz.x < (world_halfsize - 4 * scaling)) {
                input_pile_xyz.push_back(xyz);
            }
        }

        layer_z += 4.5 * scaling;
        theoritical_layer_z += scaling;
    }
    // Note: AddClumps can be called multiple times before initialization to add more clumps to the system.
    auto the_pile = DEMSim.AddClumps(my_template, input_pile_xyz);
    the_pile->SetFamily(0);
    std::cout << "Total Number of Particles: " << input_pile_xyz.size() << std::endl;

    // Excavator moves along x axis with velocity 0.1 m/s
    DEMSim.SetFamilyPrescribedLinVel(1, "0.1", "0", "0");

    // Excavator rotates about its own frame with angular velocity pi/8 rad/s
    DEMSim.SetFamilyPrescribedAngVel(2, "0", "0", "3.14 / 8");
    // At the same time, the excavator also moves at an angle of x degrees with respect to the x axis with 0.05 m/s
    // magnitude
    DEMSim.SetFamilyPrescribedLinVel(2, "0.045", "0", "0.065");

    // Go back to the original position
    DEMSim.SetFamilyPrescribedLinVel(3, "-0.1", "0", "0");
    DEMSim.SetFamilyPrescribedAngVel(3, "0", "0", "0");

    float step_size = 5e-6;
    DEMSim.InstructBoxDomainDimension({-world_halfsize, world_halfsize}, {-world_halfsize, world_halfsize},
                                      {-world_halfsize, world_halfsize});
    DEMSim.InstructBoxDomainBoundingBC("top_open", mat_type_walls);
    DEMSim.SetInitTimeStep(step_size);
    DEMSim.SetGravitationalAcceleration(make_float3(0, 0, -9.81));
    DEMSim.Initialize();

    path out_dir = current_path();
    out_dir += "/DemoOutput_Komatsu1";
    create_directory(out_dir);

    unsigned int fps = 20;
    unsigned int out_steps = (unsigned int)(1.0 / (fps * step_size));

    std::cout << "Output at " << fps << " FPS" << std::endl;
    unsigned int currframe = 0;
    unsigned int curr_step = 0;

    // Settle
    AdvanceSimulation(DEMSim, 1., step_size, out_steps, out_dir, curr_step, currframe);
    // DEMSim.ChangeFamily(10, 1);
    // DEMSim.EnableContactBetweenFamilies(0, 10);

    std::chrono::high_resolution_clock::time_point start = std::chrono::high_resolution_clock::now();

    // This is the family that moves the bucket along x axis with velocity 0.1 m/s
    // We do this for 3 seconds
    // DEMSim.EnableContactBetweenFamilies(0, 10);
    DEMSim.ChangeFamily(10, 1);

    AdvanceSimulation(DEMSim, 3, step_size, out_steps, out_dir, curr_step, currframe);

    // Rotate the bucket and move it upwards
    DEMSim.ChangeFamily(1, 2);
    AdvanceSimulation(DEMSim, 1.5, step_size, out_steps, out_dir, curr_step, currframe);

    // Go back to the original position
    DEMSim.ChangeFamily(2, 3);
    AdvanceSimulation(DEMSim, 3, step_size, out_steps, out_dir, curr_step, currframe);

    // Then rest in-place for a while
    DEMSim.ChangeFamily(3, 10);
    AdvanceSimulation(DEMSim, 0.5, step_size, out_steps, out_dir, curr_step, currframe);

    std::chrono::high_resolution_clock::time_point end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> time_sec = std::chrono::duration_cast<std::chrono::duration<double>>(end - start);
    std::cout << time_sec.count() << " seconds (wall time) to finish the 10-second plowing simulation." << std::endl;

    DEMSim.ShowTimingStats();
    DEMSim.ClearTimingStats();

    std::cout << "DEMdemo_Komatsu exiting..." << std::endl;
    return 0;
}
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

    auto mat_type_walls = DEMSim.LoadMaterial({{"E", 1e8}, {"nu", 0.3}, {"CoR", 0.3}, {"mu", 0.5}});
    auto mat_type_particles = DEMSim.LoadMaterial({{"E", 1e8}, {"nu", 0.3}, {"CoR", 0.}, {"mu", 0.58}});
    // auto mat_type_particles = DEMSim.LoadMaterial({{"E", 1e7}, {"nu", 0.3}, {"CoR", 0.3}, {"mu", 0.5}});
    // If you don't have this line, then CoR between wall material and granular material will be 0.5 (average of the
    // two).
    DEMSim.SetMaterialPropertyPair("CoR", mat_type_walls, mat_type_particles, 0.);
    DEMSim.SetMaterialPropertyPair("mu", mat_type_walls, mat_type_particles, 0.5);

    // Define the terrain particle templates
    // Calculate its mass and MOI
    float mass = 2.6e3 * 4. / 3. * math_PI * 2 * 1 * 1;
    float3 MOI = make_float3(1. / 5. * mass * (1 * 1 + 2 * 2), 1. / 5. * mass * (1 * 1 + 2 * 2),
                             1. / 5. * mass * (1 * 1 + 1 * 1));
    // We can scale this general template to make it smaller, like a DEM particle that you would actually use
    float scaling = 0.003; // Pretty much 5 mm particle
    std::shared_ptr<DEMClumpTemplate> my_template =
        DEMSim.LoadClumpType(mass, MOI, GetDEMEDataFile("clumps/ellipsoid_2_1_1.csv"), mat_type_particles);
    my_template->Scale(scaling);

    // Generate initial clumps for piling
    float spacing = 2. * scaling;
    float fill_halfwidth = world_halfsize - 4. * scaling;
    float fill_height = world_halfsize * 1.5;
    float fill_bottom = bowl_bottom + 3. * scaling;
    float max_length_at_layer;
    float max_width_at_layer;
    PDSampler sampler(spacing);
    // Use a PDSampler-based clump generation process. For PD sampler it is better to do it layer by layer.
    std::vector<float3> input_pile_xyz;
    float layer_z = 0;
    while (layer_z < fill_height) {
        max_length_at_layer = world_halfsize;
        max_width_at_layer = (fill_height - layer_z) / tan_alpha;
        float3 sample_center = make_float3(0.35, 0, fill_bottom + layer_z);
        auto layer_xyz = sampler.SampleBox(sample_center, make_float3(max_width_at_layer, max_length_at_layer, 0));

        for(auto& xyz : layer_xyz) {
            if(xyz.y <= max_width_at_layer){
                input_pile_xyz.push_back(xyz);
            }
        }

        layer_z += 4.5 * scaling;
    }
    // Note: AddClumps can be called multiple times before initialization to add more clumps to the system.
    auto the_pile = DEMSim.AddClumps(my_template, input_pile_xyz);
    the_pile->SetFamily(0);
    std::cout << "Total Number of Particles: "<<input_pile_xyz.size() << std::endl;

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
    DEMSim.DoDynamicsThenSync(1.5);

    std::chrono::high_resolution_clock::time_point start = std::chrono::high_resolution_clock::now();

    // Use 3 seconds to plow
    AdvanceSimulation(DEMSim, 3, step_size, out_steps, out_dir, curr_step, currframe);

    AdvanceSimulation(DEMSim, 2, step_size, out_steps, out_dir, curr_step, currframe);

    AdvanceSimulation(DEMSim, 2, step_size, out_steps, out_dir, curr_step, currframe);

    AdvanceSimulation(DEMSim, 3, step_size, out_steps, out_dir, curr_step, currframe);

    std::chrono::high_resolution_clock::time_point end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> time_sec = std::chrono::duration_cast<std::chrono::duration<double>>(end - start);
    std::cout << time_sec.count() << " seconds (wall time) to finish the 10-second plowing simulation." << std::endl;

    DEMSim.ShowTimingStats();
    DEMSim.ClearTimingStats();

    std::cout << "DEMdemo_Komatsu exiting..." << std::endl;
    return 0;
}

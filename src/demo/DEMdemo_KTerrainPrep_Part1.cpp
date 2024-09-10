//  Copyright (c) 2021, SBEL GPU Development Team
//  Copyright (c) 2021, University of Wisconsin - Madison
//
//	SPDX-License-Identifier: BSD-3-Clause

#include <core/ApiVersion.h>
#include <core/utils/ThreadManager.h>
#include <DEM/API.h>
#include <DEM/HostSideHelpers.hpp>
#include <DEM/utils/Samplers.hpp>

#include <cstdio>
#include <chrono>
#include <filesystem>
#include <map>
#include <random>
#include <cmath>

// =============================================================================
// In GRCPrep demo series, we try to prepare a sample of the GRC simulant, which
// are supposed to be used for extraterrestrial rover mobility simulations. It is
// made of particles of various sizes and shapes following a certain distribution.
// In Part1, it creates several batches of clumps and let them settle at the bottom
// of the domain.
// =============================================================================

using namespace deme;
using namespace std::filesystem;

int main() {
    DEMSolver DEMSim;
    DEMSim.SetVerbosity(INFO);
    DEMSim.SetOutputFormat(OUTPUT_FORMAT::CSV);
    // DEMSim.SetOutputContent(OUTPUT_CONTENT::FAMILY);
    DEMSim.SetOutputContent(OUTPUT_CONTENT::XYZ);

    srand(759);

    // Define materials
    auto mat_type_terrain = DEMSim.LoadMaterial({{"E", 1e9}, {"nu", 0.3}, {"CoR", 0.3}, {"mu", 0.5}});
    auto mat_type_wheel = DEMSim.LoadMaterial({{"E", 1e9}, {"nu", 0.3}, {"CoR", 0.3}, {"mu", 0.5}});

    // Define the simulation world
    double world_x_size = 1.0;
    double world_y_size = 1.0;
    double world_z_size = 2.0;
    float sampleheight = 0.15;      // for the generation of the random material
    double bottom = -sampleheight;
    float size_z_batch = 3 * sampleheight;
    float sample_halfheight = size_z_batch / 2;
    // DEMSim.InstructBoxDomainDimension(world_y_size, world_y_size, world_z_size);
    DEMSim.InstructBoxDomainDimension({-world_x_size / 1.5, world_x_size / 1.5},
                                      {-world_y_size / 1.5, world_y_size / 1.5}, {2.0 * bottom, world_z_size});

    auto top_plane = DEMSim.AddWavefrontMeshObject("../data/mesh/box.obj", mat_type_terrain);
    top_plane->SetInitPos(make_float3(0, 0, bottom - 0.010));
    top_plane->SetMass(1.);
    top_plane->Scale(make_float3(world_x_size / 2.0, world_y_size / 2.0, 2.0));
    top_plane->SetFamily(10);
    DEMSim.SetFamilyFixed(10);
    // Add 5 bounding planes around the simulation world, and leave the top open
    DEMSim.InstructBoxDomainBoundingBC("top_open", mat_type_terrain);

    DEMSim.AddBCPlane(make_float3(0, 0, bottom), make_float3(0, 0, 1), mat_type_terrain);

    // Shake the terrain to increase bulk density
    std::string shake_pattern_xz = " 0.002 * sin( 50 * 2 * deme::PI * t)";
    std::string shake_pattern_y = " 0.002 * sin( 200 * 2 * deme::PI * t)";
    DEMSim.SetFamilyPrescribedLinVel(9, shake_pattern_xz, shake_pattern_xz, shake_pattern_y);

    // Undersized has been skipped for now - Needs to be added back again sometime
    // Define the terrain particle templates
    // Calculate its mass and MOI
    float terrain_density = 2.6e3;
    float oversized_volume = 5.36207323370400;
    float oversized_mass = terrain_density * oversized_volume;
    float3 oversized_MOI = make_float3(1.444845, 5.4947616, 5.8418724) * terrain_density;
    float large_volume = 3.9052404;
    float large_mass = terrain_density * large_volume;
    float3 large_MOI = make_float3(1.0629883, 1.864263, 2.2812725) * terrain_density;
    float medium_volume = 3.570524727;
    float medium_mass = terrain_density * medium_volume;
    float3 medium_MOI = make_float3(0.9200936, 1.8532078, 2.0483853) * terrain_density;
    // float small_volume = 3.570524727;
    // float small_mass = terrain_density * small_volume;
    // float3 small_MOI = make_float3(0.9200936, 1.8532078, 2.0483853) * terrain_density;
    // float fine_volume = 3.570524727;
    // float fine_mass = terrain_density * fine_volume;
    // float3 fine_MOI = make_float3(0.9200936, 1.8532078, 2.0483853) * terrain_density;
    // Scale the template we just created
    std::vector<double> scales = {0.04040, 0.019392727, 0.007272727, 0.003636363, 0.001818182};
    // std::vector<double> scales = {0.007272727};
    // Then load it to system
    std::shared_ptr<DEMClumpTemplate> oversized_template = DEMSim.LoadClumpType(
        oversized_mass, oversized_MOI, GetDEMEDataFile("clumps/komatsu_oversized.csv"), mat_type_terrain);
    std::shared_ptr<DEMClumpTemplate> large_template =
        DEMSim.LoadClumpType(large_mass, large_MOI, GetDEMEDataFile("clumps/komatsu_large.csv"), mat_type_terrain);
    std::shared_ptr<DEMClumpTemplate> medium_template =
        DEMSim.LoadClumpType(medium_mass, medium_MOI, GetDEMEDataFile("clumps/komatsu_medium.csv"), mat_type_terrain);
    std::vector<std::shared_ptr<DEMClumpTemplate>> ground_particle_templates = {
        oversized_template, large_template, medium_template, DEMSim.Duplicate(medium_template),
        DEMSim.Duplicate(medium_template)};

    // std::vector<std::shared_ptr<DEMClumpTemplate>> ground_particle_templates = {medium_template};

    // Now scale those templates
    for (int i = 0; i < scales.size(); i++) {
        std::shared_ptr<DEMClumpTemplate>& my_template = ground_particle_templates.at(i);
        // Note the mass and MOI are also scaled in the process, automatically. But if you are not happy with this, you
        // can always manually change mass and MOI afterwards.
        my_template->Scale(scales.at(i));
        // Give these templates names, 0000, 0001 etc.
        char t_name[20];
        sprintf(t_name, "%04d", i);
        my_template->AssignName(std::string(t_name));
    }

    // Instatiate particles with a probability that is in line with their weight distribution.
    // For now, Undersized percentage has been added to fine
    std::vector<double> weight_perc = {0.05, 0.15, 0.25, 0.25, 0.3};
    // std::vector<double> weight_perc = {1};
    std::vector<double> grain_perc;
    for (int i = 0; i < scales.size(); i++) {
        grain_perc.push_back(weight_perc.at(i) / std::pow(scales.at(i), 3));
    }
    {
        double tmp = vector_sum(grain_perc);
        std::for_each(grain_perc.begin(), grain_perc.end(), [tmp](double& p) { p /= tmp; });
        std::cout << "Percentage of grains add up to " << vector_sum(grain_perc) << std::endl;
    }
    std::random_device r;
    std::default_random_engine e1(r());
    // Distribution that defines different weights (17, 10, etc.) for numbers.
    std::discrete_distribution<int> discrete_dist(grain_perc.begin(), grain_perc.end());

    // Sampler to use
    // HCPSampler sampler(scales.at(0) * 5);
    // These dimensions are according to the largest particle unfortunately
    GridSampler sampler(make_float3(0.18, 0.08, 0.07));
    // Make ready for simulation
    float step_size = 1e-6;
    DEMSim.SetInitTimeStep(step_size);
    DEMSim.SetGravitationalAcceleration(make_float3(0, 0, -9.81));
    // Max velocity info is generally just for the solver's reference and the user do not have to set it. The solver
    // wouldn't take into account a vel larger than this when doing async-ed contact detection: but this vel won't
    // happen anyway and if it does, something already went wrong.
    DEMSim.SetMaxVelocity(15.);
    // Error out vel is used to force the simulation to abort when something goes wrong.
    DEMSim.SetErrorOutVelocity(15.);
    DEMSim.SetExpandSafetyMultiplier(1.2);
    DEMSim.SetInitBinNumTarget(1e7);
    DEMSim.Initialize();

    float time_end = 30.0;

    path out_dir = current_path();
    out_dir += "/DemoOutput_KTerrainPrep_Part1";
    create_directory(out_dir);
    unsigned int currframe = 0;
    unsigned int curr_step = 0;

    float settletoleranceFactor = 0.18;
    
    // float targetMass = 0.60 * world_x_size * world_y_size * sampleheight * terrain_density;
    float sample_halfwidth_x = (world_y_size * 0.8) / 2;
    float sample_halfwidth_y = (world_y_size * 0.9) / 2;
    float offset_z = size_z_batch / 2 + settletoleranceFactor;
    float settle_frame_time = 0.1;
    float settle_batch_time = 0.1;

    auto max_z_finder = DEMSim.CreateInspector("clump_max_z");
    auto max_vel_finder = DEMSim.CreateInspector("clump_max_absv");
    auto totalMass = DEMSim.CreateInspector("clump_mass");
    float height = 0;
    float current_z = bottom;
    float maxvel = max_vel_finder->GetValue();
    float mass = totalMass->GetValue();

    std::chrono::high_resolution_clock::time_point start = std::chrono::high_resolution_clock::now();
    bool consolidationEnds = false;
    while (!consolidationEnds) {
        if (current_z < (world_z_size - 2 * offset_z)) {
            // DEMSim.ClearCache(); // Clearing cache is no longer needed
            float3 sample_center = make_float3(0, 0, offset_z + current_z);
            std::cout << "Sample center: " << sample_center.x << ", " << sample_center.y << ", " << sample_center.z
                      << std::endl;
            std::vector<std::shared_ptr<DEMClumpTemplate>> heap_template_in_use;
            std::vector<unsigned int> heap_family;
            // Sample and add heap particles
            auto heap_particles_xyz = sampler.SampleBox(
                sample_center, make_float3(sample_halfwidth_x, sample_halfwidth_y, sample_halfheight));
            for (unsigned int i = 0; i < heap_particles_xyz.size(); i++) {
                int ind = std::round(discrete_dist(e1));
                heap_template_in_use.push_back(ground_particle_templates.at(ind));
                heap_family.push_back(ind);
            }
            auto heap_particles = DEMSim.AddClumps(heap_template_in_use, heap_particles_xyz);
            // Give ground particles a small initial velocity so they `collapse' at the start of the simulation
            heap_particles->SetVel(make_float3(0.00, 0, -0.05));
            heap_particles->SetFamilies(heap_family);
            DEMSim.UpdateClumps();
            std::cout << "Current number of clumps: " << DEMSim.GetNumClumps() << std::endl;
        }

        // Allow for some settling
        // Must DoDynamicsThenSync (not DoDynamics), as adding entities to the simulation is only allowed at a sync-ed
        // point of time.

        std::cout << "Frame: " << currframe << std::endl;
        char filename[200];
        sprintf(filename, "%s/DEMdemo_output_%04d.csv", out_dir.c_str(), currframe);
        DEMSim.WriteSphereFile(std::string(filename));
        DEMSim.DoDynamicsThenSync(settle_frame_time);
        current_z = max_z_finder->GetValue();
        maxvel = max_vel_finder->GetValue();
        mass = totalMass->GetValue();

        std::cout << "Total mass: " << mass << std::endl;

        if (DEMSim.GetNumClumps() > 0.25e6) {
            consolidationEnds = true;
            DEMSim.ChangeFamily(10, 9);
            std::cout << "Consolidating for one second" << std::endl;
            auto max_z_finder = DEMSim.CreateInspector("clump_max_z");
            auto min_z_finder = DEMSim.CreateInspector("clump_min_z");
            for (float t = 0; t < 1.00; t += 0.025) {
                std::cout << "Frame: " << currframe << std::endl;
                char filename[200];
                sprintf(filename, "%s/DEMdemo_output_%04d.csv", out_dir.c_str(), currframe);
                DEMSim.WriteSphereFile(std::string(filename));

                currframe++;
                float terrain_max_z = max_z_finder->GetValue();
                float terrain_min_z = min_z_finder->GetValue();
                std::cout << "Consolidation: " << terrain_max_z - terrain_min_z << std::endl;
                DEMSim.DoDynamics(0.025);
            }
            break;
        }
        currframe++;
        std::cout << "Current z: " << current_z << " with max vel: " << maxvel << std::endl;

        DEMSim.ShowThreadCollaborationStats();
    }

    // Settle for some time more
    DEMSim.DoDynamicsThenSync(1.0);
    std::chrono::high_resolution_clock::time_point end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> time_sec = std::chrono::duration_cast<std::chrono::duration<double>>(end - start);
    std::cout << time_sec.count() << " seconds (wall time) to finish the simulation" << std::endl;

    char cp_filename[200];
    sprintf(cp_filename, "%s/KTerrain_250k.csv", out_dir.c_str());
    DEMSim.WriteClumpFile(std::string(cp_filename));

    DEMSim.ShowThreadCollaborationStats();
    DEMSim.ClearThreadCollaborationStats();

    char cnt_filename[200];
    sprintf(cnt_filename, "%s/KTerrain_Contact_pairs_250k.csv", out_dir.c_str());
    DEMSim.WriteContactFile(std::string(cnt_filename));

    std::cout << "DEMdemo_KTerrainPrep_Part1 exiting..." << std::endl;
    return 0;
}

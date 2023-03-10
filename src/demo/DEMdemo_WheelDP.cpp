//  Copyright (c) 2021, SBEL GPU Development Team
//  Copyright (c) 2021, University of Wisconsin - Madison
//
//	SPDX-License-Identifier: BSD-3-Clause

// =============================================================================
// NOTE!!! You have to first finish ALL 3 GRCPrep demos to obtain the GRC_20e6.csv
// file, put it in the current directory, then run this demo.
// A wheel drawbar-pull test, featuring Curiosity wheel geometry and GRC-1 simulant.
// =============================================================================

#include <DEM/API.h>
#include <DEM/HostSideHelpers.hpp>
#include <DEM/utils/Samplers.hpp>

#include <chrono>
#include <cmath>
#include <cstdio>
#include <filesystem>
#include <map>
#include <random>

using namespace deme;

const double math_PI = 3.1415927;

int main() {
    std::filesystem::path out_dir = std::filesystem::current_path();
    out_dir += "/DemoOutput_WheelDP";
    std::filesystem::create_directory(out_dir);

    float TRs[] = {0.05, 0.1, 0.2, 0.3, 0.5, 0.7, 0.8};
    unsigned int run_mode = 0;
    unsigned int currframe = 0;

    for (float TR : TRs) {
        DEMSolver DEMSim;
        DEMSim.SetVerbosity(INFO);
        DEMSim.SetOutputFormat(OUTPUT_FORMAT::CSV);
        DEMSim.SetOutputContent(OUTPUT_CONTENT::ABSV);
        DEMSim.SetMeshOutputFormat(MESH_FORMAT::VTK);
        DEMSim.SetContactOutputContent(OWNER | FORCE | POINT);

        // E, nu, CoR, mu, Crr...
        auto mat_type_wheel = DEMSim.LoadMaterial({{"E", 1e9}, {"nu", 0.3}, {"CoR", 0.3}, {"mu", 0.5}, {"Crr", 0.00}});
        auto mat_type_terrain =
            DEMSim.LoadMaterial({{"E", 1e9}, {"nu", 0.3}, {"CoR", 0.3}, {"mu", 0.5}, {"Crr", 0.00}});

        // `World'
        float G_mag = 9.81;
        float step_size = 1e-6;
        double world_size_y = 0.52;
        double world_size_x = 1.53;
        double world_size_z = 4.0;
        DEMSim.InstructBoxDomainDimension(world_size_x, world_size_y, world_size_z);
        DEMSim.InstructBoxDomainBoundingBC("top_open", mat_type_terrain);
        float bottom = -0.5;
        auto bot_wall = DEMSim.AddBCPlane(make_float3(0, 0, bottom), make_float3(0, 0, 1), mat_type_terrain);
        auto bot_wall_tracker = DEMSim.Track(bot_wall);

        // Define the wheel geometry
        float wheel_rad = 0.25;
        float wheel_width = 0.25;
        float wheel_mass = 480.0 / G_mag;
        float wheel_IYY = wheel_mass * wheel_rad * wheel_rad / 2;
        float wheel_IXX = (wheel_mass / 12) * (3 * wheel_rad * wheel_rad + wheel_width * wheel_width);
        auto wheel =
            DEMSim.AddWavefrontMeshObject(GetDEMEDataFile("mesh/rover_wheels/curiosity_wheel.obj"), mat_type_wheel);
        wheel->SetMass(wheel_mass);
        wheel->SetMOI(make_float3(wheel_IXX, wheel_IYY, wheel_IXX));
        // Give the wheel a family number so we can potentially add prescription
        wheel->SetFamily(1);
        // Track it
        auto wheel_tracker = DEMSim.Track(wheel);

        // Define the terrain particle templates
        // Calculate its mass and MOI
        float terrain_density = 2.6e3;
        double clump_vol = 5.5886717;
        float mass = terrain_density * clump_vol;
        float3 MOI = make_float3(2.928, 2.6029, 3.9908) * terrain_density;
        // Then load it to system
        std::shared_ptr<DEMClumpTemplate> my_template =
            DEMSim.LoadClumpType(mass, MOI, GetDEMEDataFile("clumps/triangular_flat.csv"), mat_type_terrain);
        my_template->SetVolume(clump_vol);
        // Make 5 copies. Note we must use DEME's duplicate method to do this, because we will make changes to the
        // templates later, using the shared_ptrs as handles. If no duplications are made, then all the changes are
        // going to be enforced on the same template, and in the end we'd not be able to get 5 distinct templates.
        std::vector<std::shared_ptr<DEMClumpTemplate>> ground_particle_templates = {
            my_template, DEMSim.Duplicate(my_template), DEMSim.Duplicate(my_template), DEMSim.Duplicate(my_template),
            DEMSim.Duplicate(my_template)};
        // Decide the scalings of the templates we just created (so that they are... like particles, not rocks)
        std::vector<double> scales = {0.00063, 0.00033, 0.00022, 0.00015, 0.00009};
        std::for_each(scales.begin(), scales.end(), [](double& r) { r *= 20.; });
        // Now scale those templates
        for (int i = 0; i < scales.size(); i++) {
            std::shared_ptr<DEMClumpTemplate>& my_template = ground_particle_templates.at(i);
            // Note the mass and MOI are also scaled in the process, automatically. But if you are not happy with this,
            // you can always manually change mass and MOI afterwards.
            my_template->Scale(scales.at(i));
            // Give these templates names, 0000, 0001 etc.
            char t_name[20];
            sprintf(t_name, "%04d", i);
            my_template->AssignName(std::string(t_name));
        }

        // Now we load clump locations from a checkpointed file
        {
            std::cout << "Making terrain..." << std::endl;
            std::unordered_map<std::string, std::vector<float3>> clump_xyz;
            std::unordered_map<std::string, std::vector<float4>> clump_quaternion;
            try {
                clump_xyz = DEMSim.ReadClumpXyzFromCsv("./GRC_20e6.csv");
                clump_quaternion = DEMSim.ReadClumpQuatFromCsv("./GRC_20e6.csv");
            } catch (...) {
                std::cout << "You will need to finish the GRCPrep demos first to obtain the checkpoint file "
                             "GRC_20e6.csv, in order to run this demo. This file is needed to generate the terrain bed."
                          << std::endl;
                return 1;
            }
            std::vector<float3> in_xyz;
            std::vector<float4> in_quat;
            std::vector<std::shared_ptr<DEMClumpTemplate>> in_types;
            unsigned int t_num = 0;
            for (int i = 0; i < scales.size(); i++) {
                char t_name[20];
                sprintf(t_name, "%04d", t_num);

                auto this_type_xyz = clump_xyz[std::string(t_name)];
                auto this_type_quat = clump_quaternion[std::string(t_name)];

                size_t n_clump_this_type = this_type_xyz.size();
                std::cout << "Loading clump " << std::string(t_name) << " which has particle num: " << n_clump_this_type
                          << std::endl;
                // Prepare clump type identification vector for loading into the system (don't forget type 0 in
                // ground_particle_templates is the template for rover wheel)
                std::vector<std::shared_ptr<DEMClumpTemplate>> this_type(n_clump_this_type,
                                                                         ground_particle_templates.at(t_num));

                // Add them to the big long vector
                in_xyz.insert(in_xyz.end(), this_type_xyz.begin(), this_type_xyz.end());
                in_quat.insert(in_quat.end(), this_type_quat.begin(), this_type_quat.end());
                in_types.insert(in_types.end(), this_type.begin(), this_type.end());
                std::cout << "Added clump type " << t_num << std::endl;
                // Our template names are 0000, 0001 etc.
                t_num++;
            }

            // Now, we don't need all particles loaded...
            std::vector<notStupidBool_t> elem_to_remove(in_xyz.size(), 0);
            for (size_t i = 0; i < in_xyz.size(); i++) {
                if (std::abs(in_xyz.at(i).y) > (world_size_y - 0.04) / 2 ||
                    std::abs(in_xyz.at(i).x) > (world_size_x - 0.08) / 2)
                    elem_to_remove.at(i) = 1;
            }
            in_xyz.erase(std::remove_if(in_xyz.begin(), in_xyz.end(),
                                        [&elem_to_remove, &in_xyz](const float3& i) {
                                            return elem_to_remove.at(&i - in_xyz.data());
                                        }),
                         in_xyz.end());
            in_quat.erase(std::remove_if(in_quat.begin(), in_quat.end(),
                                         [&elem_to_remove, &in_quat](const float4& i) {
                                             return elem_to_remove.at(&i - in_quat.data());
                                         }),
                          in_quat.end());
            in_types.erase(std::remove_if(in_types.begin(), in_types.end(),
                                          [&elem_to_remove, &in_types](const auto& i) {
                                              return elem_to_remove.at(&i - in_types.data());
                                          }),
                           in_types.end());
            DEMClumpBatch base_batch(in_xyz.size());
            base_batch.SetTypes(in_types);
            base_batch.SetPos(in_xyz);
            base_batch.SetOriQ(in_quat);

            /*
            std::vector<float> x_shift_dist = {-0.5, 0.5};
            std::vector<float> y_shift_dist = {0};
            // Add some patches of such graular bed
            for (float x_shift : x_shift_dist) {
                for (float y_shift : y_shift_dist) {
                    DEMClumpBatch another_batch = base_batch;
                    std::vector<float3> my_xyz = in_xyz;
                    std::for_each(my_xyz.begin(), my_xyz.end(), [x_shift, y_shift](float3& xyz) {
                        xyz.x += x_shift;
                        xyz.y += y_shift;
                    });
                    another_batch.SetPos(my_xyz);
                    DEMSim.AddClumps(another_batch);
                }
            }
            */

            DEMSim.AddClumps(base_batch);
        }

        // Now add a plane to compress the sample
        // auto compressor = DEMSim.AddExternalObject();
        // compressor->AddPlane(make_float3(0, 0, 0), make_float3(0, 0, -1), mat_type_terrain);
        // compressor->SetFamily(2);
        // auto compressor_tracker = DEMSim.Track(compressor);

        // Families' prescribed motions
        float w_r = math_PI / 12;
        float v_ref = w_r * wheel_rad;

        double sim_end = 6.;
        // Note: this wheel is not `dictated' by our prescrption of motion because it can still fall onto the ground
        // (move freely linearly)
        DEMSim.SetFamilyPrescribedAngVel(1, "0", to_string_with_precision(w_r), "0", false);

        // `Real sim' family number
        DEMSim.SetFamilyPrescribedAngVel(2, "0", to_string_with_precision(w_r), "0", false);
        // Note: this wheel is not `dictated' by our prescrption of motion (hence the false argument), because it
        // can sink into the ground (move on Z dir); but its X and Y motions are explicitly controlled.
        DEMSim.SetFamilyPrescribedLinVel(2, to_string_with_precision(v_ref * (1. - TR)), "0", "none", false);

        // Some inspectors
        auto max_z_finder = DEMSim.CreateInspector("clump_max_z");
        auto min_z_finder = DEMSim.CreateInspector("clump_min_z");
        auto total_mass_finder = DEMSim.CreateInspector("clump_mass");
        auto max_v_finder = DEMSim.CreateInspector("clump_max_absv");

        DEMSim.SetInitTimeStep(step_size);
        DEMSim.SetGravitationalAcceleration(make_float3(0, 0, -G_mag));
        DEMSim.SetCDUpdateFreq(30);
        // Max velocity info is generally just for the solver's reference and the user do not have to set it. The solver
        // wouldn't take into account a vel larger than this when doing async-ed contact detection: but this vel won't
        // happen anyway and if it does, something already went wrong.
        DEMSim.SetMaxVelocity(50.);
        DEMSim.SetExpandSafetyMultiplier(1.1);
        DEMSim.SetInitBinSize(2 * scales.at(2));
        DEMSim.Initialize();

        // Compress until dense enough
        unsigned int curr_step = 0;
        unsigned int fps = 10;
        unsigned int out_steps = (unsigned int)(1.0 / (fps * step_size));
        double frame_time = 1.0 / fps;
        unsigned int report_ps = 1000;
        unsigned int report_steps = (unsigned int)(1.0 / (report_ps * step_size));
        std::cout << "Output at " << fps << " FPS" << std::endl;

        // Put the wheel in place, then let the wheel sink in initially
        float max_z = max_z_finder->GetValue();
        wheel_tracker->SetPos(make_float3(-0.45, 0, max_z + 0.02 + wheel_rad));
        for (double t = 0; t < 0.5; t += frame_time) {
            char filename[200], meshname[200];
            std::cout << "Outputting frame: " << currframe << std::endl;
            sprintf(filename, "%s/DEMdemo_output_%04d.csv", out_dir.c_str(), currframe);
            sprintf(meshname, "%s/DEMdemo_mesh_%04d.vtk", out_dir.c_str(), currframe++);
            DEMSim.WriteSphereFile(std::string(filename));
            DEMSim.WriteMeshFile(std::string(meshname));

            // float3 pos = bot_wall_tracker->Pos();
            // float3 force = bot_wall_tracker->ContactAcc();
            // std::cout << "Bottom wall pos: " << pos.x << ", " << pos.y << ", " << pos.z << std::endl;
            // std::cout << "Bottom wall force: " << force.x << ", " << force.y << ", " << force.z << std::endl;
            DEMSim.DoDynamicsThenSync(frame_time);
        }

        // Switch wheel from free fall into DP test
        DEMSim.ChangeFamily(1, 2);
        DEMSim.UpdateStepSize(step_size * 2);

        for (double t = 0; t < sim_end; t += step_size, curr_step++) {
            if (curr_step % out_steps == 0) {
                char filename[200], meshname[200];
                std::cout << "Outputting frame: " << currframe << std::endl;
                sprintf(filename, "%s/DEMdemo_output_%04d.csv", out_dir.c_str(), currframe);
                sprintf(meshname, "%s/DEMdemo_mesh_%04d.vtk", out_dir.c_str(), currframe++);
                DEMSim.WriteSphereFile(std::string(filename));
                DEMSim.WriteMeshFile(std::string(meshname));
                DEMSim.ShowThreadCollaborationStats();
            }

            if (curr_step % report_steps == 0) {
                float3 forces = wheel_tracker->ContactAcc();
                forces *= wheel_mass;
                std::cout << "Current run mode: " << run_mode << std::endl;
                std::cout << "Time: " << t << std::endl;
                std::cout << "Force on wheel: " << forces.x << ", " << forces.y << ", " << forces.z << std::endl;
                std::cout << "Drawbar pull coeff: " << forces.x / (wheel_mass * G_mag) << std::endl;
                std::cout << "Max system velocity: " << max_v_finder->GetValue() << std::endl;
            }

            DEMSim.DoDynamics(step_size);
        }

        run_mode++;
        DEMSim.ShowTimingStats();
    }

    std::cout << "WheelDP demo exiting..." << std::endl;
    return 0;
}
//  Copyright (c) 2021, SBEL GPU Development Team
//  Copyright (c) 2021, University of Wisconsin - Madison
//
//      SPDX-License-Identifier: BSD-3-Clause

// =============================================================================
// NOTE!!! You have to first finish ALL GRCPrep demos to obtain the GRC_3e6.csv
// file, put it in the current directory, then run this demo.
// A slip v.s. slope test, featuring a simplified Viper wheel and GRC-1 simulant.
// The wheel in this simulation is represented by a collection of spheres (like
// a ``big'' clump), unlike in WheelDP where the wheel is mesh.
// WARNING: This is a huge simulation with millions of particles.
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

int main(int argc, char* argv[]) {
    float ang_velo_list[] = {0.70, 1.39, 2.09, 2.78, 4.17};
    float lin_velo_list[] = {0.05, 0.1, 0.15, 0.2, 0.3};

    if (argc != 2) {
        std::cout << "useage: ./WheelSlopeSlip [testID 1-5]" << std::endl;
        return -1;
    }

    int TestID = std::atoi(argv[1]);

    std::filesystem::path out_dir = std::filesystem::current_path();
    out_dir += "/rassor_single_drum_added_weight_Test" + std::to_string(TestID);
    std::filesystem::create_directory(out_dir);
    std::cout << "output directory: " << out_dir << std::endl;
    // `World'
    float G_mag = 9.81;
    float step_size = 5e-6;
    double world_size_y = 0.5;
    double world_size_x = 2.04;
    double world_size_z = 4.0;

    // Define the wheel geometry
    float wheel_rad = 0.2;
    float wheel_width = 0.205;
    float wheel_mass = 2.5;
    float total_pressure = 7.5 * G_mag;
    float added_pressure = (total_pressure - wheel_mass * G_mag);
    float wheel_IYY = 0.02;
    float wheel_IXX = 0.0058;

    //    float Slopes_deg[] = {0, 2.5, 5, 7.5, 10};
    float Slopes_deg[] = {0};

    unsigned int run_mode = 0;
    unsigned int currframe = 0;

    for (float Slope_deg : Slopes_deg) {
        DEMSolver DEMSim;
        DEMSim.SetVerbosity(INFO);
        DEMSim.SetOutputFormat(OUTPUT_FORMAT::CSV);
        DEMSim.SetOutputContent(OUTPUT_CONTENT::ABSV);
        DEMSim.SetMeshOutputFormat(MESH_FORMAT::VTK);
        // DEMSim.SetContactOutputContent(OWNER | FORCE | POINT);

        // E, nu, CoR, mu, Crr...
        auto mat_type_wheel = DEMSim.LoadMaterial({{"E", 5e7}, {"nu", 0.3}, {"CoR", 0.5}, {"mu", 0.5}, {"Crr", 0.00}});
        auto mat_type_terrain =
            DEMSim.LoadMaterial({{"E", 5e7}, {"nu", 0.3}, {"CoR", 0.5}, {"mu", 0.5}, {"Crr", 0.00}});
        // If you don't have this line, then mu between drum material and granular material will be the average of the
        // two.
        DEMSim.SetMaterialPropertyPair("mu", mat_type_wheel, mat_type_terrain, 0.8);

        DEMSim.InstructBoxDomainDimension(world_size_x, world_size_y, world_size_z);
        DEMSim.InstructBoxDomainBoundingBC("top_open", mat_type_terrain);
        float bottom = -0.5;
        auto bot_wall = DEMSim.AddBCPlane(make_float3(0, 0, bottom), make_float3(0, 0, 1), mat_type_terrain);
        auto bot_wall_tracker = DEMSim.Track(bot_wall);

        auto wheel_template = DEMSim.LoadClumpType(wheel_mass, make_float3(wheel_IXX, wheel_IYY, wheel_IXX),
                                                   GetDEMEDataFile("clumps/RassorWheelSimple.csv"), mat_type_wheel);

        std::cout << "successfully loaded wheel clump\n" << std::endl;
        // The file contains no wheel particles size info, so let's manually set them
        double wheel_particle_radius = 0.0025;
        wheel_template->radii = std::vector<float>(wheel_template->nComp, wheel_particle_radius);
        // Instantiate this wheel
        auto wheel = DEMSim.AddClumps(wheel_template, make_float3(0));
        // Give the wheel a family number so we can potentially add prescription
        wheel->SetFamily(10);
        // Track it
        auto wheel_tracker = DEMSim.Track(wheel);

        // Define the terrain particle templates
        // Calculate its mass and MOI
        float terrain_density = 2.6e3;
        float volume1 = 4.2520508;
        float mass1 = terrain_density * volume1;
        float3 MOI1 = make_float3(1.6850426, 1.6375114, 2.1187753) * terrain_density;
        float volume2 = 2.1670011;
        float mass2 = terrain_density * volume2;
        float3 MOI2 = make_float3(0.57402126, 0.60616378, 0.92890173) * terrain_density;
        // Scale the template we just created
        std::vector<double> scales = {0.014, 0.0075833, 0.0044, 0.003, 0.002, 0.0018333, 0.0017};
        // Then load it to system
        std::shared_ptr<DEMClumpTemplate> my_template2 =
            DEMSim.LoadClumpType(mass2, MOI2, GetDEMEDataFile("clumps/triangular_flat_6comp.csv"), mat_type_terrain);
        std::shared_ptr<DEMClumpTemplate> my_template1 =
            DEMSim.LoadClumpType(mass1, MOI1, GetDEMEDataFile("clumps/triangular_flat.csv"), mat_type_terrain);
        std::vector<std::shared_ptr<DEMClumpTemplate>> ground_particle_templates = {my_template2,
                                                                                    DEMSim.Duplicate(my_template2),
                                                                                    my_template1,
                                                                                    DEMSim.Duplicate(my_template1),
                                                                                    DEMSim.Duplicate(my_template1),
                                                                                    DEMSim.Duplicate(my_template1),
                                                                                    DEMSim.Duplicate(my_template1)};
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
                // from GRC part 2
                std::string soil_input = "GRC_3e6.csv";
                clump_xyz = DEMSim.ReadClumpXyzFromCsv(soil_input);
                clump_quaternion = DEMSim.ReadClumpQuatFromCsv(soil_input);
            } catch (...) {
                std::cout << "You will need to finish the GRCPrep demos first to obtain the checkpoint file "
                             "GRC_3e6.csv, in order to run this demo. This file is needed to generate the terrain bed."
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
                if (std::abs(in_xyz.at(i).y) > (world_size_y - 0.1) / 2.)
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

            // Make 2 copies of the batch of particles, then feed into simulation.
            std::vector<float> x_shift_dist = {-0.5, 0.5};
            std::vector<float> y_shift_dist = {0};
            // Add some patches of such graular bed
            for (float x_shift : x_shift_dist) {
                for (float y_shift : y_shift_dist) {
                    DEMClumpBatch batch_to_add = base_batch;
                    std::vector<float3> my_xyz = in_xyz;
                    std::for_each(my_xyz.begin(), my_xyz.end(), [x_shift, y_shift](float3& xyz) {
                        xyz.x += x_shift;
                        xyz.y += y_shift;
                    });
                    batch_to_add.SetPos(my_xyz);
                    DEMSim.AddClumps(batch_to_add);
                }
            }
        }

        float motor_ang_velo = ang_velo_list[TestID - 1];
        float motor_lin_velo = lin_velo_list[TestID - 1];
        double G_ang = Slope_deg * math_PI / 180.;
        std::cout << "anglar velo: " << motor_ang_velo << ", lin velo: " << motor_lin_velo << std::endl;
        float travel_dist = (world_size_x - 2. * wheel_rad);
        std::cout << "travel dist: " << travel_dist << " m" << std::endl;
        float travel_time = travel_dist / motor_lin_velo;

        double sim_end = std::ceil(travel_time);
        std::cout << "simulation time: " << sim_end << std::endl;
        // Note: this wheel is not `dictated' by our prescrption of motion because it can still fall onto the ground
        // (move freely linearly)
        DEMSim.SetFamilyPrescribedAngVel(1, "0", to_string_with_precision(motor_ang_velo), "0", false);
        DEMSim.SetFamilyPrescribedLinVel(1, to_string_with_precision(motor_lin_velo), "0", "none", false);
        DEMSim.AddFamilyPrescribedAcc(1, to_string_with_precision(-added_pressure * std::sin(G_ang) / wheel_mass),
                                      "none", to_string_with_precision(-added_pressure * std::cos(G_ang) / wheel_mass));
        DEMSim.SetFamilyFixed(10);

        // Some inspectors
        auto max_z_finder = DEMSim.CreateInspector("clump_max_z");
        auto min_z_finder = DEMSim.CreateInspector("clump_min_z");
        auto total_mass_finder = DEMSim.CreateInspector("clump_mass");
        auto max_v_finder = DEMSim.CreateInspector("clump_max_absv");

        float3 this_G = make_float3(-G_mag * std::sin(G_ang), 0, -G_mag * std::cos(G_ang));
        DEMSim.SetGravitationalAcceleration(this_G);

        DEMSim.SetInitTimeStep(step_size);
        DEMSim.SetCDUpdateFreq(15);
        // Max velocity info is generally just for the solver's reference and the user do not have to set it. The solver
        // wouldn't take into account a vel larger than this when doing async-ed contact detection: but this vel won't
        // happen anyway and if it does, something already went wrong.
        DEMSim.SetMaxVelocity(50.);
        // Error out vel is used to force the simulation to abort when something goes wrong and sim diverges.
        DEMSim.SetErrorOutVelocity(60.);
        DEMSim.SetExpandSafetyMultiplier(1.1);
        // You usually don't have to worry about initial bin size. In very rare cases, init bin size is so bad that auto
        // bin size adaption is effectless, and you should notice in that case kT runs extremely slow. Then in that case
        // setting init bin size may save the simulation.
        // DEMSim.SetInitBinSize(2 * scales.at(2));

        // 256 or 512 are common choices. Note that in cases where the force model is modified, too many registers may
        // be used in the kernel, so we have to reduce this number to use 256. In other cases (and most cases), 512 is
        // fine and may make the code run a bit faster. Usually, the user do not have to call
        // SetForceCalcThreadsPerBlock if they don't know the implication.
        DEMSim.SetForceCalcThreadsPerBlock(512);
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
        float max_z = -0.35;
        float init_x = -0.6;
        wheel_tracker->SetPos(make_float3(init_x, 0, max_z + wheel_rad));
        for (double t = 0; t < 0.2; t += frame_time) {
            DEMSim.DoDynamicsThenSync(frame_time);
        }

        DEMSim.ChangeFamily(10, 1);

        // filestream
        std::ofstream file;
        std::string filename = out_dir.string() + "/driving_torque.csv";
        file.open(filename);
        file << "time,pos_x,pos_z,driving_torque_local\n";

        bool start_measure = false;
        for (double t = 0; t < sim_end; t += step_size, curr_step++) {
            if (curr_step % out_steps == 0) {
                char filename[200], meshname[200];
                std::cout << "Outputting frame: " << currframe << std::endl;
                sprintf(filename, "%s/DEMdemo_output_%04d.csv", out_dir.c_str(), currframe);
                // sprintf(meshname, "%s/DEMdemo_mesh_%04d.vtk", out_dir.c_str(), currframe);
                DEMSim.WriteSphereFile(std::string(filename));
                // DEMSim.WriteMeshFile(std::string(meshname));
                DEMSim.ShowThreadCollaborationStats();
                currframe++;
            }

            if (t >= 2. && !start_measure) {
                start_measure = true;
            }

            if (curr_step % report_steps == 0 && start_measure) {
                float3 pos = wheel_tracker->Pos();
                float drum_torque_glb = wheel_tracker->ContactAngAccGlobal().y * wheel_IYY;
                float drum_torque_loc = wheel_tracker->ContactAngAccLocal().y * wheel_IYY;

                std::cout << "Time: " << t << std::endl;
                std::cout << "X    : " << wheel_tracker->Pos().x << std::endl;
                std::cout << "pos_Z: " << pos.z << std::endl;
                std::cout << "Max system velocity: " << max_v_finder->GetValue() << std::endl;

                file << t << "," << pos.x << ", " << pos.z << "," << drum_torque_loc << "\n";
            }

            DEMSim.DoDynamics(step_size);
        }

        file.close();
        run_mode++;
        DEMSim.ShowTimingStats();
        DEMSim.ShowAnomalies();
    }

    std::cout << "DEMdemo_WheelSlopeSlip demo exiting..." << std::endl;
    return 0;
}
//  Copyright (c) 2021, SBEL GPU Development Team
//  Copyright (c) 2021, University of Wisconsin - Madison
//
//	SPDX-License-Identifier: BSD-3-Clause

#include <DEM/API.h>
#include <DEM/HostSideHelpers.hpp>
#include <DEM/utils/Samplers.hpp>

#include <chrono>
#include <cmath>
#include <cstdio>
#include <filesystem>
#include <map>
#include <random>
#include <cstdlib>
#include <vector>
#include <array>

using namespace deme;

const double math_PI = 3.1415927;

unsigned long long comb(int n, int k) {
    if (k < 0 || k > n)
        return 0;
    if (k > n - k)
        k = n - k;
    unsigned long long result = 1;
    for (int i = 1; i <= k; ++i) {
        result *= n - i + 1;
        result /= i;
    }
    return result;
}

std::vector<double> cBernstein(int degree, double t) {
    std::vector<double> B(degree + 1);
    for (int i = 0; i <= degree; ++i) {
        B[i] = comb(degree, i) * std::pow(t, i) * std::pow(1.0 - t, degree - i);
    }
    return B;
}

double deviationCausedHighest(double cp_deviation, double rad, double width) {
    if (cp_deviation <= 0.) {
        return 0.;
    }
    const int num_CPs = 4;
    double t = 0.5;
    std::vector<std::array<double, 2>> outer_CPs(num_CPs);
    outer_CPs[0] = {0., width / 2};
    outer_CPs[1] = {rad * cp_deviation, width / 2 - width / 3};
    outer_CPs[2] = {rad * cp_deviation, width / 2 - width / 3 * 2};
    outer_CPs[3] = {0., -width / 2};

    std::vector<double> b = cBernstein(num_CPs - 1, t);

    // Matrix multiplication (b @ outer_CPs)
    std::array<double, 2> eval_pnt = {0., 0.};
    for (int i = 0; i < num_CPs; ++i) {
        eval_pnt[0] += b[i] * outer_CPs[i][0];
        eval_pnt[1] += b[i] * outer_CPs[i][1];
    }

    return eval_pnt[0];
}

int main(int argc, char* argv[]) {
    int cur_test = atoi(argv[1]);

    std::filesystem::path out_dir = std::filesystem::current_path();
    out_dir += "/" + std::string(argv[8]) + "_" + std::to_string(cur_test);
    std::filesystem::create_directory(out_dir);

    // `World'
    float G_mag = atof(argv[11]);
    float step_size = 1e-5;     // 5e-6;
    double world_size_y = 0.6;  // 0.52;
    double world_size_x = 4.;
    float safe_x = 1.6;
    double world_size_z = 4.0;
    float w_r = atof(argv[10]);
    double sim_end = 10.;
    float z_adv_targ = 0.2;

    // Define the wheel geometry
    float wheel_rad = atof(argv[2]);
    float eff_mass = atof(argv[3]);
    float wheel_width = atof(argv[5]);
    float wheel_mass = 5.;  // 8.7;
    float total_pressure = eff_mass * G_mag;
    float added_pressure = (total_pressure - wheel_mass * G_mag);
    float wheel_IYY = wheel_mass * wheel_rad * wheel_rad / 2;
    float wheel_IXX = (wheel_mass / 12) * (3 * wheel_rad * wheel_rad + wheel_width * wheel_width);
    float grouser_height = atof(argv[4]);

    float Slope_deg = atof(argv[6]);
    float cp_dev;
    cp_dev = deviationCausedHighest(atof(argv[7]), wheel_rad, wheel_width);

    std::cout << "cur_test: " << cur_test << std::endl;
    std::cout << "wheel_rad: " << wheel_rad << std::endl;
    std::cout << "eff_mass: " << eff_mass << std::endl;
    std::cout << "grouser_height: " << grouser_height << std::endl;
    std::cout << "wheel_width: " << wheel_width << std::endl;
    std::cout << "Slope_deg: " << Slope_deg << std::endl;
    std::cout << "cp_dev: " << cp_dev << std::endl;
    std::cout << "Out dir: " << out_dir << std::endl;
    std::cout <<"mesh file: " << argv[9] << std::endl;
    std::cout << "Ang vel: " << w_r << std::endl;
    std::cout << "G_mag: " << G_mag << std::endl;
    std::cout <<"mu wheel" << atof(argv[12]) << std::endl;

    {
        DEMSolver DEMSim;
        DEMSim.SetVerbosity(VERBOSITY::INFO);
        DEMSim.SetOutputFormat(OUTPUT_FORMAT::CSV);
        DEMSim.SetOutputContent(OUTPUT_CONTENT::ABSV);
        DEMSim.SetMeshOutputFormat(MESH_FORMAT::VTK);
        DEMSim.SetCollectAccRightAfterForceCalc(true);

        // E, nu, CoR, mu, Crr...
        float mu = 0.4;
        float mu_wheel = atof(argv[12]);
        float mu_wall = 1.;
        auto mat_type_wall =
            DEMSim.LoadMaterial({{"E", 1e9}, {"nu", 0.3}, {"CoR", 0.4}, {"mu", mu_wall}, {"Crr", 0.00}});
        auto mat_type_wheel =
            DEMSim.LoadMaterial({{"E", 1e9}, {"nu", 0.3}, {"CoR", 0.4}, {"mu", mu_wheel}, {"Crr", 0.00}});
        auto mat_type_terrain = DEMSim.LoadMaterial({{"E", 1e9}, {"nu", 0.3}, {"CoR", 0.4}, {"mu", mu}, {"Crr", 0.00}});
        DEMSim.SetMaterialPropertyPair("mu", mat_type_wheel, mat_type_terrain, mu_wheel);
        DEMSim.SetMaterialPropertyPair("mu", mat_type_wall, mat_type_terrain, mu_wall);

        DEMSim.InstructBoxDomainDimension(world_size_x, world_size_y, world_size_z);
        DEMSim.InstructBoxDomainBoundingBC("top_open", mat_type_wall);

        float bottom = -0.5;
        auto bot_wall = DEMSim.AddBCPlane(make_float3(0, 0, bottom), make_float3(0, 0, 1), mat_type_wall);
        auto bot_wall_tracker = DEMSim.Track(bot_wall);

        auto wheel = DEMSim.AddWavefrontMeshObject(std::string(argv[9]), mat_type_wheel);
        wheel->SetMass(wheel_mass);
        wheel->SetMOI(make_float3(wheel_IXX, wheel_IYY, wheel_IXX));
        // Give the wheel a family number so we can potentially add prescription
        wheel->SetFamily(11);
        DEMSim.SetFamilyFixed(11);
        DEMSim.DisableContactBetweenFamilies(11, 0);
        // Track it
        auto wheel_tracker = DEMSim.Track(wheel);

        // Define the terrain particle templates
        // Calculate its mass and MOI
        float terrain_density = 2.6e3;
        float volume1 = 4.2520508;
        float mass1 = terrain_density * volume1;
        float3 MOI1 = make_float3(1.6850426, 1.6375114, 2.1187753) * terrain_density;
        // Scale the template we just created
        std::vector<double> scales = {0.008};
        // Then load it to system
        std::shared_ptr<DEMClumpTemplate> my_template1 =
            DEMSim.LoadClumpType(mass1, MOI1, GetDEMEDataFile("clumps/triangular_flat.csv"), mat_type_terrain);
        std::vector<std::shared_ptr<DEMClumpTemplate>> ground_particle_templates = {my_template1,
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
            auto clump_xyz = DEMSim.ReadClumpXyzFromCsv("./GRC_3e6.csv");
            auto clump_quaternion = DEMSim.ReadClumpQuatFromCsv("./GRC_3e6.csv");
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
                // Prepare clump type identification vector for loading into the system (don't forget type 0 in
                // ground_particle_templates is the template for rover wheel)
                std::vector<std::shared_ptr<DEMClumpTemplate>> this_type(n_clump_this_type,
                                                                         ground_particle_templates.at(t_num));

                // Add them to the big long vector
                in_xyz.insert(in_xyz.end(), this_type_xyz.begin(), this_type_xyz.end());
                in_quat.insert(in_quat.end(), this_type_quat.begin(), this_type_quat.end());
                in_types.insert(in_types.end(), this_type.begin(), this_type.end());
                // Our template names are 0000, 0001 etc.
                t_num++;
            }

            // Now, we don't need all particles loaded...
            std::vector<notStupidBool_t> elem_to_remove(in_xyz.size(), 0);
            for (size_t i = 0; i < in_xyz.size(); i++) {
                if (std::abs(in_xyz.at(i).y) > (world_size_y - 0.05) / 2 ||
                    std::abs(in_xyz.at(i).x) > (world_size_x - 0.05) / 2)
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
            DEMSim.AddClumps(base_batch);
        }

        // Now add a plane to compress the sample
        // auto compressor = DEMSim.AddExternalObject();
        // compressor->AddPlane(make_float3(0, 0, 0), make_float3(0, 0, -1), mat_type_terrain);
        // compressor->SetFamily(2);
        // auto compressor_tracker = DEMSim.Track(compressor);

        // Families' prescribed motions (Earth)
        float v_ref = w_r * (wheel_rad + cp_dev + grouser_height);
        double G_ang = Slope_deg * math_PI / 180.;

        // Note: this wheel is not `dictated' by our prescrption of motion because it can still fall onto the ground
        // (move freely linearly)
        DEMSim.SetFamilyPrescribedAngVel(1, "0", to_string_with_precision(w_r), "0", false);
        DEMSim.AddFamilyPrescribedAcc(1, to_string_with_precision(-added_pressure * std::sin(G_ang) / wheel_mass),
                                      "none", to_string_with_precision(-added_pressure * std::cos(G_ang) / wheel_mass));
        DEMSim.SetFamilyFixed(10);
        DEMSim.DisableContactBetweenFamilies(10, 10);
        DEMSim.DisableContactBetweenFamilies(10, 255);

        // Some inspectors
        auto max_z_finder = DEMSim.CreateInspector("clump_max_z");
        // auto min_z_finder = DEMSim.CreateInspector("clump_min_z");
        // auto total_mass_finder = DEMSim.CreateInspector("clump_mass");
        // auto partial_mass_finder = DEMSim.CreateInspector("clump_mass", "return (Z <= -0.41);");
        // auto max_v_finder = DEMSim.CreateInspector("clump_max_absv");

        float3 this_G = make_float3(-G_mag * std::sin(G_ang), 0, -G_mag * std::cos(G_ang));
        DEMSim.SetGravitationalAcceleration(this_G);

        DEMSim.SetInitTimeStep(step_size);
        DEMSim.SetCDUpdateFreq(10);
        DEMSim.SetExpandSafetyAdder(v_ref);
        DEMSim.SetCDNumStepsMaxDriftMultipleOfAvg(1);
        DEMSim.SetCDNumStepsMaxDriftAheadOfAvg(5);
        DEMSim.SetErrorOutVelocity(50.);
        DEMSim.Initialize();

        // Put the wheel in place, then let the wheel sink in initially
        float init_x = -0.6;
        if (Slope_deg < 21) {
            init_x = -1.4;
        }

        float settle_time = 0.4;
        { DEMSim.DoDynamicsThenSync(settle_time); }

        // Put the wheel in place, then let the wheel sink in initially
        float max_z = max_z_finder->GetValue();
        wheel_tracker->SetPos(make_float3(init_x, 0, max_z + 0.03 + cp_dev + grouser_height + wheel_rad));


        DEMSim.ChangeFamily(11, 1);

        {
            std::chrono::high_resolution_clock::time_point start = std::chrono::high_resolution_clock::now();
            // DEMSim.DoDynamicsThenSync(2.);

            float3 V = wheel_tracker->Vel();
            float x1 = wheel_tracker->Pos().x;
            float z1 = x1 * std::sin(G_ang);
            float x2, z2, z_adv = 0.;

            float t;
            int report_fps = 100;
            int report_steps = (unsigned int)(1.0 / (report_fps * step_size));;
            float report_time = report_steps * step_size;
            unsigned int cur_step = 0;
            double energy = 0.;
            unsigned int report_num = 0;

            const int output_fps = 10; // This is just for vtk files -> This fps is actually applied as report_fps/output_fps
            unsigned int out_steps = (unsigned int)(1.0 / (output_fps * step_size));
            unsigned int curr_frame = 0;
            for (t = 0; t < sim_end; t += report_time, report_num++)
            {
                if (report_num % out_steps == 0) {
                    std::cout << "Output Frame: " << curr_frame << std::endl;
                    char filename[200], meshname[200];
                    sprintf(filename, "%s/DEMdemo_output_%04d_frame_%04d.csv", out_dir.c_str(), cur_test, curr_frame);
                    sprintf(meshname, "%s/DEMdemo_mesh_%04d_frame_%04d.vtk", out_dir.c_str(), cur_test, curr_frame);
                    DEMSim.WriteSphereFile(std::string(filename));
                    DEMSim.WriteMeshFile(std::string(meshname));
                    curr_frame++;
                }

                float adv = x2 - x1;
                float eff_energy = eff_mass * z_adv * G_mag;
                std::cout << "Time: " << t << std::endl;
                std::cout << "Slip: " << 1. - adv / (v_ref * t) << std::endl;
                std::cout << "Energy: " << energy << std::endl;
                std::cout << "Power: " << energy / t << std::endl;
                std::cout << "Efficiency: " << eff_energy / energy << std::endl;
                std::cout << "Z advance: " << z_adv << std::endl;

                if (z_adv >= z_adv_targ)
                    break;
                float3 angAcc = wheel_tracker->ContactAngAccLocal();
                energy += std::abs((double)angAcc.y * wheel_IYY * w_r * (double)report_time);
                x2 = wheel_tracker->Pos().x;
                z_adv = x2 * std::sin(G_ang) - z1;
                if (x2 > safe_x)
                    break;

                DEMSim.DoDynamics(report_time);
            }
            std::chrono::high_resolution_clock::time_point end = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double> time_sec =
                std::chrono::duration_cast<std::chrono::duration<double>>(end - start);
            std::cout << "Runtime: " << (time_sec.count()) << std::endl;
        }
    }

    return 0;
}
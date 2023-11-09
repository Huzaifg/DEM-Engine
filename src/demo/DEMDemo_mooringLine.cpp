//  Copyright (c) 2021, SBEL GPU Development Team
//  Copyright (c) 2021, University of Wisconsin - Madison
//
//	SPDX-License-Identifier: BSD-3-Clause

// =============================================================================
// Fracture
// =============================================================================

#include <DEM/API.h>
#include <DEM/HostSideHelpers.hpp>
#include <DEM/utils/Samplers.hpp>

#include <chrono>
#include <cmath>
#include <cstdio>
#include <filesystem>
#include <iostream>
#include <fstream>

using namespace deme;

const double math_PI = 3.1415927;

int main() {
    DEMSolver DEMSim;
    DEMSim.SetVerbosity(INFO);
    DEMSim.SetOutputFormat(OUTPUT_FORMAT::CSV);
    DEMSim.SetOutputContent(OUTPUT_CONTENT::VEL);
    DEMSim.SetMeshOutputFormat(MESH_FORMAT::VTK);
    DEMSim.SetContactOutputContent(OWNER | FORCE | POINT | CNT_WILDCARD);

    DEMSim.SetErrorOutAvgContacts(200);
    // DEMSim.SetForceCalcThreadsPerBlock(256);
    //  E, nu, CoR, mu, Crr...
    auto mat_type_container =
        DEMSim.LoadMaterial({{"E", 20e9}, {"nu", 0.3}, {"CoR", 0.7}, {"mu", 0.80}, {"Crr", 0.10}});
    auto mat_type_particle = DEMSim.LoadMaterial({{"E", 60e9}, {"nu", 0.20}, {"CoR", 0.5}, {"mu", 0.5}, {"Crr", 0.05}});
    // If you don't have this line, then values will take average between 2 materials, when they are in contact
    DEMSim.SetMaterialPropertyPair("CoR", mat_type_container, mat_type_particle, 0.2);
    DEMSim.SetMaterialPropertyPair("Crr", mat_type_container, mat_type_particle, 0.5);
    DEMSim.SetMaterialPropertyPair("mu", mat_type_container, mat_type_particle, 0.5);
    // We can specify the force model using a file.
    auto my_force_model = DEMSim.ReadContactForceModel("ForceModelMooring.cu");

    // Those following lines are needed. We must let the solver know that those var names are history variable etc.
    my_force_model->SetMustHaveMatProp({"E", "nu", "CoR", "mu"});
    my_force_model->SetMustPairwiseMatProp({"CoR", "mu", "Crr"});
    // Pay attention to the extra per-contact wildcard `unbroken' here.
    my_force_model->SetPerContactWildcards({"initialLength", "innerInteraction"});

    float world_size = 10.0;
    float container_diameter = 0.05;
    float terrain_density = 2.0e3;
    float sphere_rad = 0.100;

    float step_size = 5e-5;


    DEMSim.InstructBoxDomainDimension(world_size, world_size, world_size);
    // No need to add simulation `world' boundaries, b/c we'll add a cylinderical container manually
    DEMSim.InstructBoxDomainBoundingBC("all", mat_type_container);
    // DEMSim.SetInitBinSize(sphere_rad * 5);
    //  Now add a cylinderical boundary along with a bottom plane
    double bottom = -0;
    double top = 0.10;

    // auto walls = DEMSim.AddExternalObject();
    // std::string shake_pattern_xz = " 0.00 * sin( 500 * 2 * deme::PI * t)";
    // walls->AddPlane(make_float3(0, 0, bottom), make_float3(0, 0, 1), mat_type_container);
    // walls->AddPlane(make_float3(0, 0, world_size / 2. - world_size / 20.), make_float3(0, 0, -1),
    // mat_type_container); walls->SetFamily(2);

    auto walls = DEMSim.AddWavefrontMeshObject("../data/granularFlow/funnel_left.obj", mat_type_container);
    float3 move = make_float3(0.05, 0.00, 0 - 2 * sphere_rad);  // z
    float4 rot = make_float4(0.7071, 0, 0, 0.7071);
    walls->Scale(make_float3(0.8, 0.1, 2.0));
    walls->Move(move, rot);
    walls->SetFamily(2);
    DEMSim.SetFamilyFixed(2);

    auto cylinder = DEMSim.AddExternalObject();
    cylinder->AddCylinder(make_float3(0), make_float3(0, 0, 1), 1.8 * container_diameter / 2., mat_type_container, 0);
    cylinder->SetFamily(10);
    DEMSim.SetFamilyFixed(10);
    // DEMSim.SetFamilyPrescribedLinVel(10, shake_pattern_xz, shake_pattern_xz, shake_pattern_xz);

    auto plate = DEMSim.AddWavefrontMeshObject("../data/granularFlow/funnel_left.obj", mat_type_container);
    move = make_float3(0.05, 0.00, top + 2.0 * sphere_rad);  // z
    rot = make_float4(0.7071, 0, 0, 0.7071);
    plate->Scale(make_float3(0.8, 0.1, 2.0));
    plate->Move(move, rot);
    plate->SetFamily(20);
    DEMSim.SetFamilyFixed(20);

    // Track it

    auto plate_tracker = DEMSim.Track(plate);

    DEMSim.SetFamilyPrescribedLinVel(21, "0", "0", to_string_with_precision(-0.010));
    DEMSim.SetFamilyPrescribedLinVel(22, "0", "0", to_string_with_precision(-0.0050));

    // auto total_mass_finder = DEMSim.CreateInspector("clump_mass");
    // auto max_v_finder = DEMSim.CreateInspector("clump_max_absv");

    // Define the terrain particle templates

    // Calculate its mass and MOI

    float sphere_vol = 4. / 3. * math_PI * sphere_rad * sphere_rad * sphere_rad;
    float mass = terrain_density * sphere_vol;
    // Then load it to system
    std::shared_ptr<DEMClumpTemplate> my_template = DEMSim.LoadSphereType(mass, sphere_rad, mat_type_particle);

    // Sampler to sample
    HCPSampler sampler(2.0 * sphere_rad);

    float fill_height = top;
    float3 fill_center = make_float3(0, 0, bottom + fill_height / 2);
    const float fill_radius = container_diameter / 2. - sphere_rad * 0.;
    float3 fill = make_float3(fill_radius, fill_radius, bottom + fill_height / 2);
    auto input_xyz = sampler.SampleBox(fill_center, fill);
    // auto input_xyz = sampler(fill_center, fill_radius, fill_height / 2 - sphere_rad * 0.);
    auto particles = DEMSim.AddClumps(my_template, input_xyz);
    particles->SetFamily(1);
    std::cout << "Total num of particles: " << particles->GetNumClumps() << std::endl;

    std::filesystem::path out_dir = std::filesystem::current_path();
    std::string nameOutFolder = "R" + std::to_string(sphere_rad) + "_Int" + std::to_string(fact_radius) + "";
    out_dir += "/DemoOutput_Fracture_" + nameOutFolder;
    remove_all(out_dir);
    create_directory(out_dir);

    // Some inspectors

    auto max_z_finder = DEMSim.CreateInspector("clump_max_z");
    auto min_z_finder = DEMSim.CreateInspector("clump_min_z");

    DEMSim.SetFamilyExtraMargin(1, fact_radius * sphere_rad);

    DEMSim.SetInitTimeStep(step_size);
    DEMSim.SetGravitationalAcceleration(make_float3(0, 0.00, 1 * -9.81));
    DEMSim.Initialize();
    // DEMSim.DisableContactBetweenFamilies(20, 1);
    std::cout << "Initial number of contacts: " << DEMSim.GetNumContacts() << std::endl;

    float sim_end = 5;
    unsigned int fps = 1000;
    unsigned int datafps = 1000;
    unsigned int modfpsGeo = datafps / fps;
    float frame_time = 1.0 / datafps;
    std::cout << "Output at " << fps << " FPS" << std::endl;
    unsigned int out_steps = (unsigned int)(1.0 / (datafps * step_size));
    unsigned int frame_count = 0;
    unsigned int step_count = 0;

    bool status_1 = true;
    bool status_2 = true;

    // DEMSim.DisableContactBetweenFamilies(10, 1);

    double L0;
    double stress;
    std::string nameOutFile = "data_R" + std::to_string(sphere_rad) + "_Int" + std::to_string(fact_radius) + ".csv";
    std::ofstream csvFile(nameOutFile);

    DEMSim.SetFamilyContactWildcardValueAll(1, "initialLength", 0.0);
    // DEMSim.SetFamilyContactWildcardValueAll(1, "damage", 0.0);
    DEMSim.SetFamilyContactWildcardValueAll(1, "unbroken", 0.0);

    // Simulation loop
    for (float t = 0; t < sim_end; t += frame_time) {
        // DEMSim.ShowThreadCollaborationStats();

        std::cout << "Initial number of contacts: " << DEMSim.GetNumContacts() << std::endl;

        if (t >= 0.0 && status_1) {
            status_1 = false;
            DEMSim.DoDynamicsThenSync(0);
            DEMSim.SetFamilyContactWildcardValueAll(1, "unbroken", 2.0);
            DEMSim.ChangeFamily(20, 21);  // start compression

            L0 = max_z_finder->GetValue() - min_z_finder->GetValue() + 2 * sphere_rad;

            std::cout << "Establishing inner forces: " << frame_count << std::endl;
        }

        if (stress > 1e6 && status_2) {
            status_2 = false;
            DEMSim.DoDynamicsThenSync(0);
            DEMSim.ChangeFamily(21, 22);  // start compression
            std::cout << "Stress in now [Pa]: " << stress << std::endl;
            std::cout << "Adapting motion velocity: " << frame_count << std::endl;
        }

        if (frame_count % 1 == 0) {
            float3 forces = plate_tracker->ContactAcc();
            float3 pos = plate_tracker->Pos();
            stress = forces.z / (container_diameter * container_diameter);
            double L = max_z_finder->GetValue() - min_z_finder->GetValue() + 2 * sphere_rad;
            std::cout << "Time: " << t << std::endl;
            std::cout << "Pos of plate: " << pos.z << std::endl;
            std::cout << "Stress [Pa]: " << stress << std::endl;
            std::cout << "Strain [-]: " << (L - L0) / L0 << std::endl;
            csvFile << (L - L0) / L0 << "; " << stress << std::endl;

            if (frame_count % modfpsGeo == 0) {
                char filename[200];
                char meshname[200];
                char cnt_filename[200];

                std::cout << "Outputting frame: " << frame_count / modfpsGeo << std::endl;
                sprintf(filename, "%s/DEMdemo_output_%04d.csv", out_dir.c_str(), frame_count / modfpsGeo);
                sprintf(meshname, "%s/DEMdemo_mesh_%04d.vtk", out_dir.c_str(), frame_count / modfpsGeo);
                sprintf(cnt_filename, "%s/DEMdemo_contact_%04d.csv", out_dir.c_str(), frame_count / modfpsGeo);

                DEMSim.WriteSphereFile(std::string(filename));
                DEMSim.WriteMeshFile(std::string(meshname));
                DEMSim.WriteContactFile(std::string(cnt_filename));
            }
        }

        DEMSim.DoDynamics(frame_time);
        frame_count++;
    }
    csvFile.close();
    DEMSim.ShowTimingStats();
    std::cout << "Fracture demo exiting..." << std::endl;
    return 0;
}

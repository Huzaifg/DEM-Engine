#include <DEM/API.h>
#include <core/ApiVersion.h>
#include <core/utils/ThreadManager.h>
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
#include "document.h"
#include "filereadstream.h"

using namespace deme;
using namespace std::filesystem;

const double math_PI = 3.1415927;

// Define a structure for baffle positions
struct BafflePosition {
  double x, y, z;
};

// Updated SimulationParams struct to include a vector of BafflePosition
struct SimulationParams {
  double bxDim;
  double byDim;
  double bzDim;
  double baffle_thickness;
  double baffle_height;
  double baffle_width;
  std::vector<BafflePosition>
      baffle_positions; // Vector to store baffle positions
  double granular_thickness;
  double granular_height;
  double granular_width;
  double granular_x; // Assuming this is the initial position of the granular
                     // material box
  double granular_y;
  double granular_z;
  double t_end;
  float time_step;
};

void readSimulationParams(SimulationParams &params,
                          const std::string &filename) {
  FILE *fp = fopen(filename.c_str(), "r");
  if (!fp) {
    std::cerr << "Failed to open JSON file: " << filename << std::endl;
    return;
  }

  char readBuffer[65536];
  rapidjson::FileReadStream is(fp, readBuffer, sizeof(readBuffer));

  rapidjson::Document doc;
  doc.ParseStream(is);
  fclose(fp);

  if (!doc.IsObject()) {
    std::cerr << "Invalid JSON format in: " << filename << std::endl;
    return;
  }

  // Assuming the JSON structure matches the names in the SimulationParams
  // struct
  params.bxDim = doc["bxDim"].GetDouble();
  params.byDim = doc["byDim"].GetDouble();
  params.bzDim = doc["bzDim"].GetDouble();
  params.baffle_thickness = doc["baffle_thickness"].GetDouble();
  params.baffle_height = doc["baffle_height"].GetDouble();
  params.baffle_width = doc["baffle_width"].GetDouble();
  // Handling baffle_positions as nested
  const rapidjson::Value &bafflePositions = doc["baffle_positions"];
  if (bafflePositions.IsArray()) {
    for (auto &pos : bafflePositions.GetArray()) {
      BafflePosition position;
      position.x = pos["x"].GetDouble();
      position.y = pos["y"].GetDouble();
      position.z = pos["z"].GetDouble();
      params.baffle_positions.push_back(position);
    }
  }
  params.granular_thickness = doc["granular_thickness"].GetDouble();
  params.granular_height = doc["granular_height"].GetDouble();
  params.granular_width = doc["granular_width"].GetDouble();
}

int main(int argc, char* argv[]) {
     std::string paramsFileNumber = argv[1];
    // Set output directory
    path out_dir = current_path();
    out_dir += "/BAFFLE_FLOW_BENCH/benchmark_" + paramsFileNumber + "/";
    create_directory(out_dir);

    float step_size = 5e-6;     // 5e-6;

    // Read the simulation domain parameters
    SimulationParams simParams;
    if (argc < 2) {
    std::cerr << "Usage: " << argv[0] << " <simulation_params.json>"
              << std::endl;
        return 1;
    }

    
    std::string params = "input_json/rtf_benchmark_jsons/demo_FSI_Baffle_RTF_Benchmark_" + paramsFileNumber + ".json";
    std::string paramsFile = (GET_DATA_PATH() / params).string();
    readSimulationParams(simParams, paramsFile);

    // Transalate baffle_positions and granular_x, granular_y, granular_z by -bxDim / 2 along x, -byDim / 2 along y and -bzDim / 2 along z
    for (auto &pos : simParams.baffle_positions) {
        pos.x -= simParams.bxDim / 2;
        pos.y -= simParams.byDim / 2;
        pos.z -= simParams.bzDim / 2;
    }

    // Granular particle radius
    // Defined here so that they have enough clearance above the ground
    float particle_radius = 0.006 / pow(2.8, 1.0/3.0);

    // Particles box left bottom corner coordinates
    // The constant is just to ensure the SPH particles don't intersect with the
    // baffle BCE markers
    simParams.granular_x = simParams.baffle_positions[0].x -
                            simParams.baffle_thickness / 2 -
                            simParams.granular_thickness - particle_radius * 2.01;
    simParams.granular_y =
        simParams.baffle_positions[0].y - simParams.baffle_width / 2;
    // simParams.granular_z = (-simParams.bzDim / 2.) + particle_radius * 2.01;
    simParams.granular_z = (-simParams.bzDim / 2.);


    // Final simulation time
    simParams.t_end = 1.0;
    simParams.time_step = 5e-6;

    
    DEMSolver DEMSim;
    DEMSim.SetVerbosity(VERBOSITY::INFO);
    DEMSim.SetOutputFormat(OUTPUT_FORMAT::CSV);
    DEMSim.SetOutputContent(OUTPUT_CONTENT::ABSV);
    DEMSim.SetMeshOutputFormat(MESH_FORMAT::VTK);
    DEMSim.SetCollectAccRightAfterForceCalc(true);
    float mu_terrain = 0.3;
    float mu_baffle = 0.2;  
    float mu_wall = 1.;
    auto mat_type_terrain =
    DEMSim.LoadMaterial({{"E", 2e6}, {"nu", 0.3}, {"CoR", 0.4}, {"mu", mu_terrain}, {"Crr", 0.3}, {"Cohesion", 40}});

    auto mat_type_baffle = DEMSim.LoadMaterial({{"E", 1e8}, {"nu", 0.24}, {"CoR", 0.05}, {"mu", mu_baffle}, {"Crr", 0.0}, {"Cohesion", 0}});

    auto mat_type_wall =
        DEMSim.LoadMaterial({{"E", 1e9}, {"nu", 0.3}, {"CoR", 0.4}, {"mu", mu_wall}, {"Crr", 0.00}, {"Cohesion", 0}});

    DEMSim.SetMaterialPropertyPair("mu", mat_type_terrain, mat_type_baffle, mu_baffle);
    DEMSim.SetMaterialPropertyPair("mu", mat_type_terrain, mat_type_wall, mu_wall);
    
    // Don't want any rolling friction between the terrain and baffles/walls
    DEMSim.SetMaterialPropertyPair("Crr", mat_type_terrain, mat_type_baffle, 0.0);
    DEMSim.SetMaterialPropertyPair("Crr", mat_type_terrain, mat_type_wall, 0.3);

    DEMSim.SetMaterialPropertyPair("Cohesion", mat_type_terrain, mat_type_baffle, 0.0);
    DEMSim.SetMaterialPropertyPair("Cohesion", mat_type_terrain, mat_type_wall, 10);

    DEMSim.SetMaterialPropertyPair("CoR", mat_type_terrain, mat_type_baffle, 0.05);
    DEMSim.SetMaterialPropertyPair("CoR", mat_type_terrain, mat_type_wall, 0.05);

    DEMSim.InstructBoxDomainDimension(simParams.bxDim, simParams.byDim, simParams.bzDim);
    DEMSim.InstructBoxDomainBoundingBC("top_open", mat_type_wall);

    // Create the 3 baffles
    // Baffle 1
    auto baffle_1 = DEMSim.AddWavefrontMeshObject((GET_DATA_PATH() / "mesh/cube.obj").string(), mat_type_baffle);
    // scale it
    baffle_1->Scale(make_float3(simParams.baffle_thickness, simParams.baffle_width, simParams.baffle_height));

    std::cout << "Total num of triangles: " << baffle_1->GetNumTriangles() << std::endl;

    // Position the baffle
    baffle_1->SetInitPos(make_float3(simParams.baffle_positions[0].x, simParams.baffle_positions[0].y, simParams.baffle_positions[0].z + 0.5 * simParams.baffle_height));
    float density = 7850; //Kg/m^3
    float volume = simParams.baffle_thickness * simParams.baffle_width * simParams.baffle_height; // m^3
    float mass = density * volume; //Kg
    baffle_1->SetMass(mass);
    // MOI is given by 1/12 * m * (a^2+b^2) where a and b are the two sides of the cuboid that intersect at the body diagonal axis
    baffle_1->SetMOI(make_float3(1.0/12.0 * mass * (simParams.baffle_width * simParams.baffle_width + simParams.baffle_height * simParams.baffle_height), 1.0/12.0 * mass * (simParams.baffle_thickness * simParams.baffle_thickness + simParams.baffle_height * simParams.baffle_height), 1.0/12.0 * mass * (simParams.baffle_thickness * simParams.baffle_thickness + simParams.baffle_width * simParams.baffle_width)));
    baffle_1->SetFamily(2);
    DEMSim.SetFamilyFixed(2);
    DEMSim.DisableContactBetweenFamilies(255, 2);
    auto baffle_1_tracker = DEMSim.Track(baffle_1);


    // Baffle 2
    auto baffle_2 = DEMSim.AddWavefrontMeshObject((GET_DATA_PATH() / "mesh/cube.obj").string(), mat_type_baffle);
    // scale it
    baffle_2->Scale(make_float3(simParams.baffle_thickness, simParams.baffle_width, simParams.baffle_height));

    std::cout << "Total num of triangles: " << baffle_2->GetNumTriangles() << std::endl;

    // Position the baffle
    baffle_2->SetInitPos(make_float3(simParams.baffle_positions[1].x, simParams.baffle_positions[1].y, simParams.baffle_positions[1].z + 0.5 * simParams.baffle_height));
    baffle_2->SetMass(mass);
    baffle_2->SetMOI(make_float3(1.0/12.0 * mass * (simParams.baffle_width * simParams.baffle_width + simParams.baffle_height * simParams.baffle_height), 1.0/12.0 * mass * (simParams.baffle_thickness * simParams.baffle_thickness + simParams.baffle_height * simParams.baffle_height), 1.0/12.0 * mass * (simParams.baffle_thickness * simParams.baffle_thickness + simParams.baffle_width * simParams.baffle_width)));
    baffle_2->SetFamily(3);
    DEMSim.SetFamilyFixed(3);
    DEMSim.DisableContactBetweenFamilies(255, 3);
    auto baffle_2_tracker = DEMSim.Track(baffle_2);

    // Baffle 3
    auto baffle_3 = DEMSim.AddWavefrontMeshObject((GET_DATA_PATH() / "mesh/cube.obj").string(), mat_type_baffle);
    // scale it
    baffle_3->Scale(make_float3(simParams.baffle_thickness, simParams.baffle_width, simParams.baffle_height));

    std::cout << "Total num of triangles: " << baffle_3->GetNumTriangles() << std::endl;

    // Position the baffle
    baffle_3->SetInitPos(make_float3(simParams.baffle_positions[2].x, simParams.baffle_positions[2].y, simParams.baffle_positions[2].z + 0.5 * simParams.baffle_height));
    baffle_3->SetMass(mass);
    baffle_3->SetMOI(make_float3(1.0/12.0 * mass * (simParams.baffle_width * simParams.baffle_width + simParams.baffle_height * simParams.baffle_height), 1.0/12.0 * mass * (simParams.baffle_thickness * simParams.baffle_thickness + simParams.baffle_height * simParams.baffle_height), 1.0/12.0 * mass * (simParams.baffle_thickness * simParams.baffle_thickness + simParams.baffle_width * simParams.baffle_width)));
    baffle_3->SetFamily(4);
    DEMSim.SetFamilyFixed(4);
    DEMSim.DisableContactBetweenFamilies(255, 4);
    auto baffle_3_tracker = DEMSim.Track(baffle_3);


    // Create the granular domain
    // This is adjusted to meet approx the same number of "particles" in FSI

    auto terrain_template = DEMSim.LoadSphereType(particle_radius * particle_radius * particle_radius * 1.8e3 * 4 / 3 * math_PI, particle_radius, mat_type_terrain);
    unsigned int num_particle = 0;
    
    float3 box_center = make_float3(simParams.granular_x + 0.5 * simParams.granular_thickness, simParams.granular_y + 0.5 * simParams.granular_width, simParams.granular_z + 0.5 * simParams.granular_height);
    
    float3 box_halfwidth = make_float3(0.5 * simParams.granular_thickness, 0.5 * simParams.granular_width, 0.5 * simParams.granular_height);

    // GridSampler
    // auto particles = DEMBoxGridSampler(box_center, box_halfwidth, particle_radius * 2.01);
    // PDSampler
    PDSampler sampler(2.01 * particle_radius);
    auto particles = sampler.SampleBox(box_center, box_halfwidth);
    
    // Custom Cohesion force model
    auto my_force_model = DEMSim.ReadContactForceModel("ForceModelWithCohesion.cu");
    // This custom force model still uses contact history arrays, so let's define it
    my_force_model->SetPerContactWildcards({"delta_time", "delta_tan_x", "delta_tan_y", "delta_tan_z"});
    my_force_model->SetMustPairwiseMatProp({"CoR", "mu", "Crr", "Cohesion"});

    
    // Initial velocites for particles
    float init_v_x = 1.;

    auto clumps = DEMSim.AddClumps(terrain_template, particles);
    num_particle += particles.size();
    clumps->SetVel(make_float3(init_v_x, 0., 0.));

    std::cout << "Total num of particles: " << num_particle << std::endl;



    DEMSim.SetInitTimeStep(simParams.time_step);
    DEMSim.SetMaxVelocity(30.);
    DEMSim.SetGravitationalAcceleration(make_float3(0, 0, -9.81));
    DEMSim.SetIntegrator("centered_difference");

    DEMSim.Initialize();

    unsigned int fps = 100;
    float frame_time = 1.0 / fps;
    unsigned int out_steps = (unsigned int)(1.0 / (fps * simParams.time_step));

    std::cout << "Output at " << fps << " FPS" << std::endl;
    unsigned int currframe = 0;

    std::chrono::high_resolution_clock::time_point start = std::chrono::high_resolution_clock::now();
    // DEMSim.DoDynamics(simParams.t_end);
    for (float t = 0; t < simParams.t_end; t += frame_time) {
        std::cout << "Frame: " << currframe << std::endl;
        char filename[200], meshfilename[200], cnt_filename[200];
        sprintf(filename, "%s/DEMdemo_output_%04d.csv", out_dir.c_str(), currframe);
        sprintf(meshfilename, "%s/DEMdemo_mesh_%04d.vtk", out_dir.c_str(), currframe);
        // sprintf(cnt_filename, "%s/Contact_pairs_%04d.csv", out_dir.c_str(), currframe);
        DEMSim.WriteSphereFile(std::string(filename));
        DEMSim.WriteMeshFile(std::string(meshfilename));
        DEMSim.WriteContactFile(std::string(cnt_filename));
        currframe++;

        DEMSim.DoDynamics(frame_time);
        // DEMSim.DoDynamics(simParams.t_end);
        DEMSim.ShowThreadCollaborationStats();

    }
    std::chrono::high_resolution_clock::time_point end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> time_ms = std::chrono::duration_cast<std::chrono::duration<double, std::milli>>(end - start);
    std::cout << "Simulation time: " << simParams.t_end << " seconds" << std::endl;
    std::cout << "Real Time Taken by Simulation: " << time_ms.count() << " milliseconds" << std::endl;
    DEMSim.ShowTimingStats();
    std::cout << "SIMEND" << std::endl;
    return 0;
    
  }
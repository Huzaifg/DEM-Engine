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
// using namespace std::filesystem;
using namespace rapidjson;

const double math_PI = 3.1415927;

double bxDim = 0.8;
double byDim = 0.8;
double bzDim = 0.8;

double t_end = 350 * 0.0025;
bool output = true;
float step_size = 5e-6;     
// double output_fps = 500; // Only works with timestep 5e-6 to get output at timestep 0.0025



// Size of the baffles
double baffle_thickness = 0.1;  // Dimension in the x direction
double baffle_height = 0.3;     // Dimension in the z direction
double baffle_width = 0.1;      // Dimension in the y direction

// ----------------------------------------------------------------------------
// Struct for storing the ranges of random parameters obtained from a JSON file
// ----------------------------------------------------------------------------
struct RandomParams {
    RandomParams() : sim_id(0) {}

    int sim_id;  // Simulation ID

    // Granular Material
    int no_granular_piles[2];   // Range of number of granular piles to sample from
    double pile_size_range[2];  // Range of dimensions along x,y and z of the granular pile
    double granular_x[2];       // Range of starting x coordinate of the granular pile
    double granular_y[2];       // Range of starting y coordinate of the granular pile
    double granular_z[2];       // Range of starting z coordinate of the granular pile
    double pile_velx_range[2];  // Range of x component of the velocity of the granular pile
    double pile_vely_range[2];  // Range of y component of the velocity of the granular pile
    double pile_velz_range[2];  // Range of z component of the velocity of the granular pile

    // Baffles
    double no_obstacles[2];  // Range of number of obstacles

    double baffle_x[2];  // Range of x coordinate of the baffle
    double baffle_y[2];  // Range of y coordinate of the baffle
    double baffle_z[2];  // Range of z coordinate of the baffle

    // Once the granualr particles are initialized, the following variables are set in order to ensure the baffles are
    // randomly placed at a sufficient distance
    std::vector<std::array<double, 3>> granular_pile_start;  // Vector of granular pile start coordinates
    std::vector<std::array<double, 3>> granular_pile_size;   // Vector of granular pile dimensions
    int sampled_granular_piles;                         // Number of granular piles sampled
};

// ----------------------------------------------------------------------------
// Read range of randomized variables from JSON file
// ----------------------------------------------------------------------------
void readRandomRanges(RandomParams& ranges, const std::string& ranges_file) {
    std::cout << "Reading parameters from: " << ranges_file << std::endl;

    FILE* fp = fopen(ranges_file.c_str(), "r");
    if (!fp) {
        std::cerr << "Invalid JSON file!" << std::endl;
        return;
    }

    char readBuffer[32768];
    FileReadStream is(fp, readBuffer, sizeof(readBuffer));
    fclose(fp);

    Document doc;

    doc.ParseStream<ParseFlag::kParseCommentsFlag>(is);
    if (!doc.IsObject()) {
        std::cerr << "Invalid JSON file!!" << std::endl;
        return;
    }

    if (doc.HasMember("sim_id")) {
        ranges.sim_id = doc["sim_id"].GetInt();
    }

    if (doc.HasMember("no_granular_piles")) {
        ranges.no_granular_piles[0] = doc["no_granular_piles"][0].GetInt();
        ranges.no_granular_piles[1] = doc["no_granular_piles"][1].GetInt();
    }

    if (doc.HasMember("pile_size_range")) {
        ranges.pile_size_range[0] = doc["pile_size_range"][0].GetDouble();
        ranges.pile_size_range[1] = doc["pile_size_range"][1].GetDouble();
    }

    if (doc.HasMember("pile_start_range")) {
        ranges.granular_x[0] = doc["pile_start_range"][0][0u].GetDouble();
        ranges.granular_x[1] = doc["pile_start_range"][0][1u].GetDouble();

        ranges.granular_y[0] = doc["pile_start_range"][1][0u].GetDouble();
        ranges.granular_y[1] = doc["pile_start_range"][1][1u].GetDouble();

        ranges.granular_z[0] = doc["pile_start_range"][2][0u].GetDouble();
        ranges.granular_z[1] = doc["pile_start_range"][2][1u].GetDouble();
    }

    if (doc.HasMember("pile_vel_range")) {
        ranges.pile_velx_range[0] = doc["pile_vel_range"][0][0u].GetDouble();
        ranges.pile_velx_range[1] = doc["pile_vel_range"][0][1u].GetDouble();

        ranges.pile_vely_range[0] = doc["pile_vel_range"][1][0u].GetDouble();
        ranges.pile_vely_range[1] = doc["pile_vel_range"][1][1u].GetDouble();

        ranges.pile_velz_range[0] = doc["pile_vel_range"][2][0u].GetDouble();
        ranges.pile_velz_range[1] = doc["pile_vel_range"][2][1u].GetDouble();
    }

    if (doc.HasMember("no_obstacles")) {
        ranges.no_obstacles[0] = doc["no_obstacles"][0].GetDouble();
        ranges.no_obstacles[1] = doc["no_obstacles"][1].GetDouble();
    }

    if (doc.HasMember("obstacle_cg_range")) {
        ranges.baffle_x[0] = doc["obstacle_cg_range"][0][0u].GetDouble();
        ranges.baffle_x[1] = doc["obstacle_cg_range"][0][1u].GetDouble();

        ranges.baffle_y[0] = doc["obstacle_cg_range"][1][0u].GetDouble();
        ranges.baffle_y[1] = doc["obstacle_cg_range"][1][1u].GetDouble();

        ranges.baffle_z[0] = doc["obstacle_cg_range"][2][0u].GetDouble();
        ranges.baffle_z[1] = doc["obstacle_cg_range"][2][1u].GetDouble();
    }
}

// ------------------------------------------
// Some utility functions for random numbers
// ------------------------------------------
int random_int(int min, int max) {
    std::random_device rd;                              // Seed for random number generator
    std::uniform_int_distribution<int> dist(min, max);  // Distribution over the range
    return dist(rd);
}

double random_double(double min, double max) {
    std::random_device rd;   // Seed random number generator
    std::mt19937 gen(rd());  // Use Mersenne Twister engine
    std::uniform_real_distribution<> dist(min, max);

    return dist(gen);
}


// ----------------------------------------------
// Check if the granular pile is within the box
// ----------------------------------------------
bool checkGranularPlacement(const double granular_start_x,
                            const double granular_start_y,
                            const double granular_start_z,
                            const double granular_thickness,
                            const double granular_width,
                            const double granular_height) {
    // Calculate end positions of granular pile
    double granular_end_x = granular_start_x + granular_thickness;
    double granular_end_y = granular_start_y + granular_width;
    double granular_end_z = granular_start_z + granular_height;

    // Calculate end positions of box
    double box_end_x = bxDim;
    double box_end_y = byDim;
    double box_end_z = bzDim;

    // Check if granular pile is completely within the box in all dimensions
    return (0. < granular_start_x && granular_end_x < box_end_x) &&
           (0. < granular_start_y && granular_end_y < box_end_y) &&
           (0. < granular_start_z && granular_end_z < box_end_z);
}

// -----------------------------------------------------------
// Check if the baffle is placed away from the granular piles
// Also checks if the baffles are placed away from each other
// -----------------------------------------------------------
bool checkBafflePlacement(const double baffle_x,
                          const double baffle_y,
                          const double baffle_z,
                          const RandomParams ranges,
                          const std::vector<std::array<double, 3>>& baffles,
                          const int baffle_index) {
    // Calculate baffle extents (considering half the dimensions)
    double baffle_x_start = baffle_x - baffle_thickness / 2;
    double baffle_y_start = baffle_y - baffle_width / 2;
    double baffle_z_start = baffle_z - baffle_height / 2;

    double baffle_x_end = baffle_x + baffle_thickness / 2;
    double baffle_y_end = baffle_y + baffle_width / 2;
    double baffle_z_end = baffle_z + baffle_height / 2;

    // Loop over all the granular piles
    for (int i = 0; i < ranges.sampled_granular_piles; i++) {
        // Calculate granular pile start and end coordinates
        double gp_start_x = ranges.granular_pile_start[i][0];
        double gp_start_y = ranges.granular_pile_start[i][1];
        double gp_start_z = ranges.granular_pile_start[i][2];

        double gp_end_x = gp_start_x + ranges.granular_pile_size[i][0];
        double gp_end_y = gp_start_y + ranges.granular_pile_size[i][1];
        double gp_end_z = gp_start_z + ranges.granular_pile_size[i][2];

        // Check for overlap -> If there is any overlap, return false
        if (!((gp_end_x < baffle_x_start) ||  // Granular pile completely left of baffle
              (gp_start_x > baffle_x_end) ||  // Granular pile completely right of baffle
              (gp_end_y < baffle_y_start) ||  // Granular pile completely below baffle
              (gp_start_y > baffle_y_end) ||  // Granular pile completely above baffle
              (gp_end_z < baffle_z_start) ||  // Granular pile completely in front of baffle
              (gp_start_z > baffle_z_end))    // Granular pile completely behind baffle
        ) {
            return false;
        }
    }

    if (baffle_index == 0)
        return true;  // No baffles to check against (first baffle is always valid

    // Loop over all the other baffles to ensure they are placed away from each other
    for (int i = 0; i < baffle_index; i++) {
        // Calculate baffle extents (considering half the dimensions)
        double other_baffle_x_start = baffles[i][0] - baffle_thickness / 2;
        double other_baffle_y_start = baffles[i][1] - baffle_width / 2;
        double other_baffle_z_start = baffles[i][2] - baffle_height / 2;

        double other_baffle_x_end = baffles[i][0] + baffle_thickness / 2;
        double other_baffle_y_end = baffles[i][1] + baffle_width / 2;
        double other_baffle_z_end = baffles[i][2] + baffle_height / 2;

        // Check for overlap -> If there is any overlap, return false
        if (!((other_baffle_x_end < baffle_x_start) ||  // Other baffle completely left of baffle
              (other_baffle_x_start > baffle_x_end) ||  // Other baffle completely right of baffle
              (other_baffle_y_end < baffle_y_start) ||  // Other baffle completely below baffle
              (other_baffle_y_start > baffle_y_end) ||  // Other baffle completely above baffle
              (other_baffle_z_end < baffle_z_start) ||  // Other baffle completely in front of baffle
              (other_baffle_z_start > baffle_z_end))    // Other baffle completely behind baffle
        ) {
            return false;
        }
    }
    // If we reach here, the baffle doesn't overlap any piles
    return true;
}

int main(int argc, char* argv[]) {

    // ========================
    // Initialize DEM stuff
    // ========================
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

    DEMSim.InstructBoxDomainDimension(bxDim, byDim, bzDim);
    DEMSim.InstructBoxDomainBoundingBC("top_open", mat_type_wall);


    // Set output directory
    std::filesystem::path out_dir = std::filesystem::current_path();
    std::string ranges_file;
    RandomParams ranges;

    if(argc < 3){
        out_dir += "/BAFFLE_FLOW_TRAIN";
        ranges_file = argv[2];
    } else{
        out_dir += "/" + (std::string(argv[3]) + "_BAFFLE_FLOW_TRAIN_" + std::string(argv[1]));
        ranges_file = argv[2];
    }

    std::cout<<"Output directory: "<<out_dir<<std::endl;
    std::cout<<"Ranges file: "<<ranges_file<<std::endl;
    
    create_directory(out_dir);
    // // Above one fails sometimes
    // if (!filesystem::create_directory(filesystem::path(out_dir))) {
    //     std::cerr << "Error creating directory " << out_dir << std::endl;
    //     return 1;
    // }


    // Read the random ranges
    readRandomRanges(ranges, ranges_file);

    // Keeping just one pile for now because of inter-pile spacing
    // int numPiles = random_int(ranges.no_granular_piles[0], ranges.no_granular_piles[1]);
    int numPiles = 1;
    ranges.sampled_granular_piles = numPiles;
    double max_vol_dif = 0.;

    // Calculate the particle radius that we need to meet the max limit of 17000 particles within a domain of 0.4*0.4*0.4 with 2.01 times the particle radius spacing
    double volume_of_each_sphere = (0.4 * 0.4 * 0.4) / (17100*3.);
    // This does not account for the additional spacing needed between each particle - this is okay, this means we will have always lesser than 17100 particles
    double particle_radius = std::cbrt(volume_of_each_sphere * 3 / (4 * math_PI)); // We need to use the same particle radius as physics chages with particle radius
    std::cout << "Particle radius: " << particle_radius << std::endl;

    // Compute the granular pile size
    std::array<double, 2> granular_pile_size = {0.0, 0.0};
    for (int i = 0; i < numPiles; i++) {
        granular_pile_size[i] = random_double(ranges.pile_size_range[0], ranges.pile_size_range[1]);
        auto vol = granular_pile_size[i] * granular_pile_size[i] * granular_pile_size[i];
        auto vol_dif = vol / (0.2 * 0.2 * 0.2);
        if (vol_dif > max_vol_dif) {
            max_vol_dif = vol_dif;
        }
    }

    std::cout << "Granular pile size: " << granular_pile_size[0] << std::endl;

    
    // Get the position of the graunla pile
    for (int i = 0; i < numPiles; i++){
        bool valid = false;

        double granular_x, granular_y, granular_z;


        int attempt = 0;
        int max_attempts = 10000;  // Set a maximum number of attempts to avoid infinite loop
        double granular_dim = 0;

        while(!valid && attempt < max_attempts){
            granular_x = random_double(ranges.granular_x[0], ranges.granular_x[1]);
            granular_y = random_double(ranges.granular_y[0], ranges.granular_y[1]);
            granular_z = random_double(ranges.granular_z[0], ranges.granular_z[1]);

            granular_dim = granular_pile_size[i];
            valid = checkGranularPlacement(granular_x, granular_y, granular_z, granular_pile_size[i], granular_pile_size[i], granular_pile_size[i]);
            attempt++;
        }
        if (attempt == max_attempts) {
            std::cerr << "Max attempts reached for placing granular pile" << std::endl;
        }

        // Transform the granular pile start coordinates to match DEM coordinates which has box (0,0,0) at the center
        granular_x -= bxDim / 2;
        granular_y -= byDim / 2;
        granular_z -= bzDim / 2;

    
        ranges.granular_pile_start.push_back({granular_x, granular_y, granular_z});
        ranges.granular_pile_size.push_back({granular_pile_size[i], granular_pile_size[i], granular_pile_size[i]});

        auto terrain_template = DEMSim.LoadSphereType(particle_radius * particle_radius * particle_radius * 1.8e3 * 4 / 3 * math_PI, particle_radius, mat_type_terrain);
        unsigned int num_particle = 0;


        float3 box_center = make_float3(granular_x + 0.5 * granular_dim, granular_y + 0.5 * granular_dim, granular_z + 0.5 * granular_dim);

        float3 box_halfwidth = make_float3(0.5 * granular_dim, 0.5 * granular_dim, 0.5 * granular_dim);


        PDSampler sampler(2.01 * particle_radius);
        auto particles = sampler.SampleBox(box_center, box_halfwidth);

        // Custom Cohesion force model
        auto my_force_model = DEMSim.ReadContactForceModel("ForceModelWithCohesion.cu");
        // This custom force model still uses contact history arrays, so let's define it
        my_force_model->SetPerContactWildcards({"delta_time", "delta_tan_x", "delta_tan_y", "delta_tan_z"});
        my_force_model->SetMustPairwiseMatProp({"CoR", "mu", "Crr", "Cohesion"});


        // Initial velocites for particles
        double pile_velx = random_double(ranges.pile_velx_range[0], ranges.pile_velx_range[1]);
        double pile_vely = random_double(ranges.pile_vely_range[0], ranges.pile_vely_range[1]);
        double pile_velz = random_double(ranges.pile_velz_range[0], ranges.pile_velz_range[1]);

        auto clumps = DEMSim.AddClumps(terrain_template, particles);
        num_particle += particles.size();
        clumps->SetVel(make_float3(pile_velx, pile_vely, pile_velz));

        std::cout << "Total num of particles: " << num_particle << std::endl;
    }

    // Initialize baffles and its positions
    int numBaffles = random_int(ranges.no_obstacles[0], ranges.no_obstacles[1]);
    std::vector<std::array<double, 3>> baffle_pos(numBaffles);
    std::vector<std::shared_ptr<DEMMeshConnected>> baffles(numBaffles);
    std::vector<std::shared_ptr<DEMTracker>> baffle_trackers(numBaffles);

    for(int i =0; i < numBaffles; i++){
        bool valid = false;
        double baffle_x, baffle_y, baffle_z;
        int attempt = 0;
        int max_attempts = 10000;  // Set a maximum number of attempts to avoid infinite loop
        while(!valid && attempt < max_attempts){
            baffle_x = random_double(ranges.baffle_x[0], ranges.baffle_x[1]);
            baffle_y = random_double(ranges.baffle_y[0], ranges.baffle_y[1]);
            baffle_z = random_double(ranges.baffle_z[0], ranges.baffle_z[1]);
            valid = checkBafflePlacement(baffle_x, baffle_y, baffle_z, ranges, baffle_pos, i);
            attempt++;
        }
        if (attempt == max_attempts) {
            std::cerr << "Max attempts reached for placing baffle " << i << std::endl;
        }

        auto baffle = DEMSim.AddWavefrontMeshObject((GET_DATA_PATH() / "mesh/cube.obj").string(), mat_type_baffle);
        // scale it
        baffle->Scale(make_float3(baffle_thickness, baffle_width, baffle_height));
        std::cout << "Total num of triangles: " << baffle->GetNumTriangles() << std::endl;

        // Transform baffle positon to DEM frame 
        baffle_x -= bxDim / 2;
        baffle_y -= byDim / 2;
        baffle_z -= bzDim / 2;
        baffle->SetInitPos(make_float3(baffle_x, baffle_y, baffle_z + 0.5 * baffle_height));
        float density = 7850; //Kg/m^3
        float volume = baffle_thickness * baffle_width * baffle_height; // m^3
        float mass = density * volume; //Kg
        baffle->SetMass(mass);
        // MOI is given by 1/12 * m * (a^2+b^2) where a and b are the two sides of the cuboid that intersect at the body diagonal axis
        baffle->SetMOI(make_float3(1.0/12.0 * mass * (baffle_width * baffle_width + baffle_height * baffle_height), 1.0/12.0 * mass * (baffle_thickness * baffle_thickness + baffle_height * baffle_height), 1.0/12.0 * mass * (baffle_thickness * baffle_thickness + baffle_width * baffle_width)));
        baffle->SetFamily(i+2);
        DEMSim.SetFamilyFixed(i+2);
        DEMSim.DisableContactBetweenFamilies(255, i+2);
        auto baffle_tracker = DEMSim.Track(baffle);


        baffle_pos.push_back({baffle_x, baffle_y, baffle_z});
        baffles.push_back(baffle);
        baffle_trackers.push_back(baffle_tracker);
    }

    DEMSim.SetInitTimeStep(step_size);
    DEMSim.SetMaxVelocity(30.);
    DEMSim.SetGravitationalAcceleration(make_float3(0, 0, -9.81));
    DEMSim.SetIntegrator("centered_difference");

    DEMSim.Initialize();

    float frame_time = 0.0025;

    unsigned int currframe = 0;

    std::chrono::high_resolution_clock::time_point start = std::chrono::high_resolution_clock::now();
    for (float t = 0; t < t_end; t += frame_time) {
        std::cout << "Frame: " << currframe << std::endl;
        char filename[200], meshfilename[200], cnt_filename[200];
        sprintf(filename, "%s/DEMdemo_output_%04d.csv", out_dir.c_str(), currframe);
        sprintf(meshfilename, "%s/DEMdemo_mesh_%04d.vtk", out_dir.c_str(), currframe);
        // sprintf(cnt_filename, "%s/Contact_pairs_%04d.csv", out_dir.c_str(), currframe);
        DEMSim.WriteSphereFile(std::string(filename));
        DEMSim.WriteMeshFile(std::string(meshfilename));
        // DEMSim.WriteContactFile(std::string(cnt_filename));
        currframe++;

        DEMSim.DoDynamics(frame_time);
    }
    std::chrono::high_resolution_clock::time_point end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> time_ms = std::chrono::duration_cast<std::chrono::duration<double, std::milli>>(end - start);
    std::cout << "Simulation time: " << t_end << " seconds" << std::endl;
    std::cout << "Real Time Taken by Simulation: " << time_ms.count() << " milliseconds" << std::endl;
    std::cout << "SIMEND" << std::endl;
    return 0;
    
  }
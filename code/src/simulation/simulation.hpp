#pragma once

#include <cmath>
#include "cgp/cgp.hpp"



// SPH Particle
struct particle_element
{
    cgp::vec3 p; // Position
    cgp::vec3 v; // Speed
    cgp::vec3 f; // Force

    float rho;      // density at this particle position
    float pressure; // pressure at this particle position

    particle_element() : p{0,0,0},v{0,0,0},f{0,0,0},rho(0),pressure(0) {}
};

// SPH simulation parameters
struct sph_parameters_structure
{
    // Influence distance of a particle (size of the kernel)
    static constexpr float h = 0.16f; // 3 is max divergence-free speed for h = 0.1f

    // Rest density (normalized to 1 - real unit should be 1000kg/m^2)
    float rho0 = 1;

     // Total mass of a particle (consider rho0 h^2)
    float m = rho0*h*h;

    // viscosity parameter
    float nu = 0.02f;   
     
    // Stiffness converting density to pressure
    // float stiffness = 8.0f;
    float stiffness = 8.0f;

    // coefficient for surface tension
    float sigma = 15.0; 
    //particles near surface have non-zero norm of normals and inside liquid it should be 0
    float norm_threshold = 0.01; // the bigger it is the more particles leave sphere with 0 gravity

    // gravity vector is opposite to this direction:
    cgp::vec3 scene_up_vector = {0,1,0}; 
    // gravity in prev frame
    bool prev_frame_zero_gravity = false;

    // Acceleration unirofm 3D grid 
    // with structures for mapping between cell position and particle id 
    std::vector<cgp::vec3> particle_id2grid_cell;
    static const int grid_size = (int)ceil(2.0 / h);
    std::array<std::vector<int>, grid_size*grid_size*grid_size> grid_cell2particle_ids;
}; 

void simulate(float dt, cgp::buffer<particle_element>& particles, 
    sph_parameters_structure& sph_parameters, const cgp::vec3 gravity_vector, 
    bool fix_gravity_direction, bool zero_gravity);

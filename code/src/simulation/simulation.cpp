#include "simulation.hpp"

using namespace cgp;

// Convert a density value to a pressure
float density_to_pressure(float rho, float rho0, float stiffness)
{
	return stiffness*(rho-rho0);
}

float W_laplacian_viscosity(vec3 const& p_i, vec3 const& p_j, float h)
{
    float const r = norm(p_i - p_j);
    if (r<=h)
        return 45.0 / (3.14159f * std::pow(h, 6)) * (h - r);
    else
        return 0.0f;
}

vec3 viscosity_term(particle_element const& p_i, particle_element const& p_j, float h)
{
    return 1 / p_j.rho * (p_j.v - p_i.v) * W_laplacian_viscosity(p_i.p, p_j.p, h);
}

vec3 W_gradient_pressure(vec3 const& p_i, vec3 const& p_j, float h)
{
    float const r = norm(p_i - p_j);
    if (r<=h)
        return -45.0 / (3.14159f * std::pow(h, 6)) * std::pow(h - r, 2) * (p_i - p_j) / r;
    else
        return {0,0,0};
}

vec3 pressure_term(particle_element const& p_i, particle_element const& p_j, float h)
{
    return 1 / (2 * p_j.rho) * (p_i.pressure + p_j.pressure) * 
            W_gradient_pressure(p_i.p, p_j.p, h);
}

float W_density(vec3 const& p_i, const vec3& p_j, float h)
{
	float const r = norm(p_i-p_j);
    
    if (r<=h)
	   return 315.0/(64.0*3.14159f*std::pow(h,9)) * std::pow(h*h-r*r, 3.0f);
    else 
        return 0.0f;
}

float density_term(particle_element const& p_i, particle_element const& p_j, float h)
{
    return W_density(p_i.p, p_j.p, h);
}

vec3 W_gradient_color_field(vec3 const& p_i, const vec3& p_j, float h)
{
    float const r = norm(p_i - p_j);
    
    if (r<=h)
       return -3*315.0/(32.0*3.14159f*std::pow(h,9)) * (h*h-r*r) * (h*h-r*r) * (p_i - p_j);
    else 
        return {0,0,0};
}

vec3 normal_term(particle_element const& p_i, particle_element const& p_j, float h)
{
    return 1 / p_j.rho * W_gradient_color_field(p_i.p, p_j.p, h);
}

float W_laplacian_color_field(vec3 const& p_i, const vec3& p_j, float h)
{
    float const r = norm(p_i - p_j);
    if (r<=h) 
       return -3*315.0/(32.0*3.14159f*std::pow(h,9)) * 
                ( 3*(h*h-r*r) * (h*h-r*r) - 4 * (h*h-r*r) * r * r );
    else 
        return 0.0f;
}

float curvature_term(particle_element const& p_i, particle_element const& p_j, float h)
{
    return 1 / p_j.rho * W_laplacian_color_field(p_i.p, p_j.p, h);
}


// for scalar fields
void update_scalar(int particle_id, bool use_same_id, float& zero_value, 
    float (*term)(particle_element const&, particle_element const&, float),
    sph_parameters_structure const& params, const buffer<particle_element>& particles)
{       
    vec3 c_ids = params.particle_id2grid_cell[particle_id];
    int g_size = params.grid_size;
    
    // loop over neighbour grid cells
    for(int x = -1; x < 2; x++)
        for(int y = -1; y < 2; y++)
            for(int z = -1; z < 2; z++)
                if ((c_ids[0]+x >= 0) and (c_ids[0]+x < params.grid_size) and
                    (c_ids[1]+y >= 0) and (c_ids[1]+y < params.grid_size) and
                    (c_ids[2]+z >= 0) and (c_ids[2]+z < params.grid_size))
                {
                    int id = (c_ids[0] + x) * g_size * g_size + (c_ids[1] + y) * g_size + 
                             c_ids[2] + z;
                    std::vector<int> p_ids = params.grid_cell2particle_ids[id];
                    // loop over particles in neighbour grid cells
                    for (int j = 0; j < p_ids.size(); j++)
                        if (use_same_id or (particle_id != p_ids[j]))
                            zero_value += term(particles[particle_id], 
                                particles[p_ids[j]], params.h);
                }
}

// for vector fields
void update_vector(int particle_id, bool use_same_id, vec3& zero_value, 
    vec3 (*term)(particle_element const&, particle_element const&, float),
    sph_parameters_structure const& params, const buffer<particle_element>& particles)
{       
    vec3 c_ids = params.particle_id2grid_cell[particle_id];
    int g_size = params.grid_size;
    
    // loop over neighbour grid cells
    for(int x = -1; x < 2; x++)
        for(int y = -1; y < 2; y++)
            for(int z = -1; z < 2; z++)
                if ((c_ids[0]+x >= 0) and (c_ids[0]+x < params.grid_size) and
                    (c_ids[1]+y >= 0) and (c_ids[1]+y < params.grid_size) and
                    (c_ids[2]+z >= 0) and (c_ids[2]+z < params.grid_size))
                {
                    int id = (c_ids[0] + x) * g_size * g_size + (c_ids[1] + y) * g_size + 
                             c_ids[2] + z;
                    std::vector<int> p_ids = params.grid_cell2particle_ids[id];
                    // loop over particles in neighbour grid cells
                    for (int j = 0; j < p_ids.size(); j++)
                        if (use_same_id or (particle_id != p_ids[j]))
                            zero_value += term(particles[particle_id], 
                                particles[p_ids[j]], params.h);
                }
}

void update_density(buffer<particle_element>& particles, float h, float m)
{
    int const N = particles.size();
    float new_density;

    for(int i=0; i<N; ++i)
    {
        new_density = 0;
        for(int j=0; j<N; ++j)
            new_density += m * W_density(particles[i].p, particles[j].p, h); 
        particles[i].rho = new_density;
    }
}

// NEW FUNCTION BELOW - with uniform grid
void update_density_new(buffer<particle_element>& particles, 
    sph_parameters_structure const& params)
{
    float m = params.m;
    int const N = particles.size();
    float new_density;

    for(int i=0; i<N; ++i)
    {
        new_density = 0;
        update_scalar(i, true, new_density, &density_term, params, particles);
        particles[i].rho = m * new_density;
    }
}

// Convert the particle density to pressure
void update_pressure(buffer<particle_element>& particles, float rho0, float stiffness)
{
	const int N = particles.size();
    for(int i=0; i<N; ++i)
        particles[i].pressure = density_to_pressure(particles[i].rho, rho0, stiffness);
}

// Compute the forces and update the acceleration of the particles
void update_force(buffer<particle_element>& particles, float h, float m, float nu)
{
    const int N = particles.size();

    for(int i=0; i<N; ++i)
    {
        // gravity    
        particles[i].f = m * vec3{0,-9.81f,0};
        
        // pressure
        vec3 F_pressure = {0,0,0};
        for(int j=0; j<N; ++j)
            if (i != j)
                F_pressure += m / (2 * particles[j].rho) * (particles[i].pressure + particles[j].pressure) * W_gradient_pressure(particles[i].p, particles[j].p, h);
        particles[i].f += - m / particles[i].rho * F_pressure;

        // viscosity 
        vec3 F_viscosity = {0,0,0};
        for(int j=0; j<N; ++j)
            F_viscosity += m / particles[j].rho * (particles[j].v - particles[i].v) * W_laplacian_viscosity(particles[i].p, particles[j].p, h);
        particles[i].f += m * nu * F_viscosity;
    }
}

// NEW FUNCTION BELOW - uniform grid
void update_force_new(buffer<particle_element>& particles, sph_parameters_structure const& params)
{
    const int N = particles.size();
    float m = params.m, nu = params.nu;
    vec3 up_vec = params.scene_up_vector;
    float up_vec_n = norm(up_vec);

    for(int i=0; i<N; ++i)
    {
        // gravity    
        if (up_vec_n > 0) 
            particles[i].f = m * -9.81f * up_vec / norm(up_vec);
        else
            particles[i].f = {0,0,0};

        vec3 F_pressure = {0,0,0};
        update_vector(i, false, F_pressure, &pressure_term, params, particles);
        particles[i].f += - m / particles[i].rho * m * F_pressure;

        vec3 F_viscosity = {0,0,0};
        update_vector(i, false, F_viscosity, &viscosity_term, params, particles);
        particles[i].f += m * nu * m * F_viscosity;

        vec3 inner_normal = {0,0,0};
        update_vector(i, true, inner_normal, &normal_term, params, particles);
        inner_normal = m * inner_normal;
        float norm_len = norm(inner_normal);
        if (norm_len > params.norm_threshold)
        {
            float curvature = 0;
            update_scalar(i, true, curvature, &curvature_term, params, particles);
            curvature = -1 / norm_len * m * curvature;

            vec3 tension_density = params.sigma * curvature * inner_normal;
            vec3 tension = tension_density * m / particles[i].rho;
            particles[i].f += tension;
        }
    }
}

void simulate(float dt, buffer<particle_element>& particles, 
    sph_parameters_structure& sph_parameters, const vec3 camera_up_vector,
     bool fix_gravity_direction, bool zero_gravity)
{
    static const int N = particles.size();
    int const g_size = sph_parameters.grid_size;
    float h = sph_parameters.h;

    // Clean grid structures
    sph_parameters.particle_id2grid_cell.clear();
    for (int x = 0; x < g_size ; x++)
        for (int y = 0; y < g_size ; y++)
            for (int z = 0; z < g_size ; z++)
                sph_parameters.grid_cell2particle_ids[x * g_size * g_size + y * g_size + z].clear();

    // Update grid structures
    for (int i = 0; i < N; i++)
    {   
        vec3 axes_ids = {
            (int)floor( (particles[i].p[0] + 1) / h),
            (int)floor( (particles[i].p[1] + 1) / h),
            (int)floor( (particles[i].p[2] + 1) / h)
        };
        sph_parameters.particle_id2grid_cell.push_back(axes_ids);
        int id = axes_ids[0] * g_size * g_size + axes_ids[1] * g_size + axes_ids[2];
        sph_parameters.grid_cell2particle_ids[id].push_back(i);
    }

    // update gravity direction if necessary

    if (zero_gravity)
    {
        sph_parameters.scene_up_vector = {0,0,0};
        sph_parameters.prev_frame_zero_gravity = true;
    }
    else
    {
        if (sph_parameters.prev_frame_zero_gravity)
        {
            sph_parameters.scene_up_vector = camera_up_vector;
            sph_parameters.prev_frame_zero_gravity = false;
        }
        if (!fix_gravity_direction)
            sph_parameters.scene_up_vector = camera_up_vector;
    }
    

	// Update values with UNIFORM GRID - new
    update_density_new(particles, sph_parameters);          // same results as old
    update_pressure(particles, sph_parameters.rho0, sph_parameters.stiffness); // same func
    update_force_new(particles, sph_parameters);            // 

    // Update values - old
    // update_density(particles, sph_parameters.h, sph_parameters.m);                   // First compute updated density
    // update_pressure(particles, sph_parameters.rho0, sph_parameters.stiffness);       // Compute associated pressure
    // update_force(particles, sph_parameters.h, sph_parameters.m, sph_parameters.nu);  // Update forces

	// Numerical integration
	float const damping = 0.005f;
	float const m = sph_parameters.m;

	for(int k=0; k<N; ++k)
	{
		vec3& p = particles[k].p;
        if (std::isnan(p[0]) or std::isnan(p[1]) or std::isnan(p[2]))
            std::cout << "nan P" << std::endl;
		vec3& v = particles[k].v;
        if (std::isnan(v[0]) or std::isnan(v[1]) or std::isnan(v[2]))
            std::cout << "nan V" << std::endl;
        vec3& f = particles[k].f;
        if (std::isnan(f[0]) or std::isnan(f[1]) or std::isnan(f[2]))
        {
            std::cout << "nan F: p " << k << ' ' << f << std::endl;
        }

		v = (1-damping)*v + dt*f/m;
		p = p + dt*v;
	}


	// Collision with the box

    // step back from wall
    // float const epsilon = 1e-3f;
    float const epsilon = 0.5 * 1e-3f;

    // particle tremble scale
    float sh = epsilon;

    for(int k=0; k<N; ++k)
    {
        vec3& p = particles[k].p;
        vec3& v = particles[k].v;

        // small perturbation to avoid alignment
        // if( p.y < -1 ) {p.y = -1 + (epsilon * rand_interval());  v.y *= -0.5f;}

        if( p.y < -1 ) {p.y = -1 + (epsilon + sh * rand_interval());  v.y *= -0.5f;}
        if( p.y > 1 )  {p.y =  1 - (epsilon + sh * rand_interval());  v.y *= -0.5f;}

        if( p.x < -1 ) {p.x = -1 + (epsilon + sh * rand_interval());  v.x *= -0.5f;}
        if( p.x > 1 )  {p.x =  1 - (epsilon + sh * rand_interval());  v.x *= -0.5f;}

        if( p.z < -1 ) {p.z = -1 + (epsilon + sh * rand_interval());  v.z *= -0.5f;}
        if( p.z > 1 )  {p.z =  1 - (epsilon + sh * rand_interval());  v.z *= -0.5f;}
    }

}
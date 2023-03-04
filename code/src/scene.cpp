#include "scene.hpp"


using namespace cgp;


void scene_structure::initialize()
{
	// Basic set-up
	// ***************************************** //

	global_frame.initialize(mesh_primitive_frame(), "Frame");
	vec3 from = { 2.0f, 1.5f, 3.0f };
	vec3 to = { 0, 0, 0 };
	vec3 up = { 0, 1, 0 }; 
	float distance = 1.1f;
	environment.camera.look_at(distance * from, to, up);

	initialize_sph();

	sphere_particle.initialize(mesh_primitive_sphere(), "Sphere particle");
	// RADIUS
	sphere_particle.transform.scaling = 0.01f;
	sphere_particle.shading.phong = { 1,1,1 };
	sphere_particle.shading.color = { 0,0,1 };

	curve_visual.color = { 1,0,0 };
	curve_visual.initialize(curve_primitive_circle(), "Curve circles");

	// Edges of the cube
	buffer<vec3> cube_wireframe_data = { {-1,-1,-1},{1,-1,-1}, {1,-1,-1},{1,1,-1}, {1,1,-1},{-1,1,-1}, {-1,1,-1},{-1,-1,-1},
		{-1,-1,1} ,{1,-1,1},  {1,-1,1}, {1,1,1},  {1,1,1}, {-1,1,1},  {-1,1,1}, {-1,-1,1},
		{-1,-1,-1},{-1,-1,1}, {1,-1,-1},{1,-1,1}, {1,1,-1},{1,1,1},   {-1,1,-1},{-1,1,1} };
	cube_wireframe.initialize(cube_wireframe_data, "cube wireframe");
}

void scene_structure::initialize_sph()
{
	// Initial particle spacing (relative to h)
	float const c = 0.7; // < 1
	float const h = sph_parameters.h;

	float minX = -0.5f, maxX = 0.5f;
	float minY = -0.5f, maxY = 0.5f;
	float minZ = -0.5f, maxZ = 0.5f;

	// Add particles to the box
	particles.clear();
	for (float x = minX; x < maxX; x = x + c * h)
	{
		for (float y = minY; y < maxY; y = y + c * h)
		{
			for (float z = minZ; z < maxZ; z = z + c * h)
			{
				particle_element particle;
				particle.p = { 	x + h / 8.0 * rand_interval(),
								y + h / 8.0 * rand_interval(),
								z + h / 8.0 * rand_interval()
							 }; 
				particles.push_back(particle);
			}
		}
	}

}



void scene_structure::display()
{
	timer.update(); // update the timer to the current elapsed time
	float const dt = 0.005f * timer.scale;
	
	simulate(dt, particles, sph_parameters, 
		environment.camera.up(), gui.fix_gravity_direction, gui.zero_gravity);

	// Basics common elements
	// ***************************************** //
	environment.light = environment.camera.position();
	

	if (gui.display_particles) {
		for (int k = 0; k < particles.size(); ++k) {
			vec3 const& p = particles[k].p;
			sphere_particle.transform.translation = p;
			draw(sphere_particle, environment);
		}
	}


	// Display the box in which the particles should stay
	draw(cube_wireframe, environment);
}


void scene_structure::display_gui()
{
	ImGui::SliderFloat("Timer scale", &timer.scale, 0.01f, 4.0f, "%0.2f");

	bool const restart = ImGui::Button("Restart");
	if (restart)
		initialize_sph();

	ImGui::Checkbox("Fix Gravity Direction", &gui.fix_gravity_direction);
	ImGui::Checkbox("Zero Gravity", &gui.zero_gravity);
}


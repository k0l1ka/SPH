#pragma once

#include "cgp/cgp.hpp"
#include "simulation/simulation.hpp"


// The element of the GUI that are not already stored in other structures
struct gui_parameters {
	bool display_color = true;
	bool display_particles = true;
	// bool display_radius = false;
	bool fix_gravity_direction = true;
	bool zero_gravity = false;


};



// The structure of the custom scene
struct scene_structure {
	
	// ****************************** //
	// Elements and shapes of the scene
	// ****************************** //

	cgp::mesh_drawable global_frame;          // The standard global frame
	cgp::scene_environment_basic environment; // Standard environment controler
	gui_parameters gui;                       // Standard GUI element storage
	
	// ****************************** //
	// Elements and shapes of the scene
	// ****************************** //
	cgp::timer_basic timer;

	sph_parameters_structure sph_parameters; // Physical parameter related to SPH
	cgp::buffer<particle_element> particles;      // Storage of the particles
	
	cgp::mesh_drawable sphere_particle; // Sphere used to display a particle
	cgp::grid_2D<cgp::vec3> sphere_color;

	cgp::curve_drawable curve_visual;   // Circle used to display the radius h of influence

	cgp::segments_drawable cube_wireframe;

	// ****************************** //
	// Functions
	// ****************************** //

	void initialize();  // Standard initialization to be called before the animation loop
	void display();     // The frame display to be called within the animation loop
	void display_gui(); // The display of the GUI, also called within the animation loop

	void initialize_sph();
};






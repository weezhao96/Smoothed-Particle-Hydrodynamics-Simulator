/* SPH.h
    Header file for the SPH class definition.
*/

#pragma once

#include <iostream>
#include <exception>
#include <cmath>
#include <cstdlib>
#include <vector>

#include <time.h>

#include "boost/multi_array.hpp"

#include "MPI_Manager.h"
#include "IO_Manager.h"
#include "const_param.h"
#include "linalg.h"

using namespace std;

class SPH{

    public:

        // Struct for Interaction Data
        /*
            A struct assigned to every particle, containing information required
            for state update computations.

            Used in SPH::RunSimulation, struct members are computed by SPH::MapNeighbour.

        */

        struct interaction_set{
            
            // Particle ID
            int id; // Particle Index in Global 1D arrays

            // Neighbour
            int n_neighbour; // No. of neighbours
            vector<int> id_neighbour; // Particle Index of neighbouring
            // - used for pressure and viscous forcing computations

            // Interaction Data
            vector<double> q; // Vector of q_ij
            vector<double> dr; // Vector of dr_ij = r_i - r_j
            vector<double> dv; // Vector of dv_ij = v_i - v_j

        };

        // Class Constructor and Destructor (defined in SPH.cpp)
        SPH(int argc, char* argv[]);
        ~SPH();

        // Main Computational Algorithm (defined in SPH.cpp)
        void RunSimulation(); // Called by main() to start time-stepping simulation.
        void PerturbParticle(); // Apply perturbaution to particles to break "steady-state" arrangement.
        void BoundaryCheck(); // Checks if particles exceed domain boundary and applying damping.
        void MapNeighbour(vector<interaction_set> &inter_set); // Compute q_ij and map neighbouring particles h apart.
        void DensityPressureComputation(vector<interaction_set> &inter_set); // Computes density and pressure.
        void RescaleMassDensityPressure(); // Rescale the mass after first set of densities computed.
        void AccelComputation(vector<interaction_set> &inter_set); // Computes forcing variables.
        void FirstLeapfrogTimeStep(); // Apply leap-frong scheme time-stepping for first timestep.
        void LeapfrogTimeStep(); // Apply leap-frog scheme time-stepping algorithm.
        void EnergyComputation(); // Compute energy of system.

        // Function Member for Simulation Setup (defined in SPH_init.cpp)
        void GetManager(); // Init MPI and I/O Managers.
        void GetSimParam(int argc, char* argv[]); // Read command line arguments and set up config.
        void SetInitCond(); // Set initial condition of particle condition based of init_cond.
        void SetDamBreak(); // dam-break initial condition.
        void SetBlockDrop(); // block-drop initial condition.
        void SetDroplet(); // droplet initial condition.
        void SetTestOneParticle();
        void SetTestTwoParticle();
        void SetTestThreeParticle();
        void SetTestFourParticle();
        void RecursionBuild(int loop_depth, int* loop_lim, int* loop_index,
                            const double (*bound)[2], int& c_n, double* x);
                            // Imitates a nested loop with loop depth parameterised by no. of spatial dimensions.

    private:

        // Managers
        IO_Manager* io_manager = nullptr;
        MPI_Manager* mpi_manager = nullptr;

        // Simulation Time Parameters
        double dt; // Time Step
        double t_end; // Total Simulation Time
        double t; // Current Simulation Time
        
        // Simulation Particle Parameters
        double h; // Particle Radius of Influence
        double k; // Gas Constant
        double rho_0; // Resting Density
        double mu; // Viscosity
        double e; // Coefficient of Restitution
        double m; // Particle Mass

        // Gravity and Initial Condition
        double g; // Gravitational Forcing
        int init_cond; // Initial Condition - specifies type of initial particle distribution

        // State Variables - Global
        int n_particle_G; // No. of Particles
        int *id_G = nullptr; // Particle ID
        double *x_G = nullptr; // Particle Position
        double *v_G = nullptr; // Particle Velocity
        double *a_G = nullptr; // Particle Acceleration
        double *rho_G = nullptr; // Particle Density
        double *p_G = nullptr; // Particle Pressure

        // State Variables - Local
        int n_particle; // No. of Particles
        int *id = nullptr; // Particle ID
        double *x = nullptr; // Particle Position
        double *v = nullptr; // Particle Velocity
        double *a = nullptr; // Particle Acceleration
        double *rho = nullptr; // Particle Density
        double *p = nullptr; // Particle Pressure

        // Energy Variables - Global Variable
        double *Ek_G = nullptr; // Kinetic Energy of each particle
        double *Ep_G = nullptr; // Potential Energy of each particle
        double *E_G = nullptr; // Total Energy of each particle
        double Ek_total_G; // Kinetic Energy
        double Ep_total_G; // Potential Energy
        double E_total_G; // Total Energy

        // Energy Variables - Local Variable
        double *Ek = nullptr; // Kinetic Energy of each particle
        double *Ep = nullptr; // Potential Energy of each particle
        double *E = nullptr; // Total Energy of each particle
        double Ek_total; // Kinetic Energy
        double Ep_total; // Potential Energy
        double E_total; // Total Energy

};
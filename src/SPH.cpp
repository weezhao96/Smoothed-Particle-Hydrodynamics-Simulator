/* SPH.cpp
    Main body of class definition for SPH class.
    Contains the definition for class constructor and destructor, and main computational algorithm for SPH.
*/

#include "SPH.h"
#include "MPI_Manager.h"
#include "IO_Manager.h"

using namespace std;

SPH::SPH(int argc, char* argv[]){

    /* Class Constructor
        Calls void::GetSimParam (parser for command line arguments and define essential simulation parameters), then
        calls void::SetInitCond to set up initial condtion of problem based on initial condition selected from 
        command line input.
    */

    // Get Managers
    GetManager();

    if (mpi_manager->GetWorldRank() == 0) cout << "-- SPH Solver and Managers instantiated. --" << endl;

    // Get Simulation Parameters
    if (mpi_manager->GetWorldRank() == 0) cout << "-- Define simulation parameters." << endl;

    GetSimParam(argc, argv);
    if (mpi_manager->GetWorldRank() == 0){

        io_manager->OutputArgs(argc, argv, mpi_manager->GetWorldSize());

    }

    if (mpi_manager->GetWorldRank() == 0) cout << "Definition completed. --" << endl;

    // Setup Initial Condition
    if (mpi_manager->GetWorldRank() == 0) cout << "-- Setup initial condition of simulation." << endl;

    SetInitCond();

    if (mpi_manager->GetWorldRank() == 0) cout << "Initial condtion setup completed. --" << endl;
    if (mpi_manager->GetWorldRank() == 0) cout << endl;

}

SPH::~SPH(){

    /*
        Class destructor
        Deallocation managers and array before SPH deallocation.
    */

    // Deallocate State Arrays
    delete[] id;
    delete[] x;
    delete[] v;
    delete[] a;
    delete[] rho;
    delete[] p;

    delete[] id_G;
    delete[] x_G;
    delete[] v_G;
    delete[] a_G;
    delete[] rho_G;
    delete[] p_G;

    // Deallocate Energy Arrays
    delete[] Ek;
    delete[] Ep;
    delete[] E;

    delete[] Ek_G;
    delete[] Ep_G;
    delete[] E_G;

    // Deallocate Managers
    if (mpi_manager->GetWorldRank() == 0){

        delete io_manager;

    }

    delete mpi_manager;

}

void SPH::RunSimulation(){

    /*
        Runs the time-stepping computation, outputs the results for every timestep.
    */

    if (mpi_manager->GetWorldRank() == 0) cout << "-- Run SPH Simulation. --" << endl;

    // Track Computational Period in Real-time
    time_t time_start, time_end, time_sim;

    if (mpi_manager->GetWorldRank() == 0) time(&time_start);

    // Perturb Particle State
    PerturbParticle();
    BoundaryCheck();

    if (mpi_manager->GetWorldRank() == 0) cout << "-- Start Time Stepping." << endl;

    // Time Parameters
    int t_count = 0; // No. of Output Files Generated
    double t_out = 0.0; // Interval Tracker for I/O Output

    // -- Generate struct for neighbour interactions
    vector<interaction_set> inter_set(n_particle);

    // Set Vector Reserve Size
    int max_neighbour = kissing_number[N_dim];
    
    int index_1D_G = mpi_manager->GetIndex(1);

    for (int i = 0; i < n_particle; ++i){
        
        inter_set[i].id = id_G[index_1D_G];
        inter_set[i].id_neighbour.reserve(max_neighbour);
        inter_set[i].q.reserve(max_neighbour);
        inter_set[i].dr.reserve(max_neighbour * N_dim);
        inter_set[i].dv.reserve(max_neighbour * N_dim);

        index_1D_G += 1;

    }

    // Compute inter_set and Rescale and Sync Mass, Density and Pressure
    MapNeighbour(inter_set);

    DensityPressureComputation(inter_set);
    RescaleMassDensityPressure();

    mpi_manager->Comm_Array_L2G2L(rho_G, rho, 1, 'D'); // Sync Density State between Processes
    mpi_manager->Comm_Array_L2G2L(p_G, p, 1, 'D'); // Sync Pressure State between Processes

    AccelComputation(inter_set);

    // Output Initial State
    mpi_manager->Comm_Array_L2G(a_G, a, 2, 'D'); // Gatherv Global Acceleration
    mpi_manager->Comm_Array_L2G(x_G, x, 2, 'D'); // Gatherv Global Position

    EnergyComputation();

    if (mpi_manager->GetWorldRank() == 0){

        io_manager->OutputState(n_particle_G, t_count, t, m, rho_G, x_G, v_G, a_G, 0);
        io_manager->OutputEnergy(n_particle_G, t_count, t, Ek_total_G, Ep_total_G, E_total_G);
        t_count += 1;

        cout << "Simulation time, t = " << t << endl;
        
    }

    // First TimeStepping and Sync States
    FirstLeapfrogTimeStep();
    BoundaryCheck();

    mpi_manager->Comm_Array_L2G2L(x_G, x, 2, 'D'); // Sync Position State between Processes
    mpi_manager->Comm_Array_L2G2L(v_G, v, 2, 'D'); // Sync Velocity State between Processes

    t += dt;
    t_out += dt;

    // -- Time Marching Loop
    while (t < (t_end - dt - eps)){

        // Update Interaction Set, Compute and Sync Thermodynamic States
        MapNeighbour(inter_set);
        DensityPressureComputation(inter_set);
        
        mpi_manager->Comm_Array_L2G2L(rho_G, rho, 1, 'D'); // Sync Density State between Processes
        mpi_manager->Comm_Array_L2G2L(p_G, p, 1, 'D'); // Pressure State between Processes

        AccelComputation(inter_set);

        // Output every Timesteps of t_out = 0.01
        if ( t_out >= 0.01 - eps ){

            t_out = 0;
            mpi_manager->Comm_Array_L2G(a_G, a, 2, 'D'); // Gatherv Global Acceleration
            EnergyComputation();

            if (mpi_manager->GetWorldRank() == 0){

                io_manager->OutputState(n_particle_G, t_count, t, m, rho_G, x_G, v_G, a_G, 0);
                io_manager->OutputEnergy(n_particle_G, t_count, t, Ek_total_G, Ep_total_G, E_total_G);

                cout << "Simulation time, t = " << t << endl;

            }
        
            t_count += 1;

        }

        // Proceed Time Step and Sync
        LeapfrogTimeStep();
        BoundaryCheck();

        mpi_manager->Comm_Array_L2G2L(x_G, x, 2, 'D'); // Sync Position State between Processes
        mpi_manager->Comm_Array_L2G2L(v_G, v, 2, 'D'); // Sync Velocity State between Processes

        t += dt;
        t_out += dt;

    }

    // Output Final State at t = t_end
    EnergyComputation();

    if (mpi_manager->GetWorldRank() == 0){

        io_manager->OutputState(n_particle_G, t_count, t, m, rho_G, x_G, v_G, a_G, 0);
        io_manager->OutputState(n_particle_G, t_count, t, m, rho_G, x_G, v_G, a_G, 1);
        io_manager->OutputEnergy(n_particle_G, t_count, t, Ek_total_G, Ep_total_G, E_total_G);

        cout << "Simulation time, t = " << t << endl;

        cout << "Time Stepping Completed. --" << endl;
        cout << endl;

        time(&time_end);
        time_sim = difftime(time_end, time_start);

        io_manager->OutputSimPeriod(time_sim);

    }
    
}

void SPH::PerturbParticle(){

    /*
        Apply position perturbation to particles, scaled from (-0.05 to 0.05)*h.
    */

    // -- Apply Perturbation

    // Unique RNG Seed
    srand(time(NULL) + rand() % (1 + mpi_manager->GetWorldRank())); // Desync Different Process with world_rank

    int index = 0;

    for (int n = 0; n < n_particle; ++n){

        for (int dim = 0; dim < N_dim; ++dim){

            x[index] = x[index] + h * static_cast<double>( rand() % 21 - 10  ) / 200.0 ;

            index += 1;

        }

    }

}

void SPH::BoundaryCheck(){

    /*
        Check B.C of simulation domain.
    */
    
    // -- Check Boundary Collision for Particles

    // Loop Index for Particle
    int index = 0; // [i_particle, dim]

    for (int n = 0; n < n_particle; ++n){
        for (int dim = 0; dim < N_dim; ++dim){

            // Lower Bound of Domain
            if (x[index] < (domain_box[dim][0] + h + eps)){

                x[index] = domain_box[dim][0] + h;

                if (fabs(v[index]) > eps){

                    v[index] *= -e;

                } else{

                    v[index] *= -1.0;

                }
            
            // Upper Bound of Domain
            } else if (x[index] > (domain_box[dim][1] - h - eps)){

                x[index] = domain_box[dim][1] - h;

                if (fabs(v[index]) > eps){

                    v[index] *= -e;

                } else{

                    v[index] *= -1.0;

                }
            
            }

            index += 1;

        }
    }
    
}

void SPH::MapNeighbour(vector<interaction_set> &inter_set){

    /*
        Computes the interaction dataset for neighbouring particles.

        Vector<interaction_set> inter_set:
        A vector holding an interaction_set sturcture for every particle in the process.

        Struct interaction_set:
            interaction_set.id - int
                local ID of particle in process, corresponds to arrangement in
                state arrays x, v, a, rho, etc

            interaction_set.n_neighbour - int
                no. of neighbours of a current particle

            interaction_set.id_neighbour - vector<int> - size = {n_neighbour}
                array of ID of neighbouring particles
            
            interaction_set.q - vector<double> - size = {n_neighbour}
                array of q_ij

            interaction_set.dr - vector<double> - size = {n_neighbour * N_dim}
                array of r_ij vector as [r_i1, r_i2, r_i3 ...]

            interaction_set.dv - vector<double> - size = {n_neighbour * N_dim}
                array of v_ij vector as [v_i1, v_i2, v_i3 ...]

            The structure stores only properties of neighbouring particles when
            q_ij < 1.0, and excludes itself.

    */

    // -- Map Neighbours
    // Index and ID
    int id_i = 0;
    int index_xi = 0; // [i_particle, dim]
    int index_xj_G = 0; // [i_particle, dim]

    // Loop Variables
    double q_ij = 0.0;
    double* x_i = new double[N_dim]; // Position of Particle i
    double* dr = new double[N_dim];; // x_i - x_j
    double* dv = new double[N_dim];; // v_i - v_j

    // Reset Interaction Set
    for (int i = 0; i < n_particle; ++i){

        inter_set[i].n_neighbour = 0;
        inter_set[i].id_neighbour.clear();
        inter_set[i].q.clear();
        inter_set[i].dr.clear();
        inter_set[i].dv.clear();

    }
    
    // Collision Checks
    for (int i = 0; i < n_particle; ++i){ // Loop through Local Array i

        id_i = id[i];
        index_xi = i * N_dim;

        // Assign x_i
        for (int dim = 0; dim < N_dim; ++dim){

            x_i[dim] = x[index_xi + dim];

        }

        for (int j = 0; j < n_particle_G; ++j){ // Loop through Global Array j
            
            if (id_i != j){

                index_xj_G = j * N_dim;
                q_ij = 0.0;

                // Compute dr and q
                for (int dim = 0; dim < N_dim; ++dim){

                    dr[dim] = x_i[dim] - x_G[index_xj_G + dim];
                    q_ij += pow(dr[dim],2);

                }

                // Collision Check
                q_ij = sqrt(q_ij) / h;

                if (q_ij < 1.0){

                    // Assign dv
                    for (int dim = 0; dim < N_dim; ++dim){

                        dv[dim] = v[index_xi + dim] - v_G[index_xj_G + dim];

                    }

                    // -- Add to interaction set
                    // Set i
                    inter_set[i].n_neighbour += 1;
                    inter_set[i].id_neighbour.push_back(id_G[j]);
                    inter_set[i].q.push_back(q_ij);

                    for (int dim = 0; dim < N_dim; ++dim){
                        
                        // Add to Set i
                        inter_set[i].dr.push_back(dr[dim]);
                        inter_set[i].dv.push_back(dv[dim]);

                    }

                }

            }

        }

    }

    // Deallocate Array
    delete[] x_i;
    delete[] dr;
    delete[] dv;

}

void SPH::DensityPressureComputation(vector<interaction_set> &inter_set){

    /*
        Computes density and pressure of particles.
    */

    // -- Density Computation
    double q = 0.0;
    interaction_set loop_set;

    for (int i = 0; i < n_particle; ++i){

        rho[i] = 0.0;
        loop_set = inter_set[i];

        // Interaction from Neighbours
        for (int j = 0; j < loop_set.n_neighbour; ++j){
            
            q = loop_set.q[j];
            rho[i] += pow((1 - pow(q,2)) , 3);

        }

        // Self-Interaction
        rho[i] += 1;

    }

    // rho <- rho * m * 4 / (pi * h^2)
    F77NAME(dscal)(n_particle, (4 * m / (pi * pow(h,2))), rho, 1);

    // -- Pressure Computation
    // p <- rho
    F77NAME(dcopy)(n_particle, rho, 1, p, 1);

    // p <- p - rho_0
    F77NAME(daxpy)(n_particle, -1.0, &rho_0, 0, p, 1);

    // p <- k * p
    F77NAME(dscal)(n_particle, k, p, 1);

    /*
    // Output rho: Uncomment for debugging
    for (int n = 0; n < n_particle; ++n){

        cout << "rho[" << n << "] = " << rho[n] << endl;

    }

    // Output p: Uncomment for debugging
    for (int n = 0; n < n_particle; ++n){

        cout << "p[" << n << "] = " << p[n] << endl;

    }
    */

}

void SPH::RescaleMassDensityPressure(){

    /*
        Apply mass correction to scale density and recompute pressure.

        Communicate rho_sum to via AllReduce.
    */

    // -- Compute Mass
    double rho_sum = 0.0;
    rho_sum = F77NAME(dasum)(n_particle, rho, 1);
    
    mpi_manager->Comm_RhoSum(rho_sum); // Get total summation of density from all process.

    m = n_particle_G * rho_0 / rho_sum;

    if (mpi_manager->GetWorldRank() == 0) cout << "Mass, m = " << m << endl;

    // -- Scale Density
    F77NAME(dscal)(n_particle, m, rho, 1);

    // -- Compute Pressure
    // p <- rho
    F77NAME(dcopy)(n_particle, rho, 1, p, 1);

    // p <- p - rho_0
    F77NAME(daxpy)(n_particle, -1.0, &rho_0, 0, p, 1);

    // p <- k * p
    F77NAME(dscal)(n_particle, k, p, 1);

    // Output rho: Uncomment for debugging
    /*
    for (int n = 0; n < n_particle; ++n){

        cout << mpi_manager->GetWorldRank() << " | rho[" << n << "] = " << rho[n] << endl;

    }

    // Output p: Uncomment for debugging
    for (int n = 0; n < n_particle; ++n){

        cout << mpi_manager->GetWorldRank() << " | p[" << n << "] = " << p[n] << endl;

    }
    */

}

void SPH::AccelComputation(vector<interaction_set> &inter_set){

    /*
        Acceleration computation for Pressure, Viscous and Gravitational Force
    */

    // -- Pressure and Viscous Force
    // Acceleration Array and Index
    double* a_p = new double[n_particle * N_dim]();
    double* a_mu = new double[n_particle * N_dim]();

    int index_a = 0; // [i_particle, dim]

    // Loop Variables
    int j = 0;

    interaction_set loop_set;

    double p_scale = 0.0;
    double mu_scale = 0.0;

    double q = 0.0;
    double* dr = new double[N_dim]();
    double* dv = new double[N_dim]();

    int index_rv = 0;

    double rho_i = 0.0;
    double rho_j = 0.0;

    double p_i = 0.0;
    double p_j = 0.0;

    // Loop
    for (int i = 0; i < n_particle; ++i){

        rho_i = rho[i];
        p_i = p[i];
        loop_set = inter_set[i];

        index_a = i * N_dim;

        for (int n = 0; n < loop_set.n_neighbour; ++n){

            j = loop_set.id_neighbour[n];

            q = loop_set.q[n];
            dr = loop_set.dr.data();
            dv = loop_set.dv.data();

            rho_j = rho_G[j];
            p_j = p_G[j];

            index_rv = n * N_dim;

            p_scale = (p_i + p_j) * pow(1-q,2) / (q * rho_i * rho_j);
            mu_scale = (1 - q) / (rho_i * rho_j);

            for (int dim = 0; dim < N_dim; ++dim){

                a_p[index_a + dim] += p_scale * dr[index_rv + dim];
                a_mu[index_a + dim] += mu_scale * dv[index_rv + dim];

            }

        }

    }

    // a_p <- 30 m / (2 pi h^3) a_p
    F77NAME(dscal)(n_particle * N_dim, 30 * m / (pi * pow(h,3) * 2), a_p, 1);

    // a_mu <- -40 mu m / (pi h^4) a_mu
    F77NAME(dscal)(n_particle * N_dim, -40 * mu * m / (pi * pow(h,4)), a_mu, 1);

    // a_mu <- a_p + a_mu
    F77NAME(daxpy)(n_particle * N_dim, 1.0, a_p, 1, a_mu, 1);
    
    // -- Gravitational Acceleration
    // a_mu <- a_mu[i][N_dim-1] + g
    F77NAME(daxpy)(n_particle, -1.0, &g, 0, (a_mu + N_dim - 1), N_dim);

    // a <- a_mu
    F77NAME(dcopy)(n_particle * N_dim, a_mu, 1, a, 1);

    delete[] a_p;
    delete[] a_mu;

}

void SPH::FirstLeapfrogTimeStep(){

    /*  
        Proceed with leap-frog time-stepping based on velocity timestep at t = 1/2.
    */

    // v <- v + a dt / 2
    F77NAME(daxpy)(n_particle * N_dim, dt / 2, a, 1, v, 1);

    // x <- x + v dt
    F77NAME(daxpy)(n_particle * N_dim, dt, v, 1, x, 1);

}

void SPH::LeapfrogTimeStep(){

    /*  
        Proceed with leap-frog time-stepping.
    */

    // v <- v + a dt
    F77NAME(daxpy)(n_particle * N_dim, dt, a, 1, v, 1);

    // x <- x + v dt
    F77NAME(daxpy)(n_particle * N_dim, dt, v, 1, x, 1);

}

void SPH::EnergyComputation(){

    /*
        Compute kinetic, potential and total energy of local process,
        and passes data to root process via MPI::CommEnergyL2G.
    */

    // -- Compute Ek and Ek_total
    // Init Loop Variables
    double Ek_loop = 0;

    // Loop
    int index = 0;

    for (int n = 0; n < n_particle; ++n){

        Ek_loop = 0;

        for (int dim = 0; dim < N_dim; ++dim){
            
            Ek_loop = Ek_loop + pow(v[index],2);
            index += 1;

        }

        Ek[n] = Ek_loop;

    }

    // Ek <- 0.5 * m * Ek
    F77NAME(dscal)(n_particle, 0.5 * m, Ek, 1);

    // Ek_total <- sum(Ek)
    Ek_total = F77NAME(dasum)(n_particle, Ek, 1);

    // -- Compute Ep and Ep_total
    // Ep <- m * g * x[i_particle][N_dim-1]
    F77NAME(dcopy)(n_particle, (x + N_dim - 1), N_dim, Ep, 1);
    F77NAME(dscal)(n_particle, m * g, Ep, 1);

    // Ep_total <- sum(Ep)
    Ep_total = F77NAME(dasum)(n_particle, Ep, 1);

    // -- Compute Total Energy
    // E <- Ep
    F77NAME(dcopy)(n_particle, Ep, 1, E, 1);

    // E <- Ek + E (data of Ep)
    F77NAME(daxpy)(n_particle, 1.0, Ek, 1, E, 1);

    E_total = Ek_total + Ep_total;

    // Perform Message Passing to Root (Reduce)
    mpi_manager->Comm_Energy_L2G(Ek_total_G, Ek_total, Ep_total_G, Ep_total, E_total_G);

    // E_total_G
    if (mpi_manager->GetWorldRank() == 0){

        E_total_G = Ek_total_G + Ep_total_G;

    }

}
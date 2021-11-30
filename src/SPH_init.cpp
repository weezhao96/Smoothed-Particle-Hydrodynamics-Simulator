/* SPH_init.cpp
    Contains member function definition to setup simulation based on required initial condtion for SPH class.
    Called by constructor in SPH::SPH.
*/

#include "SPH.h"

void SPH::GetManager(){

    /*
        Instantiate MPI Process and I/O Managers.
    */

    // IO & MPI Manager
    mpi_manager = new MPI_Manager();

    if (mpi_manager->GetWorldRank() == 0){

        io_manager = new IO_Manager("output_folder/");

    }

}

void SPH::GetSimParam(int argc, char* argv[]){

    /*
        Function to acts as setter for simulation parameters and parse command line arugments.
        If command line specified undefined initial condition, it throws a runtime error.

        Initial Condition:
        int - init_cond - specifies type of initial partile distribution
        (1) ic-dam-break; (2) ic-block-drop; (3) ic-droplet
        (4) ic-one-particle; (5) ic-two-particles; (6) ic-four-particles

        Instantiated as 0, unless matched command line parameters.
    */

    // Simulation Particle Parameters
    k = 2000.0; // Gas Constant
    rho_0 = 1000.0; // Resting Density
    mu = 1.0; // Viscosity
    h = 0.01; // Radius of Influence (can be overwritten by command line arg)
    e = 0.5; // Coefficient of Restitution
    m = 1.0; // Particle Mass

    // Gravitational Parameter
    g = 9.81;

    // Initial Condition
    init_cond = 0;
    t = 0; // Current Simulation time = 0

    // Parse Command Line Arguments
    for (int i = 0; i < argc; ++i){

        if (init_cond == 0){

            if (std::string(argv[i]) == "--ic-dam-break"){

                init_cond = 1;

            } else if (std::string(argv[i]) == "--ic-block-drop"){

                init_cond = 2;

            } else if (std::string(argv[i]) == "--ic-droplet"){

                init_cond = 3;

            } else if (std::string(argv[i]) == "--ic-one-particle"){

                init_cond = 4;

            } else if (std::string(argv[i]) == "--ic-two-particles"){

                init_cond = 5;

            } else if (std::string(argv[i]) == "--ic-three-particles"){

                init_cond = 6;

            } else if (std::string(argv[i]) == "--ic-four-particles"){

                init_cond = 7;

            }

        }
    
        if (std::string(argv[i]) == "--dt"){

            dt = atof(argv[i+1]);

        } else if(std::string(argv[i]) == "--T"){

            t_end = atof(argv[i+1]);

        } else if(std::string(argv[i]) == "--h"){

            h = atof(argv[i+1]);

        }

    }
    
    // Display Simulation Properties
    if (mpi_manager->GetWorldRank() == 0){
        
        // Display Initial Condition
        switch (init_cond)
        {
        case 1:
            cout << "Simulation Mode: ic-dam-break." << endl;
            break;

        case 2:
            cout << "Simulation Mode: ic-block-drop." << endl;
            break;

        case 3:
            cout << "Simulation Mode: ic-droplet." << endl;
            break;

        case 4:
            cout << "Simulation Mode: ic-one-particle." << endl;
            break;

        case 5:
            cout << "Simulation Mode: ic-two-particles." << endl;
            break;

        case 6:
            cout << "Simulation Mode: ic-three-particles." << endl;
            break;

        case 7:
            cout << "Simulation Mode: ic-four-particles." << endl;
            break;

        default:
            break;
        }

        cout << "Time Step, dt = " << dt << "." << endl;
        cout << "Simulation Period, t_end = " << t_end << "." << endl;
        cout << "Particle Radius of Influence, h = " << h << "." << endl;

    }

    // ThroW Exception if Initial Conditions Invalid
    if (init_cond == 0){
        throw std::logic_error("Invalid Initial Condition.");
    }

}

void SPH::SetInitCond(){

    /*
        Function to setup initial condtion based on init_cond variable.

        Computes the total number of simulation particles and broadcast the
        results. Then a call to MPI_Manager::RunLoadBalance is performed to
        identify the local number of particles for each process.
    */

    // -- Compute No. of Particles and Initial State - Global
    if (mpi_manager->GetWorldRank() == 0){

        switch (init_cond)
        {
        case 1:
            SetDamBreak();
            break;

        case 2:
            SetBlockDrop();
            break;

        case 3:
            SetDroplet();
            break;

        case 4:
            SetTestOneParticle();
            break;

        case 5:
            SetTestTwoParticle();
            break;

        case 6:
            SetTestThreeParticle();       
            break;

        case 7:
            SetTestFourParticle();       
            break;

        default:
            break;
        }

        cout << "No. of Global Particles, n_particle_G = " << n_particle_G << endl;

    }

    // -- Broadcast n_particle_G and Compute n_particle (Load Balancing)
    mpi_manager->RunLoadBalance(n_particle_G, n_particle);

    // -- Broadcast x_G and Define Global Arrays
    // Broadcast x_G
    if (mpi_manager->GetWorldRank() != 0) x_G = new double[n_particle_G * N_dim];

    mpi_manager->Comm_Bcast(x_G, n_particle_G * N_dim, 0, 'D');

    // Define 2-D Arrays
    v_G = new double[n_particle_G * N_dim]();
    a_G = new double[n_particle_G * N_dim]();

    // Define 1-D Arrays
    Ek_G = new double[n_particle_G]();
    Ep_G = new double[n_particle_G]();
    E_G = new double[n_particle_G]();

    rho_G = new double[n_particle_G]();
    p_G = new double[n_particle_G]();

    // Assign Global ID of Particles
    id_G = new int[n_particle_G]();

    for (int i = 0; i < n_particle_G; ++i){

        id_G[i] = i;

    }

    // Define Energy Scalars
    Ek_total_G = 0.0;
    Ep_total_G = 0.0;
    E_total_G = 0.0;

    // -- Define Local Arrays
    // Define 2-D Arrays
    v = new double[n_particle * N_dim]();
    a = new double[n_particle * N_dim]();

    // Define 1-D Arrays
    Ek = new double[n_particle]();
    Ep = new double[n_particle]();
    E = new double[n_particle]();

    rho = new double[n_particle]();
    p = new double[n_particle]();

    id = new int[n_particle]();

    // Define Energy Scalars
    Ek_total = 0.0;
    Ep_total = 0.0;
    E_total = 0.0;

    // -- Assign x and id from x_G, id_G
    // Allocation
    x = new double[n_particle * N_dim];
    id = new int[n_particle * N_dim];

    // Loop Index - Get Start of Local Particle in Global Array,  index_v[world_rank]
    int index_1D_G = mpi_manager->GetIndex(1);
    int index_2D_G = mpi_manager->GetIndex(2);
    int index_2D = 0;

    // Loop
    for (int i = 0; i < n_particle; ++i){

        id[i] = id_G[index_1D_G];
        index_1D_G += 1;

        for (int dim = 0; dim < N_dim; ++dim){

            x[index_2D] = x_G[index_2D_G];

            index_2D += 1;
            index_2D_G += 1;

        }

    }

}

void SPH::SetDamBreak(){

    /*
        Setup initial condtion based on dam-break initial condition.
    */

    // Check Overlap between Initial Particle Domain with Domain Box
    double init_domain[N_dim][2] = { };

    for (int dim = 0; dim < N_dim; ++dim){

        if (box_dam_break[dim][0] <= domain_box[dim][0] + h){
            
            init_domain[dim][0] = box_dam_break[dim][0] + h;

        } else{

            init_domain[dim][0] = box_dam_break[dim][0];

        }

        if (box_dam_break[dim][1] >= domain_box[dim][1] - h){
            
            init_domain[dim][1] = box_dam_break[dim][1] - h;

        } else{

            init_domain[dim][1] = box_dam_break[dim][1];

        }

    }

    // No. of Particles
    int n_particle_dim[N_dim];

    for (int i = 0; i < N_dim; ++i){
        n_particle_dim[i] = static_cast<int>(floor( (init_domain[i][1] - init_domain[i][0] + eps) / h ) + 1);
    }

    n_particle_G = 1;

    for (int i = 0; i < N_dim; ++i){
        n_particle_G *= n_particle_dim[i];
    }

    // Specify Initial Position
    x_G = new double[n_particle_G * N_dim];

    int dim_index[N_dim] = {};
    int c_n = 0; // n-th particle

    RecursionBuild(0, n_particle_dim, dim_index, init_domain, c_n, x_G);

}

void SPH::SetBlockDrop(){

    /*
        Setup initial condtion based on block-drop initial condition.
    */

    // Check Overlap between Initial Particle Domain with Domain Box
    double init_domain[N_dim][2] = { };

    for (int dim = 0; dim < N_dim; ++dim){

        if (box_block_drop[dim][0] <= domain_box[dim][0] + h){
            
            init_domain[dim][0] = box_block_drop[dim][0] + h;

        } else{

            init_domain[dim][0] = box_block_drop[dim][0];

        }

        if (box_block_drop[dim][1] >= domain_box[dim][1] - h){
            
            init_domain[dim][1] = box_block_drop[dim][1] - h;

        } else{

            init_domain[dim][1] = box_block_drop[dim][1];

        }

    }
   
    // No. of Particles
    int n_particle_dim[N_dim];

    for (int i = 0; i < N_dim; ++i){
        n_particle_dim[i] = static_cast<int>(floor( (init_domain[i][1] - init_domain[i][0] + eps) / h) + 1);
    }

    n_particle_G = 1;

    for (int i = 0; i < N_dim; ++i){
        n_particle_G *= n_particle_dim[i];
    }

    // Specify Initial Position
    x_G = new double[n_particle * N_dim];

    int dim_index[N_dim] = {};
    int c_n = 0; // n-th particle

    RecursionBuild(0, n_particle_dim, dim_index, init_domain, c_n, x_G);

}

void SPH::SetDroplet(){

    /*
        Setup initial condtion based on droplet initial condtion.

        A set of concentric circles with radius h apart are filled is particles.
        On each circle, the particles are spaced at an angle dtheta apart which leads to
        || x_i - x_j || ~ h.
    */

    // Dimension Check
    if (N_dim != 2){
        throw std::logic_error("Droplet Config only supported for N_dim = 2.");
    }

    // Domain of Initial Particle: Droplet
    double rad_droplet = 0.1; // Radius of Droplet
    double center_droplet[N_dim][1] = {{0.5},{0.7}}; // Center of Droplet

    // Radius Interval - Generate radius of concentric circles
    int n_rad = static_cast<int>(floor( (rad_droplet + eps) / h)); // No. of concetric circles
    double* rad_int = new double[n_rad]; // Radius of concentric circles

    for (int r = 0; r < n_rad; ++r){

        rad_int[r] = rad_droplet - r * h;

    }

    // No. of Particle
    double* h_theta = new double[n_rad]; // Theta interval for h spacing
    double* n_int = new double[n_rad]; // No. of particles for each concentric circle
    double theta = 0.0;

    n_particle_G = 0;

    for (int r = 0; r < n_rad; ++r){

        theta = acos(1 - 0.5 * pow((h/rad_int[r]),2));
        n_int[r] = static_cast<int>( floor((2 * pi) / theta) );
        n_particle_G += n_int[r];

        h_theta[r] = (2 * pi) / static_cast<double>(n_int[r]);

    }

    n_particle_G += 1; // Add particle at droplet center

    // Specify Initial Position
    x_G = new double[n_particle_G * N_dim];

    double c_theta; // Loop variable, current theta
    double rad; // Circle radius
    int index_x = 0;

    for (int r = 0; r < n_rad; ++r){

        rad = rad_int[r];
        c_theta = -pi / 2.0;

        for (int n = 0; n < n_int[r]; ++n){
            
            // x0 position
            x_G[index_x] = rad * cos(c_theta) + center_droplet[0][0];
            index_x += 1;

            // x1 position
            x_G[index_x] = rad * sin(c_theta) + center_droplet[1][0];
            index_x += 1;

            c_theta += h_theta[r];

        }

    }

    // Particle at Droplet Center
    x_G[index_x] = center_droplet[0][0];
    index_x += 1;
    x_G[index_x] = center_droplet[1][0];

}

void SPH::SetTestOneParticle(){

    /*
        Setup initial condtion for one particle test.
    */

    // No. of Particle
    n_particle_G = 1;

    // Specify Initial Position
    x_G = new double[n_particle_G * N_dim];

    for (int dim = 0; dim < N_dim; ++dim){
        x_G[dim] = init_pos_1[dim][0];
    }

}

void SPH::SetTestTwoParticle(){

    /*
        Setup initial condtion for two particle test.
    */

    // Dimension Check
    if (N_dim != 2){
        throw std::logic_error("Droplet Config only supported for N_dim = 2.");
    }

    // No. of Particle
    n_particle_G = 2;

    // Init Position
    double init_pos_2[N_dim][n_particle_G] = {{0.5,0.5},{0.5,h}};

    // Specify Initial Position
    x_G = new double[N_dim * n_particle_G];

    int index = 0;
    for (int n = 0; n < n_particle_G; ++n){
        for (int dim = 0; dim < N_dim; ++dim){

            x_G[index] = init_pos_2[dim][n];
            index += 1;

        }
    }

}

void SPH::SetTestThreeParticle(){

    /*
        Setup initial condtion for three particle test.
    */

    // Dimension Check
    if (N_dim != 2){
        throw std::logic_error("Droplet Config only supported for N_dim = 2.");
    }

    // No. of Particle
    n_particle_G = 3;

    // Init Position
    double init_pos_3[N_dim][n_particle_G] = {{0.5,0.495,0.505},{0.5,h,h}};

    // Specify Initial Position
    x_G = new double [n_particle_G * N_dim];

    int index = 0;
    
    for (int n = 0; n < n_particle_G; ++n){
        for (int dim = 0; dim < N_dim; ++dim){
            
            x_G[index] = init_pos_3[dim][n];
            index += 1;

        }
    }

}

void SPH::SetTestFourParticle(){

    /*
        Setup initial condtion for four particle test.
    */

    // No. of Particle
    n_particle_G = 4;

    // Specify Initial Position
    x_G = new double[n_particle_G * N_dim];

    int index = 0;
    for (int n = 0; n < n_particle_G; ++n){
        for (int dim = 0; dim < N_dim; ++dim){

            x_G[index] = init_pos_4[dim][n];
            index += 1;

        }    
    }

}

void SPH::RecursionBuild(int loop_depth, int* loop_lim, int* loop_index, 
                         const double (*bound)[2], int& c_n, double* x){
    
    /*  

        loop_index = c_index
        loop_lim = n_particle_dim

        Fills the particle position array with spaced position based on the box bounded initial condition,
        with a equi-spaced length h.

        The recursion imitates a "nested-loop" but with loop depth parameterised by the no. of spatial dimensions.
        The deepest recursion call, fills the x array with particle position, while the outer calls act to permute
        the postion corresponding to its spatial dimension.

        Example of nested loop for a 3D problem, with bound = {{0.0, 0.0, 0.0} , {1.0, 1.0, 1.0}}:

        c_n = 0; // current n-th particle

        for (int i = 0; i < loop_lim[0]; ++i){
            for (int j = 0; j < loop_lim[1]; ++j){
                for (int k = 0; k < loop_lim[2]; ++k){

                    loop_index = {i,j,k};
                    index = c_n * N_dim;

                    x[index] = 0.0 + i * (1.0 - 0.0) / h;
                    x[index + 1] = 0.0 + j * (1.0 - 0.0) / h;
                    x[index + 2] = 0.0 + k * (1.0 - 0.0) / h;

                    c_n += 1;

                }
            }
        }

    */

    while (loop_index[loop_depth] < loop_lim[loop_depth]){

        if (loop_depth == (N_dim-1)){

            int index = c_n * N_dim;

            // Fill Array
            for (int dim = 0; dim < N_dim; ++dim){
                
                x[index] = bound[dim][0] + h * loop_index[dim];

                index += 1;

            }

            c_n += 1; // Increment n-th particle

        } else{

            RecursionBuild(loop_depth + 1, loop_lim, loop_index, bound, c_n, x);

        }

        loop_index[loop_depth] += 1; // Increment spatial index of loop_depth

    }

    loop_index[loop_depth] %= loop_lim[loop_depth];

    return;

}
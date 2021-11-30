/* IO_Manager.cpp
    Main body of class definition for IO_Manager class.
    Manages I/O of SPH results.
*/

#include "IO_Manager.h"
#include "const_param.h"

using namespace std;

IO_Manager::IO_Manager(string path = "output_folder/"){

    /* Class Constructor
        Output folder checks (overwrite or create).
        Create ofstream for I/O of state and energy of SPH particles.
    */

    // Define Output Paths
    output_path = path;
    subpath_state = output_path + "state/";
    subpath_energy = output_path;
    
    // Check Folder Path and Init ofstream
    Init_OutputFolder();

}

IO_Manager::~IO_Manager(){

    /* Class Destructor
        Close fstreams.
    */

    io_energy.close();

}

void IO_Manager::Init_OutputFolder(){

    /*
        Init ofstream for the energy and state filepath.
    */

    // Folder Check (Overwrite and Creation)
    if (boost::filesystem::exists(output_path)){

        boost::filesystem::remove_all(output_path);
        boost::filesystem::create_directory(output_path);
        boost::filesystem::create_directory(subpath_state);

    }else{

        boost::filesystem::create_directory(output_path);
        boost::filesystem::create_directory(subpath_state);

    }
    
    // Init Energy Output Stream
    string file_energy = subpath_energy + "energy.txt";
    io_energy.open(file_energy, ios::out);

    // Particle Energy Output File Header
    io_energy << setw(10) << "t";
    io_energy << setw(15) << "E_kinetic";
    io_energy << setw(15) << "E_potential";
    io_energy << setw(15) << "E_total" << endl;

}

void IO_Manager::OutputArgs(int argc, char* argv[], int world_size){

    /*
        Output terminal arguments.
    */

    // Init Output Stream
    string file_args;
    file_args = output_path + "sim_args.txt";

    io_args.open(file_args, ios::out);

    for (int i = 0; i < argc; ++i){

        io_args << argv[i] << " ";

    }

    io_args << "| --np " << world_size << " ";

    cout << endl;

    io_args.close();

}   

void IO_Manager::OutputState(int &n_particle, int &t_count, double &t,
                             double &m, double* rho,
                             double *x, double* v, double* a,
                             int bool_end){

    /*
        Generates formatted output for state variables (position, velocity & acceleration)
        at specified time steps.

        Filename : <subpath_state>t_step_<t_count>_state.txt if bool_end == 0
                   <output_path>output.txt if bool_end == 1 (for output of final particle state)
    */
   
    // Init Output Stream
    string file_state;

    if (bool_end == 1){

        file_state = output_path + "output.txt";

        io_state.open(file_state, ios::out);
        io_state.precision(7);

        // -- Output File Header
        // x-state
        for (int dim = 0; dim < N_dim; ++dim){
            io_state << setw(15) << "x_" + to_string(dim);
        }

        io_state << endl;

        // -- Output Array Data
        int index = 0;

        for (int n = 0; n < n_particle; ++n){

            // x-state
            index = n * N_dim;

            for (int dim = 0; dim < N_dim; ++dim){

                io_state << setw(15) << x[index];
                index += 1;

            }

            io_state << endl;

        }


    } else{

        file_state = subpath_state + "t_step_" + to_string(t_count) + "_state.txt";

        io_state.open(file_state, ios::out);
        io_state.precision(7);

        // -- Output File Header
        io_state << setw(10) << "t = " << to_string(t) << endl;
        io_state << setw(10) << "m = " << to_string(m) << endl;
        io_state << setw(10) << "n";

        io_state << setw(15) << "rho";

        // x-state
        for (int dim = 0; dim < N_dim; ++dim){
            io_state << setw(15) << "x_" + to_string(dim);
        }

        // v-state
        for (int dim = 0; dim < N_dim; ++dim){
            io_state << setw(15) << "v_" + to_string(dim);
        }

        // a-state
        for (int dim = 0; dim < N_dim; ++dim){
            io_state << setw(15) << "a_" + to_string(dim);
        }

        io_state << endl;

        // -- Output Array Data
        int index = 0;

        for (int n = 0; n < n_particle; ++n){
        
            // n-th particle
            io_state << setw(10) << n;

            // Density
            io_state << setw(15) << rho[n];

            // x-state
            index = n * N_dim;

            for (int dim = 0; dim < N_dim; ++dim){

                io_state << setw(15) << x[index];
                index += 1;

            }

            // v-state
            index = n * N_dim;

            for (int dim = 0; dim < N_dim; ++dim){

                io_state << setw(15) << v[index];
                index += 1;

            }

            // a-state
            index = n * N_dim;
            
            for (int dim = 0; dim < N_dim; ++dim){

                io_state << setw(15) << a[index];
                index += 1;

            }

            io_state << endl;

        }

        io_state.close();

    }

}

void IO_Manager::OutputEnergy(int &n_particle, int &t_count, double &t,
                              double &Ek_total, double &Ep_total, double &E_total){

    /*
        Generates formatted output for energy variables (kinetic, potential, total)
        at specified time steps.
    */

    io_energy.precision(7);

    // Output Data
    io_energy << setw(10) << t;
    io_energy << setw(15) << Ek_total;
    io_energy << setw(15) << Ep_total;
    io_energy << setw(15) << E_total;
    io_energy << endl;

}

void IO_Manager::OutputSimPeriod(time_t &sim_time){

    /*
        Output simulation period to energy file.
    */

    io_energy << "Simulation Period, T = " << sim_time << " s." << endl;

}
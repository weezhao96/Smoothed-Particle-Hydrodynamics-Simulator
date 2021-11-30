/* IO_Manager.h
    Header file for IO_Manager class definition.
*/

#pragma once

#include <fstream>
#include <iostream>
#include <iomanip>

#include <time.h>

#include "boost/filesystem.hpp"

using namespace std;

class IO_Manager{

    public:

        // Class Constructor and Destructor
        IO_Manager(string output_path);
        ~IO_Manager();

        // Init I/O Stream
        void Init_OutputFolder();

        // Output Arguments
        void OutputArgs(int argc, char* argv[], int world_size); // Output terminal arguments.

        // Output Results
        void OutputState(int& n_particle, int &t_count, double &t,
                         double& m, double* rho,
                         double* x, double* v, double* a,
                         int bool_end); // Output state for every step
        void OutputEnergy(int &n_particle, int &t_count, double &t, 
                          double &Ek, double &Ep, double &E); // Output energy time evolution.
        void OutputSimPeriod(time_t &sim_time); // Append simluation period to energy file.

    private:
        
        // Output Path
        string output_path; // Main Output Path.
        string subpath_state; // Subdirectory for state output
        string subpath_energy; // Subdirectory for energy output.

        // I/O Stream
        ofstream io_args;
        ofstream io_state;
        ofstream io_energy;

};
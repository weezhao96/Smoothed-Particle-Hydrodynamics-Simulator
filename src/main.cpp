/* SPH_main.cpp
    Main cpp function to instantiate SPH class, calls SPH::RunSimulation to start the simulation and handles
    deallocation of SPH class when complete.
*/

#include "SPH.h"

using namespace std;

int main(int argc, char* argv[]){

    MPI_Init(NULL, NULL);

    SPH *solver = new SPH(argc, argv);

    solver->RunSimulation();

    delete solver;

    MPI_Finalize();

}
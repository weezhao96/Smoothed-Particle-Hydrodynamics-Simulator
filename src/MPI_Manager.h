/* MPI_Manager.h
    Header file for MPI_Manager class.
*/

#pragma once

#include <iostream>
#include <vector>

#include <mpi.h>
#include <unistd.h>

#include "const_param.h"

using namespace std;

class MPI_Manager{

    public:

        // Class Constructor and Destructor
        MPI_Manager();
        ~MPI_Manager();

        // Getter
        int GetWorldRank();
        int GetWorldSize();
        int GetIndex(int nd_array);

        // Define Global and Local No. of Particles and Local Index (called by SPH::SetInitCond)
        void RunLoadBalance(int &n_particle_G, int &n_particle);
        
        // AllGatherv for Arrays
        void Comm_Array_L2G2L(void* array_G, void* array_L, int nd_array, char mpi_T);

        // Scatterv / Gatherv for Arrays
        void Comm_Array_G2L(void* array_G, void* array_L, int nd_array, char mpi_T); // Global -> Local Scatterv
        void Comm_Array_L2G(void* array_G, void* array_L, int nd_array, char mpi_T); // Local -> Global Gatherv

        // Message Passing for Mass Compute (called by SPH::RescaleMassDensityPressure)
        void Comm_RhoSum(double &rho_sum);

        // Message Passing for Energy (called by SPH::EnergyComputation)
        void Comm_Energy_L2G(double &Ek_total_G, double Ek_total,
                             double &Ep_total_G, double Ep_total,
                             double &E_total_G); // Local -> Global

        // Broadcast
        void Comm_Bcast(void* p_data, int count, int root_rank, char mpi_T);

        // SyncProcess - MPI_Barrier
        void SyncProcess();

    private:

        // MPI Parameters
        int world_rank;
        int world_size;

        // 1D Array Scatterv and Gatherv Count and Index (for rho, p, Ek, Ep, E)
        int* count_1D_v = nullptr; // count_1D_v[world_rank] = n_particle_L
        int* index_1D_v = nullptr; // index_1D_v[world_rank] = sum(count_1D_v[0: world_rank - 1])

        // 2D Array Scatterv and Gatherv Count and Index (for x, v, a)
        int* count_2D_v = nullptr; // count_2D_v[world_rank] = n_particle_L * N_dim
        int* index_2D_v = nullptr; // index_2D_v[world_rank] = sum(count_2D_v[0: world_rank - 1])

};
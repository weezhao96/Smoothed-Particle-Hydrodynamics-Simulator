/* MPI_Manager.cpp
    Main body of class definition for MPI_Manager class.
    Manages initialisation for message passing and communication when the time-stepping algorithm runs.
*/

#include "MPI_Manager.h"

using namespace std;

MPI_Manager::MPI_Manager(){

    /*
        Instantiate MPI and Process Parameters.
    */

    // Define Rank and Size
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    cout << "Rank = " << world_rank << " | World Size = " << world_size << endl;

}

MPI_Manager::~MPI_Manager(){

    /*
        Deallocation of arrays.
    */

    // Deallocation
    delete[] count_1D_v;
    delete[] count_2D_v;

    delete[] index_1D_v;
    delete[] index_2D_v;

}

int MPI_Manager::GetWorldRank(){

    /*
        Return world rank of process.
    */

    return world_rank;

}

int MPI_Manager::GetWorldSize(){

    /*
        Return world size of MPI_comm.
    */

    return world_size;

}

int MPI_Manager::GetIndex(int nd_array){

    /*
        Return starting index of world_rank in global array.
    */

    if (nd_array == 1){

        return index_1D_v[world_rank];

    } else if (nd_array == 2){

        return index_2D_v[world_rank];

    } else{

        return -1;

    }

}

void MPI_Manager::RunLoadBalance(int &n_particle_G, int &n_particle){

    /*
        Compute Local No. of Particles (Load Balancing) for each process
        and
        Compute count_v and index_v : count and index offset parameters for MPI_{}v operations

        count_v : int* - size = {world_size}
            count_v[world_rank] = n_particle_L

        index_v : int* - size = {world_size}
            index_v[world_rank] : index of first local particle in global 1D array (ie: p, rho)

        count_dim_v : int* - size = {world_size}
            count_dim_v[world_rank] = n_particle_L * N_dim = count_v[world_rank] * N_dim

        index_dim_v : int* - size = {world_size}
            index_dim_v[world_rank] : index of first local particle in global 2D array for (ie: x, v, a)
    */

    // Broadcast Global No. of Particles
    MPI_Bcast(&n_particle_G, 1, MPI_INT, 0, MPI_COMM_WORLD);

    // Allocate Index and Count for Scatterv and Gatherv
    count_1D_v = new int[world_size]();
    index_1D_v = new int[world_size]();

    count_2D_v = new int[world_size]();
    index_2D_v = new int[world_size]();

    // Compute Index and Count for Scatterv and Gatherv
    int index_start, index_end;
    int n_rem = n_particle_G % world_size; // Remainder
    int n_min = (n_particle_G - n_rem) / world_size; // Min No. of Particles

    for (int rank = 0; rank < world_size; ++rank){

        if (rank < n_rem){

            index_start = (n_min + 1) * rank;
            index_end = (n_min + 1) * (rank + 1);

        } else{

            index_start = (n_min + 1) * n_rem + n_min * (rank - n_rem);
            index_end = (n_min + 1) * n_rem + n_min * (rank - n_rem + 1);

        }

        // Scalar Array Count / Index
        count_1D_v[rank] = index_end - index_start;
        index_1D_v[rank] = index_start;

        // Vector Array Count / Index
        count_2D_v[rank] = count_1D_v[rank] * N_dim;
        index_2D_v[rank] = index_1D_v[rank] * N_dim;

    }

    // Assign Local No. of Particle
    n_particle = count_1D_v[world_rank];
    cout << "Rank = " << world_rank << " | No. of Local Particle = " << n_particle << endl;
    cout << "Rank = " << world_rank << " | Global Index = " << index_1D_v[world_rank] << endl;

    sleep(0.5); // Wait for cout buffer latency. Generates uniform outputs by processes.

}

void MPI_Manager::Comm_Array_L2G2L(void* array_G, void* array_L, int nd_array, char mpi_T){

    /*
        Collate local data and transmit to whole comm.
    */

    // For 1D Array
    if (nd_array == 1){

        if (mpi_T == 'D'){
            
            MPI_Allgatherv(array_L, count_1D_v[world_rank], MPI_DOUBLE,
                           array_G, count_1D_v, index_1D_v, MPI_DOUBLE,
                           MPI_COMM_WORLD);

        } else if (mpi_T == 'I'){

            MPI_Allgatherv(array_L, count_1D_v[world_rank], MPI_INT,
                           array_G, count_1D_v, index_1D_v, MPI_INT,
                           MPI_COMM_WORLD);

        }

    // For 2D Array
    } else if (nd_array == 2){

        if (mpi_T == 'D'){
            
            MPI_Allgatherv(array_L, count_2D_v[world_rank], MPI_DOUBLE,
                           array_G, count_2D_v, index_2D_v, MPI_DOUBLE,
                           MPI_COMM_WORLD);

        } else if (mpi_T == 'I'){

            MPI_Allgatherv(array_L, count_2D_v[world_rank], MPI_INT,
                           array_G, count_2D_v, index_2D_v, MPI_INT,
                           MPI_COMM_WORLD);

        }

    }



}

void MPI_Manager::Comm_Array_G2L(void* array_G, void* array_L, int nd_array, char mpi_T){

    /*
        Global to Local Scatterv for arrays.

        Calls the relevant Scatterv based on nd_array and mpi_T.

        array_G : pointer to global array
        array_L : pointer to local array
        nd_array : dimension of array_G and array_L (1 for 1D array, 2 for 2D array)
        mpi_T : type ('D' for MPI_DOUBLE, 'I' for MPI_INT)
    */

    // For 1D Array
    if (nd_array == 1){

        if (mpi_T == 'D'){

            MPI_Scatterv(array_G, count_1D_v, index_1D_v, MPI_DOUBLE,
                         array_L, count_1D_v[world_rank], MPI_DOUBLE,
                         0, MPI_COMM_WORLD);

        } else if (mpi_T == 'I'){

            MPI_Scatterv(array_G, count_1D_v, index_1D_v, MPI_INT,
                         array_L, count_1D_v[world_rank], MPI_INT,
                         0, MPI_COMM_WORLD);

        }

    // For 2D Array
    } else if (nd_array == 2){

        if (mpi_T == 'D'){

            MPI_Scatterv(array_G, count_2D_v, index_2D_v, MPI_DOUBLE,
                         array_L, count_2D_v[world_rank], MPI_DOUBLE,
                         0, MPI_COMM_WORLD);

        } else if (mpi_T == 'I'){

            MPI_Scatterv(array_G, count_2D_v, index_2D_v, MPI_INT,
                         array_L, count_2D_v[world_rank], MPI_INT,
                         0, MPI_COMM_WORLD);

        }

    }

} 

void MPI_Manager::Comm_Array_L2G(void* array_G, void* array_L, int nd_array, char mpi_T){

    /*
        Local to Global Gatherv for arrays.

        Calls the relevant Gatherv based on nd_array and mpi_T.

        array_G : pointer to global array
        array_L : pointer to local array
        nd_array : dimension of array_G and array_L (1 for 1D array, 2 for 2D array)
        mpi_T : type ('D' for MPI_DOUBLE, 'I' for MPI_INT)

    */

    // For 1D Array
    if (nd_array == 1){

        if (mpi_T == 'D'){
            
            MPI_Gatherv(array_L, count_1D_v[world_rank], MPI_DOUBLE,
                        array_G, count_1D_v, index_1D_v, MPI_DOUBLE,
                        0, MPI_COMM_WORLD);

        } else if (mpi_T == 'I'){

            MPI_Gatherv(array_L, count_1D_v[world_rank], MPI_INT,
                        array_G, count_1D_v, index_1D_v, MPI_INT,
                        0, MPI_COMM_WORLD);

        }

    // For 2D Array
    } else if (nd_array == 2){

        if (mpi_T == 'D'){
            
            MPI_Gatherv(array_L, count_2D_v[world_rank], MPI_DOUBLE,
                        array_G, count_2D_v, index_2D_v, MPI_DOUBLE,
                        0, MPI_COMM_WORLD);

        } else if (mpi_T == 'I'){

            MPI_Gatherv(array_L, count_2D_v[world_rank], MPI_INT,
                        array_G, count_2D_v, index_2D_v, MPI_INT,
                        0, MPI_COMM_WORLD);

        }

    }

}

void MPI_Manager::Comm_RhoSum(double &rho_sum){

    /*
        All-to-all communcation of rho_sum.
    */

    MPI_Allreduce(&rho_sum, &rho_sum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

}

void MPI_Manager::Comm_Energy_L2G(double &Ek_total_G, double Ek_total,
                                  double &Ep_total_G, double Ep_total,
                                  double &E_total_G){

    /*
        Reduce Communication of Energy Variables to world_rank = 0
    */

    // Kinetic Energy
    MPI_Reduce(&Ek_total, &Ek_total_G, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

    // Potential Energy
    MPI_Reduce(&Ep_total, &Ep_total_G, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

    // Total Energy
    if (world_rank == 0){

        E_total_G = Ek_total_G + Ep_total_G;

    }
    
}

void MPI_Manager::Comm_Bcast(void* p_data, int count, int root_rank, char mpi_T){

    /*
        Broadcast data.
    */

    if (mpi_T == 'D'){

        MPI_Bcast(p_data, count, MPI_DOUBLE, root_rank, MPI_COMM_WORLD);

    } else if(mpi_T == 'I'){

        MPI_Bcast(p_data, count, MPI_INT, root_rank, MPI_COMM_WORLD);

    }

}

void MPI_Manager::SyncProcess(){

    /*
        Calls MPI_Barrier()
    */

    MPI_Barrier(MPI_COMM_WORLD);

}
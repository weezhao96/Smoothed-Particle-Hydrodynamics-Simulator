/* const_definition.h
    Header file for constant variables definitions.
    The computational loops of SPH in the source code is written in a manner agnostic to the number of spatial dimensions.
    To vary the dimension of the computational domain or the initial box in the various initial condtions, change the
    variables and arrays defined in this header file.

    N_dim: no. of dimensions of the computational box domain

*/

#pragma once

#ifndef const_param
#define const_param

// Machine Epsilon
#include <limits>
const double eps = std::numeric_limits<double>::epsilon();

// pi
const double pi = 3.141592653589;

// No. of Spatial Dimensions
const int N_dim = 2;

/*
    Bounding Box Array Convetion:
    1st Dimension: Spatial dimension, (eg: x,y,z)
    2nd Dimension: Range of box in spatial dimension [n_dim]
*/

// Domain of Computational Box
const double domain_box[N_dim][2] = {{0.0,1.0},{0.0,1.0}};

// Domain of Initial Particle: Dam Break
const double box_dam_break[N_dim][2] = {{0.0,0.2},{0.0,0.2}};

// Domain of Initial Particle: Block Drop
const double box_block_drop[N_dim][2] = {{0.1,0.3},{0.3,0.6}};

/*
    Particle Position:
    1st Dimension: Spatial dimension, (eg: x,y,z)
    2nd Dimension: Postion of particle in spatial dimension [n_dim]
*/

// Initial Particle Positon: 1 Particle Test
const double init_pos_1[N_dim][1] = {{0.5},{0.5}};

// Initial Particle Positon: 4 Particle Test
const double init_pos_4[N_dim][4] = {{0.505,0.515,0.51,0.5},
                                     {0.5,0.5,0.45,0.45}};

// Kissing Number (Max No. of Neigbours of a Sphere)
// kissing_number[N] : kissing number at N-dimensional space
const int kissing_number[5] = {0, 2, 6, 12, 24};

#endif
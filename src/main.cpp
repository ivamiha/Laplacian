// File       : main.cpp    
// Created    : Tue Mar 16 2021  9:46:31 am CET (+0100)
// Author     : Ivan Mihajlovic Milin
// Description: File with main function for testing Laplacian function
// Copyright 2021 ETH Zurich. All Rights Reserved.

#include "Cubism/Block/DataLab.h"
#include "Cubism/Block/Field.h"
#include "Cubism/Common.h"
#include "Cubism/Compiler.h"

#include "Laplacian.h"
#include <iostream>

using namespace Cubism;

int main() 
{    
    // step 1 - pre-processing
    ///////////////////////////////////////////////////////////////////////////

    // welcome user & request input for simulation configurations
    int method;             // user defined sim configuration 
    printf("\n------------------------------------------------------------\n");
    printf("T E S T   L A P L A C I A N   C O M P U T E   K E R N E L\n\n");
    printf("Optimized Laplacian compute kernel utilizing CubismNova\n");
    printf("============================================================\n");
    printf("Selection menu:\n");
    printf("1) run simulation utilizing 2nd-order CDS\n");
    printf("2) run simulation utilizing 4th-order CDS\n");
    printf("3) run code verification of 2nd-order CDS implementation\n");
    printf("4) run code verification of 4th-order CDS implementation\n");
    printf("Please input desired simulation configuration:\n");
    std::cin >> method; 
    // ensure that selected method is indeed an integer in program menu
    method = round(method);
    while (method < 1 || method > 4) {
        printf("Please ensure input number is correct:\n");
        std::cin >> method;
        method = round(method);
    }
    
    // initialize simulation variables to be used throughout
    double time = 5000.0;   // simulation duration [s]
    double domain = 1.0;    // computational domain in each dimension [m]
    int N = 25;             // grid points on domain in each directiom []
    double D = 0.00002;     // diffusivity of fluid [m²/s]
    int print_freq = 10;    // frequency of residual print to terminal []
    double h = domain/N;    // stepsize in all directions (dx,dy,dz) [m]

    // identifers to be used for creating & managing 3D concentration fields
    using Field = Block::Field<double, EntityType::Cell, 3>;
    using IRange = typename Field::IndexRangeType; 
    using MIndex = typename IRange::MultiIndex;
    using DataLab = Block::DataLab<Field>;
    using Stencil = typename DataLab::StencilType; 

    // create 3D scalar fields to be used 
    MIndex elements(N);                 // N*N*N = N³ cells
    IRange element_domain(elements);    // generate element domain
    Field f(element_domain);            // scalar field for storage
    Field f_new(element_domain);        // scalar field for new solution 
    // define function which will return scalar field
    auto field_f = [&](const MIndex &) -> Field & { return f; };

    // generate DataLab object for easy & efficient ghost cell treatment
    DataLab dlab;
    const Stencil s2(-1,2);                 // stencil for 2nd-order CDS
    const Stencil s4(-2,3);                 // stencil for 4th-order CDS
    dlab.allocate(s2, f.getIndexRange());   // allocate memory
    dlab.loadData(MIndex(0), field_f);      // load data w/ periodic BCs
    // define function which will return data lab



    // step 2 - running the simulation/code verification
    ///////////////////////////////////////////////////////////////////////////
    
    // timestep based on 2nd-order CDS stability conditions
    // TODO: incorporate also dt calculation for 4th-order stability
    double dt = h*h / (2*D); 

    // loop through time 
    for (double t = 0.0; t < time; t += dt) {
        // 1. compute Laplacian for given time step & store in field
        Laplacian(dlab, f);
        // 2. advance the results (or perhaps do this in Laplacian?)
        // 3. exchange results between field & dlab for next step

    }



    // step 3 - post-processing
    ///////////////////////////////////////////////////////////////////////////

    return 0; 
}

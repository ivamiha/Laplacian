// File       : Laplacian2f.h
// Created    : Tue Mar 16 2021  9:50:39 am CET (+0100)
// Author     : Ivan Mihajlovic Milin
// Description: 2nd-order flat-indexing Laplacian kernel specifications
// Copyright 2021 ETH Zurich. All Rights Reserved.

#ifndef LAPLACIAN2F_H
#define LAPLACIAN2F_H

// XXX: [fabianw@mavt.ethz.ch; 2021-03-25] You do not need that here.  The
// general rule is not to define `using namespace` in header files because it
// will propagate wherever you include that file.  In local compilation units
// (.cpp files) it can make sense, otherwise always use explicit scoping ->
// Cubism::
using namespace Cubism;

/** 
 * @brief Naive 2nd-order flat-indexing Laplacian kernel    
 * @param sol Current solution stored in data structure DataLab 
 * @param tmp Temporary storage of new solution in data structure Field
 *
 * @rst Naive Laplacian compute kernel utilizing 2nd-order central
 * discretization scheme (CDS). Implemented with CubismNova's flat indexing.  
 * @endrst 
 * */

// XXX: [fabianw@mavt.ethz.ch; 2021-03-25] In the current master branch
// (fa84c47) I have removed the `DataLab` type with a `FieldLab`.  You should
// pass a `Field` in the second argument instead of the lower level `Data` type.
// For example:
//
// template <typename Lab>
// void Laplacian2f(const Lab &sol, typename Lab::FieldType &tmp)
//
template <typename DataLab>
void Laplacian2f(DataLab &sol, 
                 typename DataLab::Data &tmp) 
{
    // XXX: [fabianw@mavt.ethz.ch; 2021-03-25] In principle you only need the
    // `MIndex` type below.  You can get that from the `Lab` directly:
    //
    // using MIndex = typename Lab::Mindex;
    //
    using DataType = typename DataLab::Data::DataType;
    using Mesh = Mesh::StructuredUniform<DataType, 
                                         DataLab::Data::IndexRangeType::Dim>;
    using MIndex = typename Mesh::MultiIndex;           
 
    // apply flat-indexing & 2nd-order CDS
    const MIndex ix{1, 0, 0}; 
    const MIndex iy{0, 1, 0};
    const MIndex iz{0, 0, 1};
    // loop over all elements in passed-in Field block  
    for (auto &i : tmp.getIndexRange()) {
        tmp[i] =   sol[i + ix] + sol[i - ix] + sol[i + iy] + sol[i - iy]
                 + sol[i + iz] + sol[i - iz] - 6*sol[i];         
    }
}

#endif // LAPLACIAN2F_H

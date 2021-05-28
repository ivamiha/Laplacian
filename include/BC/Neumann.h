// File       : Neumann.h
// Created    : Tue May 18 2021 05:46:39 PM (+0200)
// Author     : Ivan Mihajlovic Milin
// Description: Neumann boundary conditions
// Copyright 2021 ETH Zurich. All Rights Reserved.

#ifndef NEUMANN_H
#define NEUMANN_H

#include "Cubism/BC/Base.h"
#include <cassert>

NAMESPACE_BEGIN(Cubism)
NAMESPACE_BEGIN(BC)

/**
 * @brief Neumann BC
 * @tparam Lab Type of ``FieldLab``
 *
 * @rst Constant value Neumann boundary condition
 * @endrst 
 * */
template <typename Lab>
class Neumann : public BC::Base<Lab>
{
    using BaseType = BC::Base<Lab>;
    using BaseType::binfo_;

    using DataType = typename Lab::DataType;

public: 
    /** 
     * @brief Main constructor
     * @param dir Direction in which to apply the boundary
     * @param side On which side along direction ``dir``
     * @param val Value of derivative at boundary
     */
    Neumann(const size_t dir, const size_t side, const DataType &val)
        : BaseType(dir, side), value_(val)
    {
        binfo_.is_periodic = false; 
    }

    /**
     * @brief Apply boundary condition
     * @param lab FieldLab on which the boundary is applied
     */
    void operator()(Lab &lab) override { apply_(lab); }

    /** 
     * @brief Name of boundary condition
     * @return Name string
     */
    std::string name() const override { return std::string("Neumann"); }

    /** 
     * @brief Get bondary value
     * @return Reference to ``DataType``
     */
    DataType &getValue() { return value_; }

    /** 
     * @brief Get boundary value
     * @return ``const`` reference to ``DataType``
     */
    const DataType &getValue() const { return value_; }

private:
    using IndexRangeType = typename Lab::IndexRangeType;
    using MultiIndex = typename Lab::MultiIndex;
    using Index = typename MultiIndex::DataType; 
    using Stencil = typename Lab::StencilType;
    using FieldType = typename Lab::FieldType;

    DataType value_;

    void apply_(Lab &lab) const
    {
        assert(binfo_.dir < IndexRangeType::Dim);
        assert(0 == binfo_.side || 1 == binfo_.side);

        const Stencil &stencil = lab.getActiveStencil();
        if (!this->isValidStencil_(stencil)) {
            return; // nothing to do; zero stencil width for binfo_.dir
        }

        MultiIndex extent;
        if (stencil.isTensorial()) {
            extent = lab.getActiveLabRange().getExtent();
        } else {
            extent = lab.getActiveRange().getExtent();
        }

        const MultiIndex sbegin = stencil.getBegin();
        const MultiIndex send = stencil.getEnd();
        MultiIndex start(0);
        if (0 == binfo_.side) {
            extent[binfo_.dir] = -sbegin[binfo_.dir];
            if (stencil.isTensorial()) {
                start = sbegin;
            } else {
                start[binfo_.dir] = sbegin[binfo_.dir];
            }
        } else {
            extent[binfo_.dir] = send[binfo_.dir] - 1;
            if (stencil.isTensorial()) {
                start = sbegin;
            }
            start[binfo_.dir] = lab.getActiveRange().getExtent()[binfo_.dir];
        }

        // obtain active field and extract associated mesh
        FieldType &field = lab.getActiveField();
        const auto &m = field.getState().mesh;
        // save extents of the field to be used for accessing desired cell
        const size_t Nx = lab.getActiveRange().getExtent()[0];
        const size_t Ny = lab.getActiveRange().getExtent()[1];
        const size_t Nz = lab.getActiveRange().getExtent()[2]; 
        // extract grid size at desired point in desired direction ``dir``
        
        DataType h;
        if (0 == binfo_.side) {
            // first cell in field contains desired grid size
            h = m->getCellSize(0)[binfo_.dir];
        }
        else {
            // access requried cell based on specified ``dir``
            if (0 == binfo_.dir) {
                h = m->getCellSize(Nx - 1)[binfo_.dir];
            }
            else if (1 == binfo_.dir) {
                h = m->getCellSize(Nx * Ny - 1)[binfo_.dir];
            } else {
                h = m->getCellSize(Nx * Ny * Nz - 1)[binfo_.dir];
            }
        }

        // generate array of size ``extent`` and loop through all elements 
        const IndexRangeType slab(extent);
        for (const auto &p : slab) {
            if (0 == binfo_.side) {
                // loop right to left, reference values on right
                lab[start + p] = value_ * h + lab[start + p + 1];
            } else {
                // loop left to right, references values on left
                lab[start + p] = value_ * h + lab[start + p - 1];
            }
        }
    }
};

NAMESPACE_END(BC)
NAMESPACE_END(Cubism)

#endif /* NEUMANN_H */

// File       : NeumannTest.cpp
// Created    : Thu May 27 2021 02:35:02 PM (+0200)
// Author     : Ivan Mihajlovic Milin
// Description: Neumann boundary tests
// Copyright 2021 ETH Zurich. All Rights Reserved.

#include "Cubism/Block/FieldLab.h"
#include "Cubism/Block/Field.h"
#include "Cubism/Grid/Cartesian.h"
#include "Cubism/Mesh/StructuredUniform.h"

#include "Neumann.h"

#include <gtest/gtest.h>
#include <string>

namespace
{
using namespace Cubism;

// template & test for non-tensorial Neumann BCs
template <size_t dir, typename FieldLab, typename IRange, typename MIndex>
void check(const FieldLab &lab,
           const IRange &range,
           MIndex start,
           const MIndex end,
           const typename FieldLab::DataType v0,
           const typename FieldLab::DataType v1)
{
    using Stencil = typename FieldLab::StencilType;
    using DataType = typename FieldLab::DataType;
    using FieldType = typename FieldLab::FieldType;
    using Vector = Core::Vector<long int, 3>;

    // test is simplified by using constant grid size throughout
    const FieldType &field = lab.getActiveField();
    const auto &m = field.getState().mesh;
    const DataType h = m->getCellSize(0)[dir];

    const Vector unitx = {1, 0, 0};
    const Vector unity = {0, 1, 0};
    const Vector unitz = {0, 0, 1};

    const Stencil &stencil = lab.getActiveStencil();
    MIndex Nslab(range.getExtent());
    Nslab[dir] = stencil.getEnd()[dir] - 1;
    const IRange halo_slab(Nslab);

    start[dir] = stencil.getBegin()[dir];   // side = 0
    for (const auto &p : halo_slab) {
        if (0 == dir) {
            EXPECT_NEAR((lab[start + p] - lab[start + p + unitx]) / h,
                                                                    v0, 0.001);
        } else if (1 == dir) {
            EXPECT_NEAR((lab[start + p] - lab[start + p + unity]) / h,
                                                                    v0, 0.001);
        } else {
            EXPECT_NEAR((lab[start + p] - lab[start + p + unitz]) / h,
                                                                    v0, 0.001);
        }
    }

    start[dir] = end[dir];                  // side = 1
    for (const auto &p : halo_slab) {
        if (0 == dir) {
            EXPECT_NEAR((lab[start + p] - lab[start + p - unitx]) / h,
                                                                    v1, 0.001);
        } else if (1 == dir) {
            EXPECT_NEAR((lab[start + p] - lab[start + p - unity]) / h, 
                                                                    v1, 0.001);
        } else {
            EXPECT_NEAR((lab[start + p] - lab[start + p - unitz]) / h,
                                                                    v1, 0.001);
        }
    }
}

TEST(BC, Neumann)
{
    using IRange = Core::IndexRange<3>;
    using MIndex = typename IRange::MultiIndex;
    using Mesh = Mesh::StructuredUniform<double, IRange::Dim>;
    using GridType = Grid::Cartesian<double, Mesh, EntityType::Cell, 0>;
    using DataType = typename GridType::DataType;
    using FieldType = typename GridType::BaseType;
    using FieldLab = Block::FieldLab<typename FieldType::FieldType>;
    using Stencil = typename FieldLab::StencilType;
    using BCVector = typename FieldType::BCVector;
    using BC = BC::Neumann<FieldLab>;

    const MIndex nblocks(4);
    const MIndex block_cells(16);
    GridType field(nblocks, block_cells);

    auto IC = [](FieldType &bf) {
        const Mesh &bm = *bf.getState().mesh;
        for (auto &ci : bm[EntityType::Cell]) {
            bf[ci] = 1.0;
        }
    };

    for (auto bf : field) {
        IC(*bf);
    }

    // value of derivative at boundary that will be applied
    const DataType gradient = 1.0;

    BCVector bcv;
    bcv.push_back(new BC(0, 0, gradient));
    bcv.push_back(new BC(0, 1, gradient));
    bcv.push_back(new BC(1, 0, gradient));
    bcv.push_back(new BC(1, 1, gradient));
    bcv.push_back(new BC(2, 0, gradient));
    bcv.push_back(new BC(2, 1, gradient));

    EXPECT_EQ(bcv[0]->name(), std::string("Neumann"));

    FieldLab lab;
    // testing methodology relies on a symmetrical stencil being used
    const Stencil s(-8, 9, false);
    lab.allocate(s, field[0].getIndexRange());
    auto findex = field.getIndexFunctor();
    lab.loadData(MIndex(0), findex, bcv);
    const auto &range = lab.getActiveRange();

    check<0>(lab, range, MIndex(0), range.getExtent(), gradient, gradient); 
    check<1>(lab, range, MIndex(0), range.getExtent(), gradient, gradient);
    check<2>(lab, range, MIndex(0), range.getExtent(), gradient, gradient);

    for (auto bc : bcv) {
        delete bc;
    }
}
} // namespace

int main(void) {
    ::testing::InitGoogleTest();
    int tests = RUN_ALL_TESTS();

    return tests;
}

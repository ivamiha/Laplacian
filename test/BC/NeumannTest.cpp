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
    MIndex Nslab(range.getExtent());
    Nslab[dir] = 3;
    const IRange halo_slab(Nslab);
    start[dir] = -3;
    for (const auto &p : halo_slab) {
        EXPECT_EQ(lab[p + start], v0);  // side = 0
    }
    start[dir] = end[dir];
    for (const auto &p : halo_slab) {
        EXPECT_EQ(lab[p + start], v1);  // side = 1
    }
}

TEST(BC, Neumann)
{
    using IRange = Core::IndexRange<3>;
    using MIndex = typename IRange::MultiIndex;
    using Mesh = Mesh::StructuredUniform<double, IRange::Dim>;
    using GridType = Grid::Cartesian<double, Mesh, EntityType::Cell, 0>;
    using FieldType = typename GridType::BaseType;
    using FieldLab = Block::FieldLab<typename FieldType::FieldType>;
    using Stencil = typename FieldLab::StencilType;
    using BCVector = typename FieldType::BCVector;
    using BC = BC::Neumann<FieldLab>;

    const MIndex nblocks(1);
    const MIndex block_cells(1);
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

    BCVector bcv;
    bcv.push_back(new BC(0, 0, 0));
    bcv.push_back(new BC(0, 1, 0));
    bcv.push_back(new BC(1, 0, 0));
    bcv.push_back(new BC(1, 1, 0));
    bcv.push_back(new BC(2, 0, 0));
    bcv.push_back(new BC(2, 1, 0));

    EXPECT_EQ(bcv[0]->name(), std::string("Neumann"));

    FieldLab lab;
    const Stencil s(-3, 4, false);
    lab.allocate(s, field[0].getIndexRange());
    auto findex = field.getIndexFunctor();
    lab.loadData(MIndex(0), findex, bcv);
    const auto &range = lab.getActiveRange();

    check<0>(lab, range, MIndex(0), range.getExtent(), 1, 1);
    check<1>(lab, range, MIndex(0), range.getExtent(), 1, 1);
    check<2>(lab, range, MIndex(0), range.getExtent(), 1, 1);

    for (auto bc : bcv) {
        delete bc;
    }

// template and test for tensorial Neumann BCs
}
} // namespace

int main(void) {
    ::testing::InitGoogleTest();
    int tests = RUN_ALL_TESTS();

    return tests;
}

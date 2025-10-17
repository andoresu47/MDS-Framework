#pragma once
#include "max_disjoint_framework.hpp"
#include <boost/graph/adjacency_list.hpp>
#include <vector>
#include <algorithm>

// --- Lattice graph with in-place “active” flags on vertices/edges ---
struct ToyVertex { bool active = true; };
struct ToyEdge   { bool active = true; }; // (unused here, but kept for symmetry)

using ToyLattice = boost::adjacency_list<
    boost::vecS, boost::vecS, boost::directedS,
    ToyVertex, ToyEdge
>;

struct ToySolution {
    std::vector<std::size_t> elems; // a subset of {0,1,2}
    bool operator==(const ToySolution& o) const { return elems == o.elems; }
};

class ToyPowersetMDS : public MaxDisjointSolutionsFramework<ToyLattice, ToySolution> {
public:
    // Build a 3-vertex “lattice” (no edges needed)
    ToyLattice build_initial_lattice() override;

    // Oracles over the current (mutated) lattice
    ToySolution O_min() override;          // smallest active vertex
    ToySolution O_max() override;          // largest  active vertex
    void        O_ds (const ToySolution& s) override; // deactivate chosen elements

    bool        are_disjoint(const ToySolution& a, const ToySolution& b) const override;

    // Define “empty” as: no active vertices left
    bool        is_empty() const override;

    // Convert internal ideals to external cut solutions
    std::vector<ToySolution>
    convert_to_solution_space(const std::vector<ToySolution>& C) const override;
};

#pragma once
#include "max_disjoint_framework.hpp"
#include <boost/graph/adjacency_list.hpp>
#include <vector>

// ---------------- Toy lattice over a 3-element “powerset” ----------------
using ToyLattice =
    boost::adjacency_list<
        boost::vecS, boost::vecS, boost::directedS,
        boost::property<boost::vertex_index_t, std::size_t>,
        boost::property<boost::edge_index_t,   std::size_t>
    >;

struct ToySolution {
    std::vector<std::size_t> elems; // subset of {0,1,2}
    bool operator==(const ToySolution& o) const { return elems == o.elems; }
};

// --------------------------------------------------------------------------
// ToyPowersetMDS: header-only declarations; definitions in .cpp
// --------------------------------------------------------------------------
class ToyPowersetMDS
: public MaxDisjointSolutionsFramework<ToyLattice, ToySolution> {
public:
    // Oracles / overrides (definitions in the .cpp)
    ToyLattice  build_initial_lattice() override;
    ToySolution O_min(const LatticeView& Lstar) override;
    ToySolution O_max(const LatticeView& Lstar) override;
    LatticeView O_ds (const LatticeView& Lstar, const ToySolution& S) override;
    bool        are_disjoint(const ToySolution& A, const ToySolution& B) const override;

    // Optional: treat view as empty if no visible vertices remain.
    std::size_t view_size(const LatticeView& Lstar) const override;
    bool        is_empty (const LatticeView& Lstar) const override;
};

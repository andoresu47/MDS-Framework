#pragma once
#include "max_disjoint_framework.hpp"

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graph_traits.hpp>
#include <vector>
#include <cstddef>

// ============================ Instance (input graph) ============================
struct MC_EdgeProps {
    long        capacity = 0;   // input capacity c_uv
    long        residual = 0;   // filled after max-flow: forward residual
    std::size_t id       = 0;   // stable original edge id [0..m)
};

using MC_G = boost::adjacency_list<
    boost::vecS, boost::vecS, boost::bidirectionalS,
    boost::no_property,
    MC_EdgeProps
>;
using Vg = boost::graph_traits<MC_G>::vertex_descriptor;
using Eg = boost::graph_traits<MC_G>::edge_descriptor;

// ============================ Compact representation of lattice ====================
// Each vertex = SCC in reverse-residual G'
// Each edge = reachability between SCCs
// Edge property includes whether it corresponds to a saturated edge in G.
struct LEdge {
    bool saturated = false;
};

using LatticeDAG = boost::adjacency_list<
    boost::vecS, boost::vecS, boost::directedS,
    boost::no_property,
    LEdge
>;

using Vh = boost::graph_traits<LatticeDAG>::vertex_descriptor;
using Eh = boost::graph_traits<LatticeDAG>::edge_descriptor;

// ============================ Solution types ===================================
// External solution = set of edges defining a minimum st cut
struct MinCutSolution {
    std::vector<std::size_t> edge_ids; // original edge ids forming δ(I)
    bool operator==(const MinCutSolution& o) const { return edge_ids == o.edge_ids; }
};

// Internal solution = Ideal over Lattice DAG (membership bitvector over |V(lattice_)|)
using Ideal = std::vector<char>;

// ============================ MinCutMaxDisjoint =======================
// Implements the three oracles directly on the SCC DAG (the lattice).
// The SCC DAG itself is used as the compact representation of the lattice of
// minimum s–t cuts; O_min / O_max / O_ds operate directly on it.
class MinCutMaxDisjoint
: public MaxDisjointSolutionsFramework<LatticeDAG, Ideal, MinCutSolution> {
public:
    explicit MinCutMaxDisjoint(const MC_G& g, Vg s, Vg t);
    explicit MinCutMaxDisjoint(MC_G&& g, Vg s, Vg t);

    // Framework hooks
    LatticeDAG  build_initial_lattice() override;       // builds compact representation of lattice
    Ideal       O_min() override;                       // minimal ideal of current sublattice
    Ideal       O_max() override;                       // largest ideal of the lattice
    void        O_ds (const Ideal& X) override;         // advance cutoff_ using saturated out-edges of X
    bool        are_disjoint(const Ideal& I1, const Ideal& I2) const override;    // compare saturated edges
    bool        is_empty() const override;              // true if no vertex remains past cutoff_

    // Convert internal ideals to external cut solutions
    std::vector<MinCutSolution>
    convert_to_solution_space(const std::vector<Ideal>& C) const override;

    // Optional accessors for testing
    const MC_G& instance() const { return G_; }
    const std::vector<int>& scc_partition() const { return scc_of_; }
    const std::vector<Vh>& topo_order() const { return topo_TS_; }

private:
    // Helpers
    void prune_to_st_core();                // O(n+m) s–t-core pruning on input graph G_
    void run_maxflow_and_fill_residual();   // Compute max-flow and fill G_[e].residual
    MinCutSolution cut_from_H_ideal(const std::vector<char>& in_I_H) const; // Cut extraction helper
    bool build_current_minimal_ideal_H(std::vector<char>& in_I_H) const;    // Ideal builder over the current sublattice
    
private:
    // Input graph (pruned)
    MC_G G_;
    Vg   s_{}, t_{};

    // Cached data for lattice & oracles
    std::vector<int> scc_of_;       // vertex -> representation id
    int scc_S_ = -1;                // SCC id of source in representation
    int scc_T_ = -1;                // SCC id of sink in representation

    // For each original edge id: create (SCC(u), SCC(v)) and mark whether it is saturated
    std::vector<std::pair<int,int>> edge_scc_pair_;
    std::vector<char> edge_is_saturated_;

    // Oracle runtime state
    std::vector<Vh> topo_TS_;   // topological order of DAG representation (lattice_)
    std::vector<int> topo_pos_; // vertex -> position in topo_TS_
    int cutoff_ = -1;           // invalid prefix endpoint in topo_TS_
};

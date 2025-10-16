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

using LatticeG = boost::adjacency_list<
    boost::vecS, boost::vecS, boost::directedS,
    boost::no_property,
    LEdge
>;

using Vh = boost::graph_traits<LatticeG>::vertex_descriptor;
using Eh = boost::graph_traits<LatticeG>::edge_descriptor;

// ============================ Solution ==========================================
struct MinCutSolution {
    std::vector<std::size_t> edge_ids; // original edge ids forming δ(I)
    bool operator==(const MinCutSolution& o) const { return edge_ids == o.edge_ids; }
};

// ============================ MinCutMaxDisjoint =======================
// Implements the three oracles directly on the SCC DAG (the lattice).
// The SCC DAG itself is used as the compact representation of the lattice of
// minimum s–t cuts; O_min / O_max / O_ds operate directly on it.
class MinCutMaxDisjoint
: public MaxDisjointSolutionsFramework<LatticeG, MinCutSolution> {
public:
    explicit MinCutMaxDisjoint(const MC_G& g, Vg s, Vg t);
    explicit MinCutMaxDisjoint(MC_G&& g, Vg s, Vg t);

    // Framework hooks
    LatticeG       build_initial_lattice() override;                        // builds compact representation of lattice
    MinCutSolution O_min(LatticeG& g) override;                             // smallest ideal
    MinCutSolution O_max(LatticeG& g) override;                             // largest ideal
    void           O_ds (LatticeG& g, const MinCutSolution& S) override;    // disjoint successors
    bool           are_disjoint(const MinCutSolution& A,
                                const MinCutSolution& B) const override;
    bool           is_empty(const LatticeG& g) const override;

     // Accessors
    const MC_G& instance() const { return G_; }
    Vg          source()   const { return s_; }
    Vg          sink()     const { return t_; }
    const std::vector<int>& scc_partition() const { return scc_of_; }
    const std::vector<Vh>& topo_order() const { return topo_TS_; }

private:
    // Max-flow helper
    void run_maxflow_and_fill_residual();

    // Cut extraction helper
    MinCutSolution cut_from_H_ideal(const std::vector<char>& in_I_H) const;

    // Ideal builder over the current lattice
    bool build_current_minimal_ideal_H(std::vector<char>& in_I_H) const;

    // Prune input graph to the s–t core (keep vertices reachable from s and that can reach t)
    void prune_to_st_core();

private:
    // ====== Stored instance ======
    MC_G G_;
    Vg   s_{}, t_{};

    // ====== SCC / lattice caches built in build_initial_lattice() ======
    // Original vertex -> SCC id in G'
    std::vector<int> scc_of_;       // vertex -> representation id
    int scc_S_ = -1;                // SCC id of source
    int scc_T_ = -1;                // SCC id of sink

    // For each original edge: map (SCC(u), SCC(v)) and whether it is saturated
    std::vector<std::pair<int,int>> edge_scc_pair_;
    std::vector<char> edge_is_saturated_;

    // ====== Oracle runtime state ======
    std::vector<Vh> topo_TS_;   // topological order of SCC DAG
    std::vector<int> topo_pos_; // vertex -> position in topo_TS_
    int cutoff_ = -1;   
};

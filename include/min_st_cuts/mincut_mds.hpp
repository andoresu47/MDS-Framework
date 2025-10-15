#pragma once
#include "max_disjoint_framework.hpp"

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graph_traits.hpp>
#include <vector>
#include <cstddef>

// ============================ Instance (input graph) ============================
struct MC_EdgeProps {           // Edge properties 
    long        capacity = 0;   // input
    long        residual = 0;   // filled after max-flow
    std::size_t id       = 0;   // stable id in [0..m)
};

using MC_G = boost::adjacency_list<
    boost::vecS, boost::vecS, boost::directedS,
    boost::no_property,
    MC_EdgeProps
>;

using Vg = boost::graph_traits<MC_G>::vertex_descriptor;
using Eg = boost::graph_traits<MC_G>::edge_descriptor;

// ============================ Lattice (compact rep) ================================
// SCC-DAG over the “tight” residual subgraph, with one lattice edge per saturated
// original edge. We require vertex_index_t and edge_index_t to drive masks.
using LatticeG = boost::adjacency_list<
    boost::vecS, boost::vecS, boost::directedS,
    boost::property<boost::vertex_index_t, std::size_t>,
    boost::property<boost::edge_index_t,   std::size_t>
>;
using Vh = boost::graph_traits<LatticeG>::vertex_descriptor;
using Eh = boost::graph_traits<LatticeG>::edge_descriptor;

// ============================ Solution ============================================
// One min s–t cut, represented as the set of original edge ids on its boundary δ(I).
struct MinCutSolution {
    std::vector<std::size_t> edge_ids;
    bool operator==(const MinCutSolution& o) const { return edge_ids == o.edge_ids; }
};

// ============================ MinCutMaxDisjoint ===================================
// Implements the three oracles on top of the framework:
//   - O_min(L*):  bottom cut within current sublattice
//   - O_max(L*):  top cut within current sublattice
//   - O_ds(L*,S): sublattice of disjoint successors of S (ban S.edge_ids)
// Disjointness = edge-disjointness on original edge ids.
class MinCutMaxDisjoint
: public MaxDisjointSolutionsFramework<LatticeG, MinCutSolution> {
public:
    // Construct from an instance (copy or move graph) and terminals s,t.
    explicit MinCutMaxDisjoint(const MC_G& g, Vg s, Vg t);
    explicit MinCutMaxDisjoint(MC_G&& g, Vg s, Vg t);

    // ---------- Oracles (override) ----------
    LatticeG        build_initial_lattice() override;
    MinCutSolution  O_min(const LatticeView& Lstar) override;
    MinCutSolution  O_max(const LatticeView& Lstar) override;
    LatticeView     O_ds (const LatticeView& Lstar, const MinCutSolution& S) override;
    bool            are_disjoint(const MinCutSolution& A,
                                 const MinCutSolution& B) const override;

    // Optional: let the framework query “size/empty” of a view (kept simple).
    std::size_t     view_size(const LatticeView& Lstar) const override;
    bool            is_empty (const LatticeView& Lstar) const override;

    // ---------- Accessors (optional) ----------
    const MC_G& instance() const { return G_; }
    Vg          source()   const { return s_; }
    Vg          sink()     const { return t_;  }

private:
    // ====== Build-time helpers (to be defined in .cpp) ======
    // Runs max-flow on a working graph and fills MC_EdgeProps::residual for forward arcs.
    void run_maxflow_and_fill_residual();

    // Build minimal/maximal ideals in the current (filtered) SCC-DAG view.
    // in_I is sized to num_vertices(*Lstar.g) and set to 0/1 membership.
    void build_minimal_ideal(const LatticeView& Lstar,
                             std::vector<char>& in_I) const;

    void build_maximal_ideal(const LatticeView& Lstar,
                             std::vector<char>& in_I) const;

    // Extract cut δ(I) from the current filtered lattice view and membership vector.
    MinCutSolution cut_from_ideal(const LatticeView& Lstar,
                                  const std::vector<char>& in_I) const;

private:
    // ====== Stored instance ======
    MC_G G_;
    Vg s_{}, t_{};

    // ====== Cached data for the lattice build ======
    // scc_of_: maps original vertices (in tight residual) -> SCC id in [0..k)
    std::vector<int> scc_of_;
    int scc_S_ = -1;   // SCC id of s
    int scc_T_ = -1;   // SCC id of t

    // Topological order of the SCC-DAG (lattice graph) for oracle sweeps.
    std::vector<Vh> topo_;

    // Map SCC id (in G') -> reduced-lattice vertex id [0..n_keep-1] or -1 if pruned
    std::vector<int> scc_to_lattice_;             

    // For each original edge id (G_[e].id), store (scc(u), scc(v)) in G'
    std::vector<std::pair<int,int>> edge_scc_pair_;

    // For each original edge id (G_[e].id), true IFF the forward arc (u->v) in G' exists
    // AND corresponds to a saturated original edge (i.e., f_ij == c_ij)
    std::vector<char> edge_is_saturated_;

    // For the reduced lattice vertex ids [0..n_keep-1], mark which belong to Succ(S) in H.
    std::vector<char> reachS_lattice_;
};

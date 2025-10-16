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
    boost::vecS, boost::vecS, boost::directedS,
    boost::no_property,
    MC_EdgeProps
>;
using Vg = boost::graph_traits<MC_G>::vertex_descriptor;
using Eg = boost::graph_traits<MC_G>::edge_descriptor;

// ============================ Lattice (SCC DAG after pruning) ====================
// We don't strictly need per-vertex/edge flags, but you can keep them if you
// later want to deactivate lattice vertices/edges too.
struct LVertex { bool active = true; };
struct LEdge   { bool active = true; };

using LatticeG = boost::adjacency_list<
    boost::vecS, boost::vecS, boost::directedS,
    LVertex,
    LEdge
>;
using Vh = boost::graph_traits<LatticeG>::vertex_descriptor;
using Eh = boost::graph_traits<LatticeG>::edge_descriptor;

// ============================ Solution ==========================================
struct MinCutSolution {
    std::vector<std::size_t> edge_ids; // original edge ids forming δ(I)
    bool operator==(const MinCutSolution& o) const { return edge_ids == o.edge_ids; }
};

// ============================ MinCutMaxDisjoint (in-place) =======================
class MinCutMaxDisjoint
: public MaxDisjointSolutionsFramework<LatticeG, MinCutSolution> {
public:
    explicit MinCutMaxDisjoint(const MC_G& g, Vg s, Vg t);
    explicit MinCutMaxDisjoint(MC_G&& g,       Vg s, Vg t);

    // Framework hooks
    LatticeG       build_initial_lattice() override;
    MinCutSolution O_min(LatticeG& g) override;
    MinCutSolution O_max(LatticeG& g) override;
    void           O_ds (LatticeG& g, const MinCutSolution& S) override;
    bool           are_disjoint(const MinCutSolution& A,
                                const MinCutSolution& B) const override;
    bool           is_empty(const LatticeG& g) const override;

    // Accessors
    const MC_G& instance() const { return G_; }
    Vg          source()   const { return s_; }
    Vg          sink()     const { return t_; }

    // Map each original vertex (0..n-1) to its reduced-lattice vertex id [0..|L|-1],
    // or -1 if that vertex’s SCC was pruned during step 4.
    std::vector<int> original_to_reduced_map() const;

private:
    // Max-flow helper
    void run_maxflow_and_fill_residual();

    // Ideal builders over the current lattice
    void build_minimal_ideal (LatticeG& g, std::vector<char>& in_I) const;
    void build_maximal_ideal (LatticeG& g, std::vector<char>& in_I) const;

    // Boundary extractor that honors edge_active_ and Succ(S)
    MinCutSolution cut_from_ideal_in_place(const std::vector<char>& in_I) const;

private:
    // ====== Stored instance ======
    MC_G G_;
    Vg   s_{}, t_{};

    // ====== SCC / lattice caches built in build_initial_lattice() ======
    // Original vertex -> SCC id in G'
    std::vector<int> scc_of_;
    int scc_S_ = -1, scc_T_ = -1;

    // SCC id in G' -> reduced lattice vertex id [0..n_keep) or -1 if pruned
    std::vector<int> scc_to_lattice_;

    // For reduced lattice vertex ids, mark those in Pred(S) in the SCC DAG
    std::vector<char> predS_scc_;

    // For each original edge id: (scc(u), scc(v)) in G'
    std::vector<std::pair<int,int>> edge_scc_pair_;

    // For each original edge id: is the forward arc saturated (f == c)?
    std::vector<char> edge_is_saturated_;

    // For each original edge id: still allowed (1) or banned (0) by O_ds
    std::vector<char> edge_active_;

    // Optional: cached topo order of the reduced lattice
    std::vector<Vh> topo_;
};

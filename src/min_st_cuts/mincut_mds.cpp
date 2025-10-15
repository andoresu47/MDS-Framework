#include "min_st_cuts/mincut_mds.hpp"

#include <boost/range/iterator_range.hpp>
#include <boost/graph/strong_components.hpp>
#include <boost/graph/topological_sort.hpp>
#include <boost/graph/push_relabel_max_flow.hpp>
#include <algorithm>
#include <queue>

// Convenience alias to the base’s view type
using View = MaxDisjointSolutionsFramework<LatticeG, MinCutSolution>::LatticeView;

// ============================ constructor =========================================

MinCutMaxDisjoint::MinCutMaxDisjoint(const MC_G& g, Vg s, Vg t)
    : G_(g), s_(s), t_(t) {} // copy version

MinCutMaxDisjoint::MinCutMaxDisjoint(MC_G&& g, Vg s, Vg t)
    : G_(std::move(g)), s_(s), t_(t) {} // move version

// ============================ oracles ======================================

LatticeG MinCutMaxDisjoint::build_initial_lattice() {
    // TODO(1): assign stable ids to original edges in G_ (0..m-1) into MC_EdgeProps::id
    // for (auto e : boost::make_iterator_range(edges(G_))) { G_[e].id = ...; }

    // TODO(2): run max flow once and fill MC_EdgeProps::residual on FORWARD arcs of G_
    run_maxflow_and_fill_residual();

    // TODO(3): build "tight residual" subgraph Z (keep saturated forward arcs where residual == 0)
    // MC_G Z(num_vertices(G_));
    // for (auto e : boost::make_iterator_range(edges(G_))) { if (G_[e].capacity>0 && G_[e].residual==0) add_edge(source(e,G_), target(e,G_), Z); }

    // TODO(4): compute strongly connected components of Z into scc_of_ (size = |V(G_)|)
    // scc_of_.assign(num_vertices(Z), -1);
    // int k = boost::strong_components(Z, make_iterator_property_map(scc_of_.begin(), get(boost::vertex_index, Z)));

    // TODO(5): record SCC ids of s and t
    // scc_S_ = scc_of_[s_];
    // scc_T_ = scc_of_[t_];

    // TODO(6): build lattice graph L with k vertices; set vertex_index = i for each i in [0..k-1]
    // LatticeG L(k);
    // for (int i=0;i<k;++i) put(boost::vertex_index, L, static_cast<Vh>(i), static_cast<std::size_t>(i));

    // TODO(7): for each saturated original edge u->v with scc(u)!=scc(v), add lattice edge (scc(u)->scc(v))
    // and set its edge_index to the original edge id (G_[e].id)
    // for (auto e : boost::make_iterator_range(edges(G_))) { ... Eh le = add_edge(cu,cv,L).first; put(boost::edge_index, L, le, G_[e].id); }

    // TODO(8): cache a topological order of L into topo_ (forward topo)
    // topo_.clear();
    // { std::vector<Vh> rev; boost::topological_sort(L, std::back_inserter(rev)); topo_.assign(rev.rbegin(), rev.rend()); }

    // TODO: return L
    return LatticeG{}; // placeholder so it compiles; replace with 'return L;'
}

MinCutSolution MinCutMaxDisjoint::O_min(const View& Lstar) {
    // TODO: compute minimal ideal in current filtered SCC-DAG view:
    // Strategy: reachability closure from scc_S_ in Lstar.filtered(), while never including scc_T_.
    // Allocate in_I = vector<char>(|V(L)|,0); mark reachable nodes (except T) as 1; force in_I[T]=0.
    std::vector<char> in_I(num_vertices(*Lstar.g), 0);
    build_minimal_ideal(Lstar, in_I);
    return cut_from_ideal(Lstar, in_I);
}

MinCutSolution MinCutMaxDisjoint::O_max(const View& Lstar) {
    // TODO: compute maximal ideal in current filtered SCC-DAG view:
    // Strategy: nodes that cannot reach T (in the visible subgraph) can be in the ideal.
    // Or equivalently, start with all=1, then drop any node that can (reach) T, preserving down-closure.
    std::vector<char> in_I(num_vertices(*Lstar.g), 1);
    build_maximal_ideal(Lstar, in_I);
    return cut_from_ideal(Lstar, in_I);
}

View MinCutMaxDisjoint::O_ds(const View& Lstar, const MinCutSolution& S) {
    // Successors inside current sublattice: ban the boundary edges used by S.
    auto next = Lstar.deep_copy(); // deep copy masks to modify independently
    // TODO: for each id in S.edge_ids, set (*next.emask)[id] = 0 (guard bounds)
    // if (next.emask) for (auto id : S.edge_ids) if (id < next.emask->size()) (*next.emask)[id] = 0;
    return next;
}

bool MinCutMaxDisjoint::are_disjoint(const MinCutSolution& A, const MinCutSolution& B) const {
    // TODO: return true iff edge_ids(A) ∩ edge_ids(B) = ∅
    // Tip: copy and sort then set_intersection into a scratch vector.
    // std::vector<std::size_t> a=A.edge_ids, b=B.edge_ids, inter;
    // std::sort(a.begin(),a.end()); std::sort(b.begin(),b.end());
    // std::set_intersection(a.begin(),a.end(), b.begin(),b.end(), std::back_inserter(inter));
    // return inter.empty();
    return true; // placeholder
}

// ============================ framework hooks ==============================

std::size_t MinCutMaxDisjoint::view_size(const View& Lstar) const {
    // Keep the default simple notion: number of visible vertices
    return Lstar.count_vertices();
}

bool MinCutMaxDisjoint::is_empty(const View& Lstar) const {
    return view_size(Lstar) == 0;
}

// ============================ helpers ======================================

void MinCutMaxDisjoint::run_maxflow_and_fill_residual() {
    // TODO:
    // Build a working residual graph WR with explicit reverse edges.
    // Keep a mapping from WR forward edges back to original G_ edges (Eg).
    // Provide to push_relabel_max_flow:
    //   - capacity map (on WR edges)
    //   - residual capacity map (on WR edges)
    //   - reverse-edge map (edge -> its mate)
    //   - vertex index map (get(vertex_index, WR))
    //
    // After max-flow, copy residual capacities on WR forward edges back into G_[ge].residual.
    //
    // See the earlier fixed version we discussed for a complete outline.
}

void MinCutMaxDisjoint::build_minimal_ideal(const View& Lstar,
                                            std::vector<char>& in_I) const {
    // Preconditions: in_I sized to |V(L)|; initially all 0
    // TODO:
    // - auto FG = Lstar.filtered();
    // - BFS/DFS from scc_S_ over FG (guard: skip if scc_S_ is masked out); do not mark scc_T_
    // - mark reachable nodes (except T) in in_I[v] = 1
    // - force in_I[scc_T_] = 0; (and in_I[scc_S_] = 1 if visible)
    //
    // Note: if scc_S_ or scc_T_ are outside the current view (masked), decide policy
    // (e.g., treat as absent; typically S must remain included and T excluded in your model).
}

void MinCutMaxDisjoint::build_maximal_ideal(const View& Lstar,
                                            std::vector<char>& in_I) const {
    // Preconditions: in_I sized to |V(L)|; typically initialized to 1
    // TODO:
    // Option A (reachability-to-T):
    //   - Compute the set R of nodes that can reach scc_T_ in FG (reverse BFS over in_edges).
    //   - Set in_I[v] = 0 for v in R; keep others 1.
    //   - Ensure down-closure: in reverse topological order, if any predecessor is 0, set v=0.
    //   - Force in_I[scc_S_] = 1 and in_I[scc_T_] = 0 (if visible).
    //
    // Option B (co-reachability from S):
    //   - Start from all 1, and clear any node that violates ideal property given topo_.
}

MinCutSolution MinCutMaxDisjoint::cut_from_ideal(const View& Lstar,
                                                 const std::vector<char>& in_I) const {
    MinCutSolution S;
    // TODO:
    // - auto FG = Lstar.filtered();
    // - For each visible lattice edge e: u->v
    //     if (in_I[u] == 1 && in_I[v] == 0) then this edge crosses δ(I)
    //       -> push get(edge_index, *Lstar.g, e) into S.edge_ids
    // - sort + unique S.edge_ids
    return S;
}

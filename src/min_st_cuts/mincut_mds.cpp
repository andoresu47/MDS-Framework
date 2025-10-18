/**
 * @file mincut_mds.cpp
 * @brief Implementation of MinCutMaxDisjoint: builds the lattice of minimum s–t cuts
 *        from a max-flow instance, and implements the oracles used by the
 *        MaxDisjointSolutionsFramework.
 *
 * High-level pipeline in build_initial_lattice():
 *   1) prune_to_st_core(): keep only vertices/edges on some s→t path.
 *   2) Number original edges with stable ids; init per-edge caches.
 *   3) run_maxflow_and_fill_residual(): compute max-flow and write forward residuals.
 *   4) Build "reverse residual" graph G' according to the specific rule:
 *        - if 0 < f < c: add (v→u) with c−f and (u→v) with f
 *        - if f == c:    add (u→v) with f
 *        - if f == 0:    add (v→u) with c
 *   5) Compute SCCs of G', and contract them into an SCC-DAG H (the lattice representation).
 *   6) Mark H-edges that correspond to saturated original edges; fill caches that
 *      map each original edge to the pair (SCC(u), SCC(v)) and whether it is saturated.
 *   7) Compute a topological order of H and cache position indices.
 *   8) Initialize all nodes as valid, except for scc of source node.
 * 
 *  Note: build_initial_lattice() might be moved to its own independent implementation. 
 *
 * Oracles:
 *   - O_min():    returns the current minimal ideal (as a bitvector over SCCs).
 *   - O_max():    returns the maximal ideal in H \ {SCC containing t}.
 *   - O_ds(X):    invalidate endpoints (and their predecessors) of saturated out-edges of X 
 *                 to ensure future ideals are disjoint from X’s induced cut.
 * 
 * Helpers:
 *   - are_disjoint(): checks edge-disjointness by converting ideals to cuts and intersecting.
 *   - is_empty(): true iff the target node is marked invalid.
 *
 * Conversion:
 *   - convert_to_solution_space(): maps a list of ideals to the corresponding min-cuts
 *     (edge-id sets) using cut_from_H_ideal().
 */

#include "min_st_cuts/mincut_mds.hpp"
#include <iostream>
#include <boost/range/iterator_range.hpp>
#include <boost/graph/strong_components.hpp>
#include <boost/graph/topological_sort.hpp>
#include <boost/graph/push_relabel_max_flow.hpp>
#include <algorithm>
#include <queue>

/**
 * @brief Copy-construct the solver from an input graph and terminals.
 * @param g Input graph (copied).
 * @param s Source vertex (descriptor within g).
 * @param t Sink vertex (descriptor within g).
 */
MinCutMaxDisjoint::MinCutMaxDisjoint(const MC_G& g, Vg s, Vg t)
    : G_(g), s_(s), t_(t) {}

/**
 * @brief Move-construct the solver from an input graph and terminals.
 * @param g Input graph (moved).
 * @param s Source vertex (descriptor within g).
 * @param t Sink vertex (descriptor within g).
 */
MinCutMaxDisjoint::MinCutMaxDisjoint(MC_G&& g, Vg s, Vg t)
    : G_(std::move(g)), s_(s), t_(t) {}

/**
 * @brief Build the compact lattice representation (SCC-DAG H) inside the solver.
 *
 * Steps:
 *  1) Prune G_ to the s–t core.
 *  2) Number edges; initialize caches.
 *  3) Run max-flow and fill forward residuals G_[e].residual.
 *  4) Construct reverse-residual graph G' using the specified rules.
 *  5) Compute SCCs of G' and their condensation DAG H.
 *  6) Mark which H-edges correspond to saturated original edges and fill per-edge caches.
 *  7) Compute a topological order topo_TS_ and the index map topo_pos_.
 *  8) Initialize all nodes as valid, except for scc of source node.
 *
 * @return The built SCC-DAG H (stored in the framework’s lattice_ by the parent class' algorithm).
 */
LatticeDAG MinCutMaxDisjoint::build_initial_lattice() {
    // Keep only vertices/edges on some s-t path.
    prune_to_st_core();

    // Number original edges by id and init per-edge arrays
    std::size_t m = 0;
    for (auto e : boost::make_iterator_range(edges(G_))) G_[e].id = m++;
    edge_scc_pair_.assign(m, {-1,-1});
    edge_is_saturated_.assign(m, 0);

    // Max-flow (fills G_[e].residual on forward arcs)
    run_maxflow_and_fill_residual();

    // Build the "reverse residual" graph G'
    using Gprime = boost::adjacency_list<
        boost::vecS, boost::vecS, boost::directedS,
        boost::no_property, 
        boost::property<boost::edge_capacity_t, long>
    >;
    Gprime Gp(num_vertices(G_));
    auto add_arc = [&](Vg u, Vg v, long rc) {
        if (rc <= 0) return;
        auto ee = add_edge(u, v, Gp);
        if (ee.second) put(boost::edge_capacity, Gp, ee.first, rc);
    };

    for (auto e : boost::make_iterator_range(edges(G_))) {
        auto u = boost::source(e, G_), v = boost::target(e, G_);
        long c  = G_[e].capacity;
        long rf = G_[e].residual;      // forward residual
        long f  = c - rf;              // actual flow

        if (0 < f && f < c) { add_arc(v,u, c-f); add_arc(u,v, f); }
        else if (f == c)     { add_arc(u,v, f); }
        else /* f==0 */      { add_arc(v,u, c); }
    }

    // SCCs on G'
    scc_of_.assign(num_vertices(Gp), -1);
    int k = boost::strong_components(
        Gp, boost::make_iterator_property_map(
                scc_of_.begin(), get(boost::vertex_index, Gp))
    );

    // Contract SCC to build compact representation H
    LatticeDAG H(k);
    for (auto e : boost::make_iterator_range(edges(Gp))) {
        int cu = scc_of_[boost::source(e, Gp)];
        int cv = scc_of_[boost::target(e, Gp)];
        if (cu != cv) add_edge(static_cast<Vh>(cu), static_cast<Vh>(cv), H);
    }

    // Mark in H edges that correspond to saturated forward original edges;
    // also fill per-edge caches mapping original edge id → (SCC(u), SCC(v)) plus saturation bit. 
    for (auto ge : boost::make_iterator_range(edges(G_))) {
        const std::size_t id = G_[ge].id;
        const int cu = scc_of_[boost::source(ge, G_)];
        const int cv = scc_of_[boost::target(ge, G_)];
        const long c  = G_[ge].capacity;
        const long rf = G_[ge].residual;
        const long f  = c - rf;
        const bool is_sat_forward = (f == c);   // saturated forward (u->v)

        // Remember this per-original-edge:
        edge_scc_pair_[id]     = {cu, cv};
        edge_is_saturated_[id] = is_sat_forward ? 1 : 0;

        if (cu != cv && is_sat_forward) {
            auto [eh, ok] = add_edge(static_cast<Vh>(cu), static_cast<Vh>(cv), H);
            // Even if the edge already exists, we just set saturated=true.
            H[eh].saturated = true;
        }
    }

    // Topological order + position map over H
    topo_TS_.clear();
    {
        std::vector<Vh> rev;
        boost::topological_sort(H, std::back_inserter(rev));
        topo_TS_.assign(rev.rbegin(), rev.rend());
    }

    // Identify S0 and T0 in H
    scc_S_ = scc_of_[s_];
    scc_T_ = scc_of_[t_];

    invalid_.assign(num_vertices(H), 0);
    invalid_[scc_S_] = 1;

    return H;   // Framework will store this in lattice_
}

// -------------------------------------------------------------------------
// Framework subroutines
// -------------------------------------------------------------------------

/**
 * @brief Return the current minimal ideal as a bitvector over SCC vertices.
 *
 * Convention used here:
 *   - invalid_ marks all nodes that are no longer part of our lattice representation.
 *   - Valid nodes define a sublattice of the original lattice.
 *   - This implementation returns the ideal of only invalid nodes; 
 *      that is, the bottom/empty ideal of the sublattice. 
 *
 * @return Ideal bitvector; empty if nothing remains.
 */
Ideal MinCutMaxDisjoint::O_min() {
    const int N = static_cast<int>(topo_TS_.size());
    if (topo_TS_[N - 1]) return Ideal{}; // nothing left
    // Ideal = all invalid vertices
    return invalid_;
}

/**
 * @brief Return the maximal ideal.
 *
 * Convention used here:
 *   - Treat the last topo position as T0 (sink-side component).
 *   - Build an ideal that includes SCCs up to (but excluding) T0.
 *   - This produces a large downset (prefix before the last vertex).
 *
 * @return Ideal bitvector representing the maximal ideal.
 */
Ideal MinCutMaxDisjoint::O_max() {
    // Find T0’s position
    const int N = static_cast<int>(topo_TS_.size());
    const int posT = N - 1;

    // Ideal = all vertices with topo_pos in (0 .. posT)
    Ideal I(topo_TS_.size(), 0);
    for (int i = 0; i < posT; ++i) {
        I[ static_cast<int>(topo_TS_[i]) ] = 1;
    }
    
    return I;
}

/**
 * @brief Disjoint-successor oracle: invalidate endpoints (and their predecessors) of saturated out-edges of X.
 *
 * Intuition:
 *   - The min-cut induced by ideal X uses saturated boundary edges u→v with u in X, v not in X.
 *   - Any future ideal that does not include such v reuses those edges, violating edge-disjointness.
 *   - We scan all saturated out-edges from vertices in X, invalidate them, and backwards-propagate
 *     until no more nodes can be invalidated.
 *
 * @param X Ideal over SCC vertices (bitvector).
 */
void MinCutMaxDisjoint::O_ds(const Ideal& X) {
    LatticeDAG& H = lattice_;
    
    // Identify saturated out-neighbors of X and invalidate them
    std::queue<Vh> q;
    const std::size_t N = invalid_.size();

    for (std::size_t u = 0; u < N; ++u) {
        if (!X[u]) continue;
        for (auto eh : boost::make_iterator_range(out_edges(static_cast<Vh>(u), H))) {
            if (!H[eh].saturated) continue;
            Vh v = boost::target(eh, H);
            if (!X[v]) {
                invalid_[v] = 1;
                q.push(v);      // Keep track of newly marked nodes
            }
        }
    }

    // Propagate down-closure along predecessors of saturated out-neighbors
    while (!q.empty()) {
        Vh x = q.front(); q.pop();
        for (auto pe : boost::make_iterator_range(in_edges(x, H))) {
            Vh p = source(pe, H);
            if (!invalid_[p]) {
                invalid_[p] = 1;
                q.push(p);
            }
        }
    }
}

/**
 * @brief Return true iff the two ideals induce edge-disjoint min-cuts.
 *
 * Implementation:
 *   - Convert both ideals to min-cuts (edge-id sets) using cut_from_H_ideal().
 *   - Intersect the two sets; return false if intersection is non-empty.
 */
bool MinCutMaxDisjoint::are_disjoint(const Ideal& I1,
                                     const Ideal& I2) const {
    MinCutSolution A = cut_from_H_ideal(I1);
    MinCutSolution B = cut_from_H_ideal(I2);                                    
    if (A.edge_ids.empty() || B.edge_ids.empty()) return true;
    std::vector<std::size_t> a=A.edge_ids, b=B.edge_ids, inter;
    std::sort(a.begin(),a.end()); std::sort(b.begin(),b.end());
    std::set_intersection(a.begin(),a.end(), b.begin(),b.end(), std::back_inserter(inter));
    return inter.empty();
}

/**
 * @brief True if target node is marked invalid.
 */
bool MinCutMaxDisjoint::is_empty() const {
    const int N = static_cast<int>(topo_TS_.size());
    return (topo_TS_[N - 1]);
}

/**
 * @brief Convert each ideal to its corresponding min-cut (edge-id set).
 *
 * Simply maps each Ideal via cut_from_H_ideal().
 */
std::vector<MinCutSolution>
MinCutMaxDisjoint::convert_to_solution_space(const std::vector<Ideal>& C) const {
    std::vector<MinCutSolution> out; out.reserve(C.size());
    for (const auto& I : C) out.push_back(cut_from_H_ideal(I));
    return out;
}

// ------------------------------------------------------------------------
// Helpers
// ------------------------------------------------------------------------

/**
 * @brief Prune input graph G_ to its s–t core: keep vertices reachable from s
 *        and that can also reach t. Rebuilds G_ compactly with reindexed vertices.
 *
 * Complexity: O(n + m)
 *
 * Effects:
 *   - s_, t_ are remapped to new indices.
 *   - Edges not in the s–t core are dropped.
 */
void MinCutMaxDisjoint::prune_to_st_core() {
    const std::size_t n = num_vertices(G_);
    if (n == 0) return;

    // forward reachability from s
    std::vector<char> reachS(n, 0);
    {
        std::queue<Vg> q;
        reachS[s_] = 1; q.push(s_);
        while (!q.empty()) {
            Vg u = q.front(); q.pop();
            for (auto e : boost::make_iterator_range(out_edges(u, G_))) {
                Vg v = boost::target(e, G_);
                if (!reachS[v]) { reachS[v] = 1; q.push(v); }
            }
        }
    }

    // backward reachability to t (walk in-edges)
    std::vector<char> reachT(n, 0);
    {
        std::queue<Vg> q;
        reachT[t_] = 1; q.push(t_);
        while (!q.empty()) {
            Vg v = q.front(); q.pop();
            for (auto e : boost::make_iterator_range(in_edges(v, G_))) {
                Vg u = boost::source(e, G_);
                if (!reachT[u]) { reachT[u] = 1; q.push(u); }
            }
        }
    }

    // map old -> new ids for s–t core
    std::vector<int> old2new(n, -1);
    int nn = 0;
    for (std::size_t u = 0; u < n; ++u)
        if (reachS[u] && reachT[u]) old2new[u] = nn++;

    // already core-only
    if ((int)n == nn) return;

    // rebuild compact graph
    MC_G Gnew(nn);
    Vg s_new = static_cast<Vg>(old2new[s_]);
    Vg t_new = static_cast<Vg>(old2new[t_]);

    for (auto e : boost::make_iterator_range(edges(G_))) {
        Vg u = boost::source(e, G_);
        Vg v = boost::target(e, G_);
        int nu = old2new[u], nv = old2new[v];
        if (nu >= 0 && nv >= 0) {
            auto ee = add_edge(static_cast<Vg>(nu), static_cast<Vg>(nv), Gnew);
            if (ee.second) Gnew[ee.first].capacity = G_[e].capacity;
        }
    }

    // swap in pruned graph & terminals
    G_.swap(Gnew);
    s_ = s_new;
    t_ = t_new;
}

/**
 * @brief Convert an ideal over H into the corresponding min-cut δ(I) in the
 *        original graph, returned as a set of edge IDs.
 *
 * Method:
 *   - For each original edge id, check:
 *       (1) edge is saturated (forward),
 *       (2) tail SCC in_I_H == 1 and head SCC in_I_H == 0.
 *     If both hold, include the id in the cut.
 */
MinCutSolution MinCutMaxDisjoint::cut_from_H_ideal(const Ideal& in_I_H) const {
    MinCutSolution S;
    S.edge_ids.reserve(edge_scc_pair_.size() / 8);

    const std::size_t M = edge_scc_pair_.size();
    for (std::size_t id = 0; id < M; ++id) {
        if (!edge_is_saturated_[id]) continue;
        auto [cu, cv] = edge_scc_pair_[id];
        if (cu < 0 || cv < 0) continue;

        const bool tail_in = (cu < (int)in_I_H.size() && in_I_H[(std::size_t)cu]);
        const bool head_in = (cv < (int)in_I_H.size() && in_I_H[(std::size_t)cv]);

        if (tail_in && !head_in)
            S.edge_ids.push_back(id);
    }

    std::sort(S.edge_ids.begin(), S.edge_ids.end());
    S.edge_ids.erase(std::unique(S.edge_ids.begin(), S.edge_ids.end()), S.edge_ids.end());
    return S;
}

/**
 * @brief Build the working residual network and run push-relabel max-flow;
 *        then copy forward residual capacities back into G_[e].residual.
 *
 * Construction:
 *   - WR (working residual) has explicit reverse edges; we maintain:
 *       * capacity map (WR)
 *       * residual capacity map (WR)
 *       * edge_index map for stable indexing
 *       * reverse-edge map (edge → its reverse)
 *       * mapping from WR forward-edge indices → original edges in G_
 *
 * After max-flow:
 *   - For each WR edge that is a forward arc corresponding to some G_ edge,
 *     copy its residual capacity into G_[ge].residual.
 */
void MinCutMaxDisjoint::run_maxflow_and_fill_residual() {
    using WR = boost::adjacency_list<
        boost::vecS, boost::vecS, boost::directedS,
        boost::no_property,
        boost::property<boost::edge_index_t, std::size_t,
        boost::property<boost::edge_capacity_t, long,
        boost::property<boost::edge_residual_capacity_t, long>>> >;

    using Ew = boost::graph_traits<WR>::edge_descriptor;
    using Vw = boost::graph_traits<WR>::vertex_descriptor;

    WR wr(num_vertices(G_));

    auto cap_map = get(boost::edge_capacity,          wr);
    auto res_map = get(boost::edge_residual_capacity, wr);
    auto eidx    = get(boost::edge_index,             wr);
    auto vidx    = get(boost::vertex_index,           wr);

    std::size_t next_index = 0;
    std::vector<Ew> rev_of;     rev_of.reserve(2*num_edges(G_));
    std::vector<char> is_fwd;   is_fwd.reserve(2*num_edges(G_));
    std::vector<Eg>  wridx2ge;  wridx2ge.reserve(2*num_edges(G_));

    auto ensure = [&](std::size_t need){
        if (rev_of.size()  < need) rev_of.resize(need);
        if (is_fwd.size()  < need) is_fwd.resize(need, 0);
        if (wridx2ge.size()< need) wridx2ge.resize(need);
    };

    for (auto ge : boost::make_iterator_range(edges(G_))) {
        auto u = boost::source(ge, G_), v = boost::target(ge, G_);

        Ew f = add_edge(u, v, wr).first;
        Ew r = add_edge(v, u, wr).first;

        put(eidx, f, next_index++);
        put(eidx, r, next_index++);
        ensure(next_index);

        put(cap_map, f, G_[ge].capacity);
        put(cap_map, r, 0L);
        put(res_map, f, 0L);
        put(res_map, r, 0L);

        auto fi = get(eidx, f);
        auto ri = get(eidx, r);
        rev_of[fi] = r; rev_of[ri] = f;

        is_fwd[fi] = 1;
        wridx2ge[fi] = ge;
    }

    auto rev_map = boost::make_iterator_property_map(rev_of.begin(), eidx);

    Vw sw = static_cast<Vw>(s_);
    Vw tw = static_cast<Vw>(t_);
    (void)boost::push_relabel_max_flow(wr, sw, tw, cap_map, res_map, rev_map, vidx);

    for (auto ewr : boost::make_iterator_range(edges(wr))) {
        auto idx = get(eidx, ewr);
        if (!is_fwd[idx]) continue;
        Eg ge = wridx2ge[idx];
        G_[ge].residual = get(res_map, ewr);
    }
}

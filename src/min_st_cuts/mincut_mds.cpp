#include "min_st_cuts/mincut_mds.hpp"

#include <boost/range/iterator_range.hpp>
#include <boost/graph/strong_components.hpp>
#include <boost/graph/topological_sort.hpp>
#include <boost/graph/push_relabel_max_flow.hpp>
#include <algorithm>
#include <queue>

// ============================ ctors =============================================

MinCutMaxDisjoint::MinCutMaxDisjoint(const MC_G& g, Vg s, Vg t)
    : G_(g), s_(s), t_(t) {}

MinCutMaxDisjoint::MinCutMaxDisjoint(MC_G&& g, Vg s, Vg t)
    : G_(std::move(g)), s_(s), t_(t) {}

// ============================ build_initial_lattice ==============================

LatticeG MinCutMaxDisjoint::build_initial_lattice() {
    // Find and delete from the input graph all nodes and arcs not on a path from s to t
    prune_to_st_core();

    // (0) Number original edges by id and init per-edge arrays
    std::size_t m = 0;
    for (auto e : boost::make_iterator_range(edges(G_))) G_[e].id = m++;
    edge_scc_pair_.assign(m, {-1,-1});
    edge_is_saturated_.assign(m, 0);

    // (1) Max-flow 
    run_maxflow_and_fill_residual();

    // (2) Build the "reverse residual" G'
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

    // (3) SCCs on G'
    scc_of_.assign(num_vertices(Gp), -1);
    int k = boost::strong_components(
        Gp, boost::make_iterator_property_map(
                scc_of_.begin(), get(boost::vertex_index, Gp))
    );

    // (4) Contract SCC to build compact representation
    LatticeG H(k);
    for (auto e : boost::make_iterator_range(edges(Gp))) {
        int cu = scc_of_[boost::source(e, Gp)];
        int cv = scc_of_[boost::target(e, Gp)];
        if (cu != cv) add_edge(static_cast<Vh>(cu), static_cast<Vh>(cv), H);
    }

    // Mark arcs that correspond to *saturated forward* original edges; fill per-edge caches
    for (auto ge : boost::make_iterator_range(edges(G_))) {
        const std::size_t id = G_[ge].id;
        const int cu = scc_of_[boost::source(ge, G_)];
        const int cv = scc_of_[boost::target(ge, G_)];
        const long c  = G_[ge].capacity;
        const long rf = G_[ge].residual;
        const long f  = c - rf;
        const bool is_sat_forward = (f == c);   // saturated forward (u->v)

        // Remember this per-original-edge (you already do this elsewhere):
        edge_scc_pair_[id]     = {cu, cv};
        edge_is_saturated_[id] = is_sat_forward ? 1 : 0;

        if (cu != cv && is_sat_forward) {
            auto [eh, ok] = add_edge(static_cast<Vh>(cu), static_cast<Vh>(cv), H);
            // Even if the edge already exists, we just set saturated=true.
            H[eh].saturated = true;
        }
    }

    // (5) Topological order + position map over H
    topo_TS_.clear();
    {
        std::vector<Vh> rev;
        boost::topological_sort(H, std::back_inserter(rev));
        topo_TS_.assign(rev.rbegin(), rev.rend());
    }

    // topo position of every vertex (H’s vertices are 0..k-1)
    topo_pos_.assign(k, -1);
    for (int i = 0; i < (int)topo_TS_.size(); ++i) {
        topo_pos_[ static_cast<int>(topo_TS_[i]) ] = i;
    }

    // Reset the invalid prefix to “none removed yet”
    cutoff_ = -1;

    return H;
}

// ============================ Oracles (in-place) =================================

MinCutSolution MinCutMaxDisjoint::O_min(LatticeG& /*g*/) {
    std::vector<char> in_I_H;
    if (!build_current_minimal_ideal_H(in_I_H)) {
        return MinCutSolution{}; // nothing left
    }
    return cut_from_H_ideal(in_I_H);
}

MinCutSolution MinCutMaxDisjoint::O_max(LatticeG& g) {
    // Find T0’s position
    const int N = static_cast<int>(topo_TS_.size());
    const int posT = N - 1;

    // Ideal = all vertices with topo_pos in (cutoff_ .. posT)
    std::vector<char> in_I_H(topo_TS_.size(), 0);
    for (int i = 0; i < posT; ++i) {
        in_I_H[ static_cast<int>(topo_TS_[i]) ] = 1;
    }
    
    return cut_from_H_ideal(in_I_H);
}

void MinCutMaxDisjoint::O_ds(LatticeG& g, const MinCutSolution& /*X*/) {
    // Reconstruct current minimal ideal X = { TS[cutoff_+1] }
    std::vector<char> in_X;
    if (!build_current_minimal_ideal_H(in_X))
        return;

    int furthest = cutoff_;

    // For each u in X, examine *saturated* out-edges u->v in the lattice and push cutoff
    for (std::size_t idx = 0; idx < in_X.size(); ++idx) {
        if (!in_X[idx]) continue;
        const Vh u = static_cast<Vh>(idx);
        for (auto eh : boost::make_iterator_range(out_edges(u, g))) {
            if (!g[eh].saturated) continue;
            const Vh v = target(eh, g);
            furthest = std::max(furthest, topo_pos_[ static_cast<int>(v) ]);
        }
    }

    // Invalidate prefix up to and including 'furthest'
    cutoff_ = furthest;
}

// ============================ Disjointness / Emptiness ===========================

bool MinCutMaxDisjoint::are_disjoint(const MinCutSolution& A,
                                     const MinCutSolution& B) const {
    if (A.edge_ids.empty() || B.edge_ids.empty()) return true;
    std::vector<std::size_t> a=A.edge_ids, b=B.edge_ids, inter;
    std::sort(a.begin(),a.end()); std::sort(b.begin(),b.end());
    std::set_intersection(a.begin(),a.end(), b.begin(),b.end(), std::back_inserter(inter));
    return inter.empty();
}

bool MinCutMaxDisjoint::is_empty(const LatticeG& /*g*/) const {
    // No valid vertex remains in the TS suffix
    return (cutoff_ + 1) >= static_cast<int>(topo_TS_.size());
}

// ============================ Helpers ============================================

// Keep only vertices reachable from s AND that can reach t; rebuild G_ compactly.
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
                Vg v = target(e, G_);
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

    // already core-only?
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

// Extract saturated forward boundary of an ideal over SCC DAG using per-original-edge caches.
MinCutSolution MinCutMaxDisjoint::cut_from_H_ideal(const std::vector<char>& in_I_H) const {
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

// Build the *current* minimal ideal over H as the singleton { TS[cutoff_+1] }.
// Return false if no valid nodes remain.
bool MinCutMaxDisjoint::build_current_minimal_ideal_H(std::vector<char>& in_I_H) const {
    const int N = static_cast<int>(topo_TS_.size());
    const int first = cutoff_ + 1;
    if (first >= N) return false;
    in_I_H.assign(topo_pos_.size(), 0);
    in_I_H[ static_cast<int>(topo_TS_[first]) ] = 1;
    return true;
}

// --- Max-flow that fills G_[e].residual on forward arcs ---------------------------
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

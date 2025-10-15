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
    // (0) Number original edges by id and init per-edge arrays
    std::size_t m = 0;
    for (auto e : boost::make_iterator_range(edges(G_))) G_[e].id = m++;
    edge_scc_pair_.assign(m, {-1,-1});
    edge_is_saturated_.assign(m, 0);
    edge_active_.assign(m, 1); // all allowed initially

    // (1) Max-flow to fill forward residuals
    run_maxflow_and_fill_residual();


    // (2) Build the "reverse residual" G' as specified (we only need SCCs)
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
    scc_S_ = (s_ < scc_of_.size()) ? scc_of_[s_] : -1;
    scc_T_ = (t_ < scc_of_.size()) ? scc_of_[t_] : -1;

    // Build SCC DAG H
    using HDAG = boost::adjacency_list<
        boost::vecS, boost::vecS, boost::directedS, boost::no_property, boost::no_property>;
    HDAG H(k);
    for (auto e : boost::make_iterator_range(edges(Gp))) {
        int cu = scc_of_[boost::source(e, Gp)];
        int cv = scc_of_[boost::target(e, Gp)];
        if (cu != cv) add_edge(static_cast<Vh>(cu), static_cast<Vh>(cv), H);
    }

    // Cache per-original-edge (cu,cv) and whether saturated forward exists
    for (auto e : boost::make_iterator_range(edges(G_))) {
        std::size_t id = G_[e].id;
        int cu = scc_of_[boost::source(e, G_)];
        int cv = scc_of_[boost::target(e, G_)];
        edge_scc_pair_[id] = {cu, cv};

        long c  = G_[e].capacity;
        long rf = G_[e].residual;
        long f  = c - rf;
        edge_is_saturated_[id] = (f == c) ? 1 : 0;
    }

    // (4) Prune T and successors; S and predecessors
    auto bfs_mark = [&](int start, auto& mark, const HDAG& Gdir) {
        if (start < 0 || start >= k) return;
        std::queue<int> q; mark.assign(k, 0); mark[start] = 1; q.push(start);
        while (!q.empty()) {
            int x = q.front(); q.pop();
            for (auto ee : boost::make_iterator_range(out_edges(static_cast<Vh>(x), Gdir))) {
                int y = static_cast<int>(target(ee, Gdir));
                if (!mark[y]) { mark[y] = 1; q.push(y); }
            }
        }
    };

    std::vector<char> succT(k, 0), succS(k, 0), predS(k, 0);
    bfs_mark(scc_T_, succT, H); // successors of T
    bfs_mark(scc_S_, succS, H); // successors of S

    // reverse adjacency for predecessors of S
    std::vector<std::vector<int>> rev(k);
    for (auto ee : boost::make_iterator_range(edges(H))) {
        int u = static_cast<int>(boost::source(ee, H));
        int v = static_cast<int>(boost::target(ee, H));
        rev[v].push_back(u);
    }
    auto bfs_mark_rev = [&](int start, auto& mark) {
        if (start < 0 || start >= k) return;
        std::queue<int> q; mark.assign(k, 0); mark[start] = 1; q.push(start);
        while (!q.empty()) {
            int x = q.front(); q.pop();
            for (int y : rev[x]) if (!mark[y]) { mark[y] = 1; q.push(y); }
        }
    };
    bfs_mark_rev(scc_S_, predS);

    // Keep-only mapping
    scc_to_lattice_.assign(k, -1);
    int n_keep = 0;
    for (int v = 0; v < k; ++v) {
        bool drop = (v == scc_S_) || (v == scc_T_) || succT[v] || predS[v];
        if (!drop) scc_to_lattice_[v] = n_keep++;
    }

    // Reduced lattice
    LatticeG L(n_keep);
    topo_.clear();
    if (n_keep > 0) {
        // Project H edges
        for (auto ee : boost::make_iterator_range(edges(H))) {
            int u = static_cast<int>(boost::source(ee, H));
            int v = static_cast<int>(boost::target(ee, H));
            int nu = scc_to_lattice_[u];
            int nv = scc_to_lattice_[v];
            if (nu >= 0 && nv >= 0) add_edge(static_cast<Vh>(nu), static_cast<Vh>(nv), L);
        }
        // Topo order
        std::vector<Vh> revtopo;
        boost::topological_sort(L, std::back_inserter(revtopo));
        topo_.assign(revtopo.rbegin(), revtopo.rend());
    }

    // Project Succ(S) to reduced lattice
    reachS_lattice_.assign(n_keep, 0);
    for (int old = 0; old < k; ++old) {
        int nv = scc_to_lattice_[old];
        if (nv >= 0 && succS[old]) reachS_lattice_[nv] = 1;
    }

    return L;
}

// ============================ Oracles (in-place) =================================

MinCutSolution MinCutMaxDisjoint::O_min(LatticeG& g) {
    std::vector<char> in_I(num_vertices(g), 0);
    build_minimal_ideal(g, in_I);               // empty set in reduced lattice
    return cut_from_ideal_in_place(in_I);       // δ(I ∪ Succ(S)) using only active & saturated edges
}

MinCutSolution MinCutMaxDisjoint::O_max(LatticeG& g) {
    std::vector<char> in_I(num_vertices(g), 1);
    build_maximal_ideal(g, in_I);               // all visible reduced-lattice vertices
    return cut_from_ideal_in_place(in_I);
}

void MinCutMaxDisjoint::O_ds(LatticeG& /*g*/, const MinCutSolution& S) {
    // Ban used original edges so future cuts are edge-disjoint
    for (auto id : S.edge_ids)
        if (id < edge_active_.size()) edge_active_[id] = 0;
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
    // No available original edges => cannot form new edge-disjoint cuts
    return std::none_of(edge_active_.begin(), edge_active_.end(),
                        [](char x){ return x!=0; }) ;
}

// ============================ Helpers ============================================

void MinCutMaxDisjoint::build_minimal_ideal(LatticeG& /*g*/,
                                            std::vector<char>& in_I) const {
    // Minimal ideal is the empty set in the reduced lattice
    std::fill(in_I.begin(), in_I.end(), 0);
}

void MinCutMaxDisjoint::build_maximal_ideal(LatticeG& g,
                                            std::vector<char>& in_I) const {
    // All current reduced-lattice vertices (we aren't deactivating lattice vertices here)
    std::fill(in_I.begin(), in_I.end(), 1);
    (void)g;
}

MinCutSolution MinCutMaxDisjoint::cut_from_ideal_in_place(const std::vector<char>& in_I) const {
    MinCutSolution S;
    S.edge_ids.reserve(edge_scc_pair_.size()/8);

    for (std::size_t id = 0; id < edge_scc_pair_.size(); ++id) {
        if (!edge_active_[id])      continue;   // banned
        if (!edge_is_saturated_[id]) continue;  // must be saturated forward arc

        auto [cu, cv] = edge_scc_pair_[id];
        if (cu < 0 || cv < 0) continue;

        int nu = (cu < (int)scc_to_lattice_.size() ? scc_to_lattice_[cu] : -1);
        int nv = (cv < (int)scc_to_lattice_.size() ? scc_to_lattice_[cv] : -1);
        if (nu < 0 || nv < 0) continue;         // pruned endpoints

        const bool tail_in = ((std::size_t)nu < in_I.size() && in_I[nu]) ||
                             (nu < (int)reachS_lattice_.size() && reachS_lattice_[nu]);
        const bool head_in = ((std::size_t)nv < in_I.size() && in_I[nv]) ||
                             (nv < (int)reachS_lattice_.size() && reachS_lattice_[nv]);

        if (tail_in && !head_in) S.edge_ids.push_back(id);
    }

    std::sort(S.edge_ids.begin(), S.edge_ids.end());
    S.edge_ids.erase(std::unique(S.edge_ids.begin(), S.edge_ids.end()), S.edge_ids.end());
    return S;
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

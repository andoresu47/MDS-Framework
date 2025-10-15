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
    // (0) Assign stable original-edge ids 0..m-1
    {
        std::size_t eid = 0;
        for (auto e : boost::make_iterator_range(edges(G_))) {
            G_[e].id = eid++;
        }
        edge_scc_pair_.assign(eid, {-1,-1});
        edge_is_saturated_.assign(eid, 0);
    }

    // (1) Run max-flow and fill forward residuals in G_ (G_[e].residual)
    run_maxflow_and_fill_residual();

    // (2) Build the reverse-residual graph G' per your spec
    using Gprime = boost::adjacency_list<
        boost::vecS, boost::vecS, boost::directedS,
        boost::no_property,
        boost::property<boost::edge_capacity_t, long>   // residual cap (not used later, but kept for clarity)
    >;

    Gprime Gp(num_vertices(G_));
    auto add_arc = [&](Vg u, Vg v, long rc) {
        if (rc <= 0) return;
        auto ee = add_edge(u, v, Gp);
        if (ee.second) put(boost::edge_capacity, Gp, ee.first, rc);
    };

    for (auto e : boost::make_iterator_range(edges(G_))) {
        auto u = boost::source(e, G_);
        auto v = boost::target(e, G_);
        const long c  = G_[e].capacity;
        const long rf = G_[e].residual;           // forward residual after max-flow
        const long f  = c - rf;                   // actual flow on (u,v)

        // Construct G' (reverse of the usual residual):
        //  - if 0 < f < c : add (v->u) with c-f   AND (u->v) with f
        //  - if f == c    : add (u->v) with c
        //  - if f == 0    : add (v->u) with c
        if (0 < f && f < c) {
            add_arc(v, u, c - f);
            add_arc(u, v, f);
        } else if (f == c) {
            add_arc(u, v, f);
        } else { // f == 0
            add_arc(v, u, c);
        }
    }

    // (3) SCCs of G' → SCC ids in [0..k-1]
    scc_of_.assign(num_vertices(Gp), -1);
    int k = boost::strong_components(
        Gp,
        boost::make_iterator_property_map(
            scc_of_.begin(), get(boost::vertex_index, Gp))
    );
    scc_S_ = (s_ < scc_of_.size()) ? scc_of_[s_] : -1;
    scc_T_ = (t_ < scc_of_.size()) ? scc_of_[t_] : -1;

    // Build SCC-DAG H
    using HDAG = boost::adjacency_list<
        boost::vecS, boost::vecS, boost::directedS,
        boost::no_property, boost::no_property>;
    HDAG H(k);
    for (auto e : boost::make_iterator_range(edges(Gp))) {
        int cu = scc_of_[boost::source(e, Gp)];
        int cv = scc_of_[boost::target(e, Gp)];
        if (cu != cv) add_edge(static_cast<Vh>(cu), static_cast<Vh>(cv), H);
    }

    // Also cache, per original edge id, its (scc(u), scc(v)) in G' and whether it's saturated forward
    for (auto e : boost::make_iterator_range(edges(G_))) {
        std::size_t id = G_[e].id;
        int cu = scc_of_[boost::source(e, G_)];
        int cv = scc_of_[boost::target(e, G_)];
        edge_scc_pair_[id] = {cu, cv};

        const long c  = G_[e].capacity;
        const long rf = G_[e].residual;
        const long f  = c - rf;
        // A forward arc (u->v) exists in G' iff f > 0 (includes saturated case f==c).
        // "Saturated forward boundary" = forward arc exists AND original edge is saturated (f == c).
        edge_is_saturated_[id] = (f == c) ? 1 : 0;
    }

    // (4) Prune: remove T and all its successors; and S and all its predecessors
    std::vector<char> succT(k, 0), predS(k, 0);

    auto bfs_mark = [&](int start, auto& mark, const HDAG& Gdir) {
        if (start < 0 || start >= k) return;
        std::queue<int> q;
        mark[start] = 1; q.push(start);
        while (!q.empty()) {
            int x = q.front(); q.pop();
            for (auto ee : boost::make_iterator_range(out_edges(static_cast<Vh>(x), Gdir))) {
                int y = static_cast<int>(target(ee, Gdir));
                if (!mark[y]) { mark[y] = 1; q.push(y); }
            }
        }
    };
    bfs_mark(scc_T_, succT, H);

    // Reverse adjacency to find predecessors
    std::vector<std::vector<int>> rev(k);
    for (auto ee : boost::make_iterator_range(edges(H))) {
        int u = static_cast<int>(boost::source(ee, H));
        int v = static_cast<int>(boost::target(ee, H));
        rev[v].push_back(u);
    }
    auto bfs_mark_reverse = [&](int start, auto& mark) {
        if (start < 0 || start >= k) return;
        std::queue<int> q;
        mark[start] = 1; q.push(start);
        while (!q.empty()) {
            int x = q.front(); q.pop();
            for (int y : rev[x]) {
                if (!mark[y]) { mark[y] = 1; q.push(y); }
            }
        }
    };
    bfs_mark_reverse(scc_S_, predS);

    // Keep only SCCs that are not dropped
    scc_to_lattice_.assign(k, -1);
    int n_keep = 0;
    for (int v = 0; v < k; ++v) {
        bool drop = (v == scc_S_) || (v == scc_T_) || succT[v] || predS[v];
        if (!drop) scc_to_lattice_[v] = n_keep++;
    }

    // Forward reachability from S (its successors in H)
    std::vector<char> succS(k, 0);
    bfs_mark(scc_S_, succS, H);

    // Build reduced lattice L
    LatticeG L(n_keep);

    // Project H edges into L (assign compact edge_index 0..E-1)
    std::size_t eidx = 0;
    for (auto ee : boost::make_iterator_range(edges(H))) {
        int u = static_cast<int>(boost::source(ee, H));
        int v = static_cast<int>(boost::target(ee, H));
        int nu = scc_to_lattice_[u];
        int nv = scc_to_lattice_[v];
        if (nu >= 0 && nv >= 0) {
            Eh le = add_edge(static_cast<Vh>(nu), static_cast<Vh>(nv), L).first;
            put(boost::edge_index, L, le, eidx++);
        }
    }

    // Cache topological order of L
    topo_.clear();
    {
        std::vector<Vh> revtopo;
        boost::topological_sort(L, std::back_inserter(revtopo));
        topo_.assign(revtopo.rbegin(), revtopo.rend());
    }

    // Project Succ(S) to the reduced lattice
    reachS_lattice_.assign(n_keep, 0);
    for (int old = 0; old < k; ++old) {
        int nv = scc_to_lattice_[old];
        if (nv >= 0 && succS[old]) reachS_lattice_[nv] = 1;
    }

    return L;
}

MinCutSolution MinCutMaxDisjoint::O_min(const View& Lstar) {
    // Minimal ideal inside the current filtered lattice.
    // With S/T pruning and Succ(S) handled in cut extraction, the minimal ideal
    // is simply the empty set over the visible lattice vertices.
    std::vector<char> in_I(num_vertices(*Lstar.g), 0);
    build_minimal_ideal(Lstar, in_I);
    return cut_from_ideal(Lstar, in_I);     // uses I ∪ Succ(S) when forming δ(I)
}

MinCutSolution MinCutMaxDisjoint::O_max(const View& Lstar) {
    // Maximal ideal inside the current filtered lattice.
    // After pruning, "all currently visible vertices" form an ideal.
    std::vector<char> in_I(num_vertices(*Lstar.g), 1);
    build_maximal_ideal(Lstar, in_I);
    return cut_from_ideal(Lstar, in_I);     // uses I ∪ Succ(S) when forming δ(I)
}

View MinCutMaxDisjoint::O_ds(const View& Lstar, const MinCutSolution& S) {
    // Disjoint successors inside the current sublattice:
    // ban (hide) every saturated original edge used by S from future cuts.
    // This is done by flipping the edge mask bits corresponding to original edge IDs.
    auto next = Lstar.deep_copy();        // independent masks to mutate safely
    if (next.emask) {
        for (auto id : S.edge_ids) {
            if (id < next.emask->size())
                (*next.emask)[id] = 0;    // hide this lattice "edge id" from future views
        }
    }
    return next;
}

bool MinCutMaxDisjoint::are_disjoint(const MinCutSolution& A,
                                     const MinCutSolution& B) const {
    // Two cuts are disjoint if they share no edge IDs.
    // We check for an empty intersection using sorted set logic.

    if (A.edge_ids.empty() || B.edge_ids.empty())
        return true; // trivially disjoint

    std::vector<std::size_t> a = A.edge_ids;
    std::vector<std::size_t> b = B.edge_ids;

    std::sort(a.begin(), a.end());
    std::sort(b.begin(), b.end());

    std::vector<std::size_t> inter;
    inter.reserve(std::min(a.size(), b.size()));

    std::set_intersection(
        a.begin(), a.end(),
        b.begin(), b.end(),
        std::back_inserter(inter)
    );

    return inter.empty();
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
    // Working residual graph: needs edge_index so we can index reverse edges.
    using WR = boost::adjacency_list<
        boost::vecS, boost::vecS, boost::directedS,
        boost::no_property,
        boost::property<boost::edge_index_t, std::size_t,
        boost::property<boost::edge_capacity_t, long,
        boost::property<boost::edge_residual_capacity_t, long>>> // residual cap
    >;

    using Ew = boost::graph_traits<WR>::edge_descriptor;
    using Vw = boost::graph_traits<WR>::vertex_descriptor;

    const std::size_t n = num_vertices(G_);
    WR wr(n);

    // Property maps on WR
    auto cap_map = get(boost::edge_capacity,           wr);
    auto res_map = get(boost::edge_residual_capacity,  wr);
    auto eidx    = get(boost::edge_index,              wr);
    auto vidx    = get(boost::vertex_index,            wr);

    // Add paired residual arcs for every original edge (u->v):
    //  forward: capacity = c_uv
    //  reverse: capacity = 0
    std::vector<Ew> forward_wr_edges; forward_wr_edges.reserve(num_edges(G_));
    std::size_t next_index = 0;

    // We'll also store reverse mapping edge_index -> reverse edge descriptor
    // and mark which indices correspond to FORWARD residual arcs
    std::vector<Ew>   rev_of;     rev_of.reserve(2 * num_edges(G_));
    std::vector<char> is_forward; is_forward.reserve(2 * num_edges(G_));

    // Keep a parallel vector mapping forward WR edge index -> original G_ edge
    std::vector<Eg> wridx2ge;

    for (auto ge : boost::make_iterator_range(edges(G_))) {
        auto u = boost::source(ge, G_);
        auto v = boost::target(ge, G_);

        // Create forward and reverse residual edges
        Ew f = add_edge(u, v, wr).first;
        Ew r = add_edge(v, u, wr).first;

        // Assign indices explicitly
        put(eidx, f, next_index++);
        put(eidx, r, next_index++);

        // Ensure rev_of / is_forward / wridx2ge are large enough
        if (rev_of.size() < next_index) { rev_of.resize(next_index); }
        if (is_forward.size() < next_index) { is_forward.resize(next_index, 0); }
        if (wridx2ge.size() < next_index) { wridx2ge.resize(next_index); }

        // Set capacities; residual map is output of max-flow (init to 0)
        put(cap_map, f, G_[ge].capacity);
        put(cap_map, r, 0L);
        put(res_map, f, 0L);
        put(res_map, r, 0L);

        // Wire reverse-edge index mapping
        auto fi = get(eidx, f);
        auto ri = get(eidx, r);
        rev_of[fi] = r;
        rev_of[ri] = f;

        // Mark forward edges and remember original edge for fi
        is_forward[fi] = 1;
        wridx2ge[fi]   = ge;
    }

    // Wrap rev_of (edge_index -> reverse Ew) as a property map
    auto rev_map = boost::make_iterator_property_map(rev_of.begin(), eidx);

    // s_, t_ are vertex descriptors in G_; with vecS vertex set, indices match
    Vw sw = static_cast<Vw>(s_);
    Vw tw = static_cast<Vw>(t_);

    // Run max-flow; residual capacities written to res_map
    (void)boost::push_relabel_max_flow(wr, sw, tw, cap_map, res_map, rev_map, vidx);

    // Copy FORWARD residual capacities back to original graph
    for (auto ewr : boost::make_iterator_range(edges(wr))) {
        auto idx = get(eidx, ewr);
        if (!is_forward[idx]) continue;               // skip reverse arcs
        Eg ge = wridx2ge[idx];
        G_[ge].residual = get(res_map, ewr);          // forward residual on (u->v)
    }
}

void MinCutMaxDisjoint::build_minimal_ideal(const View& Lstar,
                                            std::vector<char>& in_I) const {
    // Minimal extension beyond Succ(S) is the empty set.
    std::fill(in_I.begin(), in_I.end(), 0);

    // If you wanted to enforce “down-closure forced by visibility,” you could
    // add logic here, but for min-cuts this is unnecessary: δ(Succ(S)) is a valid
    // minimum cut, and further additions are handled by O_min/O_ds iterations.
    (void)Lstar; // silence unused param if not using it
}


void MinCutMaxDisjoint::build_maximal_ideal(const View& Lstar,
                                            std::vector<char>& in_I) const {
    // Start with all zeros, then set 1 for visible vertices in the current filtered view.
    std::fill(in_I.begin(), in_I.end(), 0);

    auto FG = Lstar.filtered();
    for (auto v : boost::make_iterator_range(vertices(FG))) {
        auto idx = get(Lstar.vimap, v);            // vertex index in [0..|L|-1]
        if ((std::size_t)idx < in_I.size()) in_I[(std::size_t)idx] = 1;
    }

    // Down-closure holds automatically since we included all visible vertices.
}

MinCutSolution MinCutMaxDisjoint::cut_from_ideal(const View& Lstar,
                                                 const std::vector<char>& in_I) const {
    MinCutSolution S;
    const std::size_t m = edge_scc_pair_.size();
    S.edge_ids.reserve(m / 8);

    for (std::size_t id = 0; id < m; ++id) {
        if (!edge_is_saturated_[id]) continue;   // only saturated forward arcs

        auto [cu, cv] = edge_scc_pair_[id];
        if (cu < 0 || cv < 0) continue;

        int nu = (cu < (int)scc_to_lattice_.size() ? scc_to_lattice_[cu] : -1);
        int nv = (cv < (int)scc_to_lattice_.size() ? scc_to_lattice_[cv] : -1);
        if (nu < 0 || nv < 0) continue;          // endpoints pruned

        // Effective ideal side is I U Succ(S)
        const bool tail_in  = ((std::size_t)nu < in_I.size() && in_I[nu]) || (nu < (int)reachS_lattice_.size() && reachS_lattice_[nu]);
        const bool head_in  = ((std::size_t)nv < in_I.size() && in_I[nv]) || (nv < (int)reachS_lattice_.size() && reachS_lattice_[nv]);

        if (tail_in && !head_in) S.edge_ids.push_back(id);
    }

    std::sort(S.edge_ids.begin(), S.edge_ids.end());
    S.edge_ids.erase(std::unique(S.edge_ids.begin(), S.edge_ids.end()), S.edge_ids.end());
    return S;
}

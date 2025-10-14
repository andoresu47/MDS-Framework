#pragma once
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/filtered_graph.hpp>
#include <boost/range/iterator_range.hpp>
#include <memory>
#include <vector>
#include <algorithm>

// -----------------------------------------------------------------------------
// MaxDisjointSolutionsFramework
//
// Abstract base for problems whose feasible solutions form a distributive lattice
// with oracles:
//   - O_min(L*): bottom element of current sublattice L*
//   - O_max(L*): top element of current sublattice L*
//   - O_ds(L*, S): sublattice of disjoint successors of S inside L*
//   - are_disjoint(A,B): predicate on solutions
//
// Compact lattice representation = a Boost graph L (often a DAG).
// Sublattices are represented *without copying* using boost::filtered_graph
// driven by mask vectors (vmask/emask).
//
// Requirements on L:
//   - L must provide a vertex index map (vertex_index_t).
//   - L must provide an edge index map (edge_index_t) if you intend to filter edges.
//     (If you don’t filter edges, you may keep emask null and the edge predicate
//      will accept everything.)
// -----------------------------------------------------------------------------
template<class LatticeGraph, class Solution>
class MaxDisjointSolutionsFramework {
public:
    using L      = LatticeGraph;
    using V      = typename boost::graph_traits<L>::vertex_descriptor;
    using E      = typename boost::graph_traits<L>::edge_descriptor;
    using VIndex = typename boost::property_map<const L, boost::vertex_index_t>::const_type;
    using EIndex = typename boost::property_map<const L, boost::edge_index_t>::const_type;

    // ----------------------- Sublattice view (zero-copy) ----------------------
    struct LatticeView {
        const L* g = nullptr;

        // Masks: 1 = visible, 0 = hidden. If nullptr, predicate defaults to "true".
        std::shared_ptr<std::vector<char>> vmask;
        std::shared_ptr<std::vector<char>> emask;

        // Index maps (cached from the graph)
        VIndex vimap{};
        EIndex eimap{};

        // Predicates for filtered_graph
        struct VPred {
            VIndex vimap{};
            std::shared_ptr<std::vector<char>> vmask;
            bool operator()(V v) const {
                if (!vmask) return true;
                auto i = get(vimap, v);
                return (*vmask)[static_cast<std::size_t>(i)];
            }
        };
        struct EPred {
            EIndex eimap{};
            std::shared_ptr<std::vector<char>> emask;
            bool operator()(E e) const {
                if (!emask) return true;
                auto i = get(eimap, e);
                return (*emask)[static_cast<std::size_t>(i)];
            }
        };

        using FG = boost::filtered_graph<const L, EPred, VPred>;

        FG filtered() const {
            return FG(*g, EPred{eimap, emask}, VPred{vimap, vmask});
        }

        std::size_t count_vertices() const {
            std::size_t c = 0;
            for (auto v : boost::make_iterator_range(vertices(filtered()))) (void)v, ++c;
            return c;
        }
        std::size_t count_edges() const {
            std::size_t c = 0;
            for (auto e : boost::make_iterator_range(edges(filtered()))) (void)e, ++c;
            return c;
        }

        // convenience: create a deep-copied view (copy-on-write semantics)
        LatticeView deep_copy() const {
            LatticeView out = *this;
            if (vmask) out.vmask = std::make_shared<std::vector<char>>(*vmask);
            if (emask) out.emask = std::make_shared<std::vector<char>>(*emask);
            return out;
        }
    };

    virtual ~MaxDisjointSolutionsFramework() = default;

    // ---------------------------- Oracle interface ----------------------------
    // Build the initial compact lattice representation L.
    virtual L build_initial_lattice() = 0;

    // Bottom/top element of current sublattice L* (operate on view.filtered()).
    virtual Solution O_min(const LatticeView& Lstar) = 0;
    virtual Solution O_max(const LatticeView& Lstar) = 0;

    // Sublattice of disjoint successors of S *inside the given sublattice* Lstar.
    // Return a NEW view (typically by copying masks and flipping bits).
    virtual LatticeView O_ds(const LatticeView& Lstar, const Solution& S) = 0;

    // Problem-specific disjointness test
    virtual bool are_disjoint(const Solution& A, const Solution& B) const = 0;

    // Optional: size/emptiness of a view (override if you have a better notion)
    virtual std::size_t view_size(const LatticeView& Lstar) const { return Lstar.count_vertices(); }
    virtual bool is_empty(const LatticeView& Lstar) const { return view_size(Lstar) == 0; }

    // ---------------------------- Helper: full view ---------------------------
    // Creates a full-lattice view (all vertices/edges visible).
    LatticeView make_full_view() {
        LatticeView P;
        P.g      = &lattice_;
        const L& gl = lattice_;
        P.vimap  = get(boost::vertex_index, gl);
        P.eimap  = get(boost::edge_index,   gl);
        // If your L does not have edge_index, either:
        //  - add it when building L; or
        //  - leave P.emask=nullptr and let EPred accept all edges.
        const auto nv = num_vertices(lattice_);
        const auto ne = num_edges(lattice_);
        P.vmask = std::make_shared<std::vector<char>>(nv, 1);
        if (ne > 0) P.emask = std::make_shared<std::vector<char>>(ne, 1);
        return P;
    }

    // ------------------------ Algorithms ---------------------
    // Unbounded version of Algorithm 1 (returns all found solutions).
    // Semantics:
    //   - Xz := O_max(P) computed ONCE from the initial lattice P (fixed guard).
    //   - Iterate: X := O_min(current), PX := O_ds(current, X).
    //   - Stop when X intersects Xz; add that final X.
    std::vector<Solution> find_max_disjoint() {
        lattice_ = build_initial_lattice();
        auto P   = make_full_view();
        if (is_empty(P)) return {};

        const Solution Xz = O_max(P); // fixed guard from original lattice
        std::vector<Solution> C; C.reserve(8);  // allocate space for ~8 solutions upfront; we usually expect a small number.

        // Start with the minimal of P and its successors IN P
        Solution X = O_min(P);
        auto PX     = O_ds(P, X);

        while (are_disjoint(X, Xz)) {
            C.push_back(X);
            if (is_empty(PX)) break;      // no disjoint successors left
            X  = O_min(PX);               // minimal in current sublattice
            PX = O_ds(PX, X);             // successors INSIDE current sublattice
        }

        // Add the final X that intersects Xz
        C.push_back(X);
        return C;
    }

    // k-bounded version (returns up to k solutions).
    std::vector<Solution> find_max_disjoint(int k) {
        std::vector<Solution> C;
        if (k <= 0) return C;

        lattice_ = build_initial_lattice();
        auto P   = make_full_view();
        if (is_empty(P)) return C;

        const Solution Xz = O_max(P); // fixed guard
        Solution X = O_min(P);
        auto PX     = O_ds(P, X);

        while (are_disjoint(X, Xz)) {
            C.push_back(X);
            if ((int)C.size() >= k) return C;  // cap at k
            if (is_empty(PX)) break;
            X  = O_min(PX);
            PX = O_ds(PX, X);
        }

        if ((int)C.size() < k) C.push_back(X);
        return C;
    }

protected:
    L lattice_; // Owned compact lattice representation
};

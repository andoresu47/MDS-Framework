#pragma once
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/graph_concepts.hpp>
#include <boost/range/iterator_range.hpp>
#include <concepts>
#include <vector>

// --------------------- Minimal concept: we only need a Boost graph ---------------------
template<class L>
concept LatticeGraphModel =
    requires(L& g,
             const L& cg,
             typename boost::graph_traits<L>::vertex_descriptor v,
             typename boost::graph_traits<L>::edge_descriptor   e)
{
    vertices(cg); edges(cg); out_edges(v, cg); // basic traversal
    // No property-map requirements here; subclasses can use any internal properties they like.
};

// --------------------------------------------------------------------------------------
// MaxDisjointSolutionsFramework
// --------------------------------------------------------------------------------------
// - Subclasses implement the oracles that *mutate* the lattice graph in place:
//     * O_min(L&):   compute the "bottom" solution in the current state of L
//     * O_max(L&):   compute the "top"    solution in the current state of L
//     * O_ds (L&, S):shrink L in-place so only disjoint successors (w.r.t. S) remain
// - The framework only orchestrates Algorithm 1 of the paper.
// - Disjointness is decided by are_disjoint(A,B).
// --------------------------------------------------------------------------------------
template<LatticeGraphModel LatticeGraph, class Solution>
class MaxDisjointSolutionsFramework {
public:
    using L = LatticeGraph;
    using V = typename boost::graph_traits<L>::vertex_descriptor;
    using E = typename boost::graph_traits<L>::edge_descriptor;

    virtual ~MaxDisjointSolutionsFramework() = default;

    // Build the initial lattice graph (already reduced/compact if you wish).
    // You can attach any vertex/edge properties (e.g., { bool active; }) and
    // use them inside the oracles to ignore “scrapped” parts.
    virtual L build_initial_lattice() = 0;

    // Oracles operating on the *current* lattice, in-place:
    virtual Solution O_min(L& g) = 0;                     // bottom element of current sublattice
    virtual Solution O_max(L& g) = 0;                     // top    element of current sublattice
    virtual void     O_ds (L& g, const Solution& s) = 0;  // shrink g to disjoint successors of s
    virtual bool     are_disjoint(const Solution& a, const Solution& b) const = 0;

    // Optional: subclasses can override if they track a more meaningful emptiness notion
    virtual bool is_empty(const L& g) const {
        // Default: empty if the graph has no vertices (override if you “deactivate” instead of deleting)
        return num_vertices(g) == 0;
    }

    // ----------------------------- Unbounded Algorithm 1 -----------------------------
    // Iterate minimal solutions with O_ds shrink, stop when intersecting Xz := O_max(P); include final X.
    std::vector<Solution> find_max_disjoint() {
        lattice_ = build_initial_lattice();
        if (is_empty(lattice_)) return {};

        const Solution Xz = O_max(lattice_);  // fixed guard from the initial state
        std::vector<Solution> C; C.reserve(8);

        Solution X = O_min(lattice_);

        while ( are_disjoint(X, Xz) ) {
            C.push_back(X);
            // shrink to disjoint successors *in place*
            O_ds(lattice_, X);
            if (is_empty(lattice_)) break;     // no successors left
            X = O_min(lattice_);
        }
        C.push_back(X); // final (possibly intersecting) solution
        return C;
    }

    // ----------------------------- k-bounded variant --------------------------------
    std::vector<Solution> find_max_disjoint(int k) {
        std::vector<Solution> C;
        if (k <= 0) return C;

        lattice_ = build_initial_lattice();
        if (is_empty(lattice_)) return C;

        const Solution Xz = O_max(lattice_);
        Solution X = O_min(lattice_);

        while ( are_disjoint(X, Xz) ) {
            C.push_back(X);
            if ((int)C.size() >= k) return C;
            O_ds(lattice_, X);
            if (is_empty(lattice_)) break;
            X = O_min(lattice_);
        }
        if ((int)C.size() < k) C.push_back(X);
        return C;
    }

protected:
    L lattice_;  // owned, mutates in place through oracles
};

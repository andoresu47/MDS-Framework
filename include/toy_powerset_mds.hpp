#pragma once
#include "max_disjoint_framework.hpp"
#include <boost/graph/adjacency_list.hpp>
#include <boost/range/iterator_range.hpp>
#include <algorithm>
#include <vector>

// ---------------- Toy subclass over a 3-element “powerset” ----------------
using ToyLattice =
    boost::adjacency_list<
        boost::vecS, boost::vecS, boost::directedS,
        boost::property<boost::vertex_index_t, std::size_t>,
        boost::property<boost::edge_index_t,   std::size_t>
    >;

struct ToySolution {
    std::vector<std::size_t> elems; // subset of {0,1,2}
    bool operator==(const ToySolution& o) const { return elems == o.elems; }
};

// --------------------------------------------------------------------------
// ToyPowersetMDS
//   A toy implementation of the MaxDisjointSolutionsFramework where
//   the "lattice" is just the powerset of {0,1,2} and solutions are
//   singletons {0}, {1}, {2}.
// --------------------------------------------------------------------------
class ToyPowersetMDS
: public MaxDisjointSolutionsFramework<ToyLattice, ToySolution> {
public:
    // Build a 3-node graph with vertex indices {0,1,2}.
    ToyLattice build_initial_lattice() override {
        return ToyLattice(3);
    }

    // Return the smallest-index visible vertex as the minimal singleton.
    ToySolution O_min(const LatticeView& Lstar) override {
        ToySolution s;
        auto FG = Lstar.filtered();
        for (auto v : boost::make_iterator_range(vertices(FG))) {
            auto id = get(boost::vertex_index, *Lstar.g, v);
            s.elems.push_back(static_cast<std::size_t>(id));
            break;
        }
        return s;
    }

    // Return the largest-index visible vertex as the maximal singleton.
    ToySolution O_max(const LatticeView& Lstar) override {
        ToySolution s;
        auto FG = Lstar.filtered();
        std::size_t last = static_cast<std::size_t>(-1);
        for (auto v : boost::make_iterator_range(vertices(FG))) {
            last = get(boost::vertex_index, *Lstar.g, v);
        }
        if (last != static_cast<std::size_t>(-1))
            s.elems.push_back(last);
        return s;
    }

    // Disjoint successors oracle: hide all vertices in S.
    LatticeView O_ds(const LatticeView& Lstar, const ToySolution& S) override {
        auto next = Lstar.deep_copy();
        if (next.vmask) {
            for (auto id : S.elems)
                if (id < next.vmask->size())
                    (*next.vmask)[id] = 0;
        }
        return next;
    }

    // Solutions are disjoint if they share no elements.
    bool are_disjoint(const ToySolution& A, const ToySolution& B) const override {
        std::vector<std::size_t> a = A.elems, b = B.elems, inter;
        std::sort(a.begin(), a.end());
        std::sort(b.begin(), b.end());
        std::set_intersection(a.begin(), a.end(), b.begin(), b.end(),
                              std::back_inserter(inter));
        return inter.empty();
    }

    // Optional: treat view as empty if no visible vertices remain.
    std::size_t view_size(const LatticeView& Lstar) const override {
        return Lstar.count_vertices();
    }
    bool is_empty(const LatticeView& Lstar) const override {
        return view_size(Lstar) == 0;
    }
};

#include "toy_example/toy_powerset_mds.hpp"
#include <boost/range/iterator_range.hpp>
#include <algorithm>
#include <vector>

// Build a 3-node graph with vertex indices {0,1,2}.
ToyLattice ToyPowersetMDS::build_initial_lattice() {
    return ToyLattice(3);
}

// Return the smallest-index visible vertex as the minimal singleton.
ToySolution ToyPowersetMDS::O_min(const LatticeView& Lstar) {
    ToySolution s;
    auto FG = Lstar.filtered();
    for (auto v : boost::make_iterator_range(vertices(FG))) {
        auto id = get(boost::vertex_index, *Lstar.g, v);
        s.elems.push_back(static_cast<std::size_t>(id));
        break; // singleton minimal
    }
    return s;
}

// Return the largest-index visible vertex as the maximal singleton.
ToySolution ToyPowersetMDS::O_max(const LatticeView& Lstar) {
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
MaxDisjointSolutionsFramework<ToyLattice, ToySolution>::LatticeView
ToyPowersetMDS::O_ds(const LatticeView& Lstar, const ToySolution& S) {
    auto next = Lstar.deep_copy();
    if (next.vmask) {
        for (auto id : S.elems)
            if (id < next.vmask->size())
                (*next.vmask)[id] = 0;
    }
    return next;
}

// Solutions are disjoint if they share no elements.
bool ToyPowersetMDS::are_disjoint(const ToySolution& A, const ToySolution& B) const {
    std::vector<std::size_t> a = A.elems, b = B.elems, inter;
    std::sort(a.begin(), a.end());
    std::sort(b.begin(), b.end());
    std::set_intersection(a.begin(), a.end(), b.begin(), b.end(),
                          std::back_inserter(inter));
    return inter.empty();
}

std::size_t ToyPowersetMDS::view_size(const LatticeView& Lstar) const {
    return Lstar.count_vertices();
}

bool ToyPowersetMDS::is_empty(const LatticeView& Lstar) const {
    return view_size(Lstar) == 0;
}

#include "toy_example/toy_powerset_mds.hpp"
#include <boost/range/iterator_range.hpp>

ToyLattice ToyPowersetMDS::build_initial_lattice() {
    ToyLattice g(3);
    for (auto v : boost::make_iterator_range(vertices(g))) g[v].active = true;
    return g;
}

ToySolution ToyPowersetMDS::O_min() {
    ToyLattice& g = lattice_;
    ToySolution s;
    for (auto v : boost::make_iterator_range(vertices(g))) {
        if (g[v].active) { s.elems.push_back(static_cast<std::size_t>(v)); break; }
    }
    return s; // empty if none active
}

ToySolution ToyPowersetMDS::O_max() {
    ToyLattice& g = lattice_;
    ToySolution s;
    std::size_t last = static_cast<std::size_t>(-1);
    for (auto v : boost::make_iterator_range(vertices(g))) {
        if (g[v].active) last = static_cast<std::size_t>(v);
    }
    if (last != static_cast<std::size_t>(-1)) s.elems.push_back(last);
    return s;
}

void ToyPowersetMDS::O_ds(const ToySolution& s) {
    ToyLattice& g = lattice_;
    // Shrink in place: mark chosen elements inactive so future solutions can’t reuse them
    for (auto id : s.elems) {
        if (id < num_vertices(g)) g[static_cast<ToyLattice::vertex_descriptor>(id)].active = false;
    }
}

bool ToyPowersetMDS::are_disjoint(const ToySolution& a, const ToySolution& b) const {
    if (a.elems.empty() || b.elems.empty()) return true;
    std::vector<std::size_t> x = a.elems, y = b.elems, inter;
    std::sort(x.begin(), x.end());
    std::sort(y.begin(), y.end());
    std::set_intersection(x.begin(), x.end(), y.begin(), y.end(), std::back_inserter(inter));
    return inter.empty();
}

bool ToyPowersetMDS::is_empty() const {
    const ToyLattice& g = lattice_;
    for (auto v : boost::make_iterator_range(vertices(g))) if (g[v].active) return false;
    return true;
}

std::vector<ToySolution>
    ToyPowersetMDS::convert_to_solution_space(const std::vector<ToySolution>& C) const {
        return std::vector<ToySolution>(C.begin(), C.end());
    }

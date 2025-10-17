#pragma once
#include <vector>
#include <cstddef>

// -----------------------------------------------------------------------------
// MaxDisjointSolutionsFramework
//
// Templated over:
//   - LatticeRepresentation: a (usually compact) data structure representing the lattice. 
//   - InternalSolution: what the oracles manipulate (e.g., an ideal)
//   - ExternalSolution: what you want returned to the caller (e.g., a cut)
//       defaults to InternalSolution (identity)
// -----------------------------------------------------------------------------
template<class LatticeRepresentation, class InternalSolution, class ExternalSolution = InternalSolution>
class MaxDisjointSolutionsFramework {
public:
    using L = LatticeRepresentation;
    using ISol = InternalSolution;
    using ESol = ExternalSolution;

    virtual ~MaxDisjointSolutionsFramework() = default;

    // Build initial lattice representation and store it in lattice_.
    virtual L build_initial_lattice() = 0;

    // Oracles act on the framework's internal state (including lattice_).
    virtual ISol O_min() = 0;                      // bottom element of current sublattice
    virtual ISol O_max() = 0;                      // top element of current sublattice
    virtual void O_ds(const ISol& s) = 0;          // shrink to disjoint successors of s

    // Disjointness test in the oracle's internal space.
    virtual bool are_disjoint(const ISol& a, const ISol& b) const = 0;

    // Progress notion for the current sublattice (e.g., "no valid vertices left").
    virtual bool is_empty() const = 0;

    // Conversion from internal to external solution space, if not used, just implement identity.
    virtual std::vector<ESol>
    convert_to_solution_space(const std::vector<ISol>& C) const = 0;

    // Algorithm
    std::vector<ESol> find_max_disjoint(int k = -1) {
        const bool bounded = (k >= 0);
        lattice_ = build_initial_lattice();
        if (is_empty()) return {};

        ISol Xz = O_max();
        std::vector<ISol> Cint; Cint.reserve(8);

        ISol X = O_min();
        O_ds(X);                        // restrict to disjoint successors of X inside current state

        while (are_disjoint(X, Xz)) {
            Cint.push_back(X);
            if (bounded && static_cast<int>(Cint.size()) == k)
                return convert_to_solution_space(Cint);
            if (is_empty()) break;
            X = O_min();
            O_ds(X);
        }

        // Add the final X that intersects Xz
        if (!bounded || static_cast<int>(Cint.size()) < k)
            Cint.push_back(X);

        return convert_to_solution_space(Cint);
    }

protected:
    L lattice_;   // owned lattice representation
};

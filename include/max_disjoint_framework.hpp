/**
 * @file max_disjoint_framework.hpp
 * @brief Generic framework for problems whose feasible solutions form a distributive lattice,
 *        supporting a bottom-up search for a maximum set of disjoint solutions.
 *
 * @tparam LatticeRepresentation  Compact in-memory representation of the lattice (owned in @c lattice_).
 * @tparam InternalSolution       The type manipulated by the oracles (e.g., an ideal/bitset over SCCs).
 * @tparam ExternalSolution       The type returned to callers (e.g., a cut as edge IDs).
 *                                Defaults to @c InternalSolution (identity conversion).
 *
 * ## Overview
 * Subclasses implement the following subroutines:
 *  - @c build_initial_lattice(): construct the initial lattice. 
 *  - @c O_min() : return the current minimal element of the (shrinking) sublattice.
 *  - @c O_max() : return the fixed maximal element of the *initial* lattice (guard).
 *  - @c O_ds(s) : shrink the current sublattice to “disjoint successors” of @p s.
 *  - @c are_disjoint(a,b) : predicate deciding disjointness in internal space.
 *  - @c is_empty() : indicate whether the current sublattice has been exhausted. 
 *
 * The algorithm @c find_max_disjoint() then iteratively:
 *  1) builds the initial lattice and assigns it to @c lattice_,
 *  2) computes the guard Xz = O_max(),
 *  3) repeatedly picks X = O_min(), records it, and calls O_ds(X),
 *  4) stops when X intersects Xz (per the framework semantics),
 *  5) returns all collected solutions (optionally capped at k) after converting
 *     them to @c ExternalSolution via @c convert_to_solution_space().
 *
 * The intent is to decouple internal mechanics (e.g., ideals) from the public
 * output type (e.g., cuts). If @c InternalSolution == @c ExternalSolution,
 * @c convert_to_solution_space can simply return its input.
 *
 * ## Minimal subclass contract
 *  - Store any state you need (e.g., indices, topological order, cursors) as members.
 *  - Implement:
 *      - build_initial_lattice()
 *      - O_min(), O_max(), O_ds()
 *      - are_disjoint(), is_empty()
 *      - convert_to_solution_space()
 *
 * ## Thread-safety
 * Instances are not thread-safe (mutable shared state in @c lattice_ and
 * subclass members). Use one instance per thread or guard externally.
 *
 * ## Complexity
 * The framework’s control flow is O(#solutions) oracle calls. Total cost is driven
 * by your oracle implementations.
 */

#pragma once
#include <vector>
#include <cstddef>

template<class LatticeRepresentation, class InternalSolution, class ExternalSolution = InternalSolution>
class MaxDisjointSolutionsFramework {
public:
    /// @name Type aliases
    /// @{
    using L = LatticeRepresentation;
    using ISol = InternalSolution;
    using ESol = ExternalSolution;
    /// @}

    virtual ~MaxDisjointSolutionsFramework() = default;

    // -------------------------------------------------------------------------
    // Subroutines (must override)
    // -------------------------------------------------------------------------

    /**
     * @brief Build the initial lattice.
     * @return The built lattice (usually the same object stored in @c lattice_).
     */
    virtual L build_initial_lattice() = 0;

    /**
     * @brief Return the current minimal element of the (possibly already shrunk) (sub)lattice stored in @c lattice_.
     */
    virtual ISol O_min() = 0;

    /**
     * @brief Return the maximal element of the lattice stored in @c lattice_.
     *
     * Typically computed once after @c build_initial_lattice and treated as a guard.
     */
    virtual ISol O_max() = 0;                     

    /**
     * @brief Shrink the current search space to the disjoint successors of @p X.
     *
     * This usually updates the subclass’s internal state so that subsequent calls 
     * to @c O_min() produce new, disjoint solutions.
     */
    virtual void O_ds(const ISol& X) = 0;          

    /**
     * @brief Disjointness predicate on internal solutions.
     */
    virtual bool are_disjoint(const ISol& a, const ISol& b) const = 0;

     /**
     * @brief Whether the current sublattice has been exhausted.
     */
    virtual bool is_empty() const = 0;

    /**
     * @brief Convert a batch of internal solutions to the external solution space.
     *
     * If @c InternalSolution == @c ExternalSolution, simply return @p C unchanged.
     */
    virtual std::vector<ESol>
    convert_to_solution_space(const std::vector<ISol>& C) const = 0;

    // -------------------------------------------------------------------------
    // Algorithm
    // -------------------------------------------------------------------------

    /**
     * @brief Find a maximum (or up to @p k) cardinality set of pairwise-disjoint solutions.
     *
     * @param k  If negative (default), return the full set. If non-negative, cap the output size at @p k.
     *
     * @return A vector of @c ExternalSolution produced by converting the internally
     *         collected solutions via @c convert_to_solution_space().
     *
     * @note Framework semantics (high-level):
     *   - Xz := O_max() is a fixed guard from the initial lattice.
     *   - Loop picks X := O_min() from the current state, then calls O_ds(X)
     *     to restrict to disjoint successors.
     *   - The loop stops when X is no longer disjoint from Xz; the final X is
     *     still appended.
     */
    std::vector<ESol> find_max_disjoint(int k = -1) {
        const bool bounded = (k >= 0);

        // 1) Build initial lattice/state
        lattice_ = build_initial_lattice();
        if (is_empty()) return {};

        // 2) Guard (fixed max of initial lattice)
        ISol Xz = O_max();

        // 3) Accumulate internal solutions
        std::vector<ISol> Cint; Cint.reserve(8);

        // 4) First pick + shrink
        ISol X = O_min();
        O_ds(X);                        

        // 5) Main loop
        while (are_disjoint(X, Xz)) {
            Cint.push_back(X);
            if (bounded && static_cast<int>(Cint.size()) == k)
                return convert_to_solution_space(Cint);
            if (is_empty()) break;
            X = O_min();
            O_ds(X);
        }

        // 6) Append the final final X that intersects Xz (if room remains)
        if (!bounded || static_cast<int>(Cint.size()) < k)
            Cint.push_back(X);

        return convert_to_solution_space(Cint);
    }

protected:
    /** 
     * @brief Owned lattice representation (subclass defines its layout and invariants). 
    */
    L lattice_;   
};

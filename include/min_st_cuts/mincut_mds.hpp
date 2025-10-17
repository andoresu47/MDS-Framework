/**
 * @file mincut_mds.hpp
 * @brief Definition of the MinCutMaxDisjoint solver — a specialization of
 *        MaxDisjointSolutionsFramework for computing multiple disjoint
 *        minimum s–t cuts in a directed, weighted graph.
 *
 * The lattice of all minimum s–t cuts is represented compactly as a SCC-DAG
 * (H) derived from the reversed residual graph of a maximum flow.
 * Each vertex of H represents a strongly connected component (SCC),
 * and each edge corresponds to at least one saturated edge in the
 * original graph’s residual structure.
 *
 * The solver implements the standard three oracles in the space of ideals of H :
 *   - O_min(): smallest ideal (minimal/leftmost s-t mincut)
 *   - O_max(): largest ideal (maximal/rightmost s-t mincut)
 *   - O_ds(): restricts the current lattice (representation) to disjoint successors of an ideal
 *
 * These are used by the base framework to iteratively extract disjoint minimum s–t cuts.
 */

#pragma once
#include "max_disjoint_framework.hpp"

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graph_traits.hpp>
#include <vector>
#include <cstddef>

// -------------------------------------------------------------------------
// Input graph representation
// -------------------------------------------------------------------------

/**
 * @struct MC_EdgeProps
 * @brief Per-edge data in the input graph.
 *
 * Each edge stores:
 *  - @c capacity : input weight (capacity c_uv)
 *  - @c residual : residual capacity (filled after running max-flow)
 *  - @c id       : edge identifier in [0..m)
 */
struct MC_EdgeProps {
    long        capacity = 0;   
    long        residual = 0;   
    std::size_t id       = 0;   
};

/**
 * @brief Input directed graph type for the minimum s–t cut instance.
 *
 * Vertices are stored as @c vecS (contiguous indexing), edges as @c vecS,
 * with bidirectional adjacency (for fast in/out traversal).
 */
using MC_G = boost::adjacency_list<
    boost::vecS, boost::vecS, boost::bidirectionalS,
    boost::no_property,
    MC_EdgeProps
>;

/// Vertex and edge descriptor types for @c MC_G.
using Vg = boost::graph_traits<MC_G>::vertex_descriptor;
using Eg = boost::graph_traits<MC_G>::edge_descriptor;

// -------------------------------------------------------------------------
// Compact lattice representation (SCC DAG)
// -------------------------------------------------------------------------

/**
 * @struct LEdge
 * @brief Edge property for the compact SCC-DAG lattice representation.
 *
 * We are interested in edges between SCCs in the residual graph 
 * that correspond to at least one saturated edge in the original graph, 
 * so we mark them as such.
 */
struct LEdge {
    bool saturated = false; ///< true if edge corresponds to a saturated original arc
};

/**
 * @brief Compact directed acyclic graph representing the lattice of minimum cuts.
 *
 * Each vertex corresponds to a strongly connected component (SCC)
 * in the (reversed) residual graph G′. Each directed edge represents 
 * reachability between SCCs. 
 */
using LatticeDAG = boost::adjacency_list<
    boost::vecS, boost::vecS, boost::directedS,
    boost::no_property,
    LEdge
>;

/// Vertex and edge descriptor types for @c LatticeDAG.
using Vh = boost::graph_traits<LatticeDAG>::vertex_descriptor;
using Eh = boost::graph_traits<LatticeDAG>::edge_descriptor;

// -------------------------------------------------------------------------
// Solution types
// -------------------------------------------------------------------------

/**
 * @struct MinCutSolution
 * @brief External solution type: a minimum s–t cut as a set of original edge IDs.
 *
 * The vector @c edge_ids stores the identifiers of edges (from the original
 * input graph) that form the boundary δ(I) of some ideal I in the lattice representation.
 */
struct MinCutSolution {
    std::vector<std::size_t> edge_ids; // original edge ids forming δ(I)
    bool operator==(const MinCutSolution& o) const { return edge_ids == o.edge_ids; }
};

/**
 * @brief Internal solution type: ideal representation over the lattice DAG.
 *
 * An ideal I is represented as a bitvector of length |V(lattice_)|,
 * where entry i = 1 if vertex i (an SCC) belongs to the ideal.
 */
using Ideal = std::vector<char>;

// -------------------------------------------------------------------------
// MinCutMaxDisjoint class
// -------------------------------------------------------------------------

/**
 * @class MinCutMaxDisjoint
 * @brief Solver implementing the maximum set of disjoint minimum s–t cuts.
 *
 * This class derives from MaxDisjointSolutionsFramework and instantiates it
 * for the lattice of minimum cuts derived from a flow network.
 *
 * The solver builds:
 *  - A maximum flow on the input graph G.
 *  - The reversed residual graph G′.
 *  - The SCC-DAG H of G′, which compactly represents the lattice of minimum cuts.
 *
 * Each oracle (O_min, O_max, O_ds) operates directly on H and a topological
 * order of its vertices, without needing a explicit lattice access.
 *
 * The algorithm proceeds as follows:
 *  1. Compute the initial lattice via build_initial_lattice().
 *  2. Repeatedly extract the minimal ideal (O_min), restrict to disjoint
 *     successors (O_ds), and test disjointness via are_disjoint().
 *  3. Return all disjoint cut solutions via convert_to_solution_space().
 */
class MinCutMaxDisjoint
: public MaxDisjointSolutionsFramework<LatticeDAG, Ideal, MinCutSolution> {
public:

    /**
     * @brief Construct from a constant input graph.
     * @param g  Input graph (will be copied)
     * @param s  Source vertex
     * @param t  Sink vertex
     */
    explicit MinCutMaxDisjoint(const MC_G& g, Vg s, Vg t);

    /**
     * @brief Construct from an rvalue input graph.
     * @param g  Input graph (will be moved)
     * @param s  Source vertex
     * @param t  Sink vertex
     */
    explicit MinCutMaxDisjoint(MC_G&& g, Vg s, Vg t);

    // ------------------------------------------------------------------------
    // Framework overrides
    // ------------------------------------------------------------------------

    /** @copydoc MaxDisjointSolutionsFramework::build_initial_lattice() */
    LatticeDAG  build_initial_lattice() override;       

    /** @copydoc MaxDisjointSolutionsFramework::O_min() */
    Ideal       O_min() override;                       

    /** @copydoc MaxDisjointSolutionsFramework::O_max() */
    Ideal       O_max() override;                       

    /** @copydoc MaxDisjointSolutionsFramework::O_ds() */
    void        O_ds (const Ideal& X) override;         

    /** @copydoc MaxDisjointSolutionsFramework::are_disjoint() */
    bool        are_disjoint(const Ideal& I1, const Ideal& I2) const override;    

    /** @copydoc MaxDisjointSolutionsFramework::is_empty() */
    bool        is_empty() const override;              

    /** @copydoc MaxDisjointSolutionsFramework::convert_to_solution_space() */
    std::vector<MinCutSolution>
    convert_to_solution_space(const std::vector<Ideal>& C) const override;

    // ------------------------------------------------------------------------
    // Optional accessors (for testing / analysis)
    // ------------------------------------------------------------------------

    /// @return reference to the (possibly pruned) input graph.
    const MC_G& instance() const { return G_; }

    /// @return vertex-to-SCC partition mapping.
    const std::vector<int>& scc_partition() const { return scc_of_; }

    /// @return topologically sorted order of SCC vertices in the lattice DAG.
    const std::vector<Vh>& topo_order() const { return topo_TS_; }

private:
    // ------------------------------------------------------------------------
    // Internal helpers
    // ------------------------------------------------------------------------

    /**
     * @brief Prune input graph to its s–t core.
     *
     * Removes all vertices and edges not lying on any path from s to t.
     * This guarantees that subsequent computations (max-flow, SCC-DAG) are
     * restricted to the relevant portion of the graph. Runs in O(V + E) time. 
     */
    void prune_to_st_core();                

    /**
     * @brief Run a maximum s–t flow on G_ and record residual capacities.
     *
     * After execution, every edge e in G_ satisfies:
     *  - @c G_[e].residual = residual capacity (forward direction)
     *  - @c G_[e].capacity = original capacity
     *  - @c G_[e].id       = unique ID
     */
    void run_maxflow_and_fill_residual();   

    /**
     * @brief Extract the set of saturated edges crossing δ(I) from an ideal I.
     *
     * @param in_I_H Bitvector marking which SCC nodes belong to the ideal.
     * @return MinCutSolution listing the original edge IDs that form the boundary.
     */
    MinCutSolution cut_from_H_ideal(const std::vector<char>& in_I_H) const; 
    
private:
    // ------------------------------------------------------------------------
    // Member data
    // ------------------------------------------------------------------------

    // --- Input graph (possibly pruned to s–t core) ---
    MC_G G_;        ///< Input instance graph
    Vg s_{};        ///< Source vertex
    Vg t_{};        ///< Sink vertex

    // --- Cached lattice data (SCC decomposition and mapping) ---
    std::vector<int> scc_of_;       ///< vertex → SCC id mapping
    int scc_S_ = -1;                ///< SCC id of source component
    int scc_T_ = -1;                ///< SCC id of sink component

    // --- Edge-level mapping to SCC DAG ---
    std::vector<std::pair<int,int>> edge_scc_pair_;     ///< (SCC(u), SCC(v)) for each original edge
    std::vector<char> edge_is_saturated_;               ///< whether original edge was saturated

    // --- Oracle runtime state ---
    std::vector<Vh> topo_TS_;   ///< Topological order of SCC DAG vertices
    std::vector<int> topo_pos_; ///< vertex → index in topo_TS_
    int cutoff_ = -1;           ///< Position of current cutoff (invalid prefix in topo_TS_)
};

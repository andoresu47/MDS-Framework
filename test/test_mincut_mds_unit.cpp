#include <boost/test/unit_test.hpp>

#include <vector>
#include <algorithm>
#include <cstddef>

#include "min_st_cuts/mincut_mds.hpp"

// Small helper to compare two sorted vectors of size_t
static void check_equal_set(const std::vector<std::size_t>& got,
                            const std::vector<std::size_t>& expected) {
    std::vector<std::size_t> a = got, b = expected;
    std::sort(a.begin(), a.end());
    std::sort(b.begin(), b.end());
    BOOST_TEST(a == b, boost::test_tools::per_element());
}

// Compute total cost of a cut using solver.instance() edge capacities
static long cut_cost(const MinCutMaxDisjoint& solver,
                     const std::vector<std::size_t>& edge_ids) {
    const MC_G& G = solver.instance();
    long cost = 0;
    for (auto id : edge_ids) {
        // Build a lookup table id -> edge once; here do a quick scan per id
        // (OK for small tests; keep it simple)
        for (auto e : boost::make_iterator_range(edges(G))) {
            if (G[e].id == id) {
                cost += G[e].capacity;
                break;
            }
        }
    }
    return cost;
}

BOOST_AUTO_TEST_CASE(test_on_dBMS23_graph)
{
    // Input:
    // 14 23 0 13
    // (u v w=1) in this exact order
    const std::size_t n = 14;
    const std::size_t s = 0, t = 13;

    // The 23 directed edges (u, v, cap=1), in the order provided
    const std::vector<std::pair<std::size_t, std::size_t>> edges_in_order = {
        {0,1},{0,2},{0,3},{0,4},{0,8},
        {1,5},
        {2,5},{2,8},
        {3,8},{3,9},
        {4,9},
        {5,6},{5,7},
        {6,10},{6,13},
        {7,11},{7,13},
        {8,7},
        {9,12},
        {10,13},
        {11,13},
        {12,11},{12,13}
    };

    // Build instance G with unit capacities (cap=1)
    MC_G G(n);
    for (auto [u,v] : edges_in_order) {
        auto e = add_edge(u, v, G);
        BOOST_REQUIRE_MESSAGE(e.second, "add_edge succeeded");
        G[e.first].capacity = 1;
    }

    // Construct solver and build lattice representation (also runs max-flow internally, etc.)
    MinCutMaxDisjoint solver(std::move(G), static_cast<Vg>(s), static_cast<Vg>(t));
    (void)solver.build_initial_lattice();

    // Get minimal and maximal ideals
    Ideal Imin = solver.O_min();
    Ideal Imax = solver.O_max();

    // Convert to cut space
    std::vector<Ideal> tmp; tmp.reserve(2);
    tmp.push_back(Imin);
    tmp.push_back(Imax);
    std::vector<MinCutSolution> cuts = solver.convert_to_solution_space(tmp);

    // Expected edge-id sets (based on input order -> IDs 0..22)
    const std::vector<std::size_t> expected_min = { 0, 1, 17, 18 };      // O_min
    const std::vector<std::size_t> expected_max = { 11, 16, 20, 22 };    // O_max

    // Sanity: we got two cuts; verify sizes
    BOOST_REQUIRE(cuts.size() == 2);

    // Compare as sets (order-insensitive)
    check_equal_set(cuts[0].edge_ids, expected_min);
    check_equal_set(cuts[1].edge_ids, expected_max);

    // Get max. set of disjoint st mincuts 
    cuts = solver.find_max_disjoint();

    // Compare as sets (again)
    check_equal_set(cuts[0].edge_ids, expected_min);
    check_equal_set(cuts[1].edge_ids, expected_max);
}

BOOST_AUTO_TEST_CASE(test_on_PQ82_graph_weighted)
{
    // Input:
    // 14 31 0 13
    const std::size_t n = 14;
    const std::size_t s = 0, t = 13;

    // The 23 directed edges (u, v, w), in input order
    const std::vector<std::tuple<std::size_t,std::size_t,long>> edges_in_order = {
        {0,1,3},
        {0,2,10},
        {0,3,16},
        {0,4,6},
        {1,5,5},
        {2,1,2},
        {2,3,1},
        {2,5,3},
        {2,7,2},
        {2,9,4},
        {3,7,3},
        {4,3,1},
        {4,7,5},
        {4,8,6},
        {5,9,9},
        {6,5,2},
        {6,7,4},
        {6,11,2},
        {7,3,3},
        {8,7,4},
        {8,12,3},
        {8,13,4},
        {9,13,12},
        {10,6,7},
        {10,9,2},
        {10,13,4},
        {11,8,3},
        {11,10,1},
        {11,12,2},
        {11,13,11},
        {12,13,3}
    };

    // Build instance G with unit capacities (cap=1)
    MC_G G(n);
    for (auto [u,v,w] : edges_in_order) {
        auto e = add_edge(u, v, G);
        BOOST_REQUIRE_MESSAGE(e.second, "add_edge succeeded");
        G[e.first].capacity = static_cast<long>(w);
    }

    // Construct solver and build lattice representation (also runs max-flow internally, etc.)
    MinCutMaxDisjoint solver(std::move(G), static_cast<Vg>(s), static_cast<Vg>(t));
    (void)solver.build_initial_lattice();

    // Expected cuts expressed in input-order edge indices:
    // Cut 1: { 0, 2, 4, 5, 6 }  (cost 18)
    // Cut 2: { 7, 11 }          (cost 18)
    const std::vector<std::size_t> expected_cut_1_input_ids = {0, 2, 4, 5, 6};
    const std::vector<std::size_t> expected_cut_2_input_ids = {7, 11};

    // Find maximal set of disjoint min-cuts
    auto cuts = solver.find_max_disjoint();

    // We expect exactly two cuts
    BOOST_REQUIRE(cuts.size() == 2);

    // Compare as sets (order-insensitive)
    check_equal_set(cuts[0].edge_ids, expected_cut_1_input_ids);
    check_equal_set(cuts[1].edge_ids, expected_cut_2_input_ids);

    // Costs must both be 18
    BOOST_CHECK_EQUAL(cut_cost(solver, cuts[0].edge_ids), 18);
    BOOST_CHECK_EQUAL(cut_cost(solver, cuts[1].edge_ids), 18);
}

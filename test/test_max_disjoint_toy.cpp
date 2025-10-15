// To run this test only: g++ -std=c++17 -O2 test_main.cpp test_max_disjoint_toy.cpp -I /path/to/boost -I /path/to/include/folder -o mds_toy_tests
#include <boost/test/unit_test.hpp>
#include "toy_example/toy_powerset_mds.hpp"

// Helper: flatten solutions for easy comparison.
static std::vector<std::vector<std::size_t>>
flatten(const std::vector<ToySolution>& v) {
    std::vector<std::vector<std::size_t>> out;
    out.reserve(v.size());
    for (auto& s : v) out.push_back(s.elems);
    return out;
}

// --------------------------- Tests ---------------------------

BOOST_AUTO_TEST_CASE(unbounded_returns_three_singletons) {
    ToyPowersetMDS mds;
    auto sols = mds.find_max_disjoint();
    auto flat = flatten(sols);

    // Expect { {0}, {1}, {2} }
    BOOST_REQUIRE_EQUAL(flat.size(), 3u);
    BOOST_TEST(flat[0] == std::vector<std::size_t>{0});
    BOOST_TEST(flat[1] == std::vector<std::size_t>{1});
    BOOST_TEST(flat[2] == std::vector<std::size_t>{2});
}

BOOST_AUTO_TEST_CASE(k_bounded_caps_at_two) {
    ToyPowersetMDS mds;
    auto sols = mds.find_max_disjoint(2);
    auto flat = flatten(sols);

    // Expect { {0}, {1} }
    BOOST_REQUIRE_EQUAL(flat.size(), 2u);
    BOOST_TEST(flat[0] == std::vector<std::size_t>{0});
    BOOST_TEST(flat[1] == std::vector<std::size_t>{1});
}

BOOST_AUTO_TEST_CASE(empty_lattice_yields_empty_result) {
    struct EmptyMDS : ToyPowersetMDS {
        ToyLattice build_initial_lattice() override { return ToyLattice(0); }
    } mds;

    BOOST_TEST(mds.find_max_disjoint().empty());
    BOOST_TEST(mds.find_max_disjoint(5).empty());
}

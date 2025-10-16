#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <cmath>
#include <cstdlib>

#include "min_st_cuts/mincut_mds.hpp"

static void die(const std::string& msg) {
    std::cerr << "Error: " << msg << "\n";
    std::exit(1);
}

int main(int argc, char** argv) {
    std::ios::sync_with_stdio(false);
    std::cin.tie(nullptr);

    // Usage: ./demo_min_cut <input_file or -> [scale]
    std::string path = "-";
    double scale = 1.0;
    if (argc >= 2) path = argv[1];
    if (argc >= 3) {
        std::istringstream iss(argv[2]);
        if (!(iss >> scale) || !(scale > 0.0)) die("invalid scale");
    }

    std::ifstream fin;
    std::istream* in = &std::cin;
    if (path != "-") {
        fin.open(path);
        if (!fin) die("cannot open file: " + path);
        in = &fin;
    }

    // Read header: n m s t
    std::size_t n, m, s, t;
    if (!(*in >> n >> m >> s >> t)) die("failed to read n m s t");

    // Build original graph G
    MC_G G(n);
    for (std::size_t i = 0; i < m; ++i) {
        std::size_t u, v; double w;
        if (!(*in >> u >> v >> w)) die("failed to read edge line " + std::to_string(i+2));
        if (u >= n || v >= n) die("vertex index out of range");

        long cap;
        if (scale != 1.0) {
            long long sc = std::llround(w * scale);
            if (sc < 0) die("negative capacity after scaling");
            cap = static_cast<long>(sc);
        } else {
            long long sc = std::llround(w);
            if (std::fabs(w - static_cast<double>(sc)) > 1e-12) {
                std::cerr << "Warning: non-integer weight without scale; rounding " << w
                          << " -> " << sc << "\n";
            }
            if (sc < 0) die("negative capacity");
            cap = static_cast<long>(sc);
        }

        auto e = add_edge(u, v, G);
        if (!e.second) die("failed to add edge");
        G[e.first].capacity = cap;
    }
    if (s >= n || t >= n) die("s or t out of range");

    // Build compact lattice
    MinCutMaxDisjoint solver(std::move(G), static_cast<Vg>(s), static_cast<Vg>(t));
    LatticeG L = solver.build_initial_lattice();

    // -------- Print reduced lattice in "n' m'\n u v..." format (0-based) --------
    std::size_t LV = num_vertices(L);
    std::size_t LE = num_edges(L);
    std::cout << LV << " " << LE << "\n";
    for (auto e : boost::make_iterator_range(edges(L))) {
        auto u = static_cast<std::size_t>(source(e, L));
        auto v = static_cast<std::size_t>(target(e, L));
        std::cout << u << " " << v << "\n";
    }

    // -------- Print supernode contents: original vertices grouped by reduced node ----
    auto orig2reduced = solver.original_to_reduced_map();   // size n, values in [-1..LV-1]
    std::vector<std::vector<std::size_t>> groups(LV);
    for (std::size_t v = 0; v < n; ++v) {
        int r = orig2reduced[v];
        if (r >= 0) groups[static_cast<std::size_t>(r)].push_back(v);
    }

    std::cout << "# Supernodes (reduced vertex id : original vertices)\n";
    for (std::size_t r = 0; r < LV; ++r) {
        std::cout << r << " :";
        for (std::size_t i = 0; i < groups[r].size(); ++i) {
            std::cout << (i==0 ? " {" : " ");
            std::cout << groups[r][i];
        }
        std::cout << (groups[r].empty() ? " {}\n" : "}\n");
    }

    // -------- (Optional) show one min and one max min-cut for a quick sanity check ----
    MinCutSolution Xmin = solver.O_min(L);
    MinCutSolution Xmax = solver.O_max(L);

    auto print_cut = [](const char* name, const MinCutSolution& S) {
        std::cout << "# " << name << " cut edge_ids = ";
        if (S.edge_ids.empty()) { std::cout << "{}\n"; return; }
        std::cout << "{";
        for (std::size_t i = 0; i < S.edge_ids.size(); ++i) {
            std::cout << S.edge_ids[i];
            if (i + 1 != S.edge_ids.size()) std::cout << ", ";
        }
        std::cout << "}\n";
    };
    print_cut("O_min", Xmin);
    print_cut("O_max", Xmax);

    return 0;
}

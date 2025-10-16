#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <cmath>
#include <cstdlib>
#include <iomanip> 
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graph_traits.hpp>

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
                std::cerr << "Warning: non-integer weight without scale; rounding "
                          << w << " -> " << sc << "\n";
            }
            if (sc < 0) die("negative capacity");
            cap = static_cast<long>(sc);
        }

        auto e = add_edge(u, v, G);
        if (!e.second) die("failed to add edge");
        G[e.first].capacity = cap;
        // residual & id will be filled by the solver
    }
    if (s >= n || t >= n) die("s or t out of range");

    // Build compact lattice (this runs max-flow internally)
    MinCutMaxDisjoint solver(std::move(G), static_cast<Vg>(s), static_cast<Vg>(t));
    LatticeG L = solver.build_initial_lattice();  

    std::cout << "# Finding max-disjoint minimum s–t cuts...\n";
    auto cuts = solver.find_max_disjoint();
    std::cout << "Number of disjoint min-cuts found: " << cuts.size() << "\n";

    // Build an edge-descriptor lookup by id for cost lookup and pretty printing
    const MC_G& Gsol = solver.instance();
    std::vector<Eg> edge_by_id(num_edges(Gsol));
    for (auto e : boost::make_iterator_range(edges(Gsol))) {
        std::size_t id = Gsol[e].id;
        if (id < edge_by_id.size())
            edge_by_id[id] = e;
    }

    int idx = 1;
    for (const auto& sol : cuts) {
        long total_cost_int = 0;  // integer cost (scaled)
        std::cout << "Cut " << idx++ << ": ";

        // Print both edge ids and (u→v)[cap] form
        if (sol.edge_ids.empty()) {
            std::cout << "{}  (cost = 0)\n";
            continue;
        }

        // Pretty list
        std::cout << "{ ";
        for (std::size_t i = 0; i < sol.edge_ids.size(); ++i) {
            auto eid = sol.edge_ids[i];
            std::cout << eid;
            if (i + 1 != sol.edge_ids.size()) std::cout << " ";
        }
        std::cout << " }   ";

        // Expanded edges and cost
        std::cout << "[ ";
        for (auto eid : sol.edge_ids) {
            if (eid < edge_by_id.size()) {
                auto e = edge_by_id[eid];
                auto uu = static_cast<std::size_t>(source(e, Gsol));
                auto vv = static_cast<std::size_t>(target(e, Gsol));
                long cap = Gsol[e].capacity;
                total_cost_int += cap;
                std::cout << "(" << uu << "->" << vv << ")[" << cap << "] ";
            }
        }
        std::cout << "]";

        // If you scaled inputs, divide to print real cost; else print integer
        if (scale != 1.0) {
            double total_cost = static_cast<double>(total_cost_int) / scale;
            std::cout << "  (cost = " << std::fixed << std::setprecision(6) << total_cost << ")\n";
        } else {
            std::cout << "  (cost = " << total_cost_int << ")\n";
        }
    }

    return 0;
}

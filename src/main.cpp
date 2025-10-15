#include "toy_example/toy_powerset_mds.hpp"
#include <iostream>

int main() {
    ToyPowersetMDS mds;

    // Unbounded (per your algorithm; final X added even if intersects top)
    auto all = mds.find_max_disjoint();
    std::cout << "Unbounded result: \n";
    for (std::size_t i=0;i<all.size();++i) {
        std::cout << "  S" << i+1 << " = {";
        for (std::size_t j=0;j<all[i].elems.size();++j) {
            std::cout << all[i].elems[j] << (j+1<all[i].elems.size() ? "," : "");
        }
        std::cout << "}\n";
    }

    // k-bounded
    auto up_to_2 = mds.find_max_disjoint(2);
    std::cout << "\nUp to k=2:\n";
    for (std::size_t i=0;i<up_to_2.size();++i) {
        std::cout << "  S" << i+1 << " = {";
        for (std::size_t j=0;j<up_to_2[i].elems.size();++j) {
            std::cout << up_to_2[i].elems[j] << (j+1<up_to_2[i].elems.size() ? "," : "");
        }
        std::cout << "}\n";
    }

    return 0; 
}

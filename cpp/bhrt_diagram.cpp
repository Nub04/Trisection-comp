// bhrt_diagram.cpp -- compute invariants of a CUSTOM trisection diagram.
//
// Reads a trisection diagram given by the homology classes of its cut
// curves in H_1(Sigma;Z) = Z^{2g}, validates it, and prints the homology,
// intersection form, signature and parity of the closed 4-manifold.
//
//   build/bhrt-diagram  <file.tri>
//   build/bhrt-diagram  --demo cp2        (built-in examples)
//
// File format (.tri) -- whitespace/'#'-comment tolerant:
//
//   g 1                 # genus
//   a 1 0               # an alpha curve, as a length-2g integer vector
//   b 0 1               # a beta curve
//   c 1 1               # a gamma curve
//
// Give exactly g rows of each of a, b, c.  Basis order is (a_1..a_g | b_1..b_g)
// with symplectic form omega((p|q),(r|s)) = p.s - q.r.

#include "bhrt_trisect.hpp"

#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

namespace {

bool parseDiagram(std::istream& in, bhrt::LagrangianDiagram& d, std::string& err) {
    int g = -1;
    std::string line;
    while (std::getline(in, line)) {
        auto h = line.find('#');
        if (h != std::string::npos) line = line.substr(0, h);
        std::istringstream ss(line);
        std::string tag;
        if (!(ss >> tag)) continue;
        if (tag == "g") {
            if (!(ss >> g) || g < 0) { err = "bad genus"; return false; }
            d.genus = g;
        } else if (tag == "a" || tag == "b" || tag == "c") {
            if (g < 0) { err = "row before genus 'g'"; return false; }
            std::vector<std::int64_t> row;
            std::int64_t x;
            while (ss >> x) row.push_back(x);
            if ((int)row.size() != 2 * g) {
                err = "row '" + tag + "' has " + std::to_string(row.size()) +
                      " entries, expected 2g=" + std::to_string(2 * g);
                return false;
            }
            if (tag == "a") d.alpha.push_back(std::move(row));
            else if (tag == "b") d.beta.push_back(std::move(row));
            else d.gamma.push_back(std::move(row));
        } else {
            err = "unknown line tag '" + tag + "'";
            return false;
        }
    }
    if (g < 0) { err = "no genus 'g' line"; return false; }
    return true;
}

bhrt::LagrangianDiagram demo(const std::string& name) {
    bhrt::LagrangianDiagram d;
    if (name == "s4") { d.genus = 0; }
    else if (name == "cp2") {
        d.genus = 1; d.alpha = {{1,0}}; d.beta = {{0,1}}; d.gamma = {{1,1}};
    } else if (name == "s2xs2") {
        d.genus = 2;
        d.alpha = {{1,0,0,0},{0,1,0,0}};
        d.beta  = {{0,0,1,0},{0,0,0,1}};
        d.gamma = {{0,1,1,0},{1,0,0,1}};
    } else if (name == "s1xs3") {
        d.genus = 1; d.alpha = {{1,0}}; d.beta = {{1,0}}; d.gamma = {{1,0}};
    }
    return d;
}

void printVec(const char* label, std::int32_t freeRank,
              const std::vector<std::int64_t>& tor) {
    std::cout << "  " << label << " = Z^" << freeRank;
    for (auto t : tor) std::cout << " (+) Z/" << t;
    std::cout << "\n";
}

int report(const bhrt::LagrangianDiagram& d) {
    auto val = bhrt::validateLagrangianDiagram(d);
    std::cout << "genus g          = " << d.genus << "\n";
    std::cout << "valid diagram    = " << (val.ok ? "yes" : "NO") << "\n";
    for (const auto& m : val.messages) std::cout << "    ! " << m << "\n";

    auto inv = bhrt::computeInvariantsFromClasses(d);
    std::cout << "homology:\n";
    std::cout << "  H_0 = Z\n";
    printVec("H_1", inv.h1_free_rank, inv.h1_torsion);
    printVec("H_2", inv.h2_free_rank, inv.h2_torsion);
    printVec("H_3", inv.h3_free_rank, inv.h3_torsion);
    std::cout << "  H_4 = Z\n";
    std::cout << "b_2              = " << inv.h2_free_rank << "\n";
    std::cout << "signature sigma  = " << inv.signature << "\n";
    std::cout << "parity           = " << inv.parity << "\n";
    if (!inv.intersection_form.empty()) {
        std::cout << "intersection form Q_X =\n";
        for (const auto& row : inv.intersection_form) {
            std::cout << "    [";
            for (std::size_t j = 0; j < row.size(); ++j)
                std::cout << (j ? " " : "") << row[j];
            std::cout << "]\n";
        }
    }
    for (const auto& a : inv.audit) std::cout << "  # " << a << "\n";
    return val.ok ? 0 : 1;
}

}  // namespace

int main(int argc, char** argv) {
    if (argc < 2) {
        std::cerr << "usage: " << argv[0] << " <file.tri>\n"
                  << "       " << argv[0] << " --demo {s4|cp2|s2xs2|s1xs3}\n";
        return 2;
    }
    bhrt::LagrangianDiagram d;
    if (std::string(argv[1]) == "--demo") {
        if (argc < 3) { std::cerr << "--demo needs a name\n"; return 2; }
        d = demo(argv[2]);
        std::cout << "# built-in demo: " << argv[2] << "\n";
    } else {
        std::ifstream fh(argv[1]);
        if (!fh) { std::cerr << "cannot open " << argv[1] << "\n"; return 1; }
        std::string err;
        if (!parseDiagram(fh, d, err)) {
            std::cerr << "parse error: " << err << "\n";
            return 1;
        }
    }
    return report(d);
}

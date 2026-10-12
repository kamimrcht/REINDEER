// Independent FASTA oracle for benchmark queries: no REINDEER/index libraries.
// g++ -O3 -march=native -fopenmp -std=c++17 test/oracle_real.cpp -o oracle-real
#include <algorithm>
#include <cmath>
#include <cctype>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>
#include <omp.h>

template<class F> void fasta(const std::string& path, F callback) {
    std::ifstream in(path);
    if (!in) throw std::runtime_error("Cannot read " + path);
    std::string line, header, sequence;
    while (std::getline(in, line)) {
        if (!line.empty() && line.back() == '\r') line.pop_back();
        if (!line.empty() && line[0] == '>') {
            if (!header.empty()) callback(header, sequence);
            header = line.substr(1);
            sequence.clear();
        } else {
            for (unsigned char base : line) sequence += std::toupper(base);
        }
    }
    if (!header.empty()) callback(header, sequence);
}

template<class F> void kmers(const std::string& sequence, F callback) {
    constexpr unsigned k = 31;
    constexpr uint64_t mask = (uint64_t{1} << (2 * k)) - 1;
    if (sequence.size() < k || sequence.find_first_not_of("ACGT") != std::string::npos) return;
    uint64_t forward = 0, reverse = 0;
    for (size_t i = 0; i < sequence.size(); ++i) {
        unsigned base = sequence[i] == 'A' ? 0 : sequence[i] == 'C' ? 1 : sequence[i] == 'G' ? 2 : 3;
        forward = ((forward << 2) | base) & mask;
        reverse = (reverse >> 2) | (uint64_t(base ^ 3) << (2 * (k - 1)));
        if (i + 1 >= k) callback(std::min(forward, reverse));
    }
}

struct Query { std::string name; std::vector<uint32_t> positions; };

int main(int argc, char** argv) {
    if (argc != 5) {
        std::cerr << "Usage: oracle-real query.fa datasets.fof output-directory threads\n";
        return 1;
    }
    std::unordered_map<uint64_t, uint32_t> targets;
    std::vector<Query> queries;
    fasta(argv[1], [&](const auto& header, const auto& sequence) {
        Query query{header.substr(0, header.find_first_of(" \t")), {}};
        kmers(sequence, [&](uint64_t key) {
            auto [it, inserted] = targets.emplace(key, targets.size());
            query.positions.push_back(it->second);
        });
        if (!query.positions.empty()) queries.push_back(std::move(query));
    });
    const auto& lookup = targets;
    std::ifstream fof(argv[2]);
    if (!fof) throw std::runtime_error("Cannot read file list");
    std::vector<std::string> files;
    std::string line;
    while (std::getline(fof, line)) if (!line.empty()) files.push_back(line);
    std::filesystem::create_directories(argv[3]);
    std::vector<std::string> errors(files.size());
    #pragma omp parallel for num_threads(std::stoi(argv[4])) schedule(dynamic)
    for (size_t color = 0; color < files.size(); ++color) {
        try {
            auto output_path = std::filesystem::path(argv[3]) / (std::to_string(color) + ".tsv");
            if (std::filesystem::exists(output_path)) continue;
            std::vector<uint16_t> values(targets.size(), 0);
            uint64_t checked = 0;
            fasta(files[color], [&](const std::string& header, const auto& sequence) {
                double abundance = 1;
                for (const auto* tag : {"km:f:", "KM:f:", "ka:f:", "KA:f:"}) {
                    auto position = header.find(tag);
                    if (position != std::string::npos) {
                        abundance = std::stod(header.substr(position + 5));
                        break;
                    }
                }
                if (!std::isfinite(abundance) || abundance < 0) throw std::runtime_error("Invalid abundance");
                auto value = static_cast<uint16_t>(std::min(abundance, 65535.0));
                kmers(sequence, [&](uint64_t key) {
                    ++checked;
                    auto found = lookup.find(key);
                    if (found != lookup.end()) values[found->second] = std::max(values[found->second], value);
                });
            });
            auto partial = output_path.string() + ".partial";
            std::ofstream out(partial);
            for (const auto& query : queries) {
                out << query.name << '\t';
                size_t begin = 0;
                while (begin < query.positions.size()) {
                    auto value = values[query.positions[begin]];
                    size_t end = begin + 1;
                    while (end < query.positions.size() && values[query.positions[end]] == value) ++end;
                    if (begin) out << ',';
                    out << begin << '-' << end - 1 << ':';
                    if (value) out << value; else out << '*';
                    begin = end;
                }
                out << '\n';
            }
            out.close();
            if (!out) throw std::runtime_error("Cannot write oracle result");
            std::filesystem::rename(partial, output_path);
            #pragma omp critical
            std::cout << "Dataset " << color << ": " << checked << " k-mer occurrences checked\n" << std::flush;
        } catch (const std::exception& error) {
            errors[color] = error.what();
        }
    }
    for (const auto& error : errors) if (!error.empty()) { std::cerr << error << '\n'; return 1; }
    return 0;
}

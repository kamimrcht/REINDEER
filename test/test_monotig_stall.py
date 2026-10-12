#!/usr/bin/env python3
"""Replay an oversized minimizer without a 100-million-k-mer fixture.

Lower the former cutoff in a temporary compilation, if it is present. This
fixture also makes an early partition succeed before a later one overflows,
so both termination and output/counter integrity are checked.
Run after make; --source can reproduce the failure with the old source.
"""
import argparse
import re
import subprocess
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
HARNESS = r'''
#include "monotig.hpp"
#include <stdexcept>
#define check(condition) do { if (!(condition)) throw std::runtime_error("Monotig replay mismatch at line " + std::to_string(__LINE__)); } while (false)

int main(int argc, char** argv) {
    const std::string input = std::string(argv[1]) + "/input.gz";
    std::string output = std::string(argv[1]) + "/output.gz";
    const uint64_t passes = std::stoull(argv[2]);
    kmer_Set_Light index(31, 10, 10, 5, 1, 0);
    index.color_mode = 1;
    index.nuc_minimizer = new uint32_t[index.minimizer_number.value()]();
    for (uint64_t i=0; i<index.mphf_number.value(); ++i)
        index.all_mphf[i].kmer_MPHF = nullptr;
    std::map<kmer, std::vector<uint16_t>> expected;
    auto visit = [&](const std::string& sequence, auto callback) {
        kmer forward = str2num(sequence.substr(0,31));
        kmer reverse = index.rcb(forward);
        callback(std::min(forward, reverse));
        for (size_t i=31; i<sequence.size(); ++i) {
            index.updateK(forward,sequence[i]);
            index.updateRCK(reverse,sequence[i]);
            callback(std::min(forward, reverse));
        }
    };
    {
        zstr::ofstream stream(input);
        auto record = [&](uint32_t minimizer, uint32_t color, uint16_t coverage, const std::string& sequence) {
            uint8_t length=sequence.size();
            stream.write(reinterpret_cast<char*>(&minimizer),4);
            stream.write(reinterpret_cast<char*>(&color),4);
            stream.write(reinterpret_cast<char*>(&coverage),2);
            stream.write(reinterpret_cast<char*>(&length),1);
            stream.write(sequence.data(),length);
            visit(sequence,[&](kmer key) {
                auto& row=expected[key]; row.resize(3);
                row[color]=std::max(row[color],coverage);
            });
        };
        uint64_t state=427;
        auto sequence = [&]() {
            std::string dna;
            for (int j=0; j<31; ++j) {
                state ^= state << 13; state ^= state >> 7; state ^= state << 17;
                dna += "ACGT"[state & 3];
            }
            return dna;
        };
        for (uint32_t mini : {0U,1U})
            for (int i=0; i<10; ++i) record(mini,0,4,sequence());
        // All these distinct k-mers share one minimizer ID. Splitting that
        // ID by modulus cannot make its partition smaller.
        for (int i=0; i<2050; ++i) {
            std::string dna=sequence();
            record(254,0,7,dna); record(254,0,2,dna); record(254,2,19,dna);
        }
        record(0,1,12,std::string(100,'A'));
    }
    index.merge_super_buckets_mem(input,3,output,passes,1);
    uint64_t dna=0, kmers=0, strings=0;
    {
        zstr::ifstream stream(output);
        while (stream.peek()!=EOF) {
            check(stream.get()=='>');
            uint32_t minimizer, size;
            stream.read(reinterpret_cast<char*>(&minimizer),4);
            stream.read(reinterpret_cast<char*>(&size),4);
            std::vector<unsigned char> compressed(size), decoded(4096);
            stream.read(reinterpret_cast<char*>(compressed.data()),size);
            check(trled(compressed.data(),size,decoded.data(),6)==size);
            std::vector<uint16_t> counts(3);
            std::memcpy(counts.data(),decoded.data(),6);
            check(stream.get()=='\n');
            std::string sequence; std::getline(stream,sequence);
            check(sequence.size()>=31);
            dna+=sequence.size(); kmers+=sequence.size()-30; ++strings;
            visit(sequence,[&](kmer key) {
                auto found=expected.find(key);
                check(found!=expected.end() && found->second==counts);
                expected.erase(found);
            });
        }
    }
    check(expected.empty());
    uint64_t recorded_dna=0, recorded_kmers=0, recorded_strings=0;
    for (uint64_t i=0; i<index.minimizer_number.value(); ++i) {
        recorded_dna+=index.nuc_minimizer[i];
        recorded_kmers+=index.all_mphf[i].mphf_size;
        recorded_strings+=index.all_buckets[i].skmer_number;
    }
    check(dna==recorded_dna && dna==index.total_nuc_number);
    check(kmers==recorded_kmers && strings==recorded_strings);
    delete[] index.nuc_minimizer;
    std::cout << "Replay passed: " << kmers << " k-mers; passes=" << passes << std::endl;
}
'''


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, default=ROOT / "src/monotig.cpp")
    args = parser.parse_args()
    objects = ["blight", "utils_b", "trlec", "trled", "utils", "reindeer", "query",
               "build_index", "eq_classes", "launch_bcalm", "matrix_operation"]
    with tempfile.TemporaryDirectory(prefix="reindeer-monotig-stall-") as directory:
        tmp = Path(directory)
        source = tmp / "monotig.cpp"
        source.write_text(re.sub(r"(?<!\d)100000000(?!\d)", "1024", args.source.read_text()))
        (tmp / "replay.cpp").write_text(HARNESS)
        command = ["g++", "-O2", "-g", "-flto", "-fopenmp", "-msse4", "-std=c++17",
                   "-Wno-deprecated-declarations", "-I" + str(ROOT / "src"),
                   str(tmp / "replay.cpp"), str(source)]
        command += [str(ROOT / (name + ".o")) for name in objects]
        command += [str(ROOT / "blight/lz4" / (name + ".o"))
                    for name in ["lz4frame", "lz4", "lz4hc", "xxhash"]]
        subprocess.run([*command, "-lz", "-o", str(tmp / "replay")], check=True, timeout=120)
        for passes in [1, 2, 4]:
            subprocess.run([str(tmp / "replay"), str(tmp), str(passes)], check=True, timeout=20)


if __name__ == "__main__":
    main()

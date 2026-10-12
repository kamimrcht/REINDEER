#include "monotig.hpp"
#include <malloc.h>
#include <unistd.h>
// Linux/glibc diagnostic replay for a binary super-bucket with up to 1,024
// datasets. Link with the normal REINDEER objects except main.o.
// Usage: replay <input.gz> <existing output directory> <workers> <trim after>
// Repeat the same real input per worker to compare arena retention and
// compressed output hashes before/after a change to merge_super_buckets_mem.
static void memory(const char* label) {
    std::ifstream in("/proc/self/statm"); uint64_t pages, rss;
    in >> pages >> rss;
    const auto m = mallinfo2();
    std::cerr << label << " rss_bytes=" << rss*sysconf(_SC_PAGESIZE)
              << " allocator_live_bytes=" << m.uordblks
              << " allocator_free_bytes=" << m.fordblks << std::endl;
}
int main(int argc, char** argv) {
    if (argc != 5) return 2;
    const int threads = std::stoi(argv[3]);
    kmer_Set_Light index(31,10,10,5,threads,0);
    index.color_mode=1;
    index.nuc_minimizer = new uint32_t[index.minimizer_number.value()]();
    for (uint64_t i=0;i<index.mphf_number.value();++i)
        index.all_mphf[i].kmer_MPHF=nullptr;
    memory("BEFORE");
    #pragma omp parallel for num_threads(threads)
    for (int i=0;i<threads;++i) {
        std::string out=std::string(argv[2])+"/out-"+std::to_string(i)+".gz";
        index.merge_super_buckets_mem(argv[1],1024,out,1,1);
    }
    memory("AFTER_MERGE");
    if (std::stoi(argv[4])) {
        malloc_trim(0);
        memory("AFTER_TRIM");
    }
    delete[] index.nuc_minimizer;
}

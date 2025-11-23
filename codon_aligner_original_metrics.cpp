// Instrumented original aligner for evaluation
#include <iostream>
#include <fstream>
#include <string>
#include <unordered_map>
#include <vector>
#include <chrono>

using namespace std;
using namespace chrono;

struct AlignmentMetrics {
    int frameshift_attempts = 0;
    int frameshift_success = 0;
    int exact_matches = 0;
    int synonymous_snps = 0;
    int nonsynonymous_snps = 0;
};

AlignmentMetrics global_metrics;
vector<int> per_read_penalties;

char translate_codon(const string& codon) {
    static const unordered_map<string, char> codon_table = {
        {"TTT",'F'},{"TTC",'F'},{"TTA",'L'},{"TTG",'L'},{"CTT",'L'},{"CTC",'L'},{"CTA",'L'},{"CTG",'L'},
        {"ATT",'I'},{"ATC",'I'},{"ATA",'I'},{"ATG",'M'},{"GTT",'V'},{"GTC",'V'},{"GTA",'V'},{"GTG",'V'},
        {"TCT",'S'},{"TCC",'S'},{"TCA",'S'},{"TCG",'S'},{"CCT",'P'},{"CCC",'P'},{"CCA",'P'},{"CCG",'P'},
        {"ACT",'T'},{"ACC",'T'},{"ACA",'T'},{"ACG",'T'},{"GCT",'A'},{"GCC",'A'},{"GCA",'A'},{"GCG",'A'},
        {"TAT",'Y'},{"TAC",'Y'},{"TAA",'*'},{"TAG",'*'},{"CAT",'H'},{"CAC",'H'},{"CAA",'Q'},{"CAG",'Q'},
        {"AAT",'N'},{"AAC",'N'},{"AAA",'K'},{"AAG",'K'},{"GAT",'D'},{"GAC",'D'},{"GAA",'E'},{"GAG",'E'},
        {"TGT",'C'},{"TGC",'C'},{"TGA",'*'},{"TGG",'W'},{"CGT",'R'},{"CGC",'R'},{"CGA",'R'},{"CGG",'R'},
        {"AGT",'S'},{"AGC",'S'},{"AGA",'R'},{"AGG",'R'},{"GGT",'G'},{"GGC",'G'},{"GGA",'G'},{"GGG",'G'}
    };
    auto it = codon_table.find(codon);
    return (it != codon_table.end()) ? it->second : '?';
}

int score_codons(const string& ref_codon, const string& query_codon) {
    if (ref_codon == query_codon) {
        global_metrics.exact_matches++;
        return 0;
    }
    if (translate_codon(ref_codon) == translate_codon(query_codon)) {
        global_metrics.synonymous_snps++;
        return 1;
    }
    global_metrics.nonsynonymous_snps++;
    return 3;
}

int align_sequences(const string& ref, const string& query) {
    size_t i = 0, j = 0;
    int total_penalty = 0;
    while (i + 2 < ref.size() && j + 2 < query.size()) {
        string rc = ref.substr(i, 3);
        string qc = query.substr(j, 3);
        if (rc == qc) {
            i += 3; j += 3;
            continue;
        }
        int penalty = score_codons(rc, qc);
        bool matched = false;
        
        global_metrics.frameshift_attempts++;
        
        for (int slide = 1; slide <= 6 && !matched; ++slide) {
            for (auto [di, dj] : {make_pair(-slide, 0), make_pair(slide, 0), make_pair(0, -slide), make_pair(0, slide)}) {
                int ni = i + 3 + di;
                int nj = j + 3 + dj;
                if (ni >= 0 && nj >= 0 && ni + 2 < ref.size() && nj + 2 < query.size()) {
                    if (ref.substr(ni, 3) == query.substr(nj, 3)) {
                        global_metrics.frameshift_success++;
                        total_penalty += 5;
                        i = ni;
                        j = nj;
                        matched = true;
                        break;
                    }
                }
            }
        }
        if (!matched) {
            total_penalty += penalty;
            i += 3; j += 3;
        }
    }
    return total_penalty;
}

vector<string> load_reads_fastq(const string& path, vector<string>& ids) {
    ifstream f(path);
    string line;
    vector<string> reads;
    int line_count = 0;
    
    while (getline(f, line)) {
        line_count++;
        if (line_count % 4 == 1) {
            ids.push_back(line.substr(1));
        }
        if (line_count % 4 == 2) {
            reads.push_back(line);
        }
    }
    return reads;
}

vector<string> load_reads_fasta(const string& path, vector<string>& ids) {
    ifstream f(path);
    string line, seq;
    vector<string> reads;
    
    while (getline(f, line)) {
        if (line.empty() || line[0] == '>') {
            if (!seq.empty()) {
                reads.push_back(seq);
                seq.clear();
            }
            if (line[0] == '>') {
                ids.push_back(line.substr(1));
            }
        } else {
            seq += line;
        }
    }
    if (!seq.empty()) reads.push_back(seq);
    return reads;
}

int main(int argc, char* argv[]) {
    if (argc != 3) {
        cerr << "Usage: ./codon_aligner_original_metrics <reference.fasta> <reads.fasta|fastq>\n";
        return 1;
    }
    
    auto start = high_resolution_clock::now();
    
    ifstream ref_file(argv[1]);
    string ref, line;
    while (getline(ref_file, line)) {
        if (line[0] != '>') ref += line;
    }
    
    vector<string> ids;
    ifstream test(argv[2]);
    char first_char = test.peek();
    test.close();
    
    vector<string> reads;
    if (first_char == '@') {
        reads = load_reads_fastq(argv[2], ids);
    } else {
        reads = load_reads_fasta(argv[2], ids);
    }
    
    int total_penalty = 0;
    
    ofstream detail_file("original_detailed_results.txt");
    detail_file << "ReadID\tPenalty\n";
    
    for (size_t idx = 0; idx < reads.size(); ++idx) {
        int penalty = align_sequences(ref, reads[idx]);
        detail_file << ids[idx] << "\t" << penalty << "\n";
        total_penalty += penalty;
        per_read_penalties.push_back(penalty);
    }
    
    detail_file.close();
    
    auto end = high_resolution_clock::now();
    double time_sec = duration_cast<duration<double>>(end - start).count();
    
    cout << "\n=== ORIGINAL ALIGNER METRICS ===\n";
    cout << "Total reads: " << reads.size() << "\n";
    cout << "Time: " << time_sec << " seconds\n";
    cout << "Speed: " << (reads.size() / time_sec) << " reads/sec\n";
    cout << "Total penalty: " << total_penalty << "\n";
    cout << "Avg penalty/read: " << (total_penalty / reads.size()) << "\n\n";
    
    cout << "=== FRAMESHIFT METRICS ===\n";
    cout << "Frameshift attempts: " << global_metrics.frameshift_attempts << "\n";
    cout << "Frameshift successes: " << global_metrics.frameshift_success << "\n";
    cout << "Exact matches: " << global_metrics.exact_matches << "\n";
    cout << "Synonymous SNPs: " << global_metrics.synonymous_snps << "\n";
    cout << "Non-synonymous SNPs: " << global_metrics.nonsynonymous_snps << "\n";
    
    return 0;
}

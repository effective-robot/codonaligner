// codon_aligner_batch.cpp
// Production-scale codon-level aligner (optimized for speed & performance)
// Batch aligns a reference against FASTA reads (stdin or file)

#include <iostream>
#include <fstream>
#include <string>
#include <unordered_map>
#include <vector>
#include <chrono>

using namespace std;
using namespace chrono;

// --- Translate codon to amino acid ---
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

// --- Score codons ---
int score_codons(const string& ref_codon, const string& query_codon) {
    if (ref_codon == query_codon) return 0;
    if (translate_codon(ref_codon) == translate_codon(query_codon)) return 1;
    return 3;
}

// --- Core alignment logic ---
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
        for (int slide = 1; slide <= 6 && !matched; ++slide) {
            for (auto [di, dj] : {make_pair(-slide, 0), make_pair(slide, 0), make_pair(0, -slide), make_pair(0, slide)}) {
                int ni = i + 3 + di;
                int nj = j + 3 + dj;
                if (ni >= 0 && nj >= 0 && ni + 2 < ref.size() && nj + 2 < query.size()) {
                    if (ref.substr(ni, 3) == query.substr(nj, 3)) {
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

// --- Load reads from FASTA ---
vector<string> load_reads(const string& path) {
    ifstream f(path);
    string line, seq;
    vector<string> reads;
    while (getline(f, line)) {
        if (line.empty() || line[0] == '>') {
            if (!seq.empty()) { reads.push_back(seq); seq.clear(); }
        } else {
            seq += line;
        }
    }
    if (!seq.empty()) reads.push_back(seq);
    return reads;
}

int main(int argc, char* argv[]) {
    if (argc != 3) {
        cerr << "Usage: ./codon_aligner_batch <reference.fasta> <reads.fasta>\n";
        return 1;
    }

    auto start = high_resolution_clock::now();

    // Load reference
    ifstream ref_file(argv[1]);
    string ref, line;
    while (getline(ref_file, line)) {
        if (line[0] != '>') ref += line;
    }

    // Load reads
    vector<string> reads = load_reads(argv[2]);

    // Process alignments
    int total_penalty = 0;
    for (size_t i = 0; i < reads.size(); ++i) {
        total_penalty += align_sequences(ref, reads[i]);
        if ((i+1) % 100000 == 0) {
            cout << "Processed " << (i+1) << " reads...\n";
        }
    }

    auto end = high_resolution_clock::now();
    double time_sec = duration_cast<duration<double>>(end - start).count();

    cout << "\n✅ Total reads aligned: " << reads.size() << "\n";
    cout << "🕒 Total time: " << time_sec << " seconds\n";
    cout << "📉 Total penalty score: " << total_penalty << "\n";
    cout << "⚡ Speed: " << (reads.size() / time_sec) << " reads/sec\n";
    return 0;
}

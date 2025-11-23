#!/bin/bash
# Three-way aligner comparison on ORF test data

echo "═══════════════════════════════════════════════════════════"
echo "THREE-WAY ALIGNER COMPARISON"
echo "═══════════════════════════════════════════════════════════"
echo ""

# Compile all aligners
echo "Compiling aligners..."
g++ -O3 -std=c++17 -o codon_aligner codon_aligner.cpp 2>&1 | grep -v "warning:" || true
g++ -O3 -std=c++17 -o codon_aligner_ont codon_aligner_ont.cpp 2>&1 | grep -v "warning:" || true
g++ -O3 -std=c++17 -o codon_aligner_orf codon_aligner_orf.cpp 2>&1 | grep -v "warning:" || true

echo "Compilation complete."
echo ""

# Create a simple version of original for ORF data
cat > codon_aligner_original_orf_test.cpp << 'CPP'
#include <iostream>
#include <fstream>
#include <string>
#include <unordered_map>
#include <vector>
#include <chrono>
#include <sstream>

using namespace std;
using namespace chrono;

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
    if (ref_codon == query_codon) return 0;
    if (translate_codon(ref_codon) == translate_codon(query_codon)) return 1;
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

vector<pair<string, string>> load_transcripts_fasta(const string& path) {
    ifstream f(path);
    string line, seq, id;
    vector<pair<string, string>> transcripts;
    while (getline(f, line)) {
        if (line.empty()) continue;
        if (line[0] == '>') {
            if (!seq.empty() && !id.empty()) {
                transcripts.push_back({id, seq});
                seq.clear();
            }
            istringstream iss(line.substr(1));
            iss >> id;
        } else {
            seq += line;
        }
    }
    if (!seq.empty() && !id.empty()) {
        transcripts.push_back({id, seq});
    }
    return transcripts;
}

vector<pair<string, string>> load_reads_fastq(const string& path) {
    ifstream f(path);
    string line;
    vector<pair<string, string>> reads;
    int line_count = 0;
    string read_id, read_seq;
    while (getline(f, line)) {
        line_count++;
        int pos = line_count % 4;
        if (pos == 1) {
            read_id = line.substr(1);
        } else if (pos == 2) {
            read_seq = line;
            reads.push_back({read_id, read_seq});
        }
    }
    return reads;
}

int main(int argc, char* argv[]) {
    if (argc != 3) {
        cerr << "Usage: " << argv[0] << " <reference.fasta> <reads.fastq>\n";
        return 1;
    }
    
    auto start_time = high_resolution_clock::now();
    
    auto transcripts = load_transcripts_fasta(argv[1]);
    unordered_map<string, string> transcript_seqs;
    for (const auto& [id, seq] : transcripts) {
        transcript_seqs[id] = seq;
    }
    
    auto reads = load_reads_fastq(argv[2]);
    
    int total_penalty = 0;
    int successful = 0;
    unordered_map<string, vector<int>> penalties_by_transcript;
    ofstream out("original_orf_results.txt");
    out << "ReadID\tTranscriptID\tPenalty\n";
    
    for (const auto& [read_id, read_seq] : reads) {
        string trans_id = read_id.substr(0, read_id.find('_', read_id.find('_') + 1));
        auto trans_it = transcript_seqs.find(trans_id);
        if (trans_it == transcript_seqs.end()) {
            continue;
        }
        
        const string& ref_seq = trans_it->second;
        int penalty = align_sequences(ref_seq, read_seq);
        
        out << read_id << "\t" << trans_id << "\t" << penalty << "\n";
        
        if (penalty < 9999) {
            total_penalty += penalty;
            successful++;
            penalties_by_transcript[trans_id].push_back(penalty);
        }
    }
    out.close();
    
    auto end_time = high_resolution_clock::now();
    double elapsed = duration_cast<duration<double>>(end_time - start_time).count();
    
    cout << "=== ORIGINAL ALIGNER (on ORF test data) ===\n";
    cout << "Total reads: " << reads.size() << "\n";
    cout << "Successful: " << successful << " (" << (100.0 * successful / reads.size()) << "%)\n";
    cout << "Total penalty: " << total_penalty << "\n";
    cout << "Avg penalty/read: " << (successful > 0 ? total_penalty / successful : 0) << "\n";
    cout << "Runtime: " << elapsed << " seconds\n";
    cout << "Speed: " << (reads.size() / elapsed) << " reads/sec\n";
    
    return 0;
}
CPP

g++ -O3 -std=c++17 -o codon_aligner_original_orf_test codon_aligner_original_orf_test.cpp
echo "Custom original aligner compiled."
echo ""

# Run Original Aligner
echo "══════════════════════════════════════"
echo "1. ORIGINAL ALIGNER"
echo "══════════════════════════════════════"
./codon_aligner_original_orf_test reference_transcripts_orf.fasta test_reads_orf.fastq
echo ""

# Run ONT-Optimized Aligner
echo "══════════════════════════════════════"
echo "2. ONT-OPTIMIZED ALIGNER"
echo "══════════════════════════════════════"

cat > codon_aligner_ont_orf_test.cpp << 'CPP2'
#include <iostream>
#include <fstream>
#include <string>
#include <unordered_map>
#include <vector>
#include <chrono>
#include <cmath>
#include <algorithm>
#include <sstream>

using namespace std;
using namespace chrono;

struct ONTErrorModel {
    static double get_error_rate(int position, int read_length) {
        double norm_pos = static_cast<double>(position) / max(read_length, 1);
        double base_rate = 0.075;
        double dist_from_center = abs(norm_pos - 0.5) * 2.0;
        double position_penalty = 0.15 * (dist_from_center * dist_from_center);
        return min(0.25, max(0.05, base_rate + position_penalty));
    }
};

bool is_homopolymer(const string& seq, size_t pos, int min_len = 3) {
    if (pos >= seq.size()) return false;
    char base = seq[pos];
    int count = 1;
    for (size_t i = pos + 1; i < seq.size() && seq[i] == base; ++i) {
        count++;
        if (count >= min_len) return true;
    }
    for (int i = static_cast<int>(pos) - 1; i >= 0 && seq[i] == base; --i) {
        count++;
        if (count >= min_len) return true;
    }
    return count >= min_len;
}

int homopolymer_length(const string& seq, size_t pos) {
    if (pos >= seq.size()) return 0;
    char base = seq[pos];
    int count = 1;
    for (size_t i = pos + 1; i < seq.size() && seq[i] == base; ++i) count++;
    for (int i = static_cast<int>(pos) - 1; i >= 0 && seq[i] == base; --i) count++;
    return count;
}

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

int score_codons_ont(const string& ref_codon, const string& query_codon, int position, int read_length) {
    if (ref_codon == query_codon) return 0;
    int base_penalty;
    if (translate_codon(ref_codon) == translate_codon(query_codon)) {
        base_penalty = 1;
    } else {
        base_penalty = 5;
    }
    double error_rate = ONTErrorModel::get_error_rate(position, read_length);
    double position_weight = 1.0 - error_rate;
    double weighted_penalty = base_penalty * max(0.5, position_weight);
    return static_cast<int>(weighted_penalty + 0.5);
}

int align_sequences_ont(const string& ref, const string& query) {
    size_t i = 0, j = 0;
    int total_penalty = 0;
    int query_length = query.size();
    
    while (i + 2 < ref.size() && j + 2 < query.size()) {
        string rc = ref.substr(i, 3);
        string qc = query.substr(j, 3);
        
        if (rc == qc) {
            i += 3; j += 3;
            continue;
        }
        
        bool ref_hp = is_homopolymer(ref, i, 3);
        bool query_hp = is_homopolymer(query, j, 3);
        
        if (ref_hp || query_hp) {
            int hp_len = max(
                ref_hp ? homopolymer_length(ref, i) : 0,
                query_hp ? homopolymer_length(query, j) : 0
            );
            total_penalty += 2;
            int skip_distance = max(3, (hp_len / 3) * 3);
            i += skip_distance;
            j += skip_distance;
            continue;
        }
        
        bool matched = false;
        int best_di = 0, best_dj = 0;
        int min_distance = 999;
        
        for (int di = -2; di <= 2; ++di) {
            for (int dj = -2; dj <= 2; ++dj) {
                if (di == 0 && dj == 0) continue;
                
                int ni = static_cast<int>(i) + 3 + di;
                int nj = static_cast<int>(j) + 3 + dj;
                
                if (ni < 0 || nj < 0 ||
                    static_cast<size_t>(ni + 2) >= ref.size() ||
                    static_cast<size_t>(nj + 2) >= query.size()) {
                    continue;
                }
                
                if (ref.substr(ni, 3) == query.substr(nj, 3)) {
                    int distance = abs(di) + abs(dj);
                    if (distance < min_distance) {
                        min_distance = distance;
                        best_di = di;
                        best_dj = dj;
                        matched = true;
                    }
                }
            }
        }
        
        if (matched) {
            total_penalty += (min_distance + 3);
            i = static_cast<size_t>(static_cast<int>(i) + 3 + best_di);
            j = static_cast<size_t>(static_cast<int>(j) + 3 + best_dj);
        } else {
            int penalty = score_codons_ont(rc, qc, j, query_length);
            total_penalty += penalty;
            i += 3;
            j += 3;
        }
        
        int codons_processed = max(static_cast<int>(i / 3), static_cast<int>(j / 3));
        if (codons_processed > 10) {
            int max_expected_penalty = codons_processed * 2;
            if (total_penalty > max_expected_penalty * 2.5) {
                return 9999;
            }
        }
    }
    
    return total_penalty;
}

vector<pair<string, string>> load_transcripts_fasta(const string& path) {
    ifstream f(path);
    string line, seq, id;
    vector<pair<string, string>> transcripts;
    while (getline(f, line)) {
        if (line.empty()) continue;
        if (line[0] == '>') {
            if (!seq.empty() && !id.empty()) {
                transcripts.push_back({id, seq});
                seq.clear();
            }
            istringstream iss(line.substr(1));
            iss >> id;
        } else {
            seq += line;
        }
    }
    if (!seq.empty() && !id.empty()) {
        transcripts.push_back({id, seq});
    }
    return transcripts;
}

vector<pair<string, string>> load_reads_fastq(const string& path) {
    ifstream f(path);
    string line;
    vector<pair<string, string>> reads;
    int line_count = 0;
    string read_id, read_seq;
    while (getline(f, line)) {
        line_count++;
        int pos = line_count % 4;
        if (pos == 1) {
            read_id = line.substr(1);
        } else if (pos == 2) {
            read_seq = line;
            reads.push_back({read_id, read_seq});
        }
    }
    return reads;
}

int main(int argc, char* argv[]) {
    if (argc != 3) {
        cerr << "Usage: " << argv[0] << " <reference.fasta> <reads.fastq>\n";
        return 1;
    }
    
    auto start_time = high_resolution_clock::now();
    
    auto transcripts = load_transcripts_fasta(argv[1]);
    unordered_map<string, string> transcript_seqs;
    for (const auto& [id, seq] : transcripts) {
        transcript_seqs[id] = seq;
    }
    
    auto reads = load_reads_fastq(argv[2]);
    
    int total_penalty = 0;
    int successful = 0;
    unordered_map<string, vector<int>> penalties_by_transcript;
    ofstream out("ont_orf_results.txt");
    out << "ReadID\tTranscriptID\tPenalty\n";
    
    for (const auto& [read_id, read_seq] : reads) {
        string trans_id = read_id.substr(0, read_id.find('_', read_id.find('_') + 1));
        auto trans_it = transcript_seqs.find(trans_id);
        if (trans_it == transcript_seqs.end()) {
            continue;
        }
        
        const string& ref_seq = trans_it->second;
        int penalty = align_sequences_ont(ref_seq, read_seq);
        
        out << read_id << "\t" << trans_id << "\t" << penalty << "\n";
        
        if (penalty < 9999) {
            total_penalty += penalty;
            successful++;
            penalties_by_transcript[trans_id].push_back(penalty);
        }
    }
    out.close();
    
    auto end_time = high_resolution_clock::now();
    double elapsed = duration_cast<duration<double>>(end_time - start_time).count();
    
    cout << "=== ONT-OPTIMIZED ALIGNER (on ORF test data) ===\n";
    cout << "Total reads: " << reads.size() << "\n";
    cout << "Successful: " << successful << " (" << (100.0 * successful / reads.size()) << "%)\n";
    cout << "Total penalty: " << total_penalty << "\n";
    cout << "Avg penalty/read: " << (successful > 0 ? total_penalty / successful : 0) << "\n";
    cout << "Runtime: " << elapsed << " seconds\n";
    cout << "Speed: " << (reads.size() / elapsed) << " reads/sec\n";
    
    return 0;
}
CPP2

g++ -O3 -std=c++17 -o codon_aligner_ont_orf_test codon_aligner_ont_orf_test.cpp
./codon_aligner_ont_orf_test reference_transcripts_orf.fasta test_reads_orf.fastq
echo ""

# Run ORF-Aware Aligner
echo "══════════════════════════════════════"
echo "3. ORF-AWARE ALIGNER"
echo "══════════════════════════════════════"

# Modify ORF aligner to output per-read results
cat > codon_aligner_orf_detailed.cpp << 'CPP3'
// Same as codon_aligner_orf.cpp but with detailed per-read output
// (Include full ORF aligner code with added output file)
CPP3

# Just use the existing ORF aligner and capture output
./codon_aligner_orf --orf-db orf_database_test.txt reference_transcripts_orf.fasta test_reads_orf.fastq > orf_aligner_output.txt

# Extract results for comparison
grep "Total penalty:" orf_aligner_output.txt
grep "Avg penalty/read:" orf_aligner_output.txt
grep "Successful alignments:" orf_aligner_output.txt
grep "Runtime:" orf_aligner_output.txt
grep "Speed:" orf_aligner_output.txt

echo ""
echo "══════════════════════════════════════"
echo "COMPARISON COMPLETE"
echo "══════════════════════════════════════"
echo "Results saved to:"
echo "  - original_orf_results.txt"
echo "  - ont_orf_results.txt"
echo "  - orf_aligner_output.txt"

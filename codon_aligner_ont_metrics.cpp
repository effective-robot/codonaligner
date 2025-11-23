// Instrumented version for evaluation - collects detailed metrics
#include <iostream>
#include <fstream>
#include <string>
#include <unordered_map>
#include <vector>
#include <chrono>
#include <cmath>
#include <algorithm>
#include <iomanip>

using namespace std;
using namespace chrono;

// Metrics tracking
struct AlignmentMetrics {
    int homopolymer_fastpath = 0;
    int frameshift_attempts = 0;
    int frameshift_success = 0;
    int early_terminations = 0;
    int position_adjustments = 0;
    int exact_matches = 0;
    int synonymous_snps = 0;
    int nonsynonymous_snps = 0;
};

AlignmentMetrics global_metrics;
vector<int> per_read_penalties;
vector<string> read_ids;

struct ONTErrorModel {
    static double get_error_rate(int position, int read_length) {
        double norm_pos = static_cast<double>(position) / max(read_length, 1);
        double base_rate = 0.075;
        double dist_from_center = abs(norm_pos - 0.5) * 2.0;
        double position_penalty = 0.15 * (dist_from_center * dist_from_center);
        return min(0.25, max(0.05, base_rate + position_penalty));
    }
    
    static double homopolymer_indel_rate(int hp_length) {
        if (hp_length < 3) return 1.0;
        if (hp_length <= 5) return 1.5;
        return 2.5;
    }
    
    static double context_multiplier(const string& seq_window) {
        if (seq_window.empty()) return 1.0;
        int at_count = 0;
        for (char c : seq_window) {
            if (c == 'A' || c == 'T') at_count++;
        }
        double at_fraction = static_cast<double>(at_count) / seq_window.size();
        if (at_fraction > 0.65) return 1.2;
        if (at_fraction < 0.35) return 1.15;
        return 1.0;
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
    for (int i = pos - 1; i >= 0 && seq[i] == base; --i) {
        count++;
        if (count >= min_len) return true;
    }
    return count >= min_len;
}

int homopolymer_length(const string& seq, size_t pos) {
    if (pos >= seq.size()) return 0;
    char base = seq[pos];
    int count = 1;
    for (size_t i = pos + 1; i < seq.size() && seq[i] == base; ++i) {
        count++;
    }
    for (int i = pos - 1; i >= 0 && seq[i] == base; --i) {
        count++;
    }
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
    if (ref_codon == query_codon) {
        global_metrics.exact_matches++;
        return 0;
    }
    
    int base_penalty;
    if (translate_codon(ref_codon) == translate_codon(query_codon)) {
        base_penalty = 1;
        global_metrics.synonymous_snps++;
    } else {
        base_penalty = 5;
        global_metrics.nonsynonymous_snps++;
    }
    
    double error_rate = ONTErrorModel::get_error_rate(position, read_length);
    double position_weight = 1.0 - error_rate;
    
    if (position_weight < 0.95) {
        global_metrics.position_adjustments++;
    }
    
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
            i += 3;
            j += 3;
            continue;
        }
        
        bool ref_hp = is_homopolymer(ref, i, 3);
        bool query_hp = is_homopolymer(query, j, 3);
        
        if (ref_hp || query_hp) {
            global_metrics.homopolymer_fastpath++;
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
        
        global_metrics.frameshift_attempts++;
        
        for (int di = -2; di <= 2; ++di) {
            for (int dj = -2; dj <= 2; ++dj) {
                if (di == 0 && dj == 0) continue;
                
                int ni = i + 3 + di;
                int nj = j + 3 + dj;
                
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
            global_metrics.frameshift_success++;
            total_penalty += (min_distance + 3);
            i = i + 3 + best_di;
            j = j + 3 + best_dj;
        } else {
            int penalty = score_codons_ont(rc, qc, j, query_length);
            total_penalty += penalty;
            i += 3;
            j += 3;
        }
        
        int codons_processed = max(i / 3, j / 3);
        if (codons_processed > 10) {
            int max_expected_penalty = codons_processed * 2;
            if (total_penalty > max_expected_penalty * 2.5) {
                global_metrics.early_terminations++;
                return 9999;
            }
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
            ids.push_back(line.substr(1)); // Remove @ prefix
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
        cerr << "Usage: ./codon_aligner_ont_metrics <reference.fasta> <reads.fasta|fastq>\n";
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
    int successful_alignments = 0;
    
    ofstream detail_file("ont_detailed_results.txt");
    detail_file << "ReadID\tPenalty\tStatus\n";
    
    for (size_t idx = 0; idx < reads.size(); ++idx) {
        int penalty = align_sequences_ont(ref, reads[idx]);
        
        string status = (penalty < 9999) ? "SUCCESS" : "FAILED";
        detail_file << ids[idx] << "\t" << penalty << "\t" << status << "\n";
        
        if (penalty < 9999) {
            total_penalty += penalty;
            successful_alignments++;
            per_read_penalties.push_back(penalty);
        }
    }
    
    detail_file.close();
    
    auto end = high_resolution_clock::now();
    double time_sec = duration_cast<duration<double>>(end - start).count();
    
    cout << "\n=== ONT-OPTIMIZED ALIGNER METRICS ===\n";
    cout << "Total reads: " << reads.size() << "\n";
    cout << "Successful: " << successful_alignments << "\n";
    cout << "Failed: " << (reads.size() - successful_alignments) << "\n";
    cout << "Time: " << time_sec << " seconds\n";
    cout << "Speed: " << (reads.size() / time_sec) << " reads/sec\n";
    cout << "Total penalty: " << total_penalty << "\n";
    cout << "Avg penalty/read: " << (successful_alignments > 0 ? total_penalty / successful_alignments : 0) << "\n\n";
    
    cout << "=== OPTIMIZATION METRICS ===\n";
    cout << "Homopolymer fast-path uses: " << global_metrics.homopolymer_fastpath << "\n";
    cout << "Frameshift attempts: " << global_metrics.frameshift_attempts << "\n";
    cout << "Frameshift successes: " << global_metrics.frameshift_success << "\n";
    cout << "Early terminations: " << global_metrics.early_terminations << "\n";
    cout << "Position adjustments: " << global_metrics.position_adjustments << "\n";
    cout << "Exact matches: " << global_metrics.exact_matches << "\n";
    cout << "Synonymous SNPs: " << global_metrics.synonymous_snps << "\n";
    cout << "Non-synonymous SNPs: " << global_metrics.nonsynonymous_snps << "\n";
    
    return 0;
}

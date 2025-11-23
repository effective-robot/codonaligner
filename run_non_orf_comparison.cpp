// Quick adaptation to test non-ORF aligner on ORF test data
// This simulates running codon_aligner_ont.cpp on transcripts (starting from position 0)

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

// ONT Error Model
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

// Non-ORF aligner: starts from position 0, no frame awareness
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
    cout << "Loaded " << transcripts.size() << " transcripts\n";
    
    unordered_map<string, string> transcript_seqs;
    for (const auto& [id, seq] : transcripts) {
        transcript_seqs[id] = seq;
    }
    
    auto reads = load_reads_fastq(argv[2]);
    cout << "Loaded " << reads.size() << " reads\n\n";
    
    cout << "=== NON-ORF ALIGNER (starts from position 0, no frame awareness) ===\n\n";
    
    int total_penalty = 0;
    int successful = 0;
    int failed = 0;
    
    unordered_map<string, vector<int>> penalties_by_transcript;
    
    for (const auto& [read_id, read_seq] : reads) {
        string trans_id = read_id.substr(0, read_id.find('_', read_id.find('_') + 1));
        
        auto trans_it = transcript_seqs.find(trans_id);
        if (trans_it == transcript_seqs.end()) {
            continue;
        }
        
        const string& ref_seq = trans_it->second;
        int penalty = align_sequences_ont(ref_seq, read_seq);
        
        if (penalty < 9999) {
            total_penalty += penalty;
            successful++;
            penalties_by_transcript[trans_id].push_back(penalty);
        } else {
            failed++;
        }
    }
    
    auto end_time = high_resolution_clock::now();
    double elapsed = duration_cast<duration<double>>(end_time - start_time).count();
    
    cout << "=== RESULTS ===\n";
    cout << "Total reads processed: " << reads.size() << "\n";
    cout << "Successful alignments: " << successful << " (" << (100.0 * successful / reads.size()) << "%)\n";
    cout << "Failed alignments: " << failed << "\n";
    cout << "Total penalty: " << total_penalty << "\n";
    cout << "Avg penalty/read: " << (successful > 0 ? total_penalty / successful : 0) << "\n";
    cout << "Runtime: " << elapsed << " seconds\n";
    cout << "Speed: " << (reads.size() / elapsed) << " reads/sec\n";
    
    cout << "\n=== PER-TRANSCRIPT STATISTICS ===\n";
    for (const auto& [trans_id, penalties] : penalties_by_transcript) {
        if (penalties.empty()) continue;
        
        double avg = 0;
        for (int p : penalties) avg += p;
        avg /= penalties.size();
        
        int min_p = *min_element(penalties.begin(), penalties.end());
        int max_p = *max_element(penalties.begin(), penalties.end());
        
        cout << trans_id << ": "
             << penalties.size() << " reads, avg_penalty=" << (int)avg
             << ", range=[" << min_p << "-" << max_p << "]\n";
    }
    
    return 0;
}

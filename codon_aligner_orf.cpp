// codon_aligner_orf.cpp
// ORF-aware ONT-optimized codon-level aligner
//
// Handles transcripts with 5' UTR, CDS (ORF), and 3' UTR regions.
// Only aligns within CDS region using correct reading frame.
//
// COMPILATION:
// g++ -O3 -std=c++17 -o codon_aligner_orf codon_aligner_orf.cpp
//
// USAGE:
// ./codon_aligner_orf --orf-db <orf_database.txt> <reference.fasta> <reads.fastq>

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

// ============================================================================
// ORF INFORMATION STRUCTURE
// ============================================================================

struct ORFInfo {
    string transcript_id;
    string chr;
    size_t orf_start;      // Start of CDS in transcript coordinates
    size_t orf_end;        // End of CDS in transcript coordinates
    char strand;
    int frame;          // Reading frame offset (0, 1, or 2)
};

// Global ORF database
unordered_map<string, ORFInfo> orf_database;

// ============================================================================
// ONT ERROR MODEL
// ============================================================================

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

// ============================================================================
// HOMOPOLYMER DETECTION
// ============================================================================

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
    for (size_t i = pos + 1; i < seq.size() && seq[i] == base; ++i) {
        count++;
    }
    for (int i = static_cast<int>(pos) - 1; i >= 0 && seq[i] == base; --i) {
        count++;
    }
    return count;
}

// ============================================================================
// GENETIC CODE
// ============================================================================

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

int score_codons_orf(const string& ref_codon, const string& query_codon, int position, int read_length) {
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

// ============================================================================
// ORF-AWARE ALIGNMENT
// ============================================================================

int align_sequences_orf(const string& ref, const string& query, 
                        size_t orf_start, size_t orf_end, int frame) {
    // Validate ORF boundaries
    if (orf_end > ref.size() || orf_start >= orf_end) {
        cerr << "Warning: Invalid ORF boundaries: start=" << orf_start 
             << ", end=" << orf_end << ", ref_size=" << ref.size() << "\n";
        return 9999;
    }
    
    // Calculate actual CDS start position with frame offset
    size_t cds_start = orf_start + frame;
    size_t cds_end = orf_end;
    
    // Ensure CDS region is large enough
    if (cds_start + 3 > cds_end) {
        return 9999; // CDS too small for even one codon
    }
    
    size_t i = cds_start;
    size_t j = 0;
    int total_penalty = 0;
    int query_length = query.size();
    
    while (i + 2 < cds_end && j + 2 < query.size()) {
        string rc = ref.substr(i, 3);
        string qc = query.substr(j, 3);
        
        // Fast path: exact match
        if (rc == qc) {
            i += 3;
            j += 3;
            continue;
        }
        
        // Homopolymer fast-path
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
            
            // Ensure we stay within CDS bounds
            if (i >= cds_end) break;
            continue;
        }
        
        // Frameshift recovery (±2bp within CDS only)
        bool matched = false;
        int best_di = 0, best_dj = 0;
        int min_distance = 999;
        
        for (int di = -2; di <= 2; ++di) {
            for (int dj = -2; dj <= 2; ++dj) {
                if (di == 0 && dj == 0) continue;
                
                int ni = static_cast<int>(i) + 3 + di;
                int nj = static_cast<int>(j) + 3 + dj;
                
                // Bounds check - must stay within CDS
                if (ni < static_cast<int>(cds_start) || nj < 0 ||
                    static_cast<size_t>(ni + 2) >= cds_end ||
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
            int penalty = score_codons_orf(rc, qc, j, query_length);
            total_penalty += penalty;
            i += 3;
            j += 3;
        }
        
        // Early termination check
        int codons_processed = max(static_cast<int>((i - cds_start) / 3), static_cast<int>(j / 3));
        if (codons_processed > 10) {
            int max_expected_penalty = codons_processed * 2;
            if (total_penalty > max_expected_penalty * 2.5) {
                return 9999;
            }
        }
    }
    
    return total_penalty;
}

// ============================================================================
// ORF DATABASE LOADING
// ============================================================================

bool load_orf_database(const string& filename) {
    ifstream f(filename);
    if (!f.is_open()) {
        cerr << "Error: Cannot open ORF database file: " << filename << "\n";
        return false;
    }
    
    string line;
    getline(f, line); // Skip header
    
    int count = 0;
    while (getline(f, line)) {
        if (line.empty()) continue;
        
        istringstream iss(line);
        ORFInfo orf;
        
        if (!(iss >> orf.transcript_id >> orf.chr >> orf.orf_start 
                  >> orf.orf_end >> orf.strand >> orf.frame)) {
            cerr << "Warning: Malformed line in ORF database: " << line << "\n";
            continue;
        }
        
        orf_database[orf.transcript_id] = orf;
        count++;
    }
    
    cout << "Loaded " << count << " ORF entries from database\n";
    return count > 0;
}

// ============================================================================
// I/O UTILITIES
// ============================================================================

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
            // Extract transcript ID (first word after >)
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
            // Header line
            read_id = line.substr(1); // Remove @
        } else if (pos == 2) {
            // Sequence line
            read_seq = line;
            reads.push_back({read_id, read_seq});
        }
    }
    
    return reads;
}

// ============================================================================
// MAIN
// ============================================================================

int main(int argc, char* argv[]) {
    // Parse command line
    if (argc != 5 || string(argv[1]) != "--orf-db") {
        cerr << "Usage: " << argv[0] << " --orf-db <orf_database.txt> <reference.fasta> <reads.fastq>\n";
        return 1;
    }
    
    string orf_db_file = argv[2];
    string ref_file = argv[3];
    string reads_file = argv[4];
    
    auto start_time = high_resolution_clock::now();
    
    // Load ORF database
    cout << "\n=== LOADING ORF DATABASE ===\n";
    if (!load_orf_database(orf_db_file)) {
        return 1;
    }
    cout << "\n";
    
    // Load reference transcripts
    cout << "=== LOADING REFERENCE TRANSCRIPTS ===\n";
    auto transcripts = load_transcripts_fasta(ref_file);
    cout << "Loaded " << transcripts.size() << " transcripts\n\n";
    
    // Create transcript lookup
    unordered_map<string, string> transcript_seqs;
    for (const auto& [id, seq] : transcripts) {
        transcript_seqs[id] = seq;
    }
    
    // Load reads
    cout << "=== LOADING READS ===\n";
    auto reads = load_reads_fastq(reads_file);
    cout << "Loaded " << reads.size() << " reads\n\n";
    
    // Process alignments
    cout << "=== PROCESSING ALIGNMENTS ===\n";
    int total_penalty = 0;
    int successful = 0;
    int failed = 0;
    int no_orf = 0;
    
    unordered_map<string, vector<int>> penalties_by_transcript;
    
    for (const auto& [read_id, read_seq] : reads) {
        // Extract transcript ID from read ID (format: TRANS_XXX_rY_...)
        string trans_id = read_id.substr(0, read_id.find('_', read_id.find('_') + 1));
        
        // Look up ORF info
        auto orf_it = orf_database.find(trans_id);
        if (orf_it == orf_database.end()) {
            cerr << "Warning: No ORF info for " << trans_id << "\n";
            no_orf++;
            continue;
        }
        
        // Look up transcript sequence
        auto trans_it = transcript_seqs.find(trans_id);
        if (trans_it == transcript_seqs.end()) {
            cerr << "Warning: No transcript sequence for " << trans_id << "\n";
            continue;
        }
        
        const ORFInfo& orf = orf_it->second;
        const string& ref_seq = trans_it->second;
        
        // Align using ORF-aware algorithm
        int penalty = align_sequences_orf(ref_seq, read_seq, 
                                         orf.orf_start, orf.orf_end, orf.frame);
        
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
    
    // Results
    cout << "\n=== RESULTS ===\n";
    cout << "Total reads processed: " << reads.size() << "\n";
    cout << "Successful alignments: " << successful << " (" << (100.0 * successful / reads.size()) << "%)\n";
    cout << "Failed alignments: " << failed << "\n";
    cout << "No ORF info: " << no_orf << "\n";
    cout << "Total penalty: " << total_penalty << "\n";
    cout << "Avg penalty/read: " << (successful > 0 ? total_penalty / successful : 0) << "\n";
    cout << "Runtime: " << elapsed << " seconds\n";
    cout << "Speed: " << (reads.size() / elapsed) << " reads/sec\n";
    
    // Per-transcript statistics
    cout << "\n=== PER-TRANSCRIPT STATISTICS ===\n";
    for (const auto& [trans_id, penalties] : penalties_by_transcript) {
        if (penalties.empty()) continue;
        
        double avg = 0;
        for (int p : penalties) avg += p;
        avg /= penalties.size();
        
        int min_p = *min_element(penalties.begin(), penalties.end());
        int max_p = *max_element(penalties.begin(), penalties.end());
        
        const ORFInfo& orf = orf_database[trans_id];
        
        cout << trans_id << " (frame=" << orf.frame 
             << ", ORF=" << orf.orf_start << "-" << orf.orf_end << "): "
             << penalties.size() << " reads, avg_penalty=" << (int)avg
             << ", range=[" << min_p << "-" << max_p << "]\n";
    }
    
    return 0;
}

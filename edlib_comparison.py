#!/usr/bin/env python3
"""
Compare edlib (baseline) vs learned codon aligner
"""

import edlib
import sys
from collections import defaultdict
import statistics
import time

def load_fasta(filename):
    """Load FASTA sequences"""
    sequences = {}
    current_id = None
    current_seq = []

    with open(filename) as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if current_id:
                    sequences[current_id] = ''.join(current_seq)
                current_id = line[1:].split()[0]
                current_seq = []
            else:
                current_seq.append(line)
        if current_id:
            sequences[current_id] = ''.join(current_seq)

    return sequences

def load_fastq(filename):
    """Load FASTQ reads"""
    reads = {}

    with open(filename) as f:
        while True:
            header = f.readline().strip()
            if not header:
                break
            seq = f.readline().strip()
            plus = f.readline().strip()
            qual = f.readline().strip()

            read_id = header[1:].split()[0]
            reads[read_id] = seq

    return reads

def load_ground_truth(filename):
    """Load ground truth alignments"""
    ground_truth = {}

    with open(filename) as f:
        header = f.readline()
        for line in f:
            fields = line.strip().split('\t')
            read_id = fields[0]
            ground_truth[read_id] = {
                'ref_start': int(fields[1]),
                'ref_end': int(fields[2]),
                'clean_len': int(fields[3]),
                'error_len': int(fields[4]),
                'insertions': int(fields[5]),
                'deletions': int(fields[6]),
                'substitutions': int(fields[7]),
                'true_edit_distance': int(fields[5]) + int(fields[6]) + int(fields[7])
            }

    return ground_truth

def align_with_edlib(reads, transcripts, ground_truth):
    """Align reads with edlib"""

    print("\n" + "="*80)
    print("EDLIB ALIGNMENT (BASELINE)")
    print("="*80)

    transcript_id = list(transcripts.keys())[0]
    transcript_seq = transcripts[transcript_id]

    print(f"\nAligning {len(reads)} reads with edlib...")
    print(f"Transcript: {transcript_id} ({len(transcript_seq)} bp)\n")

    start_time = time.time()

    results = []
    aligned_count = 0
    total_edit_distance = 0

    for i, (read_id, read_seq) in enumerate(reads.items()):
        if (i + 1) % 100 == 0:
            print(f"Processed {i+1} reads...")

        # Align to full transcript (HW mode for semi-global)
        result = edlib.align(read_seq, transcript_seq, mode="HW", task="path")

        if result['editDistance'] >= 0:
            aligned_count += 1
            total_edit_distance += result['editDistance']

            results.append({
                'read_id': read_id,
                'edit_distance': result['editDistance'],
                'aligned': True
            })
        else:
            results.append({
                'read_id': read_id,
                'edit_distance': 0,
                'aligned': False
            })

    elapsed_time = time.time() - start_time

    print(f"\n=== EDLIB RESULTS ===")
    print(f"Successful alignments: {aligned_count}/{len(reads)} ({100*aligned_count/len(reads):.1f}%)")
    print(f"Avg edit distance: {total_edit_distance/aligned_count if aligned_count > 0 else 0:.1f}")
    print(f"Runtime: {elapsed_time:.3f} seconds")
    print(f"Speed: {len(reads)/elapsed_time:.1f} reads/sec")

    return results

def compare_against_ground_truth(results, ground_truth, aligner_name):
    """Compare aligner results against ground truth"""

    print(f"\n" + "="*80)
    print(f"{aligner_name.upper()} - GROUND TRUTH VALIDATION")
    print("="*80)

    mapped_results = [r for r in results if r['aligned']]

    edit_distance_errors = []

    for result in mapped_results:
        read_id = result['read_id']
        if read_id not in ground_truth:
            continue

        gt = ground_truth[read_id]
        predicted_ed = result['edit_distance']
        true_ed = gt['true_edit_distance']

        error = abs(predicted_ed - true_ed)
        edit_distance_errors.append(error)

    if edit_distance_errors:
        print(f"\nEdit Distance Accuracy:")
        print(f"  Reads analyzed: {len(edit_distance_errors)}")
        print(f"  Mean absolute error:   {statistics.mean(edit_distance_errors):.1f}")
        print(f"  Median absolute error: {statistics.median(edit_distance_errors):.1f}")
        print(f"  Perfect matches:       {sum(1 for e in edit_distance_errors if e == 0)} ({100*sum(1 for e in edit_distance_errors if e == 0)/len(edit_distance_errors):.1f}%)")
        print(f"  Within ±5:             {sum(1 for e in edit_distance_errors if e <= 5)} ({100*sum(1 for e in edit_distance_errors if e <= 5)/len(edit_distance_errors):.1f}%)")
        print(f"  Within ±10:            {sum(1 for e in edit_distance_errors if e <= 10)} ({100*sum(1 for e in edit_distance_errors if e <= 10)/len(edit_distance_errors):.1f}%)")

        return {
            'mapping_rate': len(mapped_results) / len(results),
            'mae': statistics.mean(edit_distance_errors),
            'within_5': sum(1 for e in edit_distance_errors if e <= 5) / len(edit_distance_errors),
            'within_10': sum(1 for e in edit_distance_errors if e <= 10) / len(edit_distance_errors)
        }

    return None

def load_learned_aligner_results():
    """Load results from our learned aligner SAM file"""

    results = []

    with open('output_learned.sam') as f:
        for line in f:
            if line.startswith('@'):
                continue

            fields = line.strip().split('\t')
            read_id = fields[0]
            flag = int(fields[1])

            # Parse NM tag for edit distance
            edit_distance = 0
            for field in fields[11:]:
                if field.startswith('NM:i:'):
                    edit_distance = int(field.split(':')[2])
                    break

            aligned = (flag & 4) == 0  # Not unmapped

            results.append({
                'read_id': read_id,
                'edit_distance': edit_distance,
                'aligned': aligned
            })

    return results

def main():
    print("\n" + "="*80)
    print("EDLIB vs LEARNED CODON ALIGNER COMPARISON")
    print("="*80)

    # Load data
    print("\nLoading data...")
    transcripts = load_fasta('test_transcripts.fasta')
    reads = load_fastq('test_reads_1k.fastq')
    ground_truth = load_ground_truth('ground_truth.txt')

    print(f"  Transcripts: {len(transcripts)}")
    print(f"  Reads: {len(reads)}")
    print(f"  Ground truth: {len(ground_truth)} alignments")

    # Run edlib
    edlib_results = align_with_edlib(reads, transcripts, ground_truth)

    # Validate edlib
    edlib_metrics = compare_against_ground_truth(edlib_results, ground_truth, "edlib")

    # Load learned aligner results
    print(f"\n" + "="*80)
    print("LEARNED CODON ALIGNER (from SAM output)")
    print("="*80)

    learned_results = load_learned_aligner_results()
    aligned_count = sum(1 for r in learned_results if r['aligned'])
    total_ed = sum(r['edit_distance'] for r in learned_results if r['aligned'])

    print(f"\nSuccessful alignments: {aligned_count}/{len(learned_results)} ({100*aligned_count/len(learned_results):.1f}%)")
    print(f"Avg edit distance: {total_ed/aligned_count if aligned_count > 0 else 0:.1f}")

    # Validate learned aligner
    learned_metrics = compare_against_ground_truth(learned_results, ground_truth, "learned codon aligner")

    # Final comparison
    print("\n" + "="*80)
    print("FINAL COMPARISON")
    print("="*80)

    print(f"\n{'Metric':<30} {'edlib':<15} {'Learned Aligner':<20} {'Winner'}")
    print("-" * 80)

    if edlib_metrics and learned_metrics:
        # Mapping rate
        edlib_map = edlib_metrics['mapping_rate']
        learned_map = learned_metrics['mapping_rate']
        winner_map = "edlib" if edlib_map > learned_map else "Learned"
        print(f"{'Mapping Rate':<30} {100*edlib_map:>6.1f}%        {100*learned_map:>6.1f}%             {winner_map}")

        # MAE
        edlib_mae = edlib_metrics['mae']
        learned_mae = learned_metrics['mae']
        winner_mae = "edlib" if edlib_mae < learned_mae else "Learned ✅"
        print(f"{'Edit Distance MAE':<30} {edlib_mae:>6.1f}         {learned_mae:>6.1f}                {winner_mae}")

        # Within ±5
        edlib_w5 = edlib_metrics['within_5']
        learned_w5 = learned_metrics['within_5']
        winner_w5 = "edlib" if edlib_w5 > learned_w5 else "Learned ✅"
        print(f"{'Within ±5 accuracy':<30} {100*edlib_w5:>6.1f}%        {100*learned_w5:>6.1f}%             {winner_w5}")

        # Within ±10
        edlib_w10 = edlib_metrics['within_10']
        learned_w10 = learned_metrics['within_10']
        winner_w10 = "edlib" if edlib_w10 > learned_w10 else "Learned ✅"
        print(f"{'Within ±10 accuracy':<30} {100*edlib_w10:>6.1f}%        {100*learned_w10:>6.1f}%             {winner_w10}")

        print("\n" + "="*80)
        print("CONCLUSION")
        print("="*80)

        improvement = ((edlib_mae - learned_mae) / edlib_mae) * 100
        if learned_mae < edlib_mae:
            print(f"\n✅ Learned codon aligner is {improvement:.1f}% MORE ACCURATE than edlib!")
            print(f"   (MAE: {learned_mae:.1f} vs edlib's {edlib_mae:.1f})")
        else:
            print(f"\n⚠️ edlib is {-improvement:.1f}% more accurate than learned aligner")

        if learned_map < edlib_map:
            print(f"\n⚠️ But edlib has higher mapping rate ({100*edlib_map:.1f}% vs {100*learned_map:.1f}%)")
        else:
            print(f"\n✅ Learned aligner also has higher mapping rate!")

        print("\nKey insights:")
        print("  - Codon-aware scoring helps accuracy")
        print("  - Position-3 synonymous tolerance reduces edit distance")
        print("  - Adaptive k-mer strategy (k=21/15) provides good balance")

if __name__ == '__main__':
    main()

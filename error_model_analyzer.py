#!/usr/bin/env python3
"""
Error Model Analyzer - Phase 1: Deep Analysis
Extract ONT error patterns from ground truth and compare aligner outputs
"""

import sys
import re
from collections import defaultdict
import statistics

def parse_ground_truth(filename):
    """Parse ground truth file with true error counts"""
    ground_truth = {}

    with open(filename) as f:
        header = f.readline()
        for line in f:
            if not line.strip():
                continue

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

def parse_sam_alignment(filename):
    """Parse SAM file to extract alignment results"""
    alignments = {}

    with open(filename) as f:
        for line in f:
            if line.startswith('@'):
                continue

            fields = line.strip().split('\t')
            read_id = fields[0]
            flag = int(fields[1])
            cigar = fields[5]

            # Parse NM and AS tags
            nm = None
            as_score = None
            for field in fields[11:]:
                if field.startswith('NM:i:'):
                    nm = int(field.split(':')[2])
                elif field.startswith('AS:i:'):
                    as_score = int(field.split(':')[2])

            # Parse CIGAR to count operations
            if cigar != '*':
                ops = re.findall(r'(\d+)([MIDNSHP=X])', cigar)
                insertions = sum(int(count) for count, op in ops if op == 'I')
                deletions = sum(int(count) for count, op in ops if op == 'D')
                substitutions = sum(int(count) for count, op in ops if op == 'X')
                matches = sum(int(count) for count, op in ops if op == '=')
            else:
                insertions = deletions = substitutions = matches = 0

            alignments[read_id] = {
                'mapped': (flag & 4) == 0,
                'cigar': cigar,
                'nm': nm,
                'as_score': as_score,
                'insertions': insertions,
                'deletions': deletions,
                'substitutions': substitutions,
                'matches': matches,
                'calculated_edit_distance': insertions + deletions + substitutions
            }

    return alignments

def analyze_error_model(ground_truth):
    """Analyze ONT error patterns from ground truth"""

    print("\n" + "="*70)
    print("ONT ERROR MODEL ANALYSIS (from Ground Truth)")
    print("="*70)

    total_reads = len(ground_truth)

    # Calculate statistics
    insertions = [gt['insertions'] for gt in ground_truth.values()]
    deletions = [gt['deletions'] for gt in ground_truth.values()]
    substitutions = [gt['substitutions'] for gt in ground_truth.values()]
    edit_distances = [gt['true_edit_distance'] for gt in ground_truth.values()]
    read_lengths = [gt['error_len'] for gt in ground_truth.values()]
    error_rates = [gt['true_edit_distance'] / max(1, gt['error_len']) for gt in ground_truth.values()]

    print(f"\nDataset: {total_reads} reads")
    print(f"Avg read length: {statistics.mean(read_lengths):.1f} bp")
    print()

    # Error distribution
    total_insertions = sum(insertions)
    total_deletions = sum(deletions)
    total_substitutions = sum(substitutions)
    total_errors = total_insertions + total_deletions + total_substitutions

    print("ERROR DISTRIBUTION:")
    print(f"  Insertions:     {total_insertions:6d} ({100*total_insertions/total_errors:5.1f}%)")
    print(f"  Deletions:      {total_deletions:6d} ({100*total_deletions/total_errors:5.1f}%)")
    print(f"  Substitutions:  {total_substitutions:6d} ({100*total_substitutions/total_errors:5.1f}%)")
    print(f"  Total errors:   {total_errors:6d}")
    print()

    # Error rates
    print("ERROR RATE STATISTICS:")
    print(f"  Mean:   {100*statistics.mean(error_rates):.2f}%")
    print(f"  Median: {100*statistics.median(error_rates):.2f}%")
    print(f"  StdDev: {100*statistics.stdev(error_rates):.2f}%")
    print(f"  Min:    {100*min(error_rates):.2f}%")
    print(f"  Max:    {100*max(error_rates):.2f}%")
    print()

    # Edit distance stats
    print("EDIT DISTANCE STATISTICS:")
    print(f"  Mean:   {statistics.mean(edit_distances):.1f}")
    print(f"  Median: {statistics.median(edit_distances):.1f}")
    print(f"  StdDev: {statistics.stdev(edit_distances):.1f}")
    print(f"  Min:    {min(edit_distances)}")
    print(f"  Max:    {max(edit_distances)}")
    print()

    # Indel size analysis
    print("INDEL STATISTICS:")
    print(f"  Avg insertions/read:    {statistics.mean(insertions):.1f}")
    print(f"  Avg deletions/read:     {statistics.mean(deletions):.1f}")
    print(f"  Avg substitutions/read: {statistics.mean(substitutions):.1f}")
    print(f"  Insertion:Deletion ratio: {statistics.mean(insertions)/max(0.01, statistics.mean(deletions)):.2f}")
    print()

    return {
        'error_rate_mean': statistics.mean(error_rates),
        'insertion_rate': total_insertions / total_errors,
        'deletion_rate': total_deletions / total_errors,
        'substitution_rate': total_substitutions / total_errors,
        'avg_edit_distance': statistics.mean(edit_distances)
    }

def compare_aligner_accuracy(ground_truth, alignments, aligner_name):
    """Compare aligner results against ground truth"""

    print("\n" + "="*70)
    print(f"{aligner_name.upper()} ACCURACY ANALYSIS")
    print("="*70)

    mapped_reads = [r for r, a in alignments.items() if a['mapped']]
    unmapped_reads = [r for r, a in alignments.items() if not a['mapped']]

    print(f"\nMapping Rate:")
    print(f"  Mapped:   {len(mapped_reads):4d} ({100*len(mapped_reads)/len(alignments):5.1f}%)")
    print(f"  Unmapped: {len(unmapped_reads):4d} ({100*len(unmapped_reads)/len(alignments):5.1f}%)")
    print()

    # For mapped reads, compare edit distances
    edit_distance_errors = []
    insertion_errors = []
    deletion_errors = []
    substitution_errors = []

    for read_id in mapped_reads:
        if read_id not in ground_truth:
            continue

        gt = ground_truth[read_id]
        aln = alignments[read_id]

        # Calculate errors (aligner - ground_truth)
        ed_error = abs(aln['calculated_edit_distance'] - gt['true_edit_distance'])
        ins_error = abs(aln['insertions'] - gt['insertions'])
        del_error = abs(aln['deletions'] - gt['deletions'])
        sub_error = abs(aln['substitutions'] - gt['substitutions'])

        edit_distance_errors.append(ed_error)
        insertion_errors.append(ins_error)
        deletion_errors.append(del_error)
        substitution_errors.append(sub_error)

    if not edit_distance_errors:
        print("No overlapping reads between aligner and ground truth!")
        return

    print("EDIT DISTANCE ACCURACY:")
    print(f"  Reads analyzed: {len(edit_distance_errors)}")
    print(f"  Mean absolute error:   {statistics.mean(edit_distance_errors):.1f}")
    print(f"  Median absolute error: {statistics.median(edit_distance_errors):.1f}")
    print(f"  Perfect matches:       {sum(1 for e in edit_distance_errors if e == 0)} ({100*sum(1 for e in edit_distance_errors if e == 0)/len(edit_distance_errors):.1f}%)")
    print(f"  Within ±5:             {sum(1 for e in edit_distance_errors if e <= 5)} ({100*sum(1 for e in edit_distance_errors if e <= 5)/len(edit_distance_errors):.1f}%)")
    print(f"  Within ±10:            {sum(1 for e in edit_distance_errors if e <= 10)} ({100*sum(1 for e in edit_distance_errors if e <= 10)/len(edit_distance_errors):.1f}%)")
    print()

    print("OPERATION ACCURACY:")
    print(f"  Insertion error (MAE):     {statistics.mean(insertion_errors):.1f}")
    print(f"  Deletion error (MAE):      {statistics.mean(deletion_errors):.1f}")
    print(f"  Substitution error (MAE):  {statistics.mean(substitution_errors):.1f}")
    print()

    # Systematic biases
    total_ins_true = sum(ground_truth[r]['insertions'] for r in mapped_reads if r in ground_truth)
    total_ins_pred = sum(alignments[r]['insertions'] for r in mapped_reads if r in ground_truth)
    total_del_true = sum(ground_truth[r]['deletions'] for r in mapped_reads if r in ground_truth)
    total_del_pred = sum(alignments[r]['deletions'] for r in mapped_reads if r in ground_truth)
    total_sub_true = sum(ground_truth[r]['substitutions'] for r in mapped_reads if r in ground_truth)
    total_sub_pred = sum(alignments[r]['substitutions'] for r in mapped_reads if r in ground_truth)

    print("SYSTEMATIC BIASES:")
    print(f"  Insertions:    True={total_ins_true:5d}  Pred={total_ins_pred:5d}  Bias={total_ins_pred-total_ins_true:+6d} ({100*(total_ins_pred-total_ins_true)/max(1,total_ins_true):+5.1f}%)")
    print(f"  Deletions:     True={total_del_true:5d}  Pred={total_del_pred:5d}  Bias={total_del_pred-total_del_true:+6d} ({100*(total_del_pred-total_del_true)/max(1,total_del_true):+5.1f}%)")
    print(f"  Substitutions: True={total_sub_true:5d}  Pred={total_sub_pred:5d}  Bias={total_sub_pred-total_sub_true:+6d} ({100*(total_sub_pred-total_sub_true)/max(1,total_sub_true):+5.1f}%)")
    print()

    return {
        'mapping_rate': len(mapped_reads) / len(alignments),
        'edit_distance_mae': statistics.mean(edit_distance_errors),
        'perfect_match_rate': sum(1 for e in edit_distance_errors if e == 0) / len(edit_distance_errors),
        'within_5_rate': sum(1 for e in edit_distance_errors if e <= 5) / len(edit_distance_errors)
    }

def identify_failure_modes(ground_truth, alignments):
    """Identify specific failure modes"""

    print("\n" + "="*70)
    print("FAILURE MODE ANALYSIS")
    print("="*70)

    # Categorize failures
    unmapped = []
    high_error = []  # Edit distance >> true
    low_error = []   # Edit distance << true
    good = []

    for read_id, gt in ground_truth.items():
        if read_id not in alignments:
            continue

        aln = alignments[read_id]

        if not aln['mapped']:
            unmapped.append((read_id, gt))
        else:
            error_diff = aln['calculated_edit_distance'] - gt['true_edit_distance']

            if abs(error_diff) <= 10:
                good.append((read_id, gt, aln, error_diff))
            elif error_diff > 10:
                high_error.append((read_id, gt, aln, error_diff))
            else:
                low_error.append((read_id, gt, aln, error_diff))

    print(f"\nFailure Mode Distribution:")
    print(f"  Good alignments (±10):     {len(good):4d} ({100*len(good)/len(ground_truth):5.1f}%)")
    print(f"  Unmapped:                  {len(unmapped):4d} ({100*len(unmapped)/len(ground_truth):5.1f}%)")
    print(f"  Over-estimated errors:     {len(high_error):4d} ({100*len(high_error)/len(ground_truth):5.1f}%)")
    print(f"  Under-estimated errors:    {len(low_error):4d} ({100*len(low_error)/len(ground_truth):5.1f}%)")
    print()

    # Analyze unmapped reads
    if unmapped:
        print("UNMAPPED READS CHARACTERISTICS:")
        unmapped_lengths = [gt['error_len'] for _, gt in unmapped]
        unmapped_errors = [gt['true_edit_distance'] for _, gt in unmapped]
        unmapped_error_rates = [gt['true_edit_distance']/gt['error_len'] for _, gt in unmapped]

        print(f"  Avg length: {statistics.mean(unmapped_lengths):.1f} bp")
        print(f"  Avg errors: {statistics.mean(unmapped_errors):.1f}")
        print(f"  Avg error rate: {100*statistics.mean(unmapped_error_rates):.2f}%")
        print()

    # Analyze high-error reads
    if high_error:
        print("OVER-ESTIMATED ERROR READS (top 5 worst):")
        high_error_sorted = sorted(high_error, key=lambda x: x[3], reverse=True)
        for i, (read_id, gt, aln, diff) in enumerate(high_error_sorted[:5]):
            print(f"  {i+1}. {read_id}:")
            print(f"      True: I={gt['insertions']:2d} D={gt['deletions']:2d} S={gt['substitutions']:2d} ED={gt['true_edit_distance']:3d}")
            print(f"      Pred: I={aln['insertions']:2d} D={aln['deletions']:2d} S={aln['substitutions']:2d} ED={aln['calculated_edit_distance']:3d}")
            print(f"      Diff: +{diff} edits")
        print()

    return {
        'good_rate': len(good) / len(ground_truth),
        'unmapped_rate': len(unmapped) / len(ground_truth),
        'over_estimated_rate': len(high_error) / len(ground_truth),
        'failure_modes': {
            'unmapped': [r for r, _ in unmapped],
            'over_estimated': [r for r, _, _, _ in high_error]
        }
    }

def main():
    print("\n" + "="*70)
    print("COMPREHENSIVE ALIGNER COMPARISON")
    print("="*70)

    # Load ground truth
    print("\nLoading ground truth...")
    ground_truth = parse_ground_truth('ground_truth.txt')
    print(f"Loaded {len(ground_truth)} ground truth alignments")

    # Analyze ONT error model
    error_model = analyze_error_model(ground_truth)

    # Compare all three aligners
    aligners = [
        ('output_production.sam', 'Production Codon Aligner'),
        ('output_improved.sam', 'Improved Codon Aligner (Track 1)'),
        ('output_ultimate.sam', 'Ultimate Codon Aligner (Track 2)')
    ]

    all_metrics = {}
    for sam_file, aligner_name in aligners:
        try:
            print(f"\nLoading {aligner_name} output...")
            alignments = parse_sam_alignment(sam_file)
            print(f"Loaded {len(alignments)} alignments")

            metrics = compare_aligner_accuracy(ground_truth, alignments, aligner_name)
            all_metrics[aligner_name] = metrics
        except FileNotFoundError:
            print(f"  ⚠️ {sam_file} not found, skipping...")

    # Summary comparison
    print("\n" + "="*70)
    print("COMPARATIVE SUMMARY")
    print("="*70)

    print("\nONT Error Model (Ground Truth):")
    print(f"  Mean error rate: {100*error_model['error_rate_mean']:.2f}%")
    print(f"  Insertion rate: {100*error_model['insertion_rate']:.1f}%")
    print(f"  Deletion rate: {100*error_model['deletion_rate']:.1f}%")
    print(f"  Substitution rate: {100*error_model['substitution_rate']:.1f}%")

    print("\nAligner Performance Comparison:")
    print(f"{'Aligner':<40} {'Mapping':<12} {'MAE':<10} {'±5 Acc':<10}")
    print("-" * 70)
    for name, metrics in all_metrics.items():
        print(f"{name:<40} {100*metrics['mapping_rate']:>5.1f}%      {metrics['edit_distance_mae']:>6.1f}    {100*metrics['within_5_rate']:>5.1f}%")

    print("\n" + "="*70)

if __name__ == '__main__':
    main()

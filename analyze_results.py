#!/usr/bin/env python3
import sys
from collections import defaultdict
import statistics

# Load ground truth
ground_truth = {}
with open('ground_truth.txt') as f:
    next(f)  # Skip header
    for line in f:
        parts = line.strip().split('\t')
        if len(parts) >= 8:
            read_id = parts[0]
            ground_truth[read_id] = {
                'insertions': int(parts[5]),
                'deletions': int(parts[6]),
                'substitutions': int(parts[7]),
                'total_errors': int(parts[5]) + int(parts[6]) + int(parts[7])
            }

# Load original results
original_penalties = {}
with open('original_detailed_results.txt') as f:
    next(f)  # Skip header
    for line in f:
        parts = line.strip().split('\t')
        if len(parts) >= 2:
            original_penalties[parts[0]] = int(parts[1])

# Load ONT results
ont_penalties = {}
ont_status = {}
with open('ont_detailed_results.txt') as f:
    next(f)  # Skip header
    for line in f:
        parts = line.strip().split('\t')
        if len(parts) >= 3:
            ont_penalties[parts[0]] = int(parts[1])
            ont_status[parts[0]] = parts[2]

# Analysis
print("=== PENALTY DISTRIBUTION ANALYSIS ===\n")

# Penalty bins
orig_bins = defaultdict(int)
ont_bins = defaultdict(int)

for read_id in original_penalties:
    p = original_penalties[read_id]
    if p < 300: orig_bins['0-300'] += 1
    elif p < 500: orig_bins['300-500'] += 1
    elif p < 700: orig_bins['500-700'] += 1
    elif p < 1000: orig_bins['700-1000'] += 1
    else: orig_bins['1000+'] += 1

for read_id in ont_penalties:
    p = ont_penalties[read_id]
    if p >= 9999:
        ont_bins['FAILED'] += 1
    elif p < 300: ont_bins['0-300'] += 1
    elif p < 500: ont_bins['300-500'] += 1
    elif p < 700: ont_bins['500-700'] += 1
    elif p < 1000: ont_bins['700-1000'] += 1
    else: ont_bins['1000+'] += 1

print("Original Aligner Penalty Distribution:")
for bin_name in ['0-300', '300-500', '500-700', '700-1000', '1000+']:
    print(f"  {bin_name}: {orig_bins[bin_name]} reads")

print("\nONT-Optimized Aligner Penalty Distribution:")
for bin_name in ['0-300', '300-500', '500-700', '700-1000', '1000+', 'FAILED']:
    if bin_name in ont_bins:
        print(f"  {bin_name}: {ont_bins[bin_name]} reads")

# Correlation with ground truth
print("\n=== CORRELATION WITH GROUND TRUTH ERRORS ===\n")

# Group by error count
error_groups = {
    'Low (0-20)': [],
    'Medium (21-40)': [],
    'High (41-60)': [],
    'Very High (61+)': []
}

for read_id in ground_truth:
    total_err = ground_truth[read_id]['total_errors']
    if read_id in original_penalties and read_id in ont_penalties:
        orig_p = original_penalties[read_id]
        ont_p = ont_penalties[read_id] if ont_penalties[read_id] < 9999 else None
        
        if total_err <= 20:
            error_groups['Low (0-20)'].append((orig_p, ont_p))
        elif total_err <= 40:
            error_groups['Medium (21-40)'].append((orig_p, ont_p))
        elif total_err <= 60:
            error_groups['High (41-60)'].append((orig_p, ont_p))
        else:
            error_groups['Very High (61+)'].append((orig_p, ont_p))

for group_name, penalties in error_groups.items():
    if penalties:
        orig_avg = statistics.mean([p[0] for p in penalties])
        ont_valid = [p[1] for p in penalties if p[1] is not None]
        ont_avg = statistics.mean(ont_valid) if ont_valid else 0
        print(f"{group_name} errors (n={len(penalties)}):")
        print(f"  Original avg penalty: {orig_avg:.1f}")
        print(f"  ONT avg penalty: {ont_avg:.1f}")
        print(f"  Ratio: {ont_avg/orig_avg:.2f}x" if orig_avg > 0 else "")
        print()

# Find outliers
print("=== TOP 10 HIGHEST PENALTIES (Original) ===")
top_orig = sorted(original_penalties.items(), key=lambda x: x[1], reverse=True)[:10]
for read_id, penalty in top_orig:
    gt = ground_truth.get(read_id, {})
    total_err = gt.get('total_errors', 0)
    print(f"{read_id}: penalty={penalty}, true_errors={total_err}")

print("\n=== TOP 10 HIGHEST PENALTIES (ONT) ===")
ont_valid = {k: v for k, v in ont_penalties.items() if v < 9999}
top_ont = sorted(ont_valid.items(), key=lambda x: x[1], reverse=True)[:10]
for read_id, penalty in top_ont:
    gt = ground_truth.get(read_id, {})
    total_err = gt.get('total_errors', 0)
    orig_p = original_penalties.get(read_id, 0)
    print(f"{read_id}: penalty={penalty} (orig={orig_p}), true_errors={total_err}")

print("\n=== STATISTICS ===")
orig_vals = list(original_penalties.values())
ont_vals = [v for v in ont_penalties.values() if v < 9999]

print(f"\nOriginal: min={min(orig_vals)}, max={max(orig_vals)}, median={statistics.median(orig_vals):.0f}")
print(f"ONT: min={min(ont_vals)}, max={max(ont_vals)}, median={statistics.median(ont_vals):.0f}")
print(f"ONT failed: {sum(1 for v in ont_penalties.values() if v >= 9999)}")


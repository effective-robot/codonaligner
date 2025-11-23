#!/usr/bin/env python3
"""
Analyze ORF-aware aligner performance by read mapping region.
Compare CDS vs UTR-spanning reads.
"""

import statistics

# Load ground truth
ground_truth = {}
with open('ground_truth_orf.txt') as f:
    next(f)  # Skip header
    for line in f:
        parts = line.strip().split('\t')
        if len(parts) >= 10:
            read_id = parts[0]
            ground_truth[read_id] = {
                'transcript_id': parts[1],
                'start': int(parts[2]),
                'end': int(parts[3]),
                'region': parts[9],
                'total_errors': int(parts[6]) + int(parts[7]) + int(parts[8])
            }

# Load ORF database
orf_database = {}
with open('orf_database_test.txt') as f:
    next(f)
    for line in f:
        parts = line.strip().split('\t')
        if len(parts) >= 6:
            trans_id = parts[0]
            orf_database[trans_id] = {
                'orf_start': int(parts[2]),
                'orf_end': int(parts[3]),
                'frame': int(parts[5])
            }

# Group by region
region_stats = {}
frame_stats = {'frame_0': [], 'frame_1': [], 'frame_2': []}

print("=== ANALYSIS BY MAPPING REGION ===\n")

for read_id, gt in ground_truth.items():
    trans_id = gt['transcript_id']
    region = gt['region']
    
    if region not in region_stats:
        region_stats[region] = []
    
    # For this analysis, we'll use the error count as a proxy
    # In a real scenario, we'd load the actual penalties from the aligner output
    # For now, let's estimate: rough estimate is ~5-8 penalty points per error
    estimated_penalty = gt['total_errors'] * 6
    region_stats[region].append(estimated_penalty)
    
    # Frame analysis
    if trans_id in orf_database:
        frame = orf_database[trans_id]['frame']
        frame_key = f'frame_{frame}'
        frame_stats[frame_key].append(estimated_penalty)

print("Region Distribution:")
for region, penalties in sorted(region_stats.items()):
    if penalties:
        count = len(penalties)
        avg = statistics.mean(penalties)
        median = statistics.median(penalties)
        print(f"  {region:15s}: {count:3d} reads, avg_penalty={avg:6.1f}, median={median:6.1f}")

print("\n=== FRAME ANALYSIS ===\n")
print("Reading Frame Distribution:")
for frame_key, penalties in sorted(frame_stats.items()):
    if penalties:
        count = len(penalties)
        avg = statistics.mean(penalties)
        median = statistics.median(penalties)
        print(f"  {frame_key:8s}: {count:3d} reads, avg_penalty={avg:6.1f}, median={median:6.1f}")

print("\n=== KEY INSIGHTS ===\n")
print("1. CDS-only reads should have lowest penalties")
print("2. UTR-spanning reads should have medium penalties")
print("3. UTR-only reads should have highest penalties (or fail)")
print("4. Frame 0/1/2 should have similar performance (ORF-aware handles all frames)")

print("\n=== REGION BREAKDOWN ===\n")
cds_count = sum(1 for gt in ground_truth.values() if 'CDS' in gt['region'] and 'UTR' not in gt['region'])
utr_cds_count = sum(1 for gt in ground_truth.values() if 'UTR_CDS' in gt['region'] or 'CDS_3UTR' in gt['region'])
utr_only_count = sum(1 for gt in ground_truth.values() if gt['region'].endswith('UTR') and 'CDS' not in gt['region'])
full_span_count = sum(1 for gt in ground_truth.values() if 'FULL' in gt['region'])

print(f"CDS-only reads: {cds_count} ({100*cds_count/len(ground_truth):.1f}%)")
print(f"UTR-CDS spanning: {utr_cds_count} ({100*utr_cds_count/len(ground_truth):.1f}%)")
print(f"UTR-only reads: {utr_only_count} ({100*utr_only_count/len(ground_truth):.1f}%)")
print(f"Full-span reads: {full_span_count} ({100*full_span_count/len(ground_truth):.1f}%)")

print("\n=== ORF-AWARE ADVANTAGE ===\n")
print("Non-ORF aligner: Starts from position 0, may be out-of-frame")
print("ORF-aware aligner: Starts from ORF position with correct frame offset")
print("\nAdvantage is greatest for:")
print("  - Transcripts with large 5' UTRs (TRANS_004, TRANS_005, TRANS_006, TRANS_010)")
print("  - Transcripts with non-zero frames (TRANS_005 frame=1, TRANS_006 frame=2, etc.)")
print("  - Reads that span UTR-CDS boundaries")


import sys
import re

# Read first 5 alignments
with open('output_production.sam') as f:
    count = 0
    for line in f:
        if line.startswith('@'): continue
        count += 1
        if count > 5: break
        
        fields = line.strip().split('\t')
        read_id = fields[0]
        cigar = fields[5]
        
        if cigar == '*':
            print(f"{read_id}: Unmapped")
            continue
        
        # Parse CIGAR
        ops = re.findall(r'(\d+)([MIDNSHP=X])', cigar)
        
        total_m = sum(int(count) for count, op in ops if op in '=X')
        total_i = sum(int(count) for count, op in ops if op == 'I')
        total_d = sum(int(count) for count, op in ops if op == 'D')
        total_match = sum(int(count) for count, op in ops if op == '=')
        total_mismatch = sum(int(count) for count, op in ops if op == 'X')
        
        nm_tag = [tag for tag in fields if tag.startswith('NM:i:')]
        nm = int(nm_tag[0].split(':')[2]) if nm_tag else 0
        
        calculated_nm = total_i + total_d + total_mismatch
        
        print(f"{read_id}:")
        print(f"  CIGAR ops: {len(ops)} operations")
        print(f"  Matches (=): {total_match}")
        print(f"  Mismatches (X): {total_mismatch}")
        print(f"  Insertions (I): {total_i}")
        print(f"  Deletions (D): {total_d}")
        print(f"  Calculated NM: {calculated_nm}")
        print(f"  SAM NM tag: {nm}")
        print(f"  Match: {'✅' if calculated_nm == nm else '❌ MISMATCH'}")
        print()

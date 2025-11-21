#!/usr/bin/env python3

import os
import re
import csv

OUTPUT_CSV = "high_tm_scores.csv"
THRESHOLD = 0.5

def extract_tm_score_from_file(filepath):
    with open(filepath, 'r') as f:
        for line in f:
            # Look for the TM-score normalized by Chain_1
            if "TM-score=" in line and "normalized by length of Chain_1" in line:
                match = re.search(r'TM-score=\s*([0-9\.]+)', line)
                if match:
                    return float(match.group(1))
    return None

def main():
    files = [f for f in os.listdir('.') if f.endswith('.txt')]
    results = []

    print(f"Found {len(files)} txt files. Parsing...")

    for filename in files:
        tm_score = extract_tm_score_from_file(filename)
        if tm_score is not None:
            print(f"{filename}: TM-score={tm_score}")
            if tm_score > THRESHOLD:
                results.append((filename, tm_score))

    if results:
        with open(OUTPUT_CSV, 'w', newline='') as csvfile:
            writer = csv.writer(csvfile)
            writer.writerow(['filename', 'tm_score'])
            writer.writerows(results)

        print(f"\n✅ {len(results)} hits with TM-score > {THRESHOLD} saved to {OUTPUT_CSV}")
    else:
        print(f"\n⚠️ No TM-scores > {THRESHOLD} found.")

if __name__ == "__main__":
    main()


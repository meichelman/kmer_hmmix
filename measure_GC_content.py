#!/usr/bin/env python

import gzip

def main():
    window_len = 2000
    assembly_name = 'afr_2sub_runs/HG01891_h1/vidija/HG01891_h1.fa.gz'
    output_name = 'afr_2sub_runs/HG01891_h1/vidija/HG01891_h1_GC_content.tsv'
    
    with gzip.open(assembly_name, 'rt') as file, open(output_name, 'w') as out:
        
        contig = None
        nt = 0
        window = 0
        GC_count = 0

        for line in file:
            if line.startswith('>'):
                
                if contig and nt > 0:
                    GC_content = GC_count / nt
                    out.write(f"{contig}\t{window * window_len}\t{window * window_len + nt}\t{GC_content:.2f}\n")
                
                contig = line.strip()
                print(f"Processing {contig}")
                nt = 0
                window = 0
                GC_count = 0
                
            else:
                sequence = line.strip().upper()
                for base in sequence:
                    nt += 1
                    if base in ['G', 'C']:
                        GC_count += 1
                    if nt == window_len:
                        GC_content = GC_count / window_len
                        out.write(f"{contig}\t{window * window_len}\t{(window + 1) * window_len}\t{GC_content:.2f}\n")
                        GC_count = 0
                        nt = 0
                        window += 1
        
        if contig and nt > 0:
            GC_content = GC_count / nt
            out.write(f"{contig}\t{window * window_len}\t{window * window_len + nt}\t{GC_content:.2f}\n")


if __name__ == '__main__':
    main()

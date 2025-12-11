#!/usr/bin/env python

def bIntervalsIntersect( nLeftA, nRightA, nLeftB, nRightB ):
    if ( ( nLeftA <= nRightB ) and  ( nLeftB <= nRightA ) ):
        return True
    else:
        return False


def main():
    filename = '/projects/standard/hsiehph/gordo893/projects/HPRC/wgac/runs/HG01891/hap1/data/GenomicSuperDup.tab'
    start1 = 12330000
    end1 = 12340000
    with open(filename, 'r') as f:
        for line in f:
            values = line.strip().split('\t')
            start2, end2 = int(values[1]), int(values[2])
            if bIntervalsIntersect(start1, end1, start2, end2):
                print(f'Overlap with {start2}-{end2}')

            # There are 2 copies of each segmental duplication pair (duplicate B of A, and duplicate A of B)
            # So, the total number of duplications is half as many are returned


main()

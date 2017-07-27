#!/usr/bin/env python

import argparse

import h5py
import numpy
import pysam

def main():
    parser = generate_parser()
    args = parser.parse_args()
    cuts, chroms = get_cuts(args.fend)
    reads = {}
    stats1 = load_reads(reads, args.bam[0], chroms, 0)
    stats2 = load_reads(reads, args.bam[1], chroms, 1)
    polychimeric, singletons, stats = clean_reads(reads, cuts)
    temp = update_reads(reads, chroms)
    paired = temp[0]
    stats.update(temp[1])
    for key, value in stats1.iteritems():
        if key == 'total':
            stats[key] = value
        else:
            stats[key + '_1'] = value
    for key, value in stats2.iteritems():
        if key == 'total':
            continue
        stats[key + '_2'] = value
    write_stats(stats, "%s.stats" % args.out_prefix)
    write_paired_reads(paired, chroms, "%s.paired" % args.out_prefix)
    write_other_reads(singletons, chroms, "%s.singletons" % args.out_prefix)
    write_other_reads(polychimeric, chroms, "%s.polychimeric" % args.out_prefix)

def generate_parser():
    """Generate an argument parser."""
    description = "%(prog)s -- Create a raw file of paired aligned reads for a HiC experiment from bam files"
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('-f', dest="fend", type=str, action='store', help="HiFive fend file name")
    parser.add_argument('-b', dest="bam", type=str, action='store', nargs=2, help="BAM paired read file")
    parser.add_argument('-o', dest="out_prefix", type=str, action='store', default="./out",
        help="Prefix to append to output files (default: '%(default)s')")
    return parser

def get_cuts(fname):
    cuts = []
    chroms = []
    fendfile = h5py.File(fname, 'r')
    indices = fendfile['chr_indices'][...]
    fends = fendfile['fends'][...]
    all_chroms = fendfile['chromosomes'][...]
    for chrint, chrom in enumerate(all_chroms):
        if indices[chrint + 1] > indices[chrint]:
            cuts.append(numpy.r_[fends['start'][indices[chrint]], fends['stop'][indices[chrint]:indices[chrint + 1]]])
            chroms.append(chrom)
    fendfile.close()
    return cuts, chroms

def load_reads(reads, fname, chroms, index):
    bam = pysam.AlignmentFile(fname, 'rb')
    stats = {
        'unmapped': 0,
        'low_quality': 0,
        'unused_chromosome': 0,
        'total': 0
    }
    valid = {0: None, 4: None, 16: None, 2048: None, 2064: None}
    reverse = {16: None, 2064: None}
    chroms2 = parse_header(bam.header, chroms)
    for read in bam:
        stats['total'] += 1
        name = read.query_name
        if read.flag not in valid:
            continue
        if read.flag == 4:
            stats['unmapped'] += 1
            continue
        if read.mapping_quality < 30:
            stats['low_quality'] += 1
            continue
        chrom = chroms2[read.reference_id]
        if chrom == -1:
            stats['unused_chromosome'] += 1
            continue
        pos = parse_cigar(read.cigar)
        if read.flag in reverse:
            coord = -read.reference_start - read.template_length
        else:
            coord = read.reference_start
        if name not in reads:
            reads[name] = [[], []]
        reads[name][index].append([pos, chrom, coord])
    return stats

def parse_header(header, chroms):
    chroms2 = numpy.zeros(len(header['SQ']), dtype=numpy.int32) - 1
    for i, seq in enumerate(header['SQ']):
        chrom = seq['SN'].strip('chr')
        if chrom in chroms:
            chroms2[i] = chroms.index(chrom)
    return chroms2

def parse_cigar(cigar):
    pos = 0
    for seg in cigar:
        if seg[0] > 0:
            pos += seg[1]
        else:
            return pos

def clean_reads(reads, cuts):
    stats = {
        'missing_outer': 0, # single side mapped to one location
        'one_sided': 0, # single side mapped to two locations
        'two_sided': 0, # both sides mapped to one location each
        'valid_chimeric': 0, # both sides mapped, one or both sides mapped to two locations, valid fend pairing
        'polychimeric': 0, # more than two locations represented in mapping
    }
    polychimeric_reads = {}
    singleton_reads = {}
    invalid = {0: None, 1: None, 4: None}
    polychimeric = {3: None, 12: None, 15: None}
    onesided = {2: None, 8: None}
    for name in reads.keys():
        side1, side2 = reads[name]
        flag = min(3, len(side1)) + min(3, len(side2)) * 4
        if flag in invalid:
            # one or both sides not mapped and no chimeric sequence
            stats['missing_outer'] += 1
            if flag == 1:
                singleton_reads[name] = [side1[0][1:]]
            else:
                singleton_reads[name] = [side2[0][1:]]
            del reads[name]
        elif flag == 5:
            # single mapping on both sides
            reads[name] = [side1[0][1:], side2[0][1:]]
            stats['two_sided'] += 1
        elif flag in polychimeric:
            # one or both sides made up of 3 mappings
            stats['polychimeric'] += 1
            new_read = []
            for pos in side1:
                new_read.append(pos[1:])
            for pos in side2:
                new_read.append(pos[1:])
            polychimeric_reads[name] = new_read
            del reads[name]
        elif flag in onesided:
            # only one side mapped but did so to two locations
            if flag == 2:
                reads[name] = [side1[0][1:], side1[1][1:]]
            else:
                reads[name] = [side2[0][1:], side2[1][1:]]
            stats['one_sided'] += 1
        else:
            # both sides mapped, one or both sides mapped to two locations
            side1.sort()
            side2.sort()
            valid = True
            if len(side2) == 2:
                coord1 = match_chimeric(side1[0], side2[1], cuts)
                if coord1 is None:
                    valid = False
            else:
                coord1 = side1[0][2]
            if len(side1) == 2:
                coord2 = match_chimeric(side2[0], side1[1], cuts)
                if coord2 is None:
                    valid = False
            else:
                coord2 = side2[0][2]
            if valid:
                reads[name] = [[side1[0][1], coord1], [side2[0][1], coord2]]
                stats['valid_chimeric'] += 1
            else:
                stats['polychimeric'] += 1
                new_read = []
                for pos in side1:
                    new_read.append(pos[1:])
                for pos in side2:
                    new_read.append(pos[1:])
                polychimeric_reads[name] = new_read
                del reads[name]
    return polychimeric_reads, singleton_reads, stats

def match_chimeric(side1, side2, cuts):
    if side1[1] != side2[1]:
        return None
    fend1 = numpy.searchsorted(cuts[side1[1]], abs(side1[2]))
    fend2 = numpy.searchsorted(cuts[side1[1]], abs(side2[2]))
    if fend1 != fend2:
        return None
    return min(side1[2], side2[2])

def update_reads(reads, chroms):
    valid_reads = []
    stats = {
        'intra-chromosomal': 0,
        'inter-chromosomal': 0,
        'pcr-duplicates': 0,
    }
    for i in range(len(chroms)):
        valid_reads.append([])
        for j in range(i + 1):
            valid_reads[i].append({})
    for name, values in reads.iteritems():
        try:
            chr1, coord1 = values[0]
            chr2, coord2 = values[1]
        except:
            print values
            continue
        if chr1 > chr2:
            chr1, chr2, coord1, coord2 = chr2, chr1, coord2, coord1
        elif chr1 == chr2 and abs(coord1) > abs(coord2):
            coord1, coord2 = coord2, coord1
        key = (coord1, coord2)
        if key not in valid_reads[chr2][chr1]:
            valid_reads[chr2][chr1][key] = None
            if chr1 == chr2:
                stats['intra-chromosomal'] += 1
            else:
                stats['inter-chromosomal'] += 1
        else:
            stats['pcr-duplicates'] += 1
    return valid_reads, stats

def write_paired_reads(reads, chroms, fname):
    output = open(fname, 'w')
    for i in range(len(chroms)):
        for j in range(i + 1):
            keys = reads[i][j].keys()
            for start, stop in keys:
                if start <= 0:
                    coord1 = str(-start)
                    strand1 = '-'
                else:
                    coord1 = str(start)
                    strand1 = '+'
                if stop <= 0:
                    coord2 = str(-stop)
                    strand2 = '-'
                else:
                    coord2 = str(stop)
                    strand2 = '+'
                print >> output, "\t".join([chroms[j], coord1, strand1, chroms[i], coord2, strand2])
    output.close()

def write_other_reads(reads, chroms, fname):
    output = open(fname, 'w')
    for name, read in reads.iteritems():
        temp = []
        for site in read:
            if site is None:
                continue
            if site[1] < 0:
                temp += [chroms[site[0]], str(-site[1]), '-']
            else:
                temp += [chroms[site[0]], str(site[1]), '+']
        print >> output, '\t'.join(temp)
    output.close()

def write_stats(stats, fname):
    output = open(fname, 'w')
    keys = stats.keys()
    keys.sort()
    for key in keys:
        print >> output, "%s\t%i" % (key, stats[key])
    output.close()
    return None


if __name__ == "__main__":
    main()

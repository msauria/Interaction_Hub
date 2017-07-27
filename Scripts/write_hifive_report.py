#!/usr/bin/env python

import sys

import yaml
import hifive
import h5py
import numpy
from pyx import *


unit.set(defaultunit="cm")
text.set(mode="latex")
text.preamble(r"\usepackage{uarial}")
text.preamble(r"\renewcommand*\familydefault{\sfdefault}")
text.preamble(r"\renewcommand*\encodingdefault{T1}")


def main():
    hcd_fname, out_fname = sys.argv[1:3]
    c = canvas.canvas()
    pos = 0
    spacer = 1.0
    GSM = yaml.load(open('../YAML/GSM/%s.yml' % hcd_fname.replace('_stats.hcd','')))
    GSE = yaml.load(open('../YAML/GSE/%s.yml' % GSM['Series']))
    c1, pos1 = plot_yaml(GSM, GSE)
    c.insert(c1, [trafo.translate(0, pos)])
    pos -= pos1 + spacer
    hcd = h5py.File(hcd_fname, 'r')
    fnames = hcd.attrs['history'].split('filelist=[')[-1].split(']')[0].split(',')
    mapping_stats = {}
    for fname in fnames:
        name = fname.split('/')[-1].split('.')[0].strip(" '")
        mapping_stats[name] = {}
        for line in open("%s.stats" % name):
            temp = line.rstrip('\n').split('\t')
            mapping_stats[name][temp[0]] = temp[1]
    c1, pos1 = plot_mapping_stats(mapping_stats)
    c.insert(c1, [trafo.translate(0, pos)])
    pos -= pos1 + spacer
    hic_stats = hcd['stats'][...]
    c1, pos1 = plot_hic_stats(hic_stats)
    c.insert(c1, [trafo.translate(0, pos)])
    pos -= pos1 + spacer
    distributions = {}
    distributions['cis_interaction_distribution'] = hcd['cis_interaction_distribution'][...]
    distributions['trans_interaction_distribution'] = hcd['trans_interaction_distribution'][...]
    distributions['insert_distribution'] = [hcd['insert_distribution'][...], hcd.attrs['maxinsert']]
    distributions['outside_insert_range'] = load_outside_insert_range_dist(hcd, hcd_fname)
    c1, pos1 = plot_distributions(distributions)
    c.insert(c1, [trafo.translate(0, pos)])
    pos -= pos1 + spacer
    c.writePDFfile(out_fname)

def plot_yaml(GSM, GSE):
    height = 0
    c = canvas.canvas()
    GSM_pairs = [
        ['Series', 'GEO Experiment ID'],
        ['iid', 'GEO Sample ID'],
        [['Status', 'Submission-Date'], 'Submission Date'],
        [['Status', 'Release-Date'], 'Release Date'],
    ]
    c1, pos1 = plot_yaml_subset(GSM_pairs, GSM)
    c.insert(c1, [trafo.translate(0, -height)])
    height += pos1
    GSE_pairs = [
        [['Series', 'Pubmed-ID'], 'Pubmed ID'],
        [['Contributor', 'Organization'], 'Organization'],
        [['Contributor', 'Person', 'Last'], 'Submitters'],
    ]
    c1, pos1 = plot_yaml_subset(GSE_pairs, GSE)
    c.insert(c1, [trafo.translate(0, -height)])
    height += pos1
    GSM_pairs = [
        [['Instrument-Model', 'Predefined'], 'Platform'],
        [['Channel', 'Organism'], 'Organism'],
        [['Channel', 'Source'], 'Source'],
        ['Partition', 'Fragmentation'],
    ]
    c1, pos1 = plot_yaml_subset(GSM_pairs, GSM)
    c.insert(c1, [trafo.translate(0, -height)])
    height += pos1
    return c, height

def plot_yaml_subset(pairs, data):
    c = canvas.canvas()
    pos = 0
    for key, label in pairs:
        result = []
        if isinstance(key, str):
            if key in data:
                c.text(4, pos, "%s :" % label, [text.halign.right, text.valign.top])
                c.text(4, pos, " %s" % data[key], [text.halign.left, text.valign.top])
                pos -= 0.4
        elif isinstance(key, list) and len(key) == 2:
            if key[0] in data and isinstance(data[key[0]], dict) and key[1] in data[key[0]]:
                c.text(4, pos, "%s :" % label, [text.halign.right, text.valign.top])
                c.text(4, pos, " %s" % data[key[0]][key[1]], [text.halign.left, text.valign.top])
                pos -= 0.4
            if key[0] in data and isinstance(data[key[0]], list):
                for entry in data[key[0]]:
                    if key[1] in entry and entry[key[1]] != 'GEO':
                        result.append(entry[key[1]])
                if len(result) > 0:
                    c.text(4, pos, "%s :" % label, [text.halign.right, text.valign.top])
                    c.text(4, pos, " %s" % ', '.join(result), [text.halign.left, text.valign.top])
                    pos -= 0.4
        elif isinstance(key, list) and len(key) == 3:
            if key[0] in data :
                for entry in data[key[0]]:
                    if key[1] in entry and key[2] in entry[key[1]]:
                        result.append(entry[key[1]][key[2]])
            if len(result) > 0:
                c.text(4, pos, "%s :" % label, [text.halign.right, text.valign.top])
                c.text(4, pos, " %s" % ', '.join(result), [text.halign.left, text.valign.top])
                pos -= 0.4
    return c, -pos

def load_outside_insert_range_dist(hcd, hcd_fname):
    if 'fendfilename' in hcd.attrs:
        fendfilename = hcd.attrs['fendfilename']
        if fendfilename[:2] == './':
            fendfilename = fendfilename[2:]
        parent_count = fendfilename.count('../')
        fendfilename = '/'.join(hcd_fname.split('/')[:-(1 + parent_count)] +
                            fendfilename.lstrip('/').split('/')[parent_count:])
    fendfile = h5py.File(fendfilename, 'r')
    chroms = fendfile['chromosomes'][...]
    fends = fendfile['fends'][...]
    chr2int = {}
    for i, chrom in enumerate(chroms):
        chr2int[chrom] = i
    dist = numpy.zeros((100, 2), dtype=numpy.float64)
    chr_indices = fendfile['chr_indices'][...]
    for chrom in chroms:
        chrint = chr2int[chrom]
        if chr_indices[chrint + 1] - chr_indices[chrint] == 0:
            continue
        if 'invalid_distribution.%s' % chrom not in hcd:
            continue
        data = hcd['invalid_distribution.%s' % chrom][...]
        start = hcd['invalid_starts'][chrint]
        bounds = numpy.r_[fends['start'][chr_indices[chrint]:chr_indices[chrint + 1]:2],
                          fends['stop'][chr_indices[chrint + 1] - 1]]
        indices = numpy.searchsorted(numpy.arange(data.shape[0]) * 100 + start + 50, bounds)
        for i in range(indices.shape[0]):
            first = max(0, indices[i] - 50)
            last = min(data.shape[0], indices[i] + 50)
            foffset = first - (indices[i] - 50)
            dist[foffset:(last - first + foffset), 0] += data[first:last, 0]
            dist[foffset:(last - first + foffset), 1] += data[first:last, 1]
    return dist

def load_outside_insert_range_dist2(hcd, hcd_fname):
    if 'fendfilename' in hcd.attrs:
        fendfilename = hcd.attrs['fendfilename']
        if fendfilename[:2] == './':
            fendfilename = fendfilename[2:]
        parent_count = fendfilename.count('../')
        fendfilename = '/'.join(hcd_fname.split('/')[:-(1 + parent_count)] +
                            fendfilename.lstrip('/').split('/')[parent_count:])
    fendfile = h5py.File(fendfilename, 'r')
    chroms = fendfile['chromosomes'][...]
    fends = fendfile['fends'][...]
    chr2int = {}
    for i, chrom in enumerate(chroms):
        chr2int[chrom] = i
    dist = numpy.zeros((100, 2), dtype=numpy.float64)
    chr_indices = fendfile['chr_indices'][...]
    for chrom in chroms:
        chrint = chr2int[chrom]
        if chr_indices[chrint + 1] - chr_indices[chrint] == 0:
            continue
        if 'invalid_distribution.%s' % chrom not in hcd:
            continue
        invalid = hcd['invalid_distribution.%s' % chrom][...]
        bounds = numpy.r_[fends['start'][chr_indices[chrint]],
                                     fends['stop'][(1 + chr_indices[chrint]):chr_indices[chrint + 1]:2]]
        positions = numpy.arange(invalid.shape[0]) * 100 + hcd['invalid_starts'][chrint] + 50
        indices = numpy.searchsorted(bounds, positions, side='right') - 1
        where = numpy.where((indices >= 0) & (indices < bounds.shape[0] - 1))[0]
        starts = bounds[:-1]
        sizes = bounds[1:] - bounds[:-1]
        for i in range(2):
            fracs = numpy.minimum(numpy.floor((positions[where] - starts[indices[where]]).astype(numpy.float64) /
                                  sizes[indices[where]] * 100.).astype(numpy.int32), 99)
            dist[:, i] += numpy.bincount(fracs, weights=invalid[where, i], minlength=dist.shape[0])
    for i in range(2):
        dist[:, i] = dist[:, i] / max(1.0, float(numpy.sum(dist[:, i])))
    return dist

def plot_mapping_stats(stats):
    c = canvas.canvas()
    height = 17 * 0.4
    names = stats.keys()
    names.sort()
    c.text(0, 1.0, "Mapping statistics", [text.halign.left, text.valign.top, text.size(2)])
    for i in range((len(names) - 1) / 3 + 1):
        c.insert(plot_mapping_table(stats, names[(i * 3):min(len(names), (i + 1) * 3)]),
            [trafo.translate(0, -height * i)])
    c1 = canvas.canvas()
    c1.insert(c, [trafo.translate(0, -1.0)])
    return c1, height * ((len(names) - 1) / 3 + 1)

def plot_mapping_table(stats, names):
    labels = [
        ["total", "Total reads"],
        ["low_quality_1", "Quality too low, side 1"],
        ["unmapped_1", "Unmapped, side 1"],
        ["unused_chromosome_1", "Non-standard chromosome, side 1"],
        ["low_quality_2", "Quality too low, side 2"],
        ["unmapped_2", "Unmapped, side 2"],
        ["unused_chromosome_2", "Non-standard chromosome, side 2"],
        ["missing_outer", "Non-hybrid read"],
        ["one_sided", "One-sided hybrid read"],
        ["two_sided", "Two-sided hybrid reads"],
        ["valid_chimeric", "Valid sequenced ligations"],
        ["polychimeric", "More than 2 fragments"],
        ["pcr-duplicates", "PCR duplicates"],
        ["intra-chromosomal", "Valid cis reads"],
        ["inter-chromosomal", "Valid trans reads"],
    ]
    c = canvas.canvas()
    width = 5.8 + 4.0 * len(names)
    height = 16 * 0.4
    for i in range(len(names)):
        c.text(7.6 + 4.0 * i, 0.75, names[i], [text.halign.center, text.valign.top])
    c.stroke(path.line(0, 0.4, width, 0.4))
    for i in range(len(labels)):
        c.text(5.75, 0.35 - 0.4 * i, labels[i][1], [text.halign.right, text.valign.top])
        for j in range(len(names)):
            num = stats[names[j]][labels[i][0]]
            total = int(stats[names[j]]['total'])
            c.text(7.8 + 4.0 * j, 0.35 - 0.4 * i, "%s (%0.1f\%%)" % (reformat_num(num), float(num) / total * 100.),
                [text.halign.center, text.valign.top])
            c.stroke(path.line(0, -0.4 * i, width, -0.4 * i))
    for i in range(len(names) + 1):
        c.stroke(path.line(5.8 + i * 4.0, 0.8, 5.8 + i * 4.0, 0.8 - height))
    return c

def plot_hic_stats(stats):
    labels = [
        ["total_reads", "Total reads"],
        ["chr_not_in_fends", "Invalid chromosome name"],
        ["out_of_bounds", "Outside coordinate range"],
        ["pcr_duplicates", "PCR duplicates"],
        ["insert_size", "Insert size too large"],
        ["same_fragment", "Circularized reads"],
        ["valid_cis_reads", "Valid cis reads"],
        ["valid_trans_reads", "Valid trans reads"],
        ["valid_cis_pairs", "Unique cis fragment pairs"],
        ["valid_trans_pairs", "Unique trans fragment pairs"],
    ]
    c = canvas.canvas()
    width = 6.0 + 4.0
    height = 10 * 0.4
    c.text(0, 1.0, "HiFive statistics", [text.halign.left, text.valign.top, text.size(2)])
    c.stroke(path.line(0, 0.4, width, 0.4))
    index = numpy.where(stats['name'] == 'total_reads')[0][0]
    total = stats['count'][index]
    for i in range(len(labels)):
        c.text(5.95, 0.35 - 0.4 * i, labels[i][1], [text.halign.right, text.valign.top])
        index = numpy.where(stats['name'] == labels[i][0])[0][0]
        num = str(stats['count'][index])
        c.text(8, 0.35 - 0.4 * i, "%s (%0.1f\%%)" % (reformat_num(num), float(num) / total * 100.),
            [text.halign.center, text.valign.top])
        c.stroke(path.line(0, -0.4 * i, width, -0.4 * i))
    for i in range(2):
        c.stroke(path.line(6.0 + i * 4.0, 0.4, 6.0 + i * 4.0, 0.4 - height))
    c1 = canvas.canvas()
    c1.insert(c, [trafo.translate(0, -0.5)])
    return c1, height + 1.0

def plot_distributions(distributions):
    height = 0
    c = canvas.canvas()
    c.insert(plot_interaction_dist(distributions['cis_interaction_distribution'], 'Cis'),
        [trafo.translate(0, 0.0)])
    height += 6.5
    c.insert(plot_interaction_dist(distributions['trans_interaction_distribution'], 'Trans'),
        [trafo.translate(0, -6.5)])
    height += 6.5
    c.insert(plot_insert_dist(distributions['insert_distribution']),
        [trafo.translate(0, -13.0)])
    height += 6.5
    if numpy.sum(distributions['outside_insert_range']) > 0:
        c.insert(plot_invalid_insert_dist(distributions['outside_insert_range']),
            [trafo.translate(0, -19.5)])
    height += 6.5
    c1 = canvas.canvas()
    c1.insert(c, [trafo.translate(1.5, -4.0)])
    return c1, height

def plot_interaction_dist(dist, label):
    c = canvas.canvas()
    height = 4.0
    width = 16.0
    maxbin = numpy.where(dist > 0)[0][-1]
    rawX = numpy.log10(numpy.maximum(0.1, numpy.arange(dist.shape[0]))).astype(numpy.float64)[:maxbin]
    Y = numpy.log10(numpy.maximum(0.1, dist[:maxbin]).astype(numpy.float64))
    minX = rawX[0]
    maxX = rawX[-1]
    X = (rawX - minX) / (maxX - minX) * width
    minY = numpy.amin(Y)
    maxY = numpy.amax(Y) * 1.05
    Y -= minY
    Y /= (maxY - minY) / height
    lpath = path.path(path.moveto(0, 0))
    for i in range(X.shape[0]):
        lpath.append(path.lineto(X[i], Y[i]))
    c.stroke(lpath)
    c.stroke(path.rect(0, 0, width, height))
    for i in range(-1, int(numpy.floor(maxX)) + 1):
        if i < 0:
            pos = 0
            val = "0"
        else:
            pos = (i - minX) / (maxX - minX) * width
            val = r"10\textsuperscript{%i}" % i
        c.stroke(path.line(pos, 0, pos, -0.15))
        c.text(pos, -0.2, val, [text.halign.center, text.valign.top])
    c.text(width * 0.5, -0.6, "Number of interactions", [text.halign.center, text.valign.top])
    for i in range(-1, int(numpy.floor(maxY)) + 1):
        if i < 0:
            pos = 0
            val = '0'
        else:
            pos = (i - minY) / (maxY - minY) * height
            val = r"10\textsuperscript{%i}" % i
        c.stroke(path.line(0, pos, -0.15, pos))
        c.text(-0.2, pos, val, [text.halign.right, text.valign.middle])
    c.text(-1.5, height * 0.5, "Number of Fragments", [text.halign.center, text.valign.top, trafo.rotate(90)])
    c.text(width * 0.5, 0.6 + height, "%s Interaction Distribution" % label, [text.halign.center, text.valign.top, text.size(2)])
    return c

def plot_insert_dist(data):
    dist, maxinsert = data
    c = canvas.canvas()
    height = 4.0
    width = 16.0
    dist[0, 1] = 1
    X = numpy.log10(dist[:, 1])
    X = (X[1:] + X[:-1]) / 2.0
    Y = dist[:-1, 0].astype(numpy.float64)
    minX = 0
    maxX = X[-1]
    X = (X - minX) / (maxX - minX) * width
    maxY = numpy.amax(Y) * 1.05
    Y /= maxY / height
    lpath = path.path(path.moveto(0, 0))
    for i in range(X.shape[0]):
        lpath.append(path.lineto(X[i], Y[i]))
    c.stroke(lpath)
    c.stroke(path.rect(0, 0, width, height))
    for i in range(int(numpy.floor(maxX) + 1)):
        pos = i / maxX * width
        c.stroke(path.line(pos, 0, pos, -0.15))
        c.text(pos, -0.2, r"10\textsuperscript{%i}" % i, [text.halign.center, text.valign.top])
    c.text(width * 0.5, -0.6, "Insert size", [text.halign.center, text.valign.top])
    N = len(str(numpy.amax(dist[:, 0]))) - 1
    for i in range(0, int(numpy.floor(maxY / 10**N) + 1), 2):
        pos = (i * 10**N) / maxY * height
        c.stroke(path.line(0, pos, -0.15, pos))
        c.text(-0.2, pos, r"%ix10\textsuperscript{%i}" % (i, N), [text.halign.right, text.valign.middle])
    c.text(-1.5, height * 0.5, "Number of Reads", [text.halign.center, text.valign.top, trafo.rotate(90)])
    pos = numpy.log10(maxinsert) / maxX * width
    c.stroke(path.line(pos, 0, pos, height), [color.rgb.red])
    c.text(width - 0.1, height - 0.1, "Max insert: %i" % maxinsert, [color.rgb.red, text.halign.right, text.valign.top])
    c.text(width * 0.5, 0.6 + height, "Insert Size Distribution", [text.halign.center, text.valign.top, text.size(2)])
    return c

def plot_invalid_insert_dist(dist):
    c = canvas.canvas()
    height = 4.0
    width = 16.0
    X = numpy.arange(-50, 51) * 100.
    Y0 = dist[:, 0].astype(numpy.float64)
    Y1 = dist[:, 1].astype(numpy.float64)
    minX = X[0]
    maxX = X[-1]
    maxY = numpy.amax(dist) * 1.05
    Y0 /= maxY / height
    Y1 /= maxY / height
    lpath0 = None
    lpath1 = None
    for i in range(X.shape[0] - 1):
        x = ((X[i] + X[i + 1]) / 2 - minX) / (maxX - minX) * width
        if lpath0 is None:
            lpath0 = path.path(path.moveto(x, Y0[i]))
        else:
            lpath0.append(path.lineto(x, Y0[i]))
        if lpath1 is None:
            lpath1 = path.path(path.moveto(x, Y1[i]))
        else:
            lpath1.append(path.lineto(x, Y1[i]))
    c.stroke(lpath0, [color.rgb.black])
    c.stroke(lpath1, [color.rgb.red])
    c.stroke(path.rect(0, 0, width, height))
    c.text(width - 0.1, height - 0.1, "Forward Strand", [text.halign.right, text.valign.top])
    c.text(width - 0.1, height - 0.5, "Reverse Strand", [text.halign.right, text.valign.top, color.rgb.red])
    for i in range(int(minX), int(maxX) + 1000, 1000):
        pos = (i - minX) / (maxX - minX) * width
        c.stroke(path.line(pos, 0, pos, -0.15))
        c.text(pos, -0.2, r"%i Kb" % (i/1000), [text.halign.center, text.valign.top])
    c.text(width * 0.5, -0.6, "Read end distance from cutsite", [text.halign.center, text.valign.top])
    N = len(str(int(numpy.floor(numpy.amax(dist))))) - 1
    for i in range(int(numpy.floor(maxY / 10**N) + 1)):
        pos = (i * 10**N) / maxY * height
        c.stroke(path.line(0, pos, -0.15, pos))
        c.text(-0.2, pos, r"%ix10\textsuperscript{%i}" % (i, N), [text.halign.right, text.valign.middle])
    c.text(-1.5, height * 0.5, "Number of Reads", [text.halign.center, text.valign.top, trafo.rotate(90)])
    c.text(width * 0.5, 0.6 + height, "Invalid Insert End Position Distribution", [text.halign.center, text.valign.top, text.size(2)])
    return c

def plot_invalid_insert_dist2(dist):
    c = canvas.canvas()
    height = 4.0
    width = 16.0
    X = numpy.arange(100).astype(numpy.float64) + 0.5
    Y0 = dist[:, 0].astype(numpy.float64)
    Y1 = dist[:, 1].astype(numpy.float64)
    minX = 0
    maxX = 100.0
    X = (X - minX) / (maxX - minX) * width
    maxY = numpy.amax(dist) * 1.05
    Y0 /= maxY / height
    Y1 /= maxY / height
    lpath0 = path.path(path.moveto(0, 0))
    lpath1 = path.path(path.moveto(0, 0))
    for i in range(X.shape[0]):
        lpath0.append(path.lineto(X[i], Y0[i]))
        lpath1.append(path.lineto(X[i], Y1[i]))
    lpath0.append(path.lineto(width, 0))
    lpath1.append(path.lineto(width, 0))
    c.stroke(lpath0, [color.rgb.black])
    c.stroke(lpath1, [color.rgb.red])
    c.stroke(path.rect(0, 0, width, height))
    for i in range(0, 101, 10):
        pos = i / maxX * width
        c.stroke(path.line(pos, 0, pos, -0.15))
        c.text(pos, -0.2, r"$%i\%%$" % i, [text.halign.center, text.valign.top])
    c.text(width * 0.5, -0.6, "Fragment position", [text.halign.center, text.valign.top])
    N = len(str(numpy.amax(dist))) - 1
    for i in range(int(numpy.floor(maxY / 10**N) + 1)):
        pos = (i * 10**N) / maxY * height
        c.stroke(path.line(0, pos, -0.15, pos))
        c.text(-0.2, pos, "$%i$x$10^{%i}$" % (i, N), [text.halign.right, text.valign.middle])
    c.text(-1.5, height * 0.5, "Number of Reads", [text.halign.center, text.valign.top, trafo.rotate(90)])
    c.text(width * 0.5, 0.6 + height, "Invalid Insert Position Distribution", [text.halign.center, text.valign.top, text.size(2)])
    return c

def reformat_num(s):
    s1 = s[-3:]
    for i in range((len(s) - 1) / 3):
        s1 = "%s,%s" % (s[max(-len(s), -3*(i + 2)):(-3 * (i + 1))], s1)
    return s1



if __name__ == "__main__":
    #try:
    main()
    #except:
    #    print sys.argv[1]

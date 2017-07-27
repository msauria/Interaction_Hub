#!/usr/bin/env python

import os

import yaml

base_dir = "/Users/msauria/projects/repos/HiC_Hosting/YAML"

def main():
    master = load_yaml_file("%s/all_GSEs.yml" % base_dir)
    GSEs = load_GSEs(master)
    GSMs = load_GSMs(GSEs)
    print GSMs.keys()
    data = {}
    for GSM in GSMs:
        genome = GSMs[GSM]['Channel']['Organism']
        if isinstance(genome, list):
            continue
        genome = genome.encode('utf-8')
        celltype = GSMs[GSM]['Channel']['Source'].encode('utf-8')
        cutter = GSMs[GSM]['Partition'].encode('utf-8')
        if genome not in data:
            data[genome] = {}
        if celltype not in data[genome]:
            data[genome][celltype] = {}
        if cutter not in data[genome][celltype]:
            data[genome][celltype][cutter] = []
        data[genome][celltype][cutter].append(GSM)
    genomes = data.keys()
    genomes.sort()
    for genome in genomes:
        print genome
        celltypes = data[genome].keys()
        celltypes.sort()
        for celltype in celltypes:
            print "  %s" % celltype
            cutters = data[genome][celltype].keys()
            cutters.sort()
            for cutter in cutters:
                data[genome][celltype][cutter].sort()
                print "    %s %i %s" % (cutter, len(data[genome][celltype][cutter]), ','.join(data[genome][celltype][cutter]))


def load_yaml_file(fname):
    if os.path.exists(fname):
        data = yaml.load(open(fname))
    else:
        data = {}
    return data

def load_GSEs(master):
    GSEs = {}
    for GSE, ID in master.iteritems():
        data = load_yaml_file("%s/GSE/%s.yml" % (base_dir, GSE))
        if len(data) > 0:
            GSEs[GSE] = data
    return GSEs

def load_GSMs(GSEs):
    GSMs = {}
    for key, GSE in GSEs.iteritems():
        for sample in GSE['Series']['Sample-Ref']:
            data = load_yaml_file("%s/GSM/%s.yml" % (base_dir, sample['ref']))
            if len(data) > 0:
                GSMs[sample['ref']] = data
    return GSMs


if __name__ == "__main__":
    main()

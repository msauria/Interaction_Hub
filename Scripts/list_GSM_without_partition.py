#!/usr/bin/env python

import sys
import os

import yaml

base_dir = "/Users/msauria/projects/repos/HiC_Hosting/YAML"

def main():
    master = load_yaml_file("%s/all_GSEs.yml" % base_dir)
    GSEs = load_GSEs(master)
    GSMs = load_GSMs(GSEs)
    missing = []
    for GSM in GSMs:
        if GSMs[GSM]['Partition'] == '?':
            missing.append(GSM)
    missing.sort()
    for line in missing:
        print line

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
        print >> sys.stderr, (key),
        for sample in GSE['Series']['Sample-Ref']:
            data = load_yaml_file("%s/GSM/%s.yml" % (base_dir, sample['ref']))
            if len(data) > 0:
                GSMs[sample['ref']] = data
    return GSMs


if __name__ == "__main__":
    main()

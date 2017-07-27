#!/usr/bin/env python

import os
import sys

import yaml

base_dir = "/Users/msauria/projects/repos/HiC_Hosting/YAML"

def main():
    start, stop, re = sys.argv[1:4]
    start_int = int(start.split('GSM')[1].split('.')[0])
    stop_int = int(stop.split('GSM')[1].split('.')[0]) + 1
    for i in range(start_int, stop_int):
        name = "GSM%i" % i
        fname = "%s/GSM/%s.yml" % (base_dir, name)
        GSM = load_yaml_file(fname)
        if len(GSM) == 0:
            continue
        GSM['Partition'] = re
        write_yaml(GSM, fname)

def write_yaml(data, fname):
    with open(fname, 'w') as f:
        yaml.dump(data, f, default_flow_style=False)

def load_yaml_file(fname):
    if os.path.exists(fname):
        data = yaml.load(open(fname))
    else:
        data = {}
    return data


if __name__ == "__main__":
    main()

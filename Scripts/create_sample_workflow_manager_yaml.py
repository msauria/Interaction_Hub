#!/usr/bin/env python

import sys
import os
import copy

import yaml

bwa_index_dir = "/cache/genomes/bwa"
output_dir = "/project/hic_database/Data/Processed"

def main():
    in_fname, out_fname = sys.argv[1:3]
    data = load_yaml_file(in_fname)
    if len(data) == 0:
        return None
    workflow = create_workflow_yaml(data)
    if workflow is not None:
        jobs_yaml = create_jobs_yaml(workflow)
        if jobs_yaml is not None:
            write_yaml(jobs_yaml, out_fname)

def create_workflow_yaml(data):
    workflow = {
        "ftp": [],
        "idxbase": {"class": "File", "path": ""},
        "partition": {"class": "File", "path": ""},
        "hcd": {"class": "File", "path": ""},
        "hcp": {"class": "File", "path": ""},
        "normalized": {"class": "File", "path": ""},
        "hdf5_file": {"class": "File", "path": ""},
        "hdf5": "",
        "name": "",
        "minbin": 0,
        "numbins": 0,
        "accessions": [],
        "reads": [],
    }
    for SRR in data['SRA']['SRR']:
        if 'accession' in SRR:
            accession = SRR['accession']
        elif 'acc' in SRR:
            accession = SRR['acc']
        workflow['ftp'].append("%s/%s/%s.sra" % (data['SRA']['ftp'], accession, accession))
        workflow['accessions'].append(accession)
        workflow['reads'].append(SRR['total_spots'])
    genomes = {
        'Homo sapiens': 'hg38',
        'Mus musculus': 'mm10',
        'Drosophila melanogaster': 'dm6',
    }
    if data['Channel']['Organism'] not in genomes:
        return None
    genome = genomes[data['Channel']['Organism']]
    workflow['idxbase']['path'] = "%s/%s/%s.fa" % (bwa_index_dir, genome, genome)
    workflow['partition']['path'] = "%s/%s_%s.fends" % (output_dir, genome, data['Partition'])
    minbins = {
        'HindIII': 1000,
        'NcoI': 1000,
        'DpnII': 250,
        'MboI': 250,
        'MluI': 1000,
        'BglII': 1000,
        'CviJI': 250,
        'MspI': 250,
        'HinfI': 250,
        'NlaIII': 250,
        'MNase': 250,
        'DNase': 1000,
    }
    workflow['minbin'] = minbins[data['Partition']]
    numbins = {
        'Homo sapiens': 75,
        'Mus musculus': 75,
        'Drosophila melanogaster': 50,
    }
    workflow['numbins'] = numbins[data['Channel']['Organism']]
    workflow['name'] = data['Accession']
    workflow['hcd']['path'] = "%s/%s.hcd" % (output_dir, data['Accession'])
    workflow['hcp']['path'] = "%s/%s.hcp" % (output_dir, data['Accession'])
    workflow['normalized']['path'] = "%s/%s_bin.hcp" % (output_dir, data['Accession'])
    workflow['hdf5'] = "%s_HM.hdf5" % (data['Accession'])
    workflow['hdf5_file']['path'] = "%s/%s_HM.hdf5" % (output_dir, data['Accession'])
    chroms = {
        'Homo sapiens': "1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X",
        'Mus musculus': "1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,X",
        'Drosophila melanogaster': "2L,2R,3L,3R,4,X",
    }
    workflow['chroms'] = chroms[data['Channel']['Organism']]
    return workflow

def create_jobs_yaml(workflow):
    job = {}
    job['prereqs'] = [hcp_job(workflow)]
    job['cpu'] = 27
    job['disk'] = 0
    job['mem'] = 4.0
    job['script'] = '/project/hic_database/Scripts/CWL/hifive-hic_normalize_binning-nodocker.cwl'
    job['id'] = "%s_norm" % workflow['name']
    job['inputs'] = ["/project/hic_database/Data/Processed/%s.hcp" % workflow['name'],
                     "/project/hic_database/Data/Processed/%s.hcd" % workflow['name']]
    job['outputs'] = ["/project/hic_database/Data/Processed/%s_bin.hcp" % workflow['name']]
    job['yaml'] = {
        'name': workflow['name'],
        'hcp': copy.deepcopy(workflow['hcp']),
        'hcd': copy.deepcopy(workflow['hcd']),
        'partition': copy.deepcopy(workflow['partition']),
        'chroms': workflow['chroms']
    }
    return job

def hcp_job(workflow):
    job = {}
    job['prereqs'] = [hcd_job(workflow)]
    job['cpu'] = 27
    job['disk'] = 0
    job['mem'] = 1.0
    job['script'] = '/project/hic_database/Scripts/CWL/hifive-hcd2hcp-nodocker.cwl'
    job['id'] = "%s_project" % workflow['name']
    job['inputs'] = ["/project/hic_database/Data/Processed/%s.hcd" % workflow['name']]
    job['outputs'] = ["/project/hic_database/Data/Processed/%s.hcp" % workflow['name']]
    job['yaml'] = {
        'minbin': workflow['minbin'],
        'numbins': workflow['numbins'],
        'name': workflow['name'],
        'hcd': copy.deepcopy(workflow['hcd']),
        'partition': copy.deepcopy(workflow['partition']),
    }
    return job

def hcd_job(workflow):
    job = {}
    job['prereqs'] = []
    for i in range(len(workflow['ftp'])):
        job['prereqs'].append(raw_job(workflow, i))
    job['cpu'] = 1
    job['disk'] = 0
    job['mem'] = 8.0
    job['script'] = '/project/hic_database/Scripts/CWL/hifive-raw2hcd.cwl'
    job['id'] = "%s_hcd" % workflow['name']
    job['inputs'] = []
    for i in range(len(workflow['accessions'])):
        job['inputs'].append("/project/hic_database/Data/Processed/%s.paired" % workflow['accessions'][i])
    job['outputs'] = ["/project/hic_database/Data/Processed/%s.hcd" % workflow['name']]
    job['yaml'] = {
        'name': workflow['name'],
        'partition': workflow['partition'],
        'raw': []
    }
    for key in workflow['accessions']:
        job['yaml']['raw'].append({
            'paired': {'class': 'File', 'path': "/project/hic_database/Data/Processed/%s.paired" % key},
            'singletons': {'class': 'File', 'path': "/project/hic_database/Data/Processed/%s.singletons" % key},
            'polychimeric': {'class': 'File', 'path': "/project/hic_database/Data/Processed/%s.polychimeric" % key},
            'stats': {'class': 'File', 'path': "/project/hic_database/Data/Processed/%s.stats" % key},
            })
    return job

def raw_job(workflow, index):
    job = {}
    job['prereqs'] = [bam_job(workflow, index, 1), bam_job(workflow, index, 2)]
    job['cpu'] = 1
    job['disk'] = 0
    job['mem'] = int(workflow['reads'][index]) / 1000000.
    job['script'] = '/project/hic_database/Scripts/CWL/python-bam2raw.cwl'
    job['id'] = "%s_raw" % workflow['accessions'][index]
    job['inputs'] = ["/project/hic_database/Data/Processed/%s_1.bam" % workflow['accessions'][index],
                     "/project/hic_database/Data/Processed/%s_2.bam" % workflow['accessions'][index]]
    job['outputs'] = ["/project/hic_database/Data/Processed/%s.paired" % workflow['accessions'][index],
                      "/project/hic_database/Data/Processed/%s.singletons" % workflow['accessions'][index],
                      "/project/hic_database/Data/Processed/%s.polychimeric" % workflow['accessions'][index],
                      "/project/hic_database/Data/Processed/%s.stats" % workflow['accessions'][index]]
    job['yaml'] = {
        'script': {"class": "File", "path": "/project/hic_database/Scripts/Python/bam2raw.py"},
        'reads': [
        {"class": "File", "path": "/project/hic_database/Data/Processed/%s_1.bam" % workflow['accessions'][index]},
        {"class": "File", "path": "/project/hic_database/Data/Processed/%s_2.bam" % workflow['accessions'][index]}],
        'partition': copy.deepcopy(workflow['partition']),
    }
    return job

def bam_job(workflow, index, side):
    job = {}
    job['prereqs'] = [sam_job(workflow, index, side)]
    job['cpu'] = 1
    job['disk'] = 1
    job['mem'] = 1.0
    job['script'] = '/project/hic_database/Scripts/CWL/samtools-sam2bam.cwl'
    job['id'] = "%s_bam_%i" % (workflow['accessions'][index], side)
    job['inputs'] = ["/project/hic_database/Data/Processed/%s_%i.sam" % (workflow['accessions'][index], side)]
    job['outputs'] = ["/project/hic_database/Data/Processed/%s_%i.bam" % (workflow['accessions'][index], side)]
    job['yaml'] = {
        'input': {"class": "File",
                  "path": "/project/hic_database/Data/Processed/%s_%i.sam" % (workflow['accessions'][index], side)},
    }
    return job

def sam_job(workflow, index, side):
    job = {}
    if side == 1:
        job['prereqs'] = [fastq_job(workflow, index)]
    else:
        job['prereqs'] = []
    job['cpu'] = 27
    job['disk'] = 0
    job['mem'] = 2.0 * 27
    job['script'] = '/project/hic_database/Scripts/CWL/bwa-mem.cwl'
    job['id'] = "%s_sam_%i" % (workflow['accessions'][index], side)
    job['inputs'] = ["/project/hic_database/Data/Processed/%s_%i.fastq" % (workflow['accessions'][index], side)]
    job['outputs'] = ["/project/hic_database/Data/Processed/%s_%i.sam" % (workflow['accessions'][index], side)]
    job['yaml'] = {
        'reads': {"class": "File",
                  "path": "/project/hic_database/Data/Processed/%s_%i.fastq" % (workflow['accessions'][index], side)},
        'idxbase': copy.deepcopy(workflow['idxbase']),
    }
    return job

def fastq_job(workflow, index):
    job = {}
    job['prereqs'] = [sra_job(workflow, index)]
    job['cpu'] = 1
    job['disk'] = 0
    job['mem'] = 1.0
    job['script'] = '/project/hic_database/Scripts/CWL/sratools-fastq_dump.cwl'
    job['id'] = "%s_fastq" % (workflow['accessions'][index])
    job['inputs'] = ["/project/hic_database/Data/Processed/%s.sra" % (workflow['accessions'][index])]
    job['outputs'] = ["/project/hic_database/Data/Processed/%s_1.fastq" % workflow['accessions'][index],
                      "/project/hic_database/Data/Processed/%s_2.fastq" % workflow['accessions'][index],]
    job['yaml'] = {
        'sra': {"class": "File",
                "path": "/project/hic_database/Data/Processed/%s.sra" % workflow['accessions'][index]},
    }
    return job

def sra_job(workflow, index):
    job = {}
    job['prereqs'] = []
    job['cpu'] = 1
    job['disk'] = 1
    job['mem'] = 1.0
    job['script'] = '/project/hic_database/Scripts/CWL/util-wget_sra.cwl'
    job['id'] = "%s_sra" % (workflow['accessions'][index])
    job['inputs'] = []
    job['outputs'] = ["/project/hic_database/Data/Processed/%s.sra" % workflow['accessions'][index]]
    job['yaml'] = {
        'ftp': workflow['ftp'][index],
    }
    return job

def load_yaml_file(fname):
    if os.path.exists(fname):
        data = yaml.load(open(fname))
    else:
        data = {}
    return data

def write_yaml(data, fname):
    with open(fname, 'w') as f:
        yaml.dump(data, f, default_flow_style=False)


if __name__ == "__main__":
    main()
#!/usr/bin/env python

import sys
import os
import subprocess
from math import ceil

import yaml

bwa_dir = "/home-3/msauria1@jhu.edu/work/cache/bwa_indices"
working_dir = "/home-3/msauria1@jhu.edu/scratch/HiCDatabase/Data"
scripts_dir = "/home-3/msauria1@jhu.edu/scratch/HiCDatabase/Scripts"
tmp_dir = "/home-3/msauria1@jhu.edu/scratch/HiCDatabase/tmp"

def main():
    in_fname, out_fname = sys.argv[1:3]
    data = load_yaml_file(in_fname)
    if len(data) == 0:
        return None
    workflow = create_workflow_yaml(data)
    log = submit_job(workflow)
    write_yaml(log, out_fname)

def create_workflow_yaml(data):
    genomes = {
        'Homo sapiens': 'hg38',
        'Mus musculus': 'mm10',
        'Drosophila melanogaster': 'dm6',
    }
    chroms = {
        'hg38': "1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X",
        'mm10': "1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,X",
        'dm6': "2L,2R,3L,3R,4,X",
    }
    numbins = {
        'hg38': 75,
        'mm10': 75,
        'dm6': 50,
    }
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
    data['genome'] = genomes[data['Channel']['Organism']]
    data['chroms'] = chroms[data['genome']]
    data['numbin'] = numbins[data['genome']]
    data['minbin'] = minbins[data['Partition']]
    data['fendfile'] = "%s/%s_%s.fends" % (working_dir, data['genome'], data['Partition'])
    workflow = create_normalization_job(data)
    return workflow

def create_normalization_job(data):
    job = {
        'fname': "%s/%s-norm.sh" % (tmp_dir, data['iid']),
        'task': "mpirun -np 24 hifive hic-normalize binning -m 500000 -v len,GC,mappability -s 20,20,10 -u even,even,fixed-const -c %s -o %s/%s_bin.hcp %s/%s.hcp" % (data['chroms'], working_dir, data['iid'], working_dir, data['iid']),
        'ntasks': 24,
        'name': "%s-norm" % (data['iid']),
        'partition': "shared,lrgmem,parallel",
        'options': [],
        'modules': [],
        'time': "1:0:0",
    }
    if data['genome'] != 'dm6':
        job['time'] = '4:0:0'
    prereqs = [create_project_job(data)]
    outputs = ['%s/%s_bin.hcp' % (working_dir, data['iid'])]
    return {'job': job, 'prereqs': prereqs, 'outputs': outputs}

def create_project_job(data):
    job = {
        'fname': "%s/%s-hcp.sh" % (tmp_dir, data['iid']),
        'task': "mpirun -np 24 hifive hic-project -f 1 -m 500000 -j %i -n %i %s/%s.hcd %s/%s.hcp" % (data['minbin'], data['numbin'], working_dir, data['iid'], working_dir, data['iid']),
        'ntasks': 24,
        'name': "%s-hcp" % (data['iid']),
        'partition': "shared,lrgmem,parallel",
        'options': [],
        'modules': [],
        'time': "1:0:0",
    }
    if data['genome'] != 'dm6':
        job['time'] = '4:0:0'
    prereqs = [create_data_job(data)]
    outputs = ['%s/%s.hcp' % (working_dir, data['iid'])]
    return {'job': job, 'prereqs': prereqs, 'outputs': outputs}

def create_data_job(data):
    prereqs = []
    raws = []
    num_reads = 0
    for SRRs in data['SRA']:
        for SRR in SRRs['SRR']:
            SRRid = SRR['acc']
            prereqs.append(create_raw_job(data, SRRid))
            raws.append("%s/%s.paired" % (working_dir, SRRid))
            num_reads += int(SRR['total_spots'])
    ntasks = min(int(ceil(float(num_reads) / 10000000.)), 48)
    time = min(72, int(ceil(float(num_reads) / 10000000.)))
    job = {
        'fname': "%s/%s-hcd.sh" % (tmp_dir, data['iid']),
        'task': "hifive hic-data --skip-duplicate-filtering -i 650 -R %s %s %s/%s.hcd" % (" -R ".join(raws), data['fendfile'], working_dir, data['iid']),
        'ntasks': ntasks,
        'name': "%s-hcd" % (data['iid']),
        'partition': "lrgmem",
        'options': [],
        'modules': [],
        'time': "%i:0:0" % time,
    }
    outputs = ['%s/%s.hcd' % (working_dir, data['iid'])]
    return {'job': job, 'prereqs': prereqs, 'outputs': outputs}

def create_raw_job(data, SRRid):
    prereqs = []
    for SRRs in data['SRA']:
        for SRR in SRRs['SRR']:
            if SRRid != SRR['acc']:
                continue
            prereqs.append(create_bam_job(data, SRRid, 1))
            prereqs.append(create_bam_job(data, SRRid, 2))
            num_reads = int(SRR['total_spots'])
    ntasks = min(int(ceil(float(num_reads) / 10000000.)), 48)
    time = min(72, int(ceil(float(num_reads) / 10000000.)))
    job = {
        'fname': "%s/%s-raw.sh" % (tmp_dir, SRRid),
        'task': "python %s/bam2raw.py -b %s/%s_1.bam %s/%s_2.bam -f %s -o %s/%s" % (scripts_dir, working_dir, SRRid, working_dir, SRRid, data['fendfile'], working_dir, SRRid),
        'ntasks': ntasks,
        'name': "%s-raw" % (SRRid),
        'partition': "lrgmem",
        'options': [],
        'modules': [],
        'time': "%i:0:0" % time,
    }
    outputs = ['%s/%s.paired' % (working_dir, SRRid)]
    return {'job': job, 'prereqs': prereqs, 'outputs': outputs}

def create_bam_job(data, SRRid, side):
    for SRRs in data['SRA']:
        for SRR in SRRs['SRR']:
            if SRRid != SRR['acc']:
                continue
            num_reads = int(SRR['total_spots'])
    prereqs = [create_sam_job(data, SRRid, side)]
    time = min(72, int(ceil(float(num_reads) / 20000000.)))
    job = {
        'fname': "%s/%s-bam%i.sh" % (tmp_dir, SRRid, side),
        'task': "samtools view -b %s/%s_%i.sam > %s/%s_%i.bam" % (working_dir, SRRid, side, working_dir, SRRid, side),
        'ntasks': 1,
        'name': "%s-bam%i" % (SRRid, side),
        'partition': "shared,lrgmem",
        'options': [],
        'modules': ['samtools'],
        'time': "%i:0:0" % time,
    }
    outputs = ['%s/%s_%i.bam' % (working_dir, SRRid, side)]
    return {'job': job, 'prereqs': prereqs, 'outputs': outputs}

def create_sam_job(data, SRRid, side):
    for SRRs in data['SRA']:
        for SRR in SRRs['SRR']:
            if SRRid != SRR['acc']:
                continue
            num_reads = int(SRR['total_spots'])
    if side == 1:
        prereqs = [create_fastq_job(data, SRRid)]
    else:
        prereqs = ['%s-fastq' % SRRid]
    time = min(72, int(ceil(float(num_reads) / 10000000.)))
    job = {
        'fname': "%s/%s-sam%i.sh" % (tmp_dir, SRRid, side),
        'task': "bwa mem -T 30 -M -t 24 %s/%s/%s %s/%s_%i.fastq > %s/%s_%i.sam" % (bwa_dir, data['genome'], data['genome'], working_dir, SRRid, side, working_dir, SRRid, side),
        'ntasks': 24,
        'name': "%s-sam%i" % (SRRid, side),
        'partition': "shared,lrgmem,parallel",
        'options': [],
        'modules': ['bwa'],
        'time': "%i:0:0" % time,
    }
    outputs = ['%s/%s_%i.sam' % (working_dir, SRRid, side)]
    return {'job': job, 'prereqs': prereqs, 'outputs': outputs}

def create_fastq_job(data, SRRid):
    for SRRs in data['SRA']:
        for SRR in SRRs['SRR']:
            if SRRid != SRR['acc']:
                continue
            num_reads = int(SRR['total_spots'])
    prereqs = [create_sra_job(data, SRRid)]
    time = min(72, int(ceil(float(num_reads) / 50000000.)))
    job = {
        'fname': "%s/%s-fastq.sh" % (tmp_dir, SRRid),
        'task': "fastq-dump --skip-technical --dumpbase --split-files --clip -O %s %s/%s.sra" % (working_dir, working_dir, SRRid),
        'ntasks': 1,
        'name': "%s-fastq" % (SRRid),
        'partition': "shared,lrgmem",
        'options': [],
        'modules': ['sratoolkit'],
        'time': "%i:0:0" % time,
    }
    outputs = ['%s/%s_1.fastq' % (working_dir, SRRid), '%s/%s_2.fastq' % (working_dir, SRRid)]
    return {'job': job, 'prereqs': prereqs, 'outputs': outputs}

def create_sra_job(data, SRRid):
    prereqs = []
    for SRRs in data['SRA']:
        for SRR in SRRs['SRR']:
            if SRRid != SRR['acc']:
                continue
            ftp = SRRs['ftp']
    job = {
        'fname': "%s/%s-sra.sh" % (tmp_dir, SRRid),
        'task': "wget -q %s/%s/%s.sra -O %s/%s.sra" % (ftp, SRRid, SRRid, working_dir, SRRid),
        'ntasks': 1,
        'name': "%s-sra" % (SRRid),
        'partition': "shared,lrgmem",
        'options': [],
        'modules': [],
        'time': "1:0:0",
    }
    outputs = ['%s/%s.sra' % (working_dir, SRRid)]
    return {'job': job, 'prereqs': prereqs, 'outputs': outputs}

def submit_job(job, log={}):
    name = job['job']['name']
    dependencies = []
    errorstate = False
    for prereq in job['prereqs']:
        if not isinstance(prereq, str):
            log = submit_job(prereq, log)
            prereq_id = prereq['job']['name']
        else:
            prereq_id = prereq
        if prereq_id in log:
            if isinstance(log[prereq_id], int):
                dependencies.append(str(log[prereq_id]))
            elif log[prereq_id] is None:
                errorstate = True
        else:
            errorstate = True
    if errorstate:
        log[name] = None
        return log
    if len(dependencies) == 0:
        outputs_fulfilled = True
        for output in job['outputs']:
            if not os.path.exists(output):
                outputs_fulfilled = False
        if outputs_fulfilled:
            log[name] = 'Done'
            return log
    else:
        job['job']['options'].append('--dependency=afterok:%s' % ':'.join(dependencies))
    log[name] = submit_sbatch(job['job'])
    return log

def submit_sbatch(job):
    script_name = job['fname']
    write_script(**job)
    temp = subprocess.Popen(["sbatch", script_name], stdout=subprocess.PIPE)
    ID = int(temp.stdout.readline().rstrip('\n').split(' ')[-1])
    return ID

def load_yaml_file(fname):
    if os.path.exists(fname):
        data = yaml.load(open(fname))
    else:
        data = {}
    return data

def write_yaml(data, fname):
    with open(fname, 'w') as f:
        yaml.dump(data, f, default_flow_style=False)

def write_script(fname, task, ntasks=1, name='', partition='shared', options=[], modules=[], time="4:0:0"):
    output = open(fname, 'w')
    print >> output, "#!/bin/bash -l"
    print >> output, ""
    print >> output, "#SBATCH"
    print >> output, "#SBATCH --job-name=%s" % name
    print >> output, "#SBATCH --time=%s" % time
    print >> output, "#SBATCH --error=/home-3/msauria1@jhu.edu/scratch/slurm_out/%s" % name
    print >> output, "#SBATCH --output=/home-3/msauria1@jhu.edu/scratch/slurm_out/%s" % name
    print >> output, "#SBATCH --nodes=1"
    print >> output, "#SBATCH --ntasks-per-node=%i" % ntasks
    print >> output, "#SBATCH --partition=%s" % partition
    for opt in options:
        print >> output, "#SBATCH %s" % opt
    for mod in modules:
        print >> output, "module load %s" % mod
    print >> output, task
    output.close()


if __name__ == "__main__":
    main()
#!/usr/bin/env python

import os
import subprocess
from ftplib import FTP
from cStringIO import StringIO
import tarfile

import urllib2
import yaml
from lxml import etree

base_dir = "/Users/msauria/projects/repos/HiC_Hosting/YAML"

def main():
    subprocess.Popen(['mkdir', '-p', '%s/GSE' % base_dir])
    subprocess.Popen(['mkdir', '-p', '%s/GSM' % base_dir])
    master = load_yaml_file("%s/all_GSEs.yml" % base_dir)
    blacklist = load_yaml_file("%s/GSE/blacklist.yml" % base_dir)
    GSEs = load_GSEs(master)
    GSMs = load_GSMs(GSEs)
    master, GSEs, GSMs = update_GSEs(master, GSEs, GSMs, blacklist)
    subprocess.Popen(['rm', '-f', '%s/GSE/temp.yml' % base_dir])
    subprocess.Popen(['rm', '-f', '%s/GSM/temp.yml' % base_dir])

def get_etree(path):
    return etree.parse(path, parser=etree.XMLParser(recover=True))

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

def update_GSEs(master, GSEs, GSMs, blacklist):
    """
    try:
        tree = get_etree(urllib2.urlopen("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=gds&term=hi-c&retmax=5000"))
        IDs = tree.xpath("//IdList")[0].getchildren()
    except:
        print "Failed to get search results"
        IDs = []
    for i in range(len(IDs))[::-1]:
        IDs[i] = IDs[i].text
        GSE_ID = "GSE%s" % IDs[i][1:].lstrip('0')
        if GSE_ID in GSEs:
            del IDs[i]
        elif GSE_ID in blacklist:
            del IDs[i]
        elif IDs[i][0]  != '2':
            del IDs[i]
    """
    IDs = ['200095476']
    print "%i Series to process" % len(IDs)
    for ID in IDs:
        GSE_ID = "GSE%s" % ID[1:].lstrip('0')
        try:
            ftp = FTP('ftp.ncbi.nlm.nih.gov')
            ftp.login()
            ftp.cwd('/geo/series/%snnn/%s/miniml' % (GSE_ID[:-3], GSE_ID))
            w = StringIO()
            ftp.retrbinary('RETR %s_family.xml.tgz' % GSE_ID, w.write)
            r = StringIO(w.getvalue())
            tar = tarfile.open(mode="r:gz", fileobj=r)
        except:
            print "Failed to download GEO data for %s" % GSE_ID
            continue
        names = tar.getnames()
        index = -1
        for i in range(len(names)):
            if names[i].count('family.xml') > 0:
                index = i
        if index == -1:
            remove = True
            print 'Family not found'
        else:
            remove = False
            tree = etree.fromstring(tar.extractfile(tar.getnames()[index]).read())
            result = tree.getchildren()
            if len(result) == 0:
                print 'No children'
                remove = True
            else:
                GSE = parse_xml(tree, ID)
                write_yaml(GSE, "%s/GSE/temp.yml" % (base_dir))
        # Filter out non-HiC series
        if not remove and 'Type' in GSE['Series'] and 'Other' in GSE['Series']['Type']:
            sra = False
            if 'Relation' in GSE['Series']:
                if isinstance(GSE['Series']['Relation'], list):
                    for entry in GSE['Series']['Relation']:
                        if entry['type'] == 'SRA':
                            sra = True
                elif GSE['Series']['Relation']['type'] == 'SRA':
                    sra = True
            if not sra:
                for sample in GSE['Sample']:
                    if 'Type' in sample and sample['Type'] == 'SRA':
                        sra = True
            if not sra:
                print 'No SRA connection'
                remove = True
        else:
            print 'Not of type: OTHER'
            remove = True
        if not remove:
            # Parse samples
            samples = GSE['Sample']
            if isinstance(samples, dict):
                samples = [samples]
            del GSE['Sample']
            for i in range(len(samples))[::-1]:
                sample = samples[i]
                if sample['Accession'] in GSMs:
                    del samples[i]
                    continue
                sample = process_sample(sample)
                if sample is None:
                    del samples[i]
                else:
                    sample['Series'] = GSE['Series']['Accession']
            if len(samples) == 0:
                print 'No valid samples'
                remove = True
        if not remove:
            GSE['Series']['Sample-Ref'] = []
            for sample in samples:
                GSE['Series']['Sample-Ref'].append({'ref': sample['Accession'], 'name': sample['Title']})
            for sample in samples:
                sample = clean_text(sample)
                GSMs[sample['Accession']] = sample
                write_yaml(sample, "%s/GSM/%s.yml" % (base_dir, sample['Accession']))
            GSE = clean_text(GSE)
            GSEs[GSE_ID] = GSE
            write_yaml(GSE, "%s/GSE/%s.yml" % (base_dir, GSE_ID))
            master[GSE['Series']['Accession']] = GSE['id']
            write_yaml(master, "%s/all_GSEs.yml" % base_dir)
        else:
            print "Blacklisted %s" % GSE_ID
            blacklist[GSE_ID] = None
            write_yaml(blacklist, "%s/GSE/blacklist.yml" % base_dir)
    return master, GSEs, GSMs

def process_sample(sample):
    write_yaml(sample, '%s/GSM/temp.yml' % base_dir)
    if sample['Library-Strategy'] not in  ['OTHER', 'Hi-C']:
        print sample['Accession'], sample['Title'], sample['Library-Strategy']
        return None
    SRAid = None
    if isinstance(sample['Relation'], dict):
        sample['Relation'] = [sample['Relation']]
    for entry in sample['Relation']:
        if entry['type'] == 'SRA':
            SRAid = entry['target'].split('term=')[1]
    if SRAid is None:
        return None
    ftp = None
    if isinstance(sample['Supplementary-Data'], str) and sample['Supplementary-Data'].count(SRAid) > 0:
        ftp = sample['Supplementary-Data']
    elif isinstance(sample['Supplementary-Data'], list):
        for entry in sample['Supplementary-Data']:
            if entry.count(SRAid) > 0:
                ftp = entry
    sample['SRA'] = {'iid': SRAid, 'ftp': ftp, 'SRR': []}
    tree = get_etree(urllib2.urlopen("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=sra&term=%s" % SRAid))
    IDs = tree.xpath("//IdList")[0].getchildren()
    for ID in IDs:
        ID = ID.text
        tree = get_etree(urllib2.urlopen("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=sra&id=%s" % ID))
        record = parse_record(tree)
        if 'Runs' in record:
            runs = record['Runs'].rstrip('/>\n').split('/><')
            for run in runs:
                run = run.split(' ')[1:]
                SRR = {'id': ID}
                for entry in run:
                    name, value = entry.split('=')
                    SRR[name] = value
                sample['SRA']['SRR'].append(SRR)
    partition = []
    for name in ['HindIII', 'NcoI', 'DpnII', 'MboI', 'MluI', 'CviJI', 'MspI', 'NlaIII', 'BglII', 'MNase', 'DNase']:
        if find_text(sample, name):
            partition.append(name)
    if len(partition) == 1:
        sample['Partition'] = partition[0]
    else:
        sample['Partition'] = '?'
    return sample

def parse_xml(tree, ID):
    record = {'id': ID}
    keys = {}
    for element in tree:
        tag = element.tag.split("MINiML}")
        if len(tag) == 2:
            tag = tag[1]
        else:
            tag = tag[0]
        if tag not in keys:
            value = parse_xml_element(element)
            if len(value) > 0:
                record[tag] = value
                keys[tag] = 1
        else:
            value = parse_xml_element(element)
            if len(value) > 0:
                if keys[tag] == 1:
                    record[tag] = [record[tag]]
                record[tag].append(value)
                keys[tag] += 1
    return record

def parse_xml_element(element):
    children = element.getchildren()
    if len(children) == 0 and element.text is not None and element.text != '':
        text = element.text.strip("\'\" \n\r\t")
        return text
    attribs = {}
    keys = {}
    if len(element.attrib) > 0:
        for key in element.attrib.keys():
            attribs[key] = element.get(key).strip("\'\" \n\r\t")
    if len(children) > 0:
        for child in children:
            tag = child.tag.split("MINiML}")
            if len(tag) == 2:
                tag = tag[1]
            else:
                tag = tag[0]
            value = parse_xml_element(child)
            if len(value) > 0:
                if tag not in keys:
                    attribs[tag] = value
                    keys[tag] = 1
                else:
                    if keys[tag] == 1:
                        attribs[tag] = [attribs[tag]]
                    attribs[tag].append(value)
                    keys[tag] += 1
    return attribs

def parse_record(tree):
    record = {}
    items = tree.xpath("//DocSum/Item")
    for item in items:
        parsed = parse_element(item)
        if parsed is not None:
            record[item.get('Name')] = parsed
    if len(record) > 0:
        return record
    return None

def write_yaml(data, fname):
    with open(fname, 'w') as f:
        yaml.dump(data, f, default_flow_style=False)

def parse_element(element):
    if 'Type' not in element.attrib:
        return None
    etype = element.get('Type')
    if etype in ['String', 'Integer', 'Float']:
        return element.text
    elif etype == 'List':
        children = []
        for child in element.getchildren():
            parsed_child = parse_element(child)
            if parsed_child is not None:
                children.append(parsed_child)
        if len(children) > 0:
            return children
    elif etype == 'Structure':
        children = {}
        for child in element.getchildren():
            parsed_child = parse_element(child)
            if parsed_child is not None:
                children[child.get('Name')] = parsed_child
        if len(children) > 0:
            return children
    return None

def find_text(entry, target):
    if isinstance(entry, str):
        return entry.count(target) > 0
    elif isinstance(entry, list):
        present = False
        for line in entry:
            present += find_text(line, target)
        return present
    elif isinstance(entry, dict):
        present = False
        for key, line in entry.iteritems():
            present += find_text(line, target)
        return present
    return False

def clean_text(data):
    if isinstance(data, str):
        data = data.strip("\'\"")
    elif isinstance(data, list):
        for i in range(len(data)):
            data[i] = clean_text(data[i])
    elif isinstance(data, dict):
        for key, value in data.iteritems():
            data[key] = clean_text(value)
    return data


if __name__ == "__main__":
    main()

#!/usr/bin/env python

'''
Each YAML file should have the following with the same name:
		1) *.hcd file (200Mb - 5Gb)
		2) *.hcp file
		3) *_bin.hcp file (should be larger than .hcp file)
		4) *_exp.hcp file (should be larger than .hcp file)
				- False if doesn't exist or same size as .hcp file
				
Number 5-10 have the same names as the SRR accessions found in each YAML file:
		5) *.paired file (must be larger than 0)
		6) *_1.bam file
		7) *_2.bam file
		8) *.stats file
		9) *.polychimeric file
		10)	 *.singletons file
'''

import sys
import glob
import os


def main():
		main_dir = os.getcwd()
		working_yaml_dir = '/scratch/groups/jtayl139/users/msauria1/HiC_Database/YAML'
		working_data_dir = '/scratch/groups/jtayl139/users/msauria1/HiC_Database/Data'
		archive_hifive_dir = '/work-zfs/jtayl139/msauria1/HiC/HiFive'
		archive_mapping_dir = '/work-zfs/jtayl139/msauria1/HiC/Mapping'
		finished_file = '../finished.txt'

		# Get yaml files from dir
		sorted_yaml_list = get_all_yaml_names(working_yaml_dir)

		# Get SRR values from each yaml file in dir gathered in previous method
		srr_dict = get_all_srr_values(sorted_yaml_list)

		# Get yaml files from finished.txt
		finished_yaml_file_list, finished_yaml_list = add_yamls_from_file(sorted_yaml_list, finished_file, main_dir)
		sorted_yaml_list.extend(finished_yaml_file_list)
		sorted_yaml_list.sort()

		# Create dictionary
		yaml_list_of_dicts = create_dictionary(sorted_yaml_list)

		# Check both file types in all possible locations (so far, 3 locations)
		yaml_list_of_dicts = check_yaml_files(yaml_list_of_dicts, working_data_dir)
		yaml_list_of_dicts = check_yaml_files(yaml_list_of_dicts, archive_hifive_dir)
		yaml_list_of_dicts = check_yaml_files(yaml_list_of_dicts, archive_mapping_dir)

		yaml_list_of_dicts = check_srr_files(yaml_list_of_dicts, srr_dict, working_data_dir)
		yaml_list_of_dicts = check_srr_files(yaml_list_of_dicts, srr_dict, archive_hifive_dir)
		yaml_list_of_dicts = check_srr_files(yaml_list_of_dicts, srr_dict, archive_mapping_dir)

		# Update 'Complete' key in the inventory dictionary
		final_yaml_list_of_dicts_inventory = complete_the_inventory(finished_yaml_list, yaml_list_of_dicts)

		# Print inventory to text file in same location as script
		os.chdir(main_dir)
		print_to_text_file(final_yaml_list_of_dicts_inventory)


# Get the names of all yaml files	   
def get_all_yaml_names(working_yaml_dir):
		os.chdir(working_yaml_dir)
		sorted_yaml_list = glob.glob('*.yml')
		sorted_yaml_list.sort()

		return sorted_yaml_list


def get_all_srr_values(sorted_yaml_list):
		srr_dict = {}

		for file in sorted_yaml_list:
				temp_srr_list = []
				with open(file) as f:
						for line in f:
								if '- acc:' in line:
										temp_srr_list.append((line.split('acc: ')[1]).strip())
				srr_dict[file] = temp_srr_list

		return srr_dict


# Get yaml files from finished.txt and add to existing yaml file list
def add_yamls_from_file(sorted_yaml_list, finished_file, main_dir):
		content = []
		os.chdir(main_dir)
		finished_yaml_file_list = []

		with open(finished_file) as f:
				for line in f:
						if (line.strip() + '.yml')  not in sorted_yaml_list:
								content.append(line.strip() + '.yml')
						finished_yaml_file_list.append(line.strip())

		return content, finished_yaml_file_list


# Create a dictionary for the statuses of all the necessary files (True or False)
def create_dictionary(sorted_yaml_list):
		yaml_list_of_dicts = []

		for yaml_file in sorted_yaml_list:
				yaml_dict = {'Name': '', 'Complete': '', '*.hcd': False, '*.hcp': False, '*_bin.hcp': False, '*_exp.hcp': False,
						'*.paired': False, '*_1.bam': False, '*_2.bam': False, '*.stats': False, '*.polychimeric': False, '*.singletons': False}
				yaml_dict['Name'] = yaml_file[:len(yaml_file)-4]
				yaml_list_of_dicts.append(yaml_dict)

		return yaml_list_of_dicts


# Update the dictionary with the statuses for each file type of each yaml file
def check_yaml_files(yaml_list_of_dicts, data_dir):
		count = 0

		# Checks for *.hcd, *.hcp, *_bin.hcp, *_exp.hcp
		for dict in yaml_list_of_dicts:
				os.chdir(data_dir)
				if os.path.isfile(dict['Name'] + '.hcd'):
						yaml_list_of_dicts[count]['*.hcd'] = True
				if os.path.isfile(dict['Name'] + '.hcp'):
						yaml_list_of_dicts[count]['*.hcp'] = True
						if os.path.isfile(dict['Name'] + '_bin.hcp'):
								if os.path.getsize(dict['Name'] + '_bin.hcp') > os.path.getsize(dict['Name'] + '.hcp'):
										yaml_list_of_dicts[count]['*_bin.hcp'] = True
						if os.path.isfile(dict['Name'] + '_exp.hcp'):
								if os.path.getsize(dict['Name'] + '_exp.hcp') > os.path.getsize(dict['Name'] + '.hcp'):
										yaml_list_of_dicts[count]['*_exp.hcp'] = True

				count += 1

		return yaml_list_of_dicts


# Check the files with the name as the srr parsed from yaml file
def check_srr_files(yaml_list_of_dicts, srr_dict, working_data_dir):
		os.chdir(working_data_dir)
		for k,v in srr_dict.items():
				# Creating lists of the statuses since there can be multiple SRR accessions in each yaml file
				paired_true_list = []
				onebam_true_list = []
				twobam_true_list = []
				stats_true_list = []
				polychimeric_true_list = []
				singleton_true_list = []

				if v:
						for srr in v:
								paired_file = srr + '.paired'
								onebam_file = srr + '_1.bam'
								twobam_file = srr + '_2.bam'
								stats_file = srr + '.stats'
								polychimeric_file = srr + '.polychimeric'
								singletons_file = srr + '.singletons'

								if os.path.isfile(paired_file) and not os.path.getsize(paired_file) == 0:
										paired_true_list.append(True)
								elif not os.path.isfile(paired_file) or os.path.getsize(paired_file) == 0:
										paired_true_list.append(False)
								if os.path.isfile(onebam_file):
										onebam_true_list.append(True)
								elif not os.path.isfile(onebam_file):
										onebam_true_list.append(False)
								if os.path.isfile(twobam_file):
										twobam_true_list.append(True)
								elif not os.path.isfile(twobam_file):
										twobam_true_list.append(False)
								if os.path.isfile(stats_file):
										stats_true_list.append(True)
								elif not os.path.isfile(stats_file):
										stats_true_list.append(False)
								if os.path.isfile(polychimeric_file):
										polychimeric_true_list.append(True)
								elif not os.path.isfile(polychimeric_file):
										polychimeric_true_list.append(False)
								if os.path.isfile(singletons_file):
										singleton_true_list.append(True)
								elif not os.path.isfile(singletons_file):
										singleton_true_list.append(False)

						for dict in yaml_list_of_dicts:
								if dict['Name'] == k[:-4]:
										if paired_true_list and not False in paired_true_list:
												dict['*.paired'] = True
										if onebam_true_list and not False in onebam_true_list:
												dict['*_1.bam'] = True
										if twobam_true_list and not False in twobam_true_list:
												dict['*_2.bam'] = True
										if stats_true_list and not False in stats_true_list:
												dict['*.stats'] = True
										if polychimeric_true_list and not False in polychimeric_true_list:
												dict['*.polychimeric'] = True
										if singleton_true_list and not False in singleton_true_list:
												dict['*.singletons'] = True

		return yaml_list_of_dicts


# Complete inventory by checking files in 'finished.txt' and looking through values in dictionary
def complete_the_inventory(finished_yaml_list, yaml_list_of_dicts):
		i = 0
		while i < len(yaml_list_of_dicts):
				if yaml_list_of_dicts[i]['Name'] in finished_yaml_list:
						yaml_list_of_dicts[i]['Complete'] = True
				elif False in yaml_list_of_dicts[i].values():
						yaml_list_of_dicts[i]['Complete'] = False
				else:
						yaml_list_of_dicts[i]['Complete'] = True
				i += 1

		return yaml_list_of_dicts

# Print to text file
def print_to_text_file(final_yaml_list_of_dicts_inventory):
		with open("Inventory.txt", "w") as file:
				for dict in final_yaml_list_of_dicts_inventory:
						file.write(str(dict) + '\n')


if __name__ == "__main__":
		main()
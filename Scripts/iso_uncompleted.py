#!/usr/bin/env python

def main():
	uncomplete_sample_list = []
	
	with open('Inventory.txt') as file:
		for line in file:
			before_com, target_com, after_com = line.partition("'Complete':")
			before_name, target_name, after_name = line.partition("'Name':")
			sample = after_name[after_name.find("'")+1:after_name.find("'", after_name.find("'") + 1)]
			
			if after_com[1:6] == 'False':
				uncomplete_sample_list.append(sample)
				
	file.close()
	
	
	with open('incomplete_samples.txt', 'w') as file:
		for sample in uncomplete_sample_list:
			file.write(sample + '\n')
			
	file.close()


if __name__ == "__main__":
		main()
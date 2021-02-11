import sys

def add_last_portion(line, portions):
	last_comma = line.rfind(',')
	if  line.find('_', last_comma) == -1:
		return
	else:
		stripped_portion = line[last_comma+1:].strip()
		last_underscore = stripped_portion.find('_')
		portions.append('{0}'.format(stripped_portion[:last_underscore]))
		portions.append('{0}'.format(stripped_portion[last_underscore + 1:]))


def reformat_line(line: str):
	portions = []
	current_comma = line.index(',')
	portions.append(line[:current_comma])
	past_comma = current_comma	
	current_comma = line.find(',', current_comma+1)
	while current_comma != -1:	
		stripped_portion = line[past_comma+1:current_comma].strip()
		underscore_index = stripped_portion.index('_')
		portions.append('{0}'.format(stripped_portion[:underscore_index]))
		portions.append('{0}'.format(stripped_portion[underscore_index+1:]))
		past_comma = current_comma
		current_comma = line.find(',', current_comma+1)
	add_last_portion(line, portions)
	return portions

with open(sys.argv[1], 'r') as original_histogram:
	all_lines = [reformat_line(line) for line in original_histogram]
	with open(sys.argv[1] + '.mot', 'w') as reformatted_histogram:
		for portions in all_lines:
			reformatted_histogram.write(', '.join(portions))
			reformatted_histogram.write('\n')

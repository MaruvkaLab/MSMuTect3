import sys, re


def get_bounded_segment(segments, start_index):
    histogram: list = []
    i = start_index
    while segments[i] != "-99999999":
        histogram.append(segments[i])
        i+=1
    return " ".join(histogram), i


def get_probabilities(segments, start_index):
    probabilities = []
    i = start_index
    while not segments[i][0] == "[":
        probabilities.append(segments[i])
        i+=1
    return " ".join(probabilities), i


def get_list(segments, start_index):
    i = start_index
    if segments[i][-1] == "]":
        return segments[i][1:-1], i
    entries = [segments[i][1:]]
    i+=1
    while segments[i][-1] != ']':
        entries.append(segments[i])
        i+=1
    entries.append(segments[i][:-1])
    return " ".join(entries), i


def get_second_freq(segments, start_index):
    i = start_index
    current_lists = 0
    while current_lists != 3:
        i+=1
        if segments[i][0] == "[":
            current_lists+=1
    return get_list(segments, i)

def unpack(current_line):
    packed_line = []
    for attribute in current_line:
        attribute_length = len(attribute)
        for i in range(attribute_length):
            packed_line.append(attribute[i])
            if i != attribute_length - 1:
                packed_line.append(" ")
        packed_line.append("\t")
    return "".join(packed_line)

def main(cleaned_mut_path: str, output_path):
    header = "Locus	Decision\tNormal_histogram\tNormal_alleles\tNormal_frequencies\tTumor_histogram\tTumor_alleles\tTumor_frequencies"
    output_lines = [header]
    with open(cleaned_mut_path, 'r') as mut_file:
        for line in mut_file:
            current_line = []
            stripped_line = line.strip()
            edited_line = re.sub("\s+", ",", stripped_line)
            all_segments = edited_line.split(",")
            current_line.append([all_segments[0]])
            current_line.append([all_segments[1]])
            histo_string, current_index = get_bounded_segment(all_segments, 2)
            current_line.append([histo_string])
            allelic_string, current_index = get_bounded_segment(all_segments, current_index+3)
            current_line.append([allelic_string])

            probabilities, current_index = get_probabilities(all_segments, current_index+1)
            current_line.append([probabilities])

            second_alleles, current_index = get_second_freq(all_segments, current_index)
            current_line.append([second_alleles])

            second_real_alleles, current_index = get_list(all_segments, current_index+1)
            current_line.append([second_real_alleles])

            second_probabilities, _ = get_list(all_segments, current_index+1)
            current_line.append([second_probabilities])

            output_lines.append(unpack(current_line))
    with open(output_path, 'w') as output_file:
        output_file.write("\n".join(output_lines))


if __name__=='__main__':
    main(sys.argv[1], sys.argv[2])

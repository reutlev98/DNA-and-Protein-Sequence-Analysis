# Reut Lev 207385741
import sys
from Bio import SeqIO
import re
import os
import subprocess
import matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt


def grafCreatingEx1Part1(template_counts, list_name_of_prosite, ax1):
    """
    The function is creating a bar chart of the occurrences of patterns in protein sequences
    :param list_name_of_prosite: the name of the pattern
    :param template_counts:dictionary the keys are patterns and the values are counts of occurrences of the patterns.
    :param ax1:  subplot axis object
    """
    patterns = list(template_counts.keys())  # Create a list of patterns
    counts = list(template_counts.values())  # Create a list of counts
    ax1.bar(list_name_of_prosite, counts)  # Create the bar chart
    # Add labels and title
    ax1.set_xlabel('PROSITE Patterns')
    ax1.set_ylabel('Number of occurrences')
    ax1.set_title('Occurrences of PROSITE Patterns in Protein Sequences')
    # Rotate x-axis labels
    ax1.set_xticks(range(len(list_name_of_prosite)))
    ax1.set_xticklabels(list_name_of_prosite, rotation=45, ha='right')


def translateToRegex(prosite_pattern):
    """
    the function is converting a PROSITE pattern to a regular expression pattern
    :param prosite_pattern: string representing a pattern (in PROSITE format)
    :return: the prosite pattern in regex pattern
    """
    """
    prosite_pattern = re.sub(r'\[([^\]]+)\]', r'[\1]', prosite_pattern)  # Convert square brackets to regex character
    prosite_pattern = re.sub(r'\{([^\}]+)\}', r'[^\1]', prosite_pattern)  # convert curly brackets to regex negated char
    prosite_pattern = re.sub(r'x', r'.', prosite_pattern)  # replace the wildcard with . (dot)
    prosite_pattern = re.sub(r'\((\d+)\)', r'{\1}', prosite_pattern)  # Replace the number in parentheses.
    # Replace the > and < with ^ or $
    prosite_pattern = re.sub(r'<', '^', prosite_pattern)
    prosite_pattern = re.sub(r'>', '$', prosite_pattern)
    prosite_pattern = prosite_pattern.replace("-", "")  # remove the "-" separator from the pattern
    """
    prosite_pattern = prosite_pattern.replace("-", "").replace("{", "[^").replace("}", "]").replace("x", ".").replace(
        ")", "}").replace("(", "{").replace("<", "^").replace(">", "$")
    return prosite_pattern


def pattern_occurrence_count(fasta_file_name, prosite_list):
    """
    count the occurrences of pattern in the fasta file.
    :param fasta_file_name: representing the name of a fasta file.
    :param prosite_list: a list of strings representing patterns in the PROSITE format.
    :return: template_counts dictionary, which contains the counts of occurrences of each pattern in the fasta file.
    """
    fasta_sequences = SeqIO.parse(open(fasta_file_name), 'fasta')  # parse the fasta file
    for index in range(len(prosite_list)):
        prosite_list[index] = translateToRegex(prosite_list[index])
    template_counts = {template: 0 for template in prosite_list}
    # iterate over fasta sequences
    for fasta in fasta_sequences:
        name, sequence = fasta.id, str(fasta.seq)
        # iterate over prosite templates
        for template in prosite_list:
            # search for template in sequence
            matches = re.finditer(template, sequence)
            # count number of occurrences
            count = len([match.start() for match in matches])
            template_counts[template] += count
    # grafCreatingEx1Part1(template_counts)
    return template_counts


def create_one_figure(template_counts, list_name_of_prosite, lengths_file_1, lengths_file_2, file_name1, file_name2):
    """
    this faction create one figure for the first part
    :param lengths_file_1: the counts in file 1
    :param lengths_file_2: the counts in file 2
    :param file_name1: fasta file name 1
    :param file_name2: fasta file name 2
    :param list_name_of_prosite: list of the prosite name
    :param template_counts: template count (for graph1)
    :return:
    """
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6))  # Create a figure with 2 subplots
    grafCreatingEx1Part1(template_counts, list_name_of_prosite, ax1)  # Plot the first graph on the first subplot
    grafCreatingEx1Part2(lengths_file_1, lengths_file_2, file_name1, file_name2, ax2)
    fig.suptitle(
        'Comparison of PROSITE Patterns and motif lengths in 2 FASTA files')  # Add a title for the entire figure
    fig.savefig('207385741.png', bbox_inches='tight')  # Save the figure as a PNG file


def grafCreatingEx1Part2(lengths1, lengths2, fasta_file1, fasta_file2, ax2):
    """
    creating the graph of Ex1Part2
    :param fasta_file2: fasta file 2
    :param fasta_file1:  fasta file1
    :param lengths2: list of the length of the pattern in the seq
    :param lengths1: list of the count of the pattern in the seq
    :param ax2: ax object
    """
    # Create a scatter plot comparing the lengths of the motifs found in the two files
    ax2.scatter(list(lengths1.keys()), list(lengths1.values()), label=fasta_file1)
    ax2.scatter(list(lengths2.keys()), list(lengths2.values()), label=fasta_file2)
    ax2.set_yscale("log")
    ax2.legend()
    ax2.set_xlabel("Length of Motif")
    ax2.set_ylabel("Number of Occurrences (log scale)")
    ax2.set_title('Numbers of appearance for each motif in fasta files')


def find_motif_lengths(fasta_file1, fasta_file2, prosite_list):
    """
    this function find the lengths of the motives
    :param fasta_file1: fasta file 1 name
    :param fasta_file2: fasta file 2 name
    :param prosite_list: list of pattern
    :return: lengths and file name
    """
    # Initialize empty dictionaries to store the lengths of the motifs found
    lengths1 = {}
    lengths2 = {}
    for index in range(len(prosite_list)):  # translate to regex format
        prosite_list[index] = translateToRegex(prosite_list[index])
    # Iterate over the two FASTA files
    for fasta, lengths in [(fasta_file1, lengths1), (fasta_file2, lengths2)]:
        for seq_record in SeqIO.parse(fasta, "fasta"):
            for template in prosite_list:
                # Iterate over the occurrences of the template in the sequence
                for match in re.finditer(template, str(seq_record.seq)):
                    start, end = match.start(), match.end()
                    length = end - start
                    if length in lengths:
                        lengths[length] += 1
                    else:
                        lengths[length] = 1

    return lengths1, lengths2


# the function with exception throw
def run_fastp1(fastp_path, fastq_file):
    """
    the function return The difference between the number of reads in beginning and after filtering
    :param fastp_path: the path for the fastp program
    :param fastq_file: path to fastq program
    :return: The difference between the number of reads in beginning and after filtering
    """
    fastq_file = os.path.abspath(fastq_file)
    fastp_path = os.path.abspath(fastp_path)
    if not os.path.exists(fastq_file):
        raise Exception(f"{fastq_file} does not exist")
    if not os.path.isfile(fastq_file):
        raise Exception(f"{fastq_file} is not a file")
    if not os.access(fastq_file, os.R_OK):
        raise Exception(f"{fastq_file} is not readable")

    # Run the fastp software using subprocess
    process = subprocess.Popen([fastp_path, '-i', fastq_file, '-o', '/dev/null'],
                               stderr=subprocess.PIPE)
    # Wait for 5 seconds
    process.wait()

    output, error = process.communicate()
    # Check if the process has returned a non-zero exit code
    if process.returncode != 0:
        raise Exception("Process returned with non-zero exit code.")
    # Search for the number of threads at the beginning
    match = re.search(r'Read1 before filtering:\ntotal reads:\s+(\d+)', error.decode())
    if match:
        read_before_filter = int(match.group(1))
    else:
        raise Exception("Unable to find the number of threads at the beginning")
    # Search for the number of threads at the end
    match = re.search(r'Read1 after filtering:\ntotal reads:\s+(\d+)', error.decode())
    if match:
        read_after_filter = int(match.group(1))
    else:
        raise Exception("Unable to find the number of threads at the end")
    # Print the difference
    print("The difference between the number of reads in beginning and after filtering is:",
          read_before_filter - read_after_filter)  #


# the correct function 1.2
def run_fastp(fastp_path, fastq_file):
    """
    the function return The difference between the number of reads in beginning and after filtering
    :param fastp_path: the path for the fastp program
    :param fastq_file: path to fastq program
    :return: The difference between the number of reads in beginning and after filtering
    """
    fastq_file = os.path.abspath(fastq_file)  # to check the path of the program
    fastp_path = os.path.abspath(fastp_path)
    if not os.path.exists(fastq_file):
        raise Exception(f"{fastq_file} does not exist")
    if not os.path.isfile(fastq_file):
        raise Exception(f"{fastq_file} is not a file")
    if not os.access(fastq_file, os.R_OK):
        raise Exception(f"{fastq_file} is not readable")

    # Run the fastp software using subprocess
    process = fastp_path + ' -i ' + fastq_file + ' -o ' + '/dev/null'
    # Wait for 5 seconds
    result = subprocess.run(process, check=True, shell=True, timeout=5, stderr=subprocess.PIPE)
    output = result.stderr.decode()
    read_before_filter = int(output.split("\n")[4].split(":")[1].strip())
    read_after_filter = int(output.split("\n")[10].split(":")[1].strip())
    print("The difference between the number of reads in beginning and after filtering is:",
          read_before_filter - read_after_filter)


def main(input_file, fastp_path):
    """
    this function read the input file and place the data in the variables and call the function for Ex1 and EX2
    :param input_file: the input file
    :param fastp_path: path to fastp program
    :return:
    """
    # read the path of fasta files from the text file
    fasta_path_list = []
    prosite_data_list = []
    with open(input_file, 'r') as f:
        for index in range(3):  # insert the fasta file into a list
            fasta_path_list.append((f.readline()).replace('\n', ''))  # remove the \n from the end
        fastq_path = f.readline().replace('\n', '')  # insert the fastq
        for line in f:  # insert the prosite data into a list
            prosite_data_list.append((line.replace('\n', '')).split(";"))  # remove the \n and split the prosite;name

        list_of_prosite = [onePrositeNAme[0] for onePrositeNAme in prosite_data_list]  # just the prosite
        list_name_of_prosite = [onePrositeNAme[1] for onePrositeNAme in prosite_data_list]  # just the name of prosite
        template_counts = pattern_occurrence_count(fasta_path_list[0], list_of_prosite)  # The first FASTA file1&prosite
        lengths_file_1, lengths_file_2 = find_motif_lengths(fasta_path_list[1], fasta_path_list[2], list_of_prosite)
        create_one_figure(template_counts, list_name_of_prosite, lengths_file_1, lengths_file_2, fasta_path_list[1],
                          fasta_path_list[2])  # create one figure for ex1

        if not os.path.exists(fasta_path_list[0]):
            raise Exception(f"{fasta_path_list[0]} does not exist")
        run_fastp(fastp_path, fastq_path)


if __name__ == "__main__":
    if len(sys.argv) != 3:  # check if the correct number of arguments were provided
        print("Error: Incorrect number of arguments. Usage: python program.py [text_file] [FASTP_software]")
        sys.exit()
    # assign the input arguments to variables
    text_file = sys.argv[1]
    FASTP_software = sys.argv[2]
    main(text_file, FASTP_software)
    #main("Ex4File/input.txt", "/home/reutlev/fastp")

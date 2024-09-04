# Reut Lev 207385741
import sys


# represent a cell
class Cell:
    # constructor
    def __init__(self, name, genome, reading_frame):
        self.name = name
        self.genome = genome
        self.reading_frame = reading_frame

    # Returns the print format of a cell
    def __str__(self):
        result = f'< {self.name} , {self.genome} >'
        return result

    # find srr in genome by its index
    def srr_find(self, genome_index):
        """
        parameter :genome_index represent index DNA sequence we search in
        return :srr_list a list of srr - simple sequence repeat
        the function gets a DNA sequence and return a list with all the repeats in seq and how many times they appear.
        #genome_index % len(self.genome) -> 6%4=2 !!!!!!
        """
        dna_seq = self.genome[genome_index % len(self.genome)]  # Finding the genome sequence by index
        srr_list = []
        for srr_len in range(1, 7):  # check every option of srr seq in the dna
            for base_index in range(len(dna_seq)):
                srr = dna_seq[base_index:base_index + srr_len]
                srr_count = 1
                search_index = base_index + srr_len
                while dna_seq[search_index:srr_len + search_index] == srr:
                    srr_count += 1
                    search_index += srr_len
                if srr_count >= 3:  # If the sequence repeats 3 times or more
                    temp = [tup[1] for tup in srr_list if srr in tup]
                    if temp:  # Duplicate treatment
                        if srr_count > temp[0]:
                            srr_list.remove(temp)
                            srr_list.append((srr, srr_count))
                    else:
                        srr_list.append((srr, srr_count))
        return srr_list

    # transcribe (from DNA to RNA)
    def transcribe(self, genome_index):
        """
        parameter :genome_index represent DNA sequence index.
        return :dna_seq a complementary DNA sequence.
        the function gets an DNA sequence and returns the complementary RNA sequence the Transcription seq from 5' to 3'
        """
        dna_seq = self.genome[genome_index % len(self.genome)]  # Finding the genome sequence by index
        if dna_seq == "":
            return None
        dna_seq = dna_seq.upper()
        dna_seq = dna_seq[::-1]
        rna_seq = ""
        for base in dna_seq:

            if base == "G":
                rna_seq += "C"
            elif base == "C":
                rna_seq += "G"
            elif base == "A":
                rna_seq += "U"
            elif base == "T":
                rna_seq += "A"
        return rna_seq

    # ex1_3
    def translate(self, genome_index):
        """
           parameter :genome_index represent DNA sequence index.
           read_frame represent a codon from which the sequence is divided into codons.
           return :protein a list with the amino acid sequence
           the function gets an RNA sequence and a number that indicates the reading frame and translates the RNA,
           from the methionine codon to the last full triplet or to the first stop codon.
            """

        read_frame = int(self.reading_frame[genome_index % len(self.reading_frame)])  # Finding the read frame by index
        rna_seq = self.transcribe(genome_index)  # Convert DNA sequence to RNA sequence

        table = {
            'AUA': 'I', 'AUC': 'I', 'AUU': 'I', 'AUG': 'M',
            'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACU': 'T',
            'AAC': 'N', 'AAU': 'N', 'AAA': 'K', 'AAG': 'K',
            'AGC': 'S', 'AGU': 'S', 'AGA': 'R', 'AGG': 'R',
            'CUA': 'L', 'CUC': 'L', 'CUG': 'L', 'CUU': 'L',
            'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCU': 'P',
            'CAC': 'H', 'CAU': 'H', 'CAA': 'Q', 'CAG': 'Q',
            'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGU': 'R',
            'GUA': 'V', 'GUC': 'V', 'GUG': 'V', 'GUU': 'V',
            'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCU': 'A',
            'GAC': 'D', 'GAU': 'D', 'GAA': 'E', 'GAG': 'E',
            'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGU': 'G',
            'UCA': 'S', 'UCC': 'S', 'UCG': 'S', 'UCU': 'S',
            'UUC': 'F', 'UUU': 'F', 'UUA': 'L', 'UUG': 'L',
            'UAC': 'Y', 'UAU': 'Y', 'UAA': '_', 'UAG': '_',
            'UGC': 'C', 'UGU': 'C', 'UGA': '_', 'UGG': 'W',
        }
        start_codon = rna_seq.find("AUG")  # if the AUG isn't in the seq return none
        if start_codon == -1:
            return None
        rna_seq = rna_seq[read_frame - 1:]  # start from the reading frame
        protein = []
        aug_list = []
        stop_list = []
        for x in range(0, len(rna_seq), 3):  # Checks where there is a start codon and an end codon
            code = rna_seq[x:x + 3]
            if code == 'AUG':
                aug_list.append(x)
            if code == 'UAA' or code == 'UGA' or code == 'UAG':
                stop_list.append(x)
        max_seq = 0
        if len(aug_list) == 0:
            return None
        for x in range(len(stop_list)):  # Checks what is the longest amino seq
            if aug_list[x] < stop_list[x]:
                if stop_list[x] - aug_list[x] > max_seq:
                    max_seq = aug_list[x]
            elif aug_list[x] > stop_list[x]:
                max_seq = aug_list[x]
        rna_seq = rna_seq[max_seq:]
        for x in range(0, len(rna_seq), 3):  # translate
            code = rna_seq[x:x + 3]
            if len(code) % 3 != 0 or table[code] == '_':
                break
            protein += table[code]
        return protein

    def repertoire(self):
        """
        :return: returns a list, for each sequence in the genome tuple with:
        repetitions that appear in that sequence
        and the translation
        """
        genome_list = []  # the list that return the repertoire
        for geneI in range(len(self.genome)):
            # srr
            srr_list = self.srr_find(geneI)
            srr_string = ""
            if len(srr_list) == 0:  # if there is no srr in the genome
                srr_string = "No simple repeats in DNA sequence"
                print(srr_string)
            else:
                # Sorting a List by the Second Element of the Tuple
                srr_list.sort(key=lambda a: a[1])
                for srrTup in srr_list:
                    srr_string += str(srrTup[0]) + "," + str(srrTup[1]) + ";"
                srr_string = srr_string[:-1]
                print(srr_string)  # print the srr of the genome

            # translate
            protein = self.translate(geneI)
            protein_string = ""
            if protein is None:  # if there is no translation in the genome
                protein_string = "Non-coding RNA"
                print(protein_string)
            else:
                for pro in protein:
                    protein_string += pro + ";"
                protein_string = protein_string[:-1]
                print("Translation: " + protein_string)  # print the translation of the genome
            tuple_repertoire = (srr_string, protein_string)
            genome_list.append(tuple_repertoire)
        return genome_list


# class that represent Stem cell
class StemCell(Cell):
    # constructor
    def __init__(self, genome, reading_frame):
        super().__init__("Stem Cell", genome, reading_frame)

    #  overrides the multiplication operator
    def __mul__(self, num_of_mul):
        """
        :return: multiplying by an integer will return a list of identical cells as the multiplier number
        """
        cell_list = []  # cell list
        new_cell = StemCell(self.genome, self.reading_frame)  # creating a new cell with same parameters
        for i in range(num_of_mul):
            cell_list.append(new_cell)  # Inserts a number of cells according to the num_of_mul
        return cell_list

    # returns a list with two exactly identical cells.
    def mitosis(self):
        """
        :return:  list of 2 identical cells using the multiplication operation we overrides
        """
        return self * 2

    # the class represent Nerve Cell
    class NerveCell(Cell):
        # constructor
        def __init__(self, genome, reading_frame, coefficient):
            super().__init__("Nerve Cell", genome, reading_frame)
            self.numeric_signal = None
            self.internal_unique_coefficient = coefficient

        # A function that receives a numeric signal
        def receive(self, numeric_signal):
            self.numeric_signal = numeric_signal

        # A function that transmits the signal after multiplying by the internal coefficient.
        def send(self):
            return self.numeric_signal * self.internal_unique_coefficient

    # the class represent Muscle Cell
    class MuscleCell(Cell):
        # constructor
        def __init__(self, genome, reading_frame, path, threshold):
            super().__init__("Muscle Cell", genome, reading_frame)
            self.threshold = threshold
            self.path = path

        # A function that receives a numeric signal and write in file.
        def receive(self, numeric_signal):
            """
           If the signal is higher than the threshold defined for this cell,
           The cell will write the received signal to the file and "I like to move it".
           """
            if numeric_signal > self.threshold:
                f = open(self.path, "a")  # open the file
                temp_str = str(numeric_signal) + ", " + "I like to move it\n"
                f.write(temp_str)  # write in the file
                f.close()  # close the file

    # this function is in the StemCell class - cause the cell to differentiate into a different type of cell.
    def differentiate(self, cell_name, list_of_parameters):
        """
        :return : new Nerve Cell or new Muscle Cell
        """
        if cell_name == "Nerve Cell":
            genome = self.genome
            reading_frame = self.reading_frame
            coefficient = list_of_parameters[0]
            return self.NerveCell(genome, reading_frame, coefficient)
        elif cell_name == "Muscle Cell":
            genome = self.genome
            reading_frame = self.reading_frame
            path = list_of_parameters[0]
            threshold = int(list_of_parameters[1])
            return self.MuscleCell(genome, reading_frame, path, threshold)


# the class represent Nerve Network
# A class that manages the transfer of the signal through the nerve cells to the target cell - a muscle cell
class NerveNetwork:
    # constractor
    def __init__(self, target_muscle, list_of_nerveCell):
        self.target_muscle = target_muscle
        self.list_of_nerveCell = list_of_nerveCell

    # receive a signal and transmit it through the neurons in the list and send it ro muscle cell.
    def signal_send(self, signal):
        """
        :param signal: receive a signal and transmit it through the neurons in the list
        The signal will pass through all the cells in the list, which will pass it from one to the other,
        and then it will be sent to the muscle cell.
        """
        for nerve_cell_i in self.list_of_nerveCell:  # passes the signal between all nerve cells
            nerve_cell_i.receive(signal)
            signal = nerve_cell_i.send()
        self.target_muscle.receive(signal)  # send the last signal to the target - muscle cell


# this function manage th reading input file and creating the cells accordingly
def main(fileName, list_of_signal):
    fileData = []  # list that contains the data we read from the file
    with open(fileName) as f:  # open file
        for line in f:  # Read data line by line
            l_temp = line.split('\t')  # split data by tab and store it in list
            fileData.append(l_temp)

    list_of_nerve = []  # list that contains the nerve cell
    target = ""  # A variable that will contain the muscle cell
    for line in fileData[1:]:  # go over the data line by line
        # Checks for the input validation using assert
        assert line[0] == "NC" or line[0] == "MC", "File illegal"
        assert line[1] and all(c in "ATGC," for c in line[1].upper()), "File illegal"
        assert line[2] and all(c in "123," for c in line[2].upper()), "File illegal"
        assert line[1].count(",") == line[2].count(","), "File illegal"
        assert (line[3].count(",") == 0 and isinstance(float(line[3]), (int, float))) or line[3].count(",") == 1, \
            "File illegal "
        l_genome = line[1].split(",")  # creating list of genome
        l_reading_frame = line[2].split(",")  # creating list of reading frame
        tempStem = StemCell(l_genome, l_reading_frame)
        if line[0] == "NC":  # create nerve cell
            list_cell_mit = tempStem.mitosis()
            for cell_i in list_cell_mit:
                nerve_cell = cell_i.differentiate("Nerve Cell", [float(line[3])])
                list_of_nerve.append(nerve_cell)  # insert the new nerveCell into the nerv cell list
        elif line[0] == "MC":  # create muscle cell
            l_param = line[3].split(",")
            target = tempStem.differentiate("Muscle Cell", l_param)
    network = NerveNetwork(target, list_of_nerve)  # creating nerve network
    for signal in list_of_signal:
        network.signal_send(int(signal))
    target.repertoire()  # print the repertoire one time regardless of the amount of signals


if __name__ == "__main__":
    signal_list = (sys.argv[2]).split(",")
    main(sys.argv[1], signal_list)  # send to the main the file name and the signal list

import sys
"""
Reut Lev , 207385741
"""
# ex1_1
def srr_find(dna_seq):
    """
    parameter :dna_seq represent DNA sequence we search in
    return :srr_list a list of srr - simple sequence repeat
    the function gets a DNA sequence and return a list with all the repeats in seq and how many times they appear in it.
    """
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


# ex1_2
def reverse_transcribe(rna_seq):
    """
    parameter :rna_seq represent RNA sequence.
    return :dna_seq a complementary DNA sequence.
    the function gets an RNA sequence and returns the complementary DNA sequence - the Transcription seq from 5' to 3'.
    """
    if rna_seq == "":
        return None
    rna_seq = rna_seq.upper()
    rna_seq = rna_seq[::-1]
    dna_seq = " "
    for base in rna_seq:

        if base == "G":
            dna_seq += "C"
        elif base == "C":
            dna_seq += "G"
        elif base == "A":
            dna_seq += "T"
        elif base == "U":
            dna_seq += "A"
    return dna_seq


# ex1_3
def translate(rna_seq, read_frame):
    """
   parameter :rna_seq represent RNA sequence we search in.
   read_frame represent a codon from which the sequence is divided into codons.
   return :protein a list with the amino acid sequence
   the function gets an RNA sequence and a number that indicates the reading frame and translates the RNA,
   from the methionine codon to the last full triplet or to the first stop codon.
    """
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
    start_codon = rna_seq.find("AUG")
    if start_codon == -1:
        return None
    rna_seq = rna_seq[read_frame - 1:]
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
        if table[code] == '_':
            break
        protein += table[code]
    return protein


def main(dna_seq, rna_seq1, rna_seq2, reading_frame):
    """
    :parameter :dna_seq for ssr_find() , rna_seq1 for reverse_transcribe() , rna_seq2 + reading_frame for translate().
    """
    # Srr
    srr_list = srr_find(dna_seq)
    if len(srr_list) == 0:
        print("No simple repeats in DNA sequence")
    else:
        # Sorting a List by the Second Element of the Tuple
        srr_list.sort(key=lambda a: a[1])
        srr_string = ""
        for srrTup in srr_list:
            srr_string += str(srrTup[0]) + "," + str(srrTup[1]) + ";"
        print(srr_string[:-1])
    # Transcribe
    print("DNA sequence:" + reverse_transcribe(rna_seq1))
    # Translate
    protein = translate(rna_seq2, int(reading_frame))
    if protein is not None:
        protein_string = ""
        for pro in protein:
            protein_string += pro + ";"
        print("Translation:" + protein_string[:-1])
    else:
        print("Non-coding RNA")


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])


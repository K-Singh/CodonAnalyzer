import Bio as bp
from Bio.Seq import Seq as sq
import numpy as np


# Fills a dictionary based on the frequency of items in the list parameter
# The list parameter is a list filled with consecutive substrings from a parent string
# that are divided such that they correspond to a key in the dictionary
# The frequency_dict parameter is an initialized dictionary whose keys correspond to
# the various substrings found in list and the parent string from which list
# was generated
def fill_frequency_dict(substring_list, frequency_dict):
    """
        Fills a dictionary based on the frequency of items in the list parameter
The list parameter is a list filled with consecutive substrings from a parent string
that are divided such that they correspond to a key in the dictionary
The frequency_dict parameter is an initialized dictionary whose keys correspond to
the various substrings found in list and the parent string from which list
was generated
    :param substring_list: A consecutive array of substrings from a parent string
    :param frequency_dict: An initialized dictionary of each substring with a frequency of 1
    """
    for i in range(0, len(substring_list)):
        item = substring_list[i]
        duplicate_bool = 0
        if item in frequency_dict:
            if frequency_dict[item] != 1:
                duplicate_bool = 1
            # print("coom" + str(pair))
        for j in range(i, len(substring_list)):
            if item == substring_list[j] and j > i:
                if item in frequency_dict and duplicate_bool == 0:
                    #  print(pair)
                    frequency_dict[item] += 1


# Returns a substring_list representing an array of consecutive substrings from the parent string
# and an initialized frequency dictionary frequency_dict whose keys are equal to each unique element in the array
# substring_list is a list of the substrings generated from parent_string
# These substrings are of length substring_length, and are iterated through using substring_step
# If substring_step is equal to substring_length, there is no overlap in substrings from the parent strings
# If substring_step is less than substring_length, there is overlap in the substrings such that the last substring_step
# characters from the previous element are found in the first substring_step characters in the next element
# For example, codon pairs would have a substring_length of 6(the length of a codon pair)
# but a substring_step of 3(the length of a single codon). This would mean an overlap of one codon
# in each pair and therefore each element in substring_list
# For parent_string = "AUGACCAUGACC", substring_list would equal ["AUGACC", "ACCAUG", "AUGACC"]
# list_iteration_num is a number representing how many times the list should be iterated through
# parent_string is the string from which the substring_list and initialized frequency list are created
def gen_frequency_dict(parent_string, list_iteration_num, substring_length, substring_step):
    """
Returns a substring_list representing an array of consecutive substrings from the parent string
and an initialized frequency dictionary frequency_dict whose keys are equal to each unique element in the array
substring_list is a list of the substrings generated from parent_string
These substrings are of length substring_length, and are iterated through using substring_step
If substring_step is equal to substring_length, there is no overlap in substrings from the parent strings
If substring_step is less than substring_length, there is overlap in the substrings such that the last substring_step
characters from the previous element are found in the first substring_step characters in the next element
For example, codon pairs would have a substring_length of 6(the length of a codon pair)
but a substring_step of 3(the length of a single codon). This would mean an overlap of one codon
in each pair and therefore each element in substring_list
For parent_string = "AUGACCAUGACC", substring_list would equal ["AUGACC", "ACCAUG", "AUGACC"]
list_iteration_num is a number representing how many times the list should be iterated through
parent_string is the string from which the substring_list and initialized frequency list are created

    :param parent_string: String from which to generate substring_list
    :param list_iteration_num: Number of times substring_list is iterated through the string
    :param substring_length: Length of each substring element in the substring list
    :param substring_step: Step value used to specify how many characters to avoid adding in next iteration
    :return: substring_list: List of substrings generated from parent_string and parameters
    :return: frequency_dict: Initialized dictionary that contains every unique element of substring_list with a
    frequency of 1
    """
    substring_list = []
    frequency_dict = {}
    counter = 0
    for i in range(0, list_iteration_num, 1):
        item = parent_string[int(counter):int(counter + substring_length)]
        substring_list.append(item)
        if item not in frequency_dict:
            frequency_dict[item] = 1
        counter += substring_step
    return substring_list, frequency_dict


# generates list of codon pairs and a dictionary assigning frequencies to each codon pair
def gen_codon_pair_freq(gene_string):
    codon_pairs = gen_frequency_dict(gene_string, int(len(gene_string) / 3 - 1), 6, 3)
    fill_frequency_dict(codon_pairs[0], codon_pairs[1])
    return codon_pairs


# returns the frequency for a specific codon pair from the frequency dictionary
def get_codon_pair_freq(codon_pair, codon_pair_freq_dict):
    if codon_pair in codon_pair_freq_dict:
        return codon_pair_freq_dict[codon_pair]
    else:
        return 0


# generates list of singular codons and a dictionary assigning frequencies to each codon
def gen_codon_freq(gene_string):
    codons = gen_frequency_dict(gene_string, int(len(gene_string) / 3), 3, 3)
    fill_frequency_dict(codons[0], codons[1])
    return [codons[0], codons[1]]


# returns the frequency for a specific codon from the frequency dictionary
def get_single_codon_freq(codon, codon_freq_dict):
    if codon in codon_freq_dict:
        return codon_freq_dict[codon]
    else:
        return 0


# generates list of amino acid pairs and a dictionary assigning frequencies to each amino acid pair
def gen_residue_pair_freq(res_string):
    residues = gen_frequency_dict(res_string, int(len(res_string) - 1), 2, 1)
    fill_frequency_dict(residues[0], residues[1])
    return [residues[0], residues[1]]


# returns the frequency of an amino acid pair from the frequency dictionary
def get_residue_pair_freq(residue_pair, residue_pair_freq_dict):
    if residue_pair in residue_pair_freq_dict:
        return residue_pair_freq_dict[residue_pair]
    else:
        return 0


# generates list of singular amino acids and a dictionary assigning frequencies to each amino acid
def gen_residue_freq(res_string):
    residues = gen_frequency_dict(res_string, int(len(res_string)), 1, 1)
    fill_frequency_dict(residues[0], residues[1])
    return [residues[0], residues[1]]


# returns the frequency of an amino acid from the frequency dictionary
def get_single_residue_freq(residue, residue_freq_dict):
    if residue in residue_freq_dict:
        return residue_freq_dict[residue]
    else:
        return 0


# translates the given string of genetic code into the amino acids it codes for
def translate_to_protein(gene_string):
    sequence = sq(gene_string)
    return sequence.translate()


# Returns the codon pair score for a given codon pair from a string of genetic code
def find_codon_pair_score(gene_string, codon_pair):
    translated_gene = str(translate_to_protein(gene_string))
    codon_pair_dict = gen_codon_pair_freq(gene_string)[1]
    codon_single_dict = gen_codon_freq(gene_string)[1]
    residue_pair_dict = gen_residue_pair_freq(translated_gene)[1]
    residue_single_dict = gen_residue_freq(translated_gene)[1]

    codon_1 = codon_pair[0:3]
    codon_2 = codon_pair[3:6]

    acid_1 = str(translate_to_protein(codon_1))
    acid_2 = str(translate_to_protein(codon_2))

    acid_pair = acid_1 + acid_2
    #print(acid_pair)
   # print(acid_1)
    #print(acid_2)
   # print(get_single_codon_freq(codon_1, codon_single_dict))
   # print(get_single_codon_freq(codon_2, codon_single_dict))

    codon_single_product = float(get_single_codon_freq(codon_1, codon_single_dict)) * \
                           float(get_single_codon_freq(codon_2, codon_single_dict))
   # print(codon_single_product)
    residue_single_product = float(get_single_residue_freq(acid_1, residue_single_dict)) * \
                             float(get_single_residue_freq(acid_2, residue_single_dict))
   # print("res" + str(residue_single_product))
    residue_pair_freq = float(get_residue_pair_freq(acid_pair, residue_pair_dict))
   # print("res"+str(residue_pair_freq))
    codon_pair_freq = float(get_codon_pair_freq(codon_pair, codon_pair_dict))
  #  print("pair_freq"+str(codon_pair_freq))
    quotient = codon_single_product / residue_single_product
  #  print(quotient)
    product = quotient * residue_pair_freq
  #  print("prod"+str(product))
    quotient = codon_pair_freq / product
   # print(quotient)
    codon_pair_score = np.log(quotient)
   # print(codon_pair_score)

    return float(codon_pair_score)


# =============DEBUG================

string = "AAATTTCCCAAAGGGAAATTTAAGTTT"

protein_seq = translate_to_protein(string)
print(protein_seq)

print(gen_codon_pair_freq(string)[0])
print(gen_codon_pair_freq(string)[1])

codonpair_freq = gen_codon_pair_freq(string)[1]
codonpair_list = gen_codon_pair_freq(string)[0]
print(string)
print(codonpair_list)
print(codonpair_freq)
print(gen_codon_freq(string)[0])
print(gen_codon_freq(string)[1])

codonsingle_freq = gen_codon_freq(string)[1]

print(get_single_codon_freq("AUG", codonsingle_freq))
print(str(protein_seq))
pairs = gen_residue_pair_freq(str(protein_seq))
print(gen_residue_freq(str(protein_seq))[0])
print(gen_residue_freq(str(protein_seq))[1])
print(pairs)

cps = find_codon_pair_score(string, "AAATTT")
print(cps)

import Bio as bp
from Bio.Seq import Seq as sq
from Bio import SeqIO as sqio
import numpy as np
from tqdm import tqdm
import warnings
import xlsxwriter

class CodonAnalyzer:

    def gen_frequency_dict(self, list_iteration_num, substring_length, substring_step, freq_dict=None, parent_string=None):
        """
        Returns a substring list representing the substrings specified by the parameters found in the parent string
        along with a dictionary representing the frequency of each unique substring in the parent string.

        :param parent_string: String from which to generate substring_list
        :param list_iteration_num: Number of times substring_list is iterated through the string
        :param substring_length: Length of each substring element in the substring list
        :param substring_step: Step value used to specify how many characters to avoid adding in next iteration
        :return: substring_list: List of substrings generated from parent_string and parameters
        :return: frequency_dict: Filled dictionary that maps substrings to their frequency in the parent string
        """
        substring_list = []
        frequency_dict = {}
        if freq_dict is not None:
            frequency_dict = freq_dict
        counter = 0
        for i in range(0, list_iteration_num, 1):
            item = parent_string[int(counter):int(counter + substring_length)]
            substring_list.append(item)
            if item not in frequency_dict:
                frequency_dict[item] = 1
            else:
                frequency_dict[item] = frequency_dict[item] + 1
            counter += substring_step
        return substring_list, frequency_dict

    # generates list of codon pairs and a dictionary assigning frequencies to each codon pair
    def gen_codon_pair_freq(self, freq_dict=None, geneStr=None, debug=False):
        """Generates frequency list for codon pairs in geneStr.
            :param freq_dict: Frequency dictionary to modify. If none is provided, a new frequency dictionary is returned
            :param geneStr: String of genetic code to search through. If none is provided the total genetic sequence
            is used.
            :param debug: Debug options
            :return: Frequency dictionary for codon pairs
        """
        if geneStr is None:
            geneStr = self.gene_string
        if debug:
            print("\n========================================================")
            print("Generating frequency dictionary for codon pairs...")
        codon_pairs = self.gen_frequency_dict(int(len(geneStr) / 3 - 1), 6, 3, freq_dict=freq_dict, parent_string=geneStr)
        if debug:
            print("Frequency dictionary generated!")

            print("Current frequency dictionary: ", codon_pairs[1])
            print("========================================================\n")
        return codon_pairs[1]

    # returns the frequency for a specific codon pair from the frequency dictionary
    def get_codon_pair_freq(self, codon_pair):
        if codon_pair in self.codon_pair_frequency:
            return self.codon_pair_frequency[codon_pair]
        else:
            return 0

    # generates list of singular codons and a dictionary assigning frequencies to each codon
    def gen_codon_freq(self, freq_dict=None, geneStr=None, debug=False):
        """Generates frequency list for codons in geneStr.
            :param freq_dict: Frequency dictionary to modify. If none is provided, a new frequency dictionary is returned
            :param geneStr: String of genetic code to search through. If none is provided the total genetic sequence
            is used.
            :param debug: Debug options
            :return: Frequency dictionary for codons
        """
        if geneStr is None:
            geneStr = self.gene_string
        if debug:
            print("\n========================================================")
            print("Generating frequency dictionary for codons...")

        codons = self.gen_frequency_dict(int(len(geneStr) / 3), 3, 3, freq_dict=freq_dict, parent_string=geneStr)
        if debug:
            print("Frequency dictionary generated!")

            print("Current frequency dictionary: ", codons[1])
            print("========================================================\n")
        return codons[1]

    # returns the frequency for a specific codon from the frequency dictionary
    def get_single_codon_freq(self, codon):
        if codon in self.codon_frequency:
            return self.codon_frequency[codon]
        else:
            return 0

    # generates list of amino acid pairs and a dictionary assigning frequencies to each amino acid pair
    def gen_residue_pair_freq(self, freq_dict=None, resStr=None, debug=False):
        """Generates frequency list for residues/amino acids in resStr.
            :param freq_dict: Frequency dictionary to modify. If none is provided, a new frequency dictionary is returned
            :param resStr: String of residues to search through. If none is provided the total residue sequence
            is used.
            :param debug: Debug options
            :return: Frequency dictionary for residue pairs
        """
        if resStr is None:
            resStr = self.residue_string
        if debug:
            print("\n========================================================")
            print("Generating frequency dictionary for residue pairs...")
        residue_pairs = self.gen_frequency_dict(int(len(resStr) - 1), 2, 1, freq_dict=freq_dict, parent_string=resStr)
        if debug:
            print("Frequency dictionary generated!")

            print("Current frequency dictionary: ", residue_pairs[1])
            print("========================================================\n")
        return residue_pairs[1]

    # returns the frequency of an amino acid pair from the frequency dictionary
    def get_residue_pair_freq(self, residue_pair):
        if residue_pair in self.residue_pair_frequency:
            return self.residue_pair_frequency[residue_pair]
        else:
            return 0

    # generates list of singular amino acids and a dictionary assigning frequencies to each amino acid
    def gen_residue_freq(self, freq_dict=None, resStr=None, debug=False):
        """Generates frequency list for residues/amino acids pairs in resStr.
            :param freq_dict: Frequency dictionary to modify. If none is provided, a new frequency dictionary is returned
            :param resStr: String of residues to search through. If none is provided the total residue sequence
            is used.
            :param debug: Debug options
            :return: Frequency dictionary for residues
        """
        if resStr is None:
            resStr = self.residue_string
        if debug:
            print("\n========================================================")
            print("Generating frequency dictionary for singular residues...")
        residues = self.gen_frequency_dict(int(len(resStr)), 1, 1, freq_dict=freq_dict, parent_string=resStr)
        if debug:
            print("Frequency dictionary generated!")

            print("Current frequency dictionary: ", residues[1])
            print("========================================================\n")
        return residues[1]

    # returns the frequency of an amino acid from the frequency dictionary
    def get_single_residue_freq(self, residue):
        if residue in self.residue_frequency:
            return self.residue_frequency[residue]
        else:
            return 0

    # translates string of genetic code into the amino acids it codes for
    def translate_to_protein(self, gene_string):
        sequence = sq(gene_string)
        return str(sequence.translate())

    # Returns the codon pair score for a given codon pair from a string of genetic code
    def find_codon_pair_score(self, codon_pair):
        """
        Returns the codon pair score for a given codon pair
        :param codon_pair: Codon pair to find score for
        :return: Calculated codon pair score for codon_pair
        """
        codon_1 = codon_pair[0:3]
        codon_2 = codon_pair[3:6]

        acid_1 = str(self.translate_to_protein(codon_1))
        acid_2 = str(self.translate_to_protein(codon_2))

        acid_pair = acid_1 + acid_2
        codon_single_product = float(self.get_single_codon_freq(codon_1)) * \
                               float(self.get_single_codon_freq(codon_2))
        residue_single_product = float(self.get_single_residue_freq(acid_1)) * \
                                 float(self.get_single_residue_freq(acid_2))
        residue_pair_freq = float(self.get_residue_pair_freq(acid_pair))
        codon_pair_freq = float(self.get_codon_pair_freq(codon_pair))
        quotient = codon_single_product / residue_single_product
        product = quotient * residue_pair_freq
        quotient = codon_pair_freq / product
        codon_pair_score = np.log(quotient)

        return float(codon_pair_score)

    def find_codon_pair_bias(self, locus_tag):
        """
        Finds the codon pair bias for a given locus/gene
        :param locus_tag: Name of gene/locus tag in gene_dictionary
        :return: Codon pair bias for given gene/locus tag
        """
        codon_pair_bias = 0
        gene_seq = self.gene_dictionary[locus_tag]
        frequency_dictionary = self.gen_codon_pair_freq(geneStr=gene_seq)
        for codon_pair in frequency_dictionary.keys():
            codon_pair_bias += self.find_codon_pair_score(codon_pair)
        return codon_pair_bias

    def gen_codon_pair_bias_dictionary(self, debug=False):
        """
        Generates a dictionary mapping genes/locus tags to their codon pair bias
        :param debug: Debug options
        :return: Dictionary with genes/locus tags as keys and their corresponding codon pair bias as values
        """
        codon_pair_dict = {}
        print("\nFinding codon pair bias for each gene...") if debug else 0
        for gene in (tqdm(self.gene_dictionary.keys()) if debug else self.gene_dictionary.keys()):
            codon_pair_bias = self.find_codon_pair_bias(gene)
            codon_pair_dict[gene] = codon_pair_bias
        print("Codon pair bias dictionary filled.") if debug else 0
        return codon_pair_dict

    def gen_codon_pair_score_dictionary(self, debug=False):
        """
        Generates a dictionary mapping codon pairs in the entire genome to their codon pair score
        :param debug: Debug options
        :return: Dictionary with codon pairs as keys and corresponding codon pair score as values
        """
        codon_pair_score_dict = {}
        print("\nFinding codon pair score for each codon pair...") if debug else 0
        for codon_pair in (tqdm(self.codon_pair_list) if debug else self.codon_pair_list):
            codon_pair_score = self.find_codon_pair_score(codon_pair)
            codon_pair_score_dict[codon_pair] = codon_pair_score
        print("Codon pair score dictionary filled.") if debug else 0
        return codon_pair_score_dict

    def __init__(self, gene_string=None, output_name=None):
        """ Create a Codon Analyzer Object """
        if gene_string != None:
            # Todo make load from string separate method
            self.gene_string = gene_string
            self.residue_string = self.translate_to_protein(gene_string)
    #
            self.codon_frequency = self.gen_codon_freq()
            self.codon_pair_frequency = self.gen_codon_pair_freq()
            self.residue_frequency = self.gen_residue_freq()
            self.residue_pair_frequency = self.gen_residue_pair_freq()
    #
            self.codon_list = list(self.codon_frequency.keys())
            self.codon_pair_list = list(self.codon_pair_frequency.keys())
            self.residue_list = list(self.residue_frequency.keys())
            self.residue_pair_list = list(self.residue_pair_frequency.keys())
            self.gene_dictionary = {"Provided Gene": self.gene_string}
            self.cps_dictionary = self.gen_codon_pair_score_dictionary()
            self.cpb_dictionary = self.gen_codon_pair_bias_dictionary()
        else:
            self.gene_string = ""
            self.residue_string = ""
            #
            self.codon_frequency = {}
            self.codon_pair_frequency = {}
            self.residue_frequency = {}
            self.residue_pair_frequency = {}
            #
            self.codon_list = []
            self.codon_pair_list = []
            self.residue_list = []
            self.residue_pair_list = []
            self.gene_dictionary = {}
            self.cps_dictionary = None
            self.cpb_dictionary = None
        if output_name is None:
            output_name = "CodonAnalyzerOutput.xlsx"
        self.output_workbook = xlsxwriter.Workbook(output_name)
        self.title_format = self.output_workbook.add_format({'bold': True})

        self.record_list = None

    # Todo Work in progress, most likely not needed unless FASTA file is used
    # def readFASTA(self, filename):
    #     sequence_list = sqio.parse(filename, "genbank")
    #     total_str = ""
    #     for seq in sequence_list:
    #         total_str = total_str + str(seq.seq)
    #         print("Sequence: ", seq.seq)
    #     print("Total Sequence: ", total_str)
    #     self.gene_string = total_str
    #     self.residue_string = self.translate_to_protein(total_str)
    #
    #     self.codon_frequency = self.gen_codon_freq()
    #     self.codon_pair_frequency = self.gen_codon_pair_freq()
    #     self.residue_frequency = self.gen_residue_freq()
    #     self.residue_pair_frequency = self.gen_residue_pair_freq()
    #
    #     self.codon_list = list(self.codon_frequency.keys())
    #     self.codon_pair_list = list(self.codon_pair_frequency.keys())
    #     self.residue_list = list(self.residue_frequency.keys())
    #     self.residue_pair_list = list(self.residue_pair_frequency.keys())

    def read_gen_bank(self, filename, fill_additional_dictionaries=False, debug=False):
        """
        Reads a genBank file and fills frequency dictionaries for the CodonAnalyzer object
        :param filename: Path to genBank file
        :param fill_additional_dictionaries: If true, fills additional dictionaries such as the codon pair score
        dictionary
        :param debug: Debug options
        """
        self.record_list = sqio.parse(filename, "genbank")

        iters = 0
        self.gene_dictionary = {}
        with warnings.catch_warnings():
            warnings.simplefilter('ignore', bp.BiopythonWarning)

            for record in sqio.parse(filename, "genbank"):
                print("Reading Features From Record", record.id) if debug else 0
                total_str = str(record.seq)
                total_res_str = self.translate_to_protein(str(record.seq))
                for feat in (tqdm(record.features) if debug else record.features):
                    if feat.type == "CDS":
                        locus_tag = feat.qualifiers['locus_tag']
                        current_seq = feat.extract(record.seq)

                        translated_current_seq = self.translate_to_protein(str(current_seq))
                        current_seq = str(current_seq)
                        self.gene_dictionary[locus_tag[0]] = current_seq
                        if iters == 0:
                            self.codon_frequency = self.gen_codon_freq(geneStr=current_seq, debug=False)
                            self.codon_pair_frequency = self.gen_codon_pair_freq(geneStr=current_seq, debug=False)
                            self.residue_frequency = self.gen_residue_freq(resStr=translated_current_seq, debug=False)
                            self.residue_pair_frequency = self.gen_residue_pair_freq(resStr=translated_current_seq, debug=False)

                        else:
                            self.codon_frequency = self.gen_codon_freq(freq_dict=self.codon_frequency, geneStr=current_seq, debug=False)
                            self.codon_pair_frequency = self.gen_codon_pair_freq(freq_dict=self.codon_pair_frequency, geneStr=current_seq, debug=False)
                            self.residue_frequency = self.gen_residue_freq(freq_dict=self.residue_frequency, resStr=translated_current_seq, debug=False)
                            self.residue_pair_frequency = self.gen_residue_pair_freq(freq_dict=self.residue_pair_frequency, resStr=translated_current_seq, debug=False)
                        iters += 1
            self.gene_string = total_str
            self.residue_string = total_res_str
            self.codon_list = list(self.codon_frequency.keys())
            self.codon_pair_list = list(self.codon_pair_frequency.keys())
            self.residue_list = list(self.residue_frequency.keys())
            self.residue_pair_list = list(self.residue_pair_frequency.keys())
            if fill_additional_dictionaries:
                self.cps_dictionary = self.gen_codon_pair_score_dictionary(debug=debug)
                self.cpb_dictionary = self.gen_codon_pair_bias_dictionary(debug=debug)

    def write_dictionary_to_xls(self, dictionary, key_title, value_title, close_workbook=False, debug=False):
        """
        Writes a given dictionary to a worksheet in the output xls workbook
        :param dictionary: Dictionary to write
        :param key_title: Title of key column
        :param value_title: Title of value column
        :param close_workbook: Set to true if this workbook is to be closed and the file written. False otherwise.
        """
        worksheet = self.output_workbook.add_worksheet()
        worksheet.set_column('A:A', 20)
        worksheet.set_column('B:B', 20)
        worksheet.write(0, 0, key_title, self.title_format)
        worksheet.write(0, 1, value_title, self.title_format)
        indx = 1
        for key in (tqdm(dictionary.keys()) if debug else dictionary.keys()):
            worksheet.write(indx, 0, key)
            worksheet.write(indx, 1, dictionary[key])
            indx += 1
        if close_workbook:
            self.output_workbook.close()

    def write_all_to_xls(self, write_additional_dictionaries=False, debug=False):
        """
        Writes all dictionaries to xls file using predetermined titles, then closes it.
        :param write_additional_dictionaries: Writes additional dictionaries to xls file
        """
        dictionary_list = []
        title_list = []
        if write_additional_dictionaries:
            dictionary_list.extend([self.cps_dictionary, self.cpb_dictionary])
            title_list.extend([["Codon Pair", "Codon Pair Score"], ["Gene Name", "Codon Pair Bias"]])
        dictionary_list.extend([self.codon_frequency, self.codon_pair_frequency, self.residue_frequency, self.residue_pair_frequency])
        title_list.extend([["Codon", "Frequency"], ["Codon Pair", "Frequency"], ["Residue", "Frequency"], ["Residue Pair", "Frequency"]])
        print("\nWriting dictionaries to xlsx file...") if debug else 0
        for frequency_dictionary_idx in range(len(dictionary_list)):
            self.write_dictionary_to_xls(dictionary_list[frequency_dictionary_idx],
                                         title_list[frequency_dictionary_idx][0], title_list[frequency_dictionary_idx][1], debug=True)
            print("Sheet #", frequency_dictionary_idx + 1, "completed\n")
        self.output_workbook.close()
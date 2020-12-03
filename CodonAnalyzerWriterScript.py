from CodonAnalyzerObject import CodonAnalyzer
import argparse
import os

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Writes CodonAnalyzer dictionaries to an xlsx file")
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("-s", "--string", type=str, default=None, help="Genetic Code String, ex: \'-s AAAGGGTTTCCC\'")
    group.add_argument("-gf", "--genbank_file", type=str, default=None, help="Genbank file, ex: \'-gf /path/to/genbank/file.gbff\'")
    parser.add_argument("-fa", "--fill_additional", action='store_true', default=False, help="Fill additional dictionaries like Codon Pair Score and Codon Pair Bias")
    parser.add_argument("-of", "--output_file", type=str, default=None, help="Writes xlsx output "
                                                                      "to filename, ex: \'-of data.xlsx\' "
                                                                      "writes dictionaries to the file \'data.xlsx\'.")
    try:
        args = parser.parse_args()
        CA = CodonAnalyzer()

        if args.string is not None:
            string = args.string
            # Todo additional dictionaries for string
            CA = CodonAnalyzer(gene_string=string)
            if args.output_file is not None:
                CA = CodonAnalyzer(gene_string=string, output_name=args.of)
                CA.write_all_to_xls(debug=True)
        elif args.genbank_file is not None:

            CA = CodonAnalyzer()
            if args.output_file is not None:
                CA = CodonAnalyzer(output_name=args.output_file)
            CA.read_gen_bank(args.genbank_file, fill_additional_dictionaries=args.fill_additional, debug=True)
            CA.write_all_to_xls(write_additional_dictionaries=args.fill_additional, debug=True)
        else:
            parser.print_help()

        if args.output_file is None:
            print("\n Finished writing dictionaries to output file")
        else:
            print("\n Finished writing dictionaries to output file", args.output_file, "")
    except:
        print("\n\nThere was an error with how you ran the script.")
        print("Look at the following for help:")
        print("-------------Help Message-------------\n")
        parser.print_help()


    # print(CA.codon_pair_frequency)
    # string_dict = CA.gen_codon_pair_bias_dictionary(debug=True)
    # print(string_dict)
    # CA.read_gen_bank("gcf.gbff", fill_additional_dictionaries=True, debug=True)
    #
    # # codon_pair_score_dict = CA.gen_codon_pair_score_dictionary(debug=True)
    # # CA.write_dictionary_to_xls(codon_pair_score_dict, "Codon Pair", "Codon Pair Score")
    # # print(codon_pair_score_dict)
    # # codon_pair_bias_dict = CA.gen_codon_pair_bias_dictionary(debug=True)
    # # CA.write_dictionary_to_xls(codon_pair_bias_dict, "Gene Name", "Codon Pair Bias", close_workbook=True)
    # CA.write_all_to_xls(write_additional_dictionaries=True)
    # # print(codon_pair_score_dict)
    # # print(codon_pair_bias_dict)



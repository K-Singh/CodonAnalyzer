import CodonAnalyzer as CodonAnalyzer


string = "AAATTTCCCAAAGGGAAATTTAAGTTT"
CBP = CodonAnalyzer.CodonAnalyzer(string)
protein_seq = CBP.translate_to_protein(string)
print(protein_seq)

print(CBP.codon_pair_frequency)

codonpair_freq = CBP.codon_pair_frequency
codonpair_list = CBP.codon_pair_list
print(string)
print(codonpair_list)
print(codonpair_freq)


codonsingle_freq = CBP.codon_frequency

print(CBP.get_single_codon_freq("AUG"))
print(str(protein_seq))
pairs = CBP.residue_pair_frequency
print(CBP.residue_frequency)
print(CBP.residue_list)
print(pairs)

cps = CBP.find_codon_pair_score("AAATTT")
cpb = CBP.find_codon_pair_bias()
print(cps)
print(cpb)
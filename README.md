# CodonAnalyzer
Returns information and provides functions regarding codons and related genetic material. This program takes either a string or a genbank file as input.

###  Frequency Dictionaries
Frequency dictionaries are the core of CodonAnalyzer. Frequency dictionaries represent python dictionaries that map codons, codon pairs, residues
and residue pairs to their respective frequency within the provided genome. Frequency dictionaries allow for quick access to information about any one piece of genetic data.

### Additional Dictionaries
Additional dictionaries refer to any datasets created through calculations made on any of the original frequency dictionaries. Currently, the only available
additional dictionaries that can be created are the codon pair bias(CPB) dictionary and the codon pair score(CPS) dictionary. The CPB dictionary maps 
genes to their calculated codon pair bias. The CPS dictionary maps codon pairs to their calculated codon pair score.

### CodonAnalyzer Functions
These functions can be accessed from the CodonAnalyzer module and object. They provide access to the various dictionaries and attributes within the CodonAnalyzer object.

### Console/Terminal Script
Allows for users who are only interested in the basic functions of CodonAnalyzer to input a file or string and get back an excel file containing the various dictionaries
generated by CodonAnalyzer.


# Lab 5: ORF Finder and Sequence Analysis Background and Guide

### Background:
Every organism has a genome that contains their entire body of DNA, becoming the blueprint for the organism. DNA (deoxyribonucleic acid) is a  long strand made up of the nucleotide bases: adenine (A), thymine (T), guanine (G), and cytosine (C). These strands contain a region called the open reading frame (ORF) that is able to be translated into an amino acid protein sequence. ORFs are certain sequences of DNA usually with a 5’ start codon and 3’ stop codon. 

For this project, I have utilized ‘ATG’ for the start codon and ‘TAG’, ‘TAA’, and ‘TGA’ for the stop codons. 


Additionally, there are some special cases such as dangling start fragments and stop fragments are used to locate ORFs. This means if a start codon was not found but a stop codon was found, this would be considered a dnagling stop ORF fragment (and vice versa for dangling start fragments). 

I have provided an example of my code should work below:

```python
>exampleguide1
AAA AAA AAA TGA CCC CCC
```
The DNA sequence is first read at frame position 1 and in sets of 3 representing codons. The program will read the sequence until a stop codon is found or it reaches the end of the DNA sequence. Here, 'TGA' is consider a stop codon, thus an ORF is located. This will repeat for each of the 6 frame positions. 

This data will be collected and output in this format:

{frame position} {start position} {end position} {ORF length}

Thus, the example should have an ORF output like this:
```python
+1 1 12 12 
```

### Project Overview:

For this project, I have designed and programmed the class OrfFinder in the file sequenceAnalysis. This class locates all the possible combinations for ORFs that meet the conditions explained above. The purpose of this program is to automate the genome annotation and locate ORFs on DNA sequences. 
With the data collected listed out in my output, this program can be helpful to provide others in the future with information about the location (start and end position) of ORFs on desired DNA sequences, as well as their lengths to be translated. 


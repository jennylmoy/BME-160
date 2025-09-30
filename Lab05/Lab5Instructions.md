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



### Summary of Program:

I have designed and programmed the class OrfFinder in the file sequenceAnalysis. This class is to locate all the possible combinations of ORFs that meet the conditions explained above. 

# Lab 6: Find Unique Substrings of Mitochondrial tRNA

### Background:
Mitochondrial transfer RNA (tRNA) are essential in synthesizing mitochondrial RNA (mRNA) into an amino acid protein chain. THey are crucial to the central dogma for the process of how genetic information is replicated and reproduced in the mitochondria.

Certain diseases can be linked to the mitochondrial genome due to mutations in the gnes of mitochdonrial tRNA (mt.tRNA). It is important to identify these mutations in tRNA genes. One method to ensure that mt.tRNA genes are not mutated is to find unique substrings in their sequence between each of the 22 distinguished tRNAs. This program will utilize sets and dictionaries to accurately sort through unique and essential substrings from the 22 tRNA sequences.

I have provided a small example to showcase how my program should work:

```python
>tRNA1
ACGU
>tRNA2
ACGA
>tRNA3
UUAC
```

First, the program should find any and all possible combinations in the sequence, maintaining the order of the sequence. For tRNA1, tRNA2, and tRNA3 respectively, the power set of all possible combination substings should include:

```python
{"A", "C", "G", "U", "AC", "CG", "GU", "ACG", "CGU", "ACGU"}
{"A", "C", "G", "AC", "CG", "GA", "ACG", "CGA", "ACGA"}
{"U", "A", "C", "UU", "UA", "AC", "UUA", "UAC", "UUAC"}
```

The program then must find unique substrings only from the substrings listed in the power sets, meaning they cannot have repeat substrings across all the tRNA power sets. For example, since all three tRNAs have "A", "C", and "U", these cannot be counted as unique substrings and thus removed. 

```python
{"GU", "CGU", "ACGU"}
{"GA", "CGA", "ACGA"}
{"UU", "UA", "UUA", "UAC", "UUAC"}
```

Lastly, the program locates only the essential substrings from the unique tRNA sets. This means if a substring within a set, say "ACGU" contains another substring listed in the same set, say "GU", the longer substrings are redundant and thus also removed. The final set of essential tRNA substrings will be formatted into this output:

```python
tRNA1
ACGU
..GU

tRNA2
ACGA
...GA

tRNA3
UUAC
UU..
...AC
```

### Project Overview:
For this project, I have designed and programmed the class tRNAFinder in the file findUnique. This class follows the structure outlined above to locate all unique and essential substrings of the 22 tRNA genes. The purpose of this program is to identify the unique substrings within the tRNA gene that differentiates them from the rest, acting as “tags” for each tRNA. This helps immensely with detecting any mutations that can disrupt the processing of genetic information. 


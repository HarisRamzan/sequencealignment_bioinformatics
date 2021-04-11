# -*- coding: utf-8 -*-
"""
Created on Wed Mar 11 01:10:15 2020

@author: haris
"""


#inporting modules and functions
from Bio.Seq import Seq
from Bio.SubsMat import MatrixInfo
from Bio.pairwise2 import format_alignment
from Bio import pairwise2
from Bio import AlignIO
###############################Question 1###########################################
sequence1 = Seq("HEAGAWGHEE")
sequence2 = Seq("PAWHEAE")

#substitution matrices
BLOSUM=MatrixInfo.blosum62
PAM= MatrixInfo.pam120

print("List of all  possible global alignments");
for globalalign in pairwise2.align.globalxx(sequence1,sequence2):   #ALL POSSIBLE GLOBAL ALIGNMENTS
  print(format_alignment(*globalalign))
print("List of all  possible local alignments");
for localalign in pairwise2.align.localxx(sequence1,sequence2):   #ALL POSSIBLE LOCAL ALIGNMENTS
  print(format_alignment(*localalign))

print("List of blosum62 (global alignments with score)");
for blosum in pairwise2.align.globaldx(sequence1,sequence2,BLOSUM):  #GLOBAL ALIGNMENT WITH BLOSUM62 SUBSTITUTION MATRIX AND SCORE
 print(format_alignment(*blosum))
print("List of pam120 (global alignments with score)");
for pam in pairwise2.align.globaldx(sequence1,sequence2,PAM):  #GLOBAL ALIGNMENT WITH PAM120 SUBSTITUTION MATRIX AND SCORE
 print(format_alignment(*pam))
 
print("List of blosum62 (local alignments with score)");
for blosum in pairwise2.align.localdx(sequence1,sequence2,BLOSUM):  #LOCAL ALIGNMENT WITH BLOSUM62 SUBSTITUTION MATRIX AND SCORE
 print(format_alignment(*blosum))
 
print("List of pam120 (local alignments with score)");
for pam in pairwise2.align.localdx(sequence1,sequence2,PAM):  #LOCAL ALIGNMENT WITH BLOSUM62 SUBSTITUTION MATRIX AND SCORE
 print(format_alignment(*pam))
 
###############################Question 2###########################################

sequences = AlignIO.read(open("PF18815_seed.txt"), "stockholm") #READING FILE
print(sequences)
print("Total sequences:",len(sequences))
for sequence in sequences:
    print("Sequence:",sequence.seq); #length of sequences
    print("Length: ",len(sequence.seq))
    
    
    
#link to download file http://pfam.xfam.org/family/PF18815#tabview=tab3 
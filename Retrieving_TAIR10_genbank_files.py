############################
from Bio import SeqIO
#from Bio.Seq import Seq
from BCBio.GFF import GFFExaminer
from BCBio import GFF
#from Bio.Alphabet import IUPAC
import pprint
import csv
import re
import sys
import string
import os

from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_dna
from Bio import SeqFeature
from Bio.SeqFeature import FeatureLocation
from Bio.SeqFeature import SeqFeature
from Bio import Entrez
from httplib import HTTPException
import subprocess
#import csv
#import fasta
#import string
#import sys
#import re
#import random
#import os
#import math
##from sets import Set
#from scipy import stats
#from collections import Counter
#from Bio.Alphabet import IUPAC
#from Bio import Entrez
#from Bio.Seq import Seq
#from Bio import SeqFeature
#from Bio import SeqRecord
#from Bio.SeqFeature import FeatureLocation
#from Bio import SeqIO
#from httplib import HTTPException

# <codecell>
#####Retreive genbank records for each chromosome of TAIR10
#gid_dict = {1:"240254421", 2:"240254678", 3:"240255695", 4:"240256243", 5:"240256493"} #GI ids for TAIR10
gid_dict = {2:"240254678"} #temporary for debugging
#
#for c in gid_dict:
#    #print "Processing records for chromosome", c
#    Entrez.email = "trudi.gulick001@umb.edu"
#    filename = "GI" + str(gid_dict[c]) + ".gbk"
#    try: #adapted from http://biopython.org/pipermail/biopython/2011-October/007555.html
#        if not os.path.isfile(filename):
#            net_handle = Entrez.efetch(db="nucleotide", rettype="gb", retmode="text", id=gid_dict[c])
#            out_handle = open(filename, "w")
#            out_handle.write(net_handle.read())
#            out_handle.close()
#            net_handle.close()
#            print "Saved gb record for chromosome", c
#                        
#                
#    except HTTPException, e:
#        print "Network problem: %s" % e
#        print "Second (and final) attempt..."
#        if not os.path.isfile(filename):
#            net_handle = Entrez.efetch(db="nucleotide", rettype="gb", retmode="text", id=gid_dict[c])
#            out_handle = open(filename, "w")
#            out_handle.write(net_handle.read())
#            out_handle.close()
#            net_handle.close()
#            print "Saved gb record for chromosome", c
#


#need to add 
cds_index = []

for c in gid_dict: 
    print "Processing chrom", c
    filename = "GI" + str(gid_dict[c]) + ".gbk"
    print filename
    for gb_record in SeqIO.parse(open(filename, "r"), "gb"): #based on http://www2.warwick.ac.uk/fac/sci/moac/people/students/peter_cock/python/genbank/
        for feature in gb_record.features:
            print feature.type
            if feature.type == "CDS":
                if "db_xref" in feature.qualifiers:
                    print feature.qualifiers["db_xref"]
                    #write code to pull out "TAIR" and write the corresponding AT id's to file
                    
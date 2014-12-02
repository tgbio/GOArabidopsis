# -*- coding: utf-8 -*-
"""
Created on Tue Dec  2 17:09:36 2014

@author: tg
"""

from Bio import SeqIO
#from Bio.Seq import Seq

from Bio import SeqFeature
#from Bio.SeqFeature import FeatureLocation
#from Bio.SeqFeature import SeqFeature
#from Bio import Entrez
#from httplib import HTTPException
#import subprocess

#####Retreive genbank records for each chromosome of TAIR10
#you'll need a good internet conection for this to work properly
#once you've downloaded the genbank files you do not need to repeat this section

gid_dict = {1:"240254421", 2:"240254678", 3:"240255695", 4:"240256243", 5:"240256493"} #GI ids for TAIR10

for c in gid_dict:
    #print "Processing records for chromosome", c
    Entrez.email = "trudi.gulick001@umb.edu"
    filename = "GI" + str(gid_dict[c]) + ".gbk"
    try: #adapted from http://biopython.org/pipermail/biopython/2011-October/007555.html
        if not os.path.isfile(filename):
            net_handle = Entrez.efetch(db="nucleotide", rettype="gb", retmode="text", id=gid_dict[c])
            out_handle = open(filename, "w")
            out_handle.write(net_handle.read())
            out_handle.close()
            net_handle.close()
            print "Saved gb record for chromosome", c
                        
                
    except HTTPException, e:
        print "Network problem: %s" % e
        print "Second (and final) attempt..."
        if not os.path.isfile(filename):
            net_handle = Entrez.efetch(db="nucleotide", rettype="gb", retmode="text", id=gid_dict[c])
            out_handle = open(filename, "w")
            out_handle.write(net_handle.read())
            out_handle.close()
            net_handle.close()
            print "Saved gb record for chromosome", c
#
for seq_record in SeqIO.parse("GI240255695.gbk", "genbank"):
    print seq_record
    print len(seq_record.features)
    #AT IDs should be in the locus tag feature
    #if locus tag in list of AT IDS
    #write gene id to a file 
        

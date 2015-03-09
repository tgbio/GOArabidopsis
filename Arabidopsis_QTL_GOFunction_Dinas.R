#source("http://bioconductor.org/biocLite.R")
#biocLite("GOFunction")
library('GOFunction')
#source("http://bioconductor.org/biocLite.R")
#biocLite("org.At.tair.db")
library('org.At.tair.db')

#retrieve mapping between TAIR ID and Entrez Gene identifiers
all.at.genes <- org.At.tairENTREZID #maps between TAIR acc. #s and NCBI Entrez Gene IDs
all.at.genes[1:5]#slices the map to just the first 5 entries (?)
mapped.at.genes = mappedkeys(all.at.genes) #limits to only those genes that have a mapping
mapped.list <- as.list(all.at.genes[mapped.at.genes]) # within all at genes, only those that are mapped, make into list

if(length(mapped.list) > 0) { #ensures that the mapped list contains values
  mapped.list[1:5] # Get the Entrez gene IDs for the first five genes
  str(mapped.list[1:5])  #the structure of mapped.list
  
  ###########FOR TESTING
  #get the gene id corresponding to AT id "AT1G01010"
  #similar to a Python dictionary
  mapped.list["AT1G01010"]######END TEST
} 

qtl_at_id_file = list("CORRECTEDTESTLOCALBLAST_QGG24G02.yg.ab1_CLS_S3_Contig10858_LettuceNCBI.csv")#names the file of at qlt

at_ids = read.csv("CORRECTEDTESTLOCALBLAST_QGG24G02.yg.ab1_CLS_S3_Contig10858_LettuceNCBI.csv",header=F)$V13 # read AT IDs of protein coding genes within region into a vector
#This cannot take a variable as an argument
gi_ids = "" #creates an empty list


#there appears to be something wrong here with the variables, confusion of names
#this section uses the at id's identified for the region and compiles a list of entrez gene id's with them via the mapping
for (i in 1:length(at_ids)){ #for variable within the length range of at_ids
  for (k in mapped.list) { #the GI from mapped list look 
    print(k) #very different from the GI from the at_ids/CSV file...
    print(at_ids[i])
    if(k==at_ids[i])#############SAY "GI" %in% mapped.list[values] #### basically going the other way so asking is this Gene ID in the values of AT ID mapped list? if so what is the key [the at id]
      gi_ids = c(gi_ids,k)}}
  
gi_char = as.character(gi_ids) #converts the list of gene ids to character vector 
#this is done because the factor values need to be characters
gi_fac = as.factor(gi_char)#converts the character list of GIs to a factor list
##########this is done because (?) ----- we never actually end up using gi_ids (or gi_fac) again

ont_type = list("MF", "BP", "CC") #defines a list of ontology types
for (q in qtl_at_id_file) { # loop over each csv filename in list qtl_at_id_file
  qtl_genes = read.csv(q,header=F)$V13 # read AT IDs of protein coding genes within QTL region into a vector
  #str(qtl_genes)      ##########Not sure why these are here
  #qtl_genes = gi_ids  ##########Not sure why these are here
  for (o in ont_type) { # for each type of ontology
    print(o) #state which ontology is being tested
    region_fragment = substr(q, 16, regexpr("csv", q, fixed=T)[1]-2) # return the region name from the csv filename
    print(region_fragment) 
    outfilename = paste(region_fragment, "_GO_", o, sep = "") # generate a unique filename for the GO output
    print(outfilename)
    num_gen = length(qtl_genes) # count the number of genes in QTL region
    print("Number of Genes:")
    print(num_gen) #state that
    print(qtl_genes) ############IS this needed? It does help show that the code is working/progressing but beyond that?
    sigTerm = GOFunction(qtl_genes, mapped.at.genes, organism="org.At.tair.db", ontology=o, fdrmethod="BY", fdrth=0.05, ppth=0.05, pcth=0.05, poth=0.05, peth=0.05, bmpSize=2000, filename=outfilename) # find GO terms enriched in region
    if(length(sigTerm) > 0){ #if the significant term is non-zero [ie if it is found to be significant]
      print("sigTerm") #print it
      print(sigTerm)
      
      #adapted from the org.At.tair.db vignette for org.At.tairGO2ALLTAIRS
      #provides mappings between a given GO identifier and all TAIR identifiers annotated at that GO term or one of its children in the GO ontology
      #mapped.list is a list with names = GO ids, values = a named list with values of AT ids and names of evidence codes
      
      csv.out = paste(outfilename, "_GO_signif_ATIDs.csv", sep = "") #write an output file
      print(csv.out) #print what is contained
      mapped.list <- as.list(org.At.tairGO2ALLTAIRS) #map from TAIRGO to AT IDs
      # Gets the tair identifiers for the significant GO term
      for (i in 1:nrow(sigTerm)){ #for each GO sigterm in the list written
        sindex = grep(as.character(sigTerm[[i,1]]), names(mapped.list), fixed=TRUE, value = FALSE)#index the significant GO term
        goid <- mapped.list[sindex] #convert to list
        stopifnot(length(goid)==1) #if there is no signif term, this stops this code
        # Gets all the tair identifiers for significant GO term
        at_id_char = as.character(goid[[1]]) #convert to character
        at_id_list = as.list(at_id_char) #convert from character to list terms
        gene_list = at_id_list[at_id_list %in% qtl_genes]# Select just those tair identifiers also in qtl_genes
        out_df = c(as.character(sigTerm[[i,1]]), gene_list) #combine the sigterms and the list of overlapping genes
        print(str(out_df))
        write.table(out_df, csv.out, append=TRUE,col.names=FALSE, row.names=FALSE, sep=",") #write the results to an output file
      }
    }
  }
} #close and complete



library('GOFunction')
library('org.At.tair.db')

#retrieve mapping between TAIR ID and Entrez Gene identifiers
all.at.genes <- org.At.tairENTREZID
all.at.genes[1:5]

#get the ATids that are mapped to gene ids
mapped.at.genes <- mappedkeys(all.at.genes)

# Convert to a list
mapped.list <- as.list(all.at.genes[mapped.at.genes])


if(length(mapped.list) > 0) {
  # Get the Entrez gene IDs for the first five genes
  mapped.list[1:5]
  
  #the structure of mapped.list
  str(mapped.list[1:5])
  
  #get the gene id corresponding to AT id "AT1G01010"
  #similar to a Python dictionary
  mapped.list["AT1G01010"]
}


qtl_at_id_file = list("AT_ids_between_markers_example.csv")

#
at_ids = read.table("AT_ids_between_markers_example.csv",header=F)$V1 # read AT IDs of protein coding genes within region into a vector

gi_ids = ""


#there appears to be something wrong here with the variables, confusion of names
for (i in 1:length(at_ids)){
  gid = mapped.list[at_ids[i]]
  print(gid)
  gi_ids = append(gi_ids, gid)
}

gi_ids = as.character(gi_ids)
gi_ids = as.factor(gi_ids)

ont_type = list("MF", "BP", "CC")

for (q in qtl_at_id_file) { # loop over each csv filename in list qtl_at_id_file
  qtl_genes = read.csv(q,header=F)$V1 # read AT IDs of protein coding genes within QTL region into a vector
  #str(qtl_genes)
  q#tl_genes = gi_ids
  
  for (o in ont_type) { # for each type of ontology
    print(o)
    region_fragment = substr(q, 16, regexpr("csv", q, fixed=T)[1]-2) # return the region name from the csv filename
    print(region_fragment)
    outfilename = paste(region_fragment, "_GO_", o, sep = "") # generate a unique filename for the GO output
    print(outfilename)
    num_gen = length(qtl_genes) # count the number of genes in QTL region
    print("Number of Genes:")
    print(num_gen)
    print(qtl_genes)
    sigTerm = GOFunction(qtl_genes, mapped.at.genes, organism="org.At.tair.db", ontology=o, fdrmethod="BY", fdrth=0.05, ppth=0.05, pcth=0.05, poth=0.05, peth=0.05, bmpSize=2000, filename=outfilename) # find GO terms enriched in region
    if(length(sigTerm) > 0){
      print("sigTerm")
      #print(sigTerm)
      
      #adapted from the org.At.tair.db vignette for org.At.tairGO2ALLTAIRS
      #provides mappings between a given GO identifier and all TAIR identifiers annotated at that GO term or one of its children in the GO ontology
      #mapped.list is a list with names = GO ids, values = a named list with values of AT ids and names of evidence codes
      
      csv.out = paste(outfilename, "_GO_signif_ATIDs.csv", sep = "")
      print(csv.out)
      mapped.list <- as.list(org.At.tairGO2ALLTAIRS)
      for (i in 1:nrow(sigTerm)){
        # Gets the tair identifiers for the significant GO term
        sindex = grep(as.character(sigTerm[[i,1]]), names(mapped.list), fixed=TRUE, value = FALSE)
        goid <- mapped.list[sindex]
        stopifnot(length(goid)==1)
        # Gets all the tair identifiers for significant GO term
        at_id_char = as.character(goid[[1]])
        at_id_list = as.list(at_id_char)
        # Select just those tair identifiers also in qtl_genes
        gene_list = at_id_list[at_id_list %in% qtl_genes]
        out_df = c(as.character(sigTerm[[i,1]]), gene_list)
        print(str(out_df))
        write.table(out_df, csv.out, append=TRUE,col.names=FALSE, row.names=FALSE, sep=",")
      }
    }
  }
}



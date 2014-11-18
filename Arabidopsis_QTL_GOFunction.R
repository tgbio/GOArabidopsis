library('GOFunction')
library('org.At.tair.db')

#retrieve mapping between TAIR ID and Entrez Gene identifiers
all.at.genes <- org.At.tairENTREZID

#converts mapped ids to a list
mapped.at.genes <- mappedkeys(all.at.genes)
#qtl_at_ids = list("AT_ids_between_CH.215L_and_HH.360L.csv", "AT_ids_between_DF.408C_and_FD.90L.csv")

qtl_at_ids = list("AT_ids_between_F4I1_and_MSAT2.22.csv", "AT_ids_between_GH.473C_and_HH.480C.csv", "AT_ids_between_PVV4_and_HH.335C.csv")
ont_type = list("MF", "BP", "CC")

for (q in qtl_at_ids) { # loop over each csv filename in list qtl_at_ids
  qtl_genes = read.csv(q,header=F)$V1 # read AT IDs of protein coding genes within QTL region into a vector
  for (o in ont_type) { # for each type of ontology
    print(o)
    region_fragment = substr(q, 16, regexpr("csv", q, fixed=T)[1]-2) # return the region name from the csv filename
    print(region_fragment)
    outfilename = paste(region_fragment, "_GO_", o, sep = "") # generate a unique filename for the GO output
    print(outfilename)
    num_gen = length(qtl_genes) # count the number of genes in QTL region
    print("Number of Genes:")
    print(num_gen)
    sigTerm = GOFunction(qtl_genes, mapped.at.genes, organism="org.At.tair.db", ontology=o, fdrmethod="BY", fdrth=0.05, ppth=0.05, pcth=0.05, poth=0.05, peth=0.05, bmpSize=2000, filename=outfilename) # find GO terms enriched in region
    if(length(sigTerm) > 0){
      print(sigTerm)
      
      #adapted from the org.At.tair.db vignette for org.At.tairGO2ALLTAIRS
      # provides mappings between a given GO identifier and all TAIR identifiers annotated at that GO term or one of its children in the GO ontology
      #xx is a list with names = GO ids, values = a named list with values of AT ids and names of evidence codes
      
      csv.out = paste(outfilename, "_GO_signif_ATIDs.csv", sep = "")
      print(csv.out)
      xx <- as.list(org.At.tairGO2ALLTAIRS)
      for (i in 1:nrow(sigTerm)){
        # Gets the tair identifiers for the significant GO term
        sindex = grep(as.character(sigTerm[[i,1]]), names(xx), fixed=TRUE, value = FALSE)
        goid <- xx[sindex]
        stopifnot(length(goid)==1)
        # Gets all the tair identifiers for significant GO term
        at_id_char = as.character(goid[[1]])
        at_id_list = as.list(at_id_char)
        # Select just those tair identifiers also in qtl_genes
        gene_list = at_id_list[at_id_list %in% qtl_genes]
        out_df = c(as.character(sigTerm[[i,1]]), gene_list)
        print(str(out_df))
        write.table(out_df, csv.out, append=TRUE,col.names=FALSE, row.names=FALSE, sep=",")
      }}
  }
}






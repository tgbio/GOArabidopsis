#Important Links: http://www.bioconductor.org/packages/release/bioc/manuals/GOFunction/man/GOFunction.pdf

library('GOFunction')
library('org.At.tair.db')

blastfile = "abc_transporter.csv" #Specifies the file for use
blastrecord = read.csv(blastfile, header=F) #Reads the file
#species = blastrecord$V18 #Picks only the column for the species name
#atonly = blastrecord[species=="Arabidopsis thaliana", ] #asserts that the species is Arabidopsis
atonlygeneIDs = as.character(blastrecord$V2) #Picks out the gene IDs of A. thaliana records
#NOTE: For testing, V1 is used - experimental code requires V13 & selection for AT

#retrieve mapping between TAIR ID and Entrez Gene identifiers then converts mapped ids to a list
mapped.at.genes <- mappedkeys(org.At.tairENTREZID)

ont_type = list("MF", "BP", "CC")

for (o in ont_type) {# for each type of ontology
  print(o)
  outfilename = "positive_control_results.csv" # generate a unique filename for the GO output
  print(outfilename)
  num_gen = length(atonlygeneIDs) # count the number of genes in QTL region
  print("Number of Genes:")
  print(num_gen)
  ####line below: do with test IDs
  sigTerm = GOFunction(atonlygeneIDs, mapped.at.genes, organism="org.At.tair.db", ontology=o, fdrmethod="BY", fdrth=0.05, ppth=0.05, pcth=0.05, poth=0.05, peth=0.05, bmpSize=2000, filename=outfilename) # find GO terms enriched in region
  print(sigTerm)
  if(length(sigTerm) > 0){
    print(sigTerm)

atonlyTAIR <- mappedkeys(org.At.tairGO2TAIR)
atonlyTGID <- mappedkeys(org.At.tairENTREZID)

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
}}}
#LOOK FOR LINKAGE BETWEEN TAIR IDENTIFIERS and GENE IDs; between GI, and GO term
#look at documentation
#use p10 in the bioonductor viniette for the code on how to pull these out
#go to NCBI, find a AT gene, take the accession#/gene id list and run as a positive control
#output will be a list of the enriched genes, that have been put into a 



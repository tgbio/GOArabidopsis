library('GOFunction')
library('org.At.tair.db')

blastfile = "CORRECTEDTESTLOCALBLAST_QGG24G02.yg.ab1_CLS_S3_Contig10858_LettuceNCBI.csv"
blastrecord = read.csv(blastfile, header=F)
species = blastrecord$V18
atonly = blastrecord[species=="Arabidopsis thaliana",]
interestGenes = as.character(atonly$V14)


#retrieve mapping between TAIR ID and Entrez Gene identifiers then converts mapped ids to a list
refGenes <- mappedkeys(org.At.tairENTREZID)

#GO function takes GeneIDs so no conversion required

ont_type = list("MF", "BP", "CC")

for (o in ont_type){# for each type of ontology
  print(o)
  outfilename = "QGG24G02.yg.ab1_CLS_S3_Contig10858_LettuceNCBI_Athaliana" # generate a unique filename for the GO output
  print(outfilename)
  num_gen = length(interestGenes) # count the number of genes in QTL region
  print("Number of Genes:")
  print(num_gen)
  ####line below: do with test IDs
  sigTerm = GOFunction(interestGenes, refGenes, organism="org.At.tair.db", ontology=o, fdrmethod="BY", fdrth=0.05, ppth=0.05, pcth=0.05, poth=0.05, peth=0.05, bmpSize=2000, filename=outfilename) # find GO terms enriched in region
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
}}
#LOOK FOR LINKAGE BETWEEN TAIR IDENTIFIERS and GENE IDs; between GI, and GO term
#look at documentation
#use p10 in the bioonductor viniette for the code on how to pull these out
#go to NCBI, find a AT gene, take the accession#/gene id list and run as a positive control
#output will be a list of the enriched genes, that have been put into a 



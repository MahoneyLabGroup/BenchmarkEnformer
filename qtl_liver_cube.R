library(qtl2)
load("~/Downloads/dataset.DO.CUBE.multissue.RData") 

data <- dataset.DO.Cube.Liver
data_annot<-as.data.frame( data$annot.samples)

get_dat_dummpy<-function(dat, name){
  dat_unique = unique(dat)
  num_column = dim(dat_unique)[1] - 1
  num_row = dim(dat)[1]
  dat_dummy = matrix(0, nrow = num_row, ncol = num_column)
  print(dat_unique)
  col_names = c()
  for (i in 1:num_column) {
    d = dat_unique[[1]][i]
    print(d)
    coln = paste0(name, d)
    dat_dummy[dat == d,i] = 1
    col_names = append(col_names, coln)
  }
  colnames(dat_dummy) = col_names
  rownames(dat_dummy) = rownames(dat)
  return(dat_dummy)
}


###Geneate covar matrix 

sex <- data_annot[, "Sex", drop=FALSE]
rownames(sex) =  data_annot[,"mouse.id"]

gen = data_annot[,"Generation", drop=FALSE]
rownames(gen) = data_annot[,"mouse.id"]

batches = data_annot[,"Batch", drop=FALSE]
rownames(batches) = data_annot[,"mouse.id"]

gen = get_dat_dummpy(gen, "gen")
batches = get_dat_dummpy(batches, "batch")

covar = data.frame((sex[, "Sex"] == "M")*1)
colnames(covar) = "sex"
covar = cbind(covar, gen)
covar = cbind(covar, batches)


##### download annotations from biomaRT
library(biomaRt)
ensembl <- useMart("ensembl", dataset="mmusculus_gene_ensembl", host="https://jul2023.archive.ensembl.org") # most recent GRCm39
annot<-getBM(c("ensembl_gene_id", "ensembl_peptide_id", "mgi_symbol", "chromosome_name", "strand", "start_position", "end_position","gene_biotype"), mart=ensembl)
gene_list <- read.csv("~/Downloads/gene_list", header=FALSE) ##  load gene list from pseudo references
colnames(gene_list) = c("gene_id")


print(paste("# of genes in pseudo references: ", length(gene_list$gene_id)))
gene_list2 <- intersect(gene_list$gene_id, dimnames(data$data$rz)[[2]] )  # find common genes between pseudo reference and liver data
print(paste("# of genes in pseudo references and data: ", length(gene_list2)))
annot2<- annot[annot$ensembl_gene_id %in% gene_list2, ]
print(paste("# of genes in pseudo references, data and annotations: ", length(annot2$ensembl_gene_id)))
annot2 = annot2[annot2$chromosome_name != "Y", ]
print(paste("# of genes after removing genes on chrY", length(annot2$ensembl_gene_id)))

#### Calculate allele effect 
res_coef = list()
for (gene_id in annot2$ensembl_gene_id){
  #print(gene_id)
  gene_info <- annot2[annot2$ensembl_gene_id == gene_id, ][1, ]
  chr <- gene_info$chr
  pos <- gene_info$start_position/1e6
  marker_id <- find_marker(map, chr, pos = pos)
  f = fit1(genoprobs[[chr]][,,marker_id], data$data$norm[,gene_id], addcovar = covar)
  ### 
  
  res_coef[[gene_id]] = append(f$lod, f$coef[c("A", "B", "C", "D", "E", "F", "G", "H")])
}
res_coef <- do.call("rbind", res_coef)
colnames(res_coef) <- c("lod", "A_J.39" ,"C57BL_6J.39",
                        "129S1_SvImJ.39", "NOD_ShiLtJ.39",
                        "NZO_HlLtJ.39", "CAST_EiJ.39",
                        "PWK_PhJ.39", "WSB_EiJ.39")

write.csv(res_coef, "~/Desktop/qtlcoef_liver_CUBE.csv") 
head(res_coef)

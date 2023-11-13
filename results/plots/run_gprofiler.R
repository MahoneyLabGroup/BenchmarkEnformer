library(gprofiler2)
data1 <- as.list(read.csv("~/Desktop/genelist_1.csv", header = FALSE))$V1
data1[1:10]

gostres1 = gost(query = data1,
               organism = "mmusculus", 
               numeric_ns="ENTREZGENE_ACC",
               domain_scope= "annotated",
               highlight = TRUE,
               user_threshold = 0.05,
               correction_method = "g_SCS"
               )
res1 = gostres1$result
res1 = res1[res1$highlighted,]
gostplot(gostres1, capped = TRUE, interactive = FALSE)
write.csv(res1,file='~/Desktop/gostres_1.csv', row.names=FALSE, quote = FALSE)


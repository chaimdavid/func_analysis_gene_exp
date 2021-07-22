#Installation and activation of the "gProfileR" package that downloads functional analysis data
#from given databases and attaches them to the appropriate genes
install.packages("gProfileR", repos = "http://cran.us.r-project.org")
library(gProfileR)
funcenr1<-gprofiler(query=as.character(cluster1.genenames),organism= "mmusculus",src_filter=c("GO:BP","KEGG"))
funcenr2<-gprofiler(query=as.character(cluster2.genenames),organism= "mmusculus",src_filter=c("GO:BP","KEGG"))
funcenr3<-gprofiler(query=as.character(cluster3.genenames),organism= "mmusculus",src_filter=c("GO:BP","KEGG"))
funcenr4<-gprofiler(query=as.character(cluster4.genenames),organism= "mmusculus",src_filter=c("GO:BP","KEGG"))
funcenr5<-gprofiler(query=as.character(cluster5.genenames),organism= "mmusculus",src_filter=c("GO:BP","KEGG"))

#Creation of character vectors that contain the "term.name" and "p.value" values of the 
#most significant genes according to given criteria
cluster1.term.name<-funcenr1$term.name[which((funcenr1$term.size<=200)&(funcenr1$p.value<=0.01))]
cluster2.term.name<-funcenr2$term.name[which((funcenr2$term.size<=200)&(funcenr2$p.value<=0.01))]
cluster3.term.name<-funcenr3$term.name[which((funcenr3$term.size<=200)&(funcenr3$p.value<=0.01))]
cluster4.term.name<-funcenr4$term.name[which((funcenr4$term.size<=200)&(funcenr4$p.value<=0.01))]
cluster5.term.name<-funcenr5$term.name[which((funcenr5$term.size<=200)&(funcenr5$p.value<=0.01))]

cluster1.pvalue<-funcenr1$p.value[which((funcenr1$term.size<=200)&(funcenr1$p.value<=0.01))]
cluster2.pvalue<-funcenr2$p.value[which((funcenr2$term.size<=200)&(funcenr2$p.value<=0.01))]
cluster3.pvalue<-funcenr3$p.value[which((funcenr3$term.size<=200)&(funcenr3$p.value<=0.01))]
cluster4.pvalue<-funcenr4$p.value[which((funcenr4$term.size<=200)&(funcenr4$p.value<=0.01))]
cluster5.pvalue<-funcenr5$p.value[which((funcenr5$term.size<=200)&(funcenr5$p.value<=0.01))]

#Creation of the dataframes that contain the gene
ssf1<-data.frame(cluster1.term.name,cluster1.pvalue)
ssf2<-data.frame(cluster2.term.name,cluster2.pvalue)
ssf3<-data.frame(cluster3.term.name,cluster3.pvalue)
ssf4<-data.frame(cluster4.term.name,cluster4.pvalue)
ssf5<-data.frame(cluster5.term.name,cluster5.pvalue)

#Arranging the dataframe with increasing order with the "arrange" function of the "tidyverse" library
library(tidyverse)
ssf1<-ssf1%>%arrange(cluster1.pvalue)
ssf2<-ssf2%>%arrange(cluster2.pvalue)
ssf3<-ssf3%>%arrange(cluster3.pvalue)
ssf4<-ssf4%>%arrange(cluster4.pvalue)
ssf5<-ssf5%>%arrange(cluster5.pvalue)
write.xlsx(ssf1,"4.Cluster 1 Functions.xlsx",col.names=T)
write.xlsx(ssf2,"4.Cluster 2 Functions.xlsx",col.names=T)
write.xlsx(ssf3,"4.Cluster 3 Functions.xlsx",col.names=T)
write.xlsx(ssf4,"4.Cluster 4 Functions.xlsx",col.names=T)
write.xlsx(ssf5,"4.Cluster 5 Functions.xlsx",col.names=T)

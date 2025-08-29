suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tibble))

args<-commandArgs(T)

counts_file<-args[1]
out_dic<-args[2]

dir.create(out_dic, recursive = T, showWarnings = F)

gene<-readRDS("reference/TIME.subtype.markers.rds")
mean_sd_v<-readRDS("reference/marker.mean_sd.rds")

counts_dat<-read.table(counts_file, head=T, row.names = 1, sep="\t", check.names = FALSE)

counts_dat<-counts_dat[intersect(gene$gene, rownames(counts_dat)),]
counts_dat<-counts_dat[rowSums(counts_dat > 0) > 0,]
counts_dat<-apply(counts_dat, 2, function(x){x/sum(x)})

for(i in 1:nrow(counts_dat)){
  gene_i<-rownames(counts_dat)[i]
  mean_i<-mean_sd_v %>% filter(gene==gene_i) %>% pull(mean_value)
  sd_i<-mean_sd_v %>% filter(gene==gene_i) %>% pull(sd_value)
  counts_dat[i,]<-(counts_dat[i,]-mean_i)/sd_i
}

our_res<-NULL
TIME.subtype_v<-c("Cps1","Cps2","Cps3","Cps4","Cps5")
for(cc in TIME.subtype_v){
  gg<-gene %>% filter((TIME.subtype==cc) & (gene %in% rownames(counts_dat))) %>% pull(gene)
  our_res<-rbind(our_res, apply(counts_dat[gg,], 2, mean))
}
rownames(our_res)<-TIME.subtype_v
TIME.subtype<-apply(our_res, 2, which.max) %>% as.data.frame() %>% mutate(TIME.subtype=rownames(our_res)[`.`]) %>% 
  select(TIME.subtype) %>% rownames_to_column(var="sample")

write.table(TIME.subtype, fs::path(out_dic, "TIME.subtype.txt"), row.names = F, quote=F, sep="\t")

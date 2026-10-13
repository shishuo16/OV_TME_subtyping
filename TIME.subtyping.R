suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(GSVA))

args<-commandArgs(T)

tpm_file<-args[1]
out_dic<-args[2]

dir.create(out_dic, recursive = T, showWarnings = F)

gene<-readRDS("reference/TIME.subtype.markers.rds")
mean_sd_v<-readRDS("reference/gene.mean_sd.rds")

tpm_dat<-read.table(tpm_file, head=T, row.names = 1, check.names = FALSE)

tpm_dat<-tpm_dat[intersect(mean_sd_v$gene, rownames(tpm_dat)),]
tpm_dat<-tpm_dat[rowSums(tpm_dat > 0) > 0,]

#scale tpm
mean_sd_v <- data.frame(gene=rownames(tpm_dat)) %>% left_join(mean_sd_v, by="gene")
tpm_dat <- sweep(sweep(as.matrix(tpm_dat), 1, mean_sd_v$mean_value, "-"), 1, mean_sd_v$sd_value, "/")

#clustering
modules <- split(gene$gene, gene$TIME.subtype)

ssgsea_par <- ssgseaParam(
  exprData = as.matrix(tpm_dat),
  geneSets = modules,
  minSize = 5,
  maxSize = 500,
  normalize = TRUE
)
ssgsea_score <- gsva(
  ssgsea_par,
  verbose = TRUE
)
ssgsea_score<-as.data.frame(ssgsea_score)
cluster<-apply(ssgsea_score, 2, which.max) %>% as.data.frame() %>% mutate(TIME.subtype=rownames(ssgsea_score)[`.`]) %>% select(TIME.subtype)
cluster$sample<-rownames(cluster)
cluster<-cluster[,c("sample", "TIME.subtype")]

write.table(cluster, fs::path(out_dic, "TIME.subtype.txt"), row.names = F, quote=F, sep="\t")

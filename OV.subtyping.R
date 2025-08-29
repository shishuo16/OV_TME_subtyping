suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tibble))
suppressPackageStartupMessages(library(GSVA))

args<-commandArgs(T)

tpm_file<-args[1]
out_dic<-args[2]

dir.create(out_dic, recursive = T, showWarnings = F)

gene<-readRDS("reference/OV.subtype.markers.rds")
tpm_ref<-readRDS("reference/TCGA-OV.tpm.rds")
colnames(tpm_ref)<-paste(colnames(tpm_ref), sample(1000:9999, length(colnames(tpm_ref))))

tpm_dat<-read.table(tpm_file, head=T, row.names = 1, sep="\t", check.names = FALSE)
tpm_dat<-tpm_dat[rowSums(tpm_dat)>0,]

tpm_dat<-tpm_dat[intersect(rownames(tpm_ref), rownames(tpm_dat)), , drop=F]
tpm_ref<-tpm_ref[intersect(rownames(tpm_ref), rownames(tpm_dat)), , drop=F]
tpm_dat_all<-cbind(tpm_dat, tpm_ref)
tpm_dat_all<-t(scale(t(tpm_dat_all)))

###GSVA
modules<-split(gene$gene, gene$OV.subtype)

gsvaPar <- gsvaParam(tpm_dat_all, modules, maxDiff = T, kcdf = "Gaussian")
gsva_res <- gsva(gsvaPar)
gsva_res<-as.data.frame(t(as.data.frame(gsva_res)))

OV.subtype<-apply(gsva_res, 2, function(x){ifelse(x<=median(x), "Low" , "High")}) %>% as.data.frame()
OV.subtype<-OV.subtype[colnames(tpm_dat), , drop=F] %>% rownames_to_column(var="sample")

write.table(OV.subtype, fs::path(out_dic, "OV.subtype.txt"), row.names = F, quote=F, sep="\t")


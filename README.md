Methods for tumor immune microenvironment (TIME) subtyping and cancer cell subtyping in ovarian cancer.

## 1. TIME subtyping

```
Rscript TIME.subtyping.R counts.txt out_dir
```

### input
counts.txt: reads counts matrix (gene x sample). Note: The matrix must include marker genes from different subtypes, totaling 175 genes, as specified in the file reference/TIME.subtype.markers.rds.

out_dir: output directory

### output
TIME.subtyping.txt: two columns: sample, TIME.subtype


## 2. Ovarian cancer cell subtyping

```
Rscript OV.subtyping.R tpm.txt out_dir
```

### input
tpm.txt: gene tpm matrix (gene x sample). Note: The matrix must include marker genes from different subtypes, totaling 55 genes, as specified in the file reference/OV.subtype.markers.rds. In addition, at least 500 other genes (as background genes for computation) must also be provided.

out_dir: output directory

### output
OV.subtyping.txt: twelve columns: sample and cell enrichment level of 11 cancer cell subtypes (High/Low)

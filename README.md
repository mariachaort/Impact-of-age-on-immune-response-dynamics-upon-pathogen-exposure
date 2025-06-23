# Impact-of-age-on-immune-response-dynamics-upon-pathogen-exposure

### Cell Type Composition Changes - Figures

This script generates figures showing how cell type composition varies across different stimulation conditions, time points, and chemistries.

**Input:** Preprocessed metadata and cell composition `.rds` files.
**Workflow:**  
  1. Loads metadata to order cell types for plotting.  
  2. Defines experimental parameters and filters valid condition combinations.  
  3. Reads cell composition data for each condition, combining into one dataset.  
  4. Uses `ggplot2` to create faceted plots highlighting significant composition changes (p < 0.05).  
  5. Saves plots as PDFs.

**Packages:** `ggplot2`, `dplyr`, `RColorBrewer`, `ggpubr`, `ggrepel`, and others.

---
### Age-DEA
### Step 1: Build Pseudobulk Data
**Script:** `Age-DEA/pseudobulk_dreamlet.R`

This script uses the `dreamlet::aggregateToPseudobulk()` function to aggregate single-cell data into pseudobulk samples, then processAssays() filters genes based on criteria such as minimum counts and minimum number of donors.

**Input:** SingleCellExperiment files for each condition. 
It saves two files:
- `output/pb.rds`: the resulting pseudobulk object.
- `output/dea_toptable.rds`: a filtered toptable used for downstream differential expression analysis.

Example usage inside the script:
```r
library(dreamlet)

# Aggregate to pseudobulk
system.time(pb <- aggregateToPseudoBulk(sce,
                                        assay = "counts",     
                                        cluster_id = "celltype", 
                                        sample_id = "assignment",
                                        fun = opt$aggr_fun,
                                        verbose = FALSE))
# Filter genes by counts and donor number
  system.time(res.proc <- processAssays(sceObj = ge_dge, 
                                        formula = form,
                                        min.cells = 5,
                                        min.count = 5,
                                        min.samples = 4,
                                        min.prop=opt$min_prop))
```
### Step 2: Linear mixed effects model
**Script:** `Age-DEA/pseudobulk_inrt_lmer.R`
This script runs gene-wise linear mixed models (LMMs) to test for associations between gene expression and a phenotype (e.g., age) using pseudobulked single-cell data.

**Input:**
- `*_pb.Rds`: Pseudobulk expression matrix (`SummarizedExperiment`), generated using `aggregateToPseudobulk()`.
- `*_dea_vp_topTable.rds`: Output from `dreamlet::dreamlet()`, containing the processed assay object with filtering info (via `processAssays()`).

**Workflow:**  
  1. Load pseudobulk expression and DE results.
  2. Extract raw counts, `voomWithDreamWeights`-transformed, and `log2CPM`-transformed matrices. Filter based on retained genes/samples from `processAssays()`. 
  3. Apply inverse-normal rank transformation per gene (optional).
  4. For each gene and normalization method (`voomWithDreamWeights`, `log2cpm`), fit an LMM using `lmerTest::lmer()`.

**Output**:
   - Results are structured in a nested list:  
     - `model`: tidy output of LMM fit.  
     - `nearZeroVar`: report on low-variance predictors.



# Impact of age on immune response dynamics upon pathogen stimulation at single-cell resolution
This work aims to uncover the age effect on the immune response when there is an exposition to external stimuli. We study both compositional and transcriptional changes in a context of pathogen challenge at different ages. Together, our results demonstrate the importance of incorporating age as a critical variable in immunological studies and in the design of more
effective, tailored therapeutic and vaccination strategies.

### Cell Type Composition Changes - Figures

This script generates figures showing how cell type composition varies across different stimulation conditions ans timepoints. 

**Input:** Preprocessed metadata and cell composition `.rds` files.

**Workflow:**  
  1. Loads metadata to order cell types for plotting.  
  2. Defines experimental parameters and filters valid condition combinations.  
  3. Reads cell composition data for each condition, combining into one dataset.  
  4. Uses `ggplot2` to create faceted plots highlighting significant composition changes (p < 0.05).  

---
### Age-DEA
### Step 1: Build Pseudobulk Data
**Script:** `Age-DEA/pseudobulk_dreamlet.R`

This script uses the `dreamlet::aggregateToPseudobulk()` function to aggregate single-cell data into pseudobulk samples, then processAssays() filters genes based on criteria such as minimum counts and minimum number of samples. 

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
**Script:** `Age_DEA/pseudobulk_inrt_lmer.R`
This script runs gene-wise linear mixed models (LMMs) to test for associations between gene expression and a phenotype (e.g., age) using pseudobulked single-cell data.

**Input:**
- `*_pb.Rds`: Pseudobulk expression matrix generated using `aggregateToPseudobulk()`.
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

### Response-DEA
**Script:** `Response_DEA/log_ratio_response.R`
This script calculates the log ratio between two conditions and fits a linear effects model. 

**Input:** 
- `*_pb.Rds`: Gene expression values in pseudobulk format. 
- Sample_metadata (data frame): Metadata including condition labels and covariates for each sample.
- Condition_levels (vector of two strings): The two conditions to compare (e.g., c("treated", "control")).

**Workflow:**  
1. Subset the pseudobulk matrix and sample_metadata to include only samples from the two specified conditions.
2. Calculate the log2 ratio of expression values between the two conditions for each gene.
3. Construct a linear model (e.g., lm(log_ratio ~ covariate1 + covariate2, data = ...)) using the log ratios as response variables.
4. Fit the model and extract relevant statistics (e.g., coefficients, p-values, adjusted p-values) for downstream interpretation.


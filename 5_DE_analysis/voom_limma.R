######################
## Voom limma - CPM ##
######################

# adapted from
# https://github.com/csoneson/conquer_comparison/blob/master/scripts/apply_voomlimma.R
suppressPackageStartupMessages(library(limma))
suppressPackageStartupMessages(library(edgeR))

run_voomlimma <- function(L) {
  message("voomlimma")
  timing = system.time({
    dge = DGEList(L$count, group = L$condt)
    dge = calcNormFactors(dge)
    design = model.matrix(~L$condt)
    vm = voom(dge, design = design, plot = TRUE)
    fit = lmFit(vm, design = design)
    fit = eBayes(fit)
    tt = topTable(fit, n = Inf, adjust.method = "BH")
  })

  tbl_results = as_tibble(tt, rownames = "geneID") %>%
                    mutate(log2fc = logFC) %>%
                    mutate(pval = tt$P.Value) %>%
                    mutate(qval = tt$adj.P.Val) %>%
                    select(geneID, log2fc, pval, qval, AveExpr, t, B) %>%
                    arrange(qval, log2fc)
  
  return(tbl_results)
}

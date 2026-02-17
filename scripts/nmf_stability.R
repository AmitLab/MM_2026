library(foreach)
library(Matrix)
library(ggrepel)
library(ggplot2)
library(ggpubr)
library(cowplot)
library(doParallel)
library(plyr)
library(MASS)
library(NMF)
library(circlize)
library(plotly)
library(reshape2)

setwd('/home/projects/amit/annaku/repos/Blueprint/scripts')
theme_set(theme_cowplot())

registerDoParallel(cores = 32)
getDoParWorkers()

zscore_path = '/home/projects/amit/annaku/repos/Blueprint/data/processed/zscore_outputs'
nmf_save_path = '/home/projects/amit/annaku/repos/Blueprint/data/processed/nmf_outputs'

ver <- '20250306'
level <- "samplelevel"

dat <- read.table(paste0(zscore_path, "/zstat_Atlas_v_", ver, "_full_", level, "_preproc.txt"), sep = '\t', header = T, row.names = 1, check.names = F)
print(dim(dat))

############################

setClass(
  "nmf_sample",
  representation(
    gpct = "numeric",
    cpct = "numeric",
    gene_clus = "list",
    cl_clus = "list"
  )
)

run_nmf_sample <- function(dat, k, gpct = 0.9, cpct = 0.9){
    row_new <- sample(rownames(dat), nrow(dat) * gpct) 
    col_new <- sample(colnames(dat), ncol(dat) * cpct) 
    dat_new <- dat[row_new, col_new]
    
    res_nmf <- nmf(t(dat_new), k, rng = 10)
    #res_nmf_map <- consensusmap(res_nmf)
    
    gene_clus <- melt(predict(res_nmf, what = "columns"))
    gene_clus <- data.frame(Gene = rownames(gene_clus), clus = gene_clus$value)
    gene_clus <- gene_clus[order(gene_clus$clus),]
    rownames(gene_clus) <- gene_clus$Gene
    
    cl_clus <- melt(predict(res_nmf, what = "rows"))
    cl_clus <- data.frame(PID = rownames(cl_clus), clus = cl_clus$value)
    cl_clus <- cl_clus[order(cl_clus$clus),]
    rownames(cl_clus) <- cl_clus$PID
    
    res <- new("nmf_sample",
               gpct = gpct,
               cpct = cpct,
               gene_clus = gene_clus,
               cl_clus = cl_clus
              )
    return(res)
}

########

cal_nmf_sim <- function(run1, run2){
    gene_both <- intersect(run1@gene_clus$Gene, run2@gene_clus$Gene)
    gene_cl_run1 <- run1@gene_clus$clus[match(gene_both, run1@gene_clus$Gene)]
    gene_cl_run2 <- run2@gene_clus$clus[match(gene_both, run2@gene_clus$Gene)]
    cf_gene <- table(gene_cl_run1, gene_cl_run2)
    sim_gene <- sum(apply(cf_gene, 1, max))/sum(cf_gene)
    
    cl_both <- intersect(run1@cl_clus$PID, run2@cl_clus$PID)
    cl_cl_run1 <- run1@cl_clus$clus[match(cl_both, run1@cl_clus$PID)]
    cl_cl_run2 <- run2@cl_clus$clus[match(cl_both, run2@cl_clus$PID)]
    cf_cl <- table(cl_cl_run1, cl_cl_run2)
    sim_cl <- sum(apply(cf_cl, 1, max))/sum(cf_cl)

    max_cols <- apply(cf_gene, 1, which.max)
    max_genes <- unlist(lapply(seq_along(max_cols), function(i) {
        gene_both[gene_cl_run1 == i & gene_cl_run2 == max_cols[i]]
    }))

    ###

    res <- structure(list(gene = sim_gene, cl = sim_cl), names = c('gene', 'cl'))
    return(list(similarity = res, intersected_genes = max_genes))
}

nruns <- 100
min_k <- 4
max_k <- 15

library(doParallel)
registerDoParallel(cores = 30)
getDoParWorkers()

gpct = 1 
cpct = 0.9 

set.seed(10)
sim <- foreach(i = seq(nruns), .combine = rbind) %dopar% {
    message(i)
    foreach(k = seq(min_k,max_k), .combine = rbind) %do% {
        run1 <- run_nmf_sample(dat, k, gpct = gpct, cpct = cpct)
        run2 <- run_nmf_sample(dat, k, gpct = gpct, cpct = cpct)
        result <- cal_nmf_sim(run1, run2)
        res <- result$similarity
        int_genes = paste(result$intersected_genes, collapse = ", ")
        res <- c(k, i, gpct, cpct, res, int_genes)
    }
}

colnames(sim) <- c('k', 'i', 'gpct', 'cpct', 'gene', 'cl', 'intersected_genes')
sim <- as.data.frame(sim)
head(sim)

sim_tosave <- as.data.frame(lapply(sim, unlist))
filename <- paste0(nmf_save_path, "/sim_", 'ver_', ver, '_', nruns, "_runs_", min_k, "_", max_k, "_withgenes.csv")
write.csv(sim_tosave, file = filename)
print(paste0('saved to ', filename))

# save plots

save_dir <- '/home/projects/amit/annaku/repos/Blueprint/figures/fig2/'

plot <- ggboxplot(sim_tosave, x = "k", y = "gene", ylim = c(0,1), add = "jitter")
filename <- paste0("genes_stab_by_", ver, '_', nruns, "_runs_", min_k, max_k,'.svg')
ggsave(paste0(save_dir,filename), plot = plot)

#

plot <- ggboxplot(sim_tosave, x = "k", y = "gene", ylim = c(0,1),)
filename <- paste0("genes_stab_by_", ver, '_', nruns, "_runs_", min_k, max_k,'_nodots.svg')
ggsave(paste0(save_dir,filename), plot = plot)

#

plot <- ggboxplot(sim_tosave, x = "k", y = "cl", ylim = c(0,1), add = "jitter")
filename <- paste0("samples_stab_by_", ver, '_', nruns, "_runs_", min_k, max_k,'.svg')
ggsave(paste0(save_dir,filename), plot = plot)

#

plot <- ggboxplot(sim_tosave, x = "k", y = "cl", ylim = c(0,1),)
filename <- paste0("samples_stab_by_", ver, '_', nruns, "_runs_", min_k, max_k,'_nodots.svg')
ggsave(paste0(save_dir,filename), plot = plot)
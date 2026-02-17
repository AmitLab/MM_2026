library(reticulate)
library(metacell)
library(Matrix)
library(foreach)
library(doParallel)
library(parallel)
library(plyr)
library(doMC)

VERSION <- '20250306'
data_path = '/home/projects/amit/annaku/repos/Blueprint/data/processed/zscore_outputs/'
repo_path = '/home/projects/amit/annaku/repos/Blueprint/'

# 1 downsample (high RAM needed)

cells <- readLines(paste0(data_path, 'cells_forz_v_', VERSION, '.txt')) 
genes <- readLines(paste0(data_path, 'genes_forz_v_', VERSION, '.txt'))
mat <- Matrix(readMM(paste0(data_path, 'matrix_forz_v_', VERSION, '.mtx')), sparse=TRUE)
colnames(mat) <- cells
rownames(mat) <- genes

message(paste("matrix dimensions:", paste(dim(mat), collapse=" x ")))

message(paste('1. starting downsampling'))
gc(full=TRUE)

target_umi <- 300
downsampled = scm_downsamp(mat, target_umi)
dim(downsampled)

missed_cells <- setdiff(colnames(mat), colnames(downsampled))
message(paste(length(missed_cells)))

if (length(missed_cells) > 0) {  
    downsampled_added <- scm_downsamp(mat[, missed_cells], target_umi)
    message(paste(dim(downsampled_added)))

    downsampled <- cbind(downsampled, downsampled_added)
    message(paste(dim(downsampled)))
}

message(paste(dim(downsampled)))

genes <-rownames(downsampled)
cells <- colnames(downsampled)
writeLines(genes, paste0(data_path, "genes_forz_down_", target_umi, "_v_", VERSION, ".mtx"))
writeLines(cells, paste0(data_path, "cells_forz_down_", target_umi, "_v_", VERSION,".mtx"))

writeMM(downsampled, paste0(data_path, "matrix_forz_down_", target_umi, "_v_", VERSION,".mtx"))

message(paste0('saved downsampled matrix to ', data_path, 'matrix_forz_down_', target_umi, '_v_', VERSION, '.mtx'))

# 2 run zscores
message(paste('2. running zscores'))

########## helper functions ##########
get_n_hip_cells_per_control <- function(cells_hip_all, tot, returnNames = T){
  n_control <- length(unique(cells_hip_all$pID))
  message(paste(n_control, "controls"))
  n <- ceiling(tot/n_control)
  message(paste("Sample", n, "cells per control"))
  res <- foreach(i = sort(unique(cells_hip_all$pID)), .combine = rbind) %do% 
    cells_hip_all[sample(rownames(cells_hip_all[cells_hip_all$pID == i,]), min(n, nrow(cells_hip_all[cells_hip_all$pID == i,]))),]
  res <- res[sample(rownames(res), tot),]
  if(returnNames){
    return(rownames(res))
  }else{
    return(res)
  }
}

cal_exp_diff <- function(grpA, grpB, mat_ds, ntime, blck = NULL){
  if(is.null(blck)){
    blck <- c("^MT-", "^MTMR", "^MTND", "^MTRN", "^MTCO", "^MTATP", 
              "^AC[0-9]", "^AL[0-9]", "^AP[0-9]", "^Gm[0-9]", "^MIR", "^TMSB", "^SNOR", "^ATP", "^DNAJ",
              "^FP[0-9]", "^FO[0-9]", "^RN7SL", "^COL[0-9]")
  }
  message(paste("grpA", length(grpA)))
  message(paste("grpB", length(grpB)))
  data <- mat_ds
  gene <- rownames(data)
  gene_bad <- foreach(i = blck, .combine = c) %do% grep(i, gene, v = T)
  gene <- setdiff(gene, gene_bad)
  message(paste("gene", length(gene)))
  A <- data[gene, grpA]
  B <- data[gene, grpB]
  AB <- cbind(A, B)
  pctA <- Matrix::rowSums(A > 0)/ncol(A)
  pctB <- Matrix::rowSums(B > 0)/ncol(B)
  avgA <- Matrix::rowMeans(A)
  avgB <- Matrix::rowMeans(B)
  wilcox.p <- apply(AB, 1, function(x) wilcox.test(x[1:ncol(A)], x[(ncol(A) + 1):ncol(AB)])[[3]])
  wilcox.p[is.na(wilcox.p)] <- 1
  zstat <- abs(qnorm(wilcox.p/2)) * sign(avgB - avgA)
  wilcox.padj <- p.adjust(wilcox.p, 'fdr')
  wilcox.padj[is.na(wilcox.padj)] <- 1
  zstatadj <- abs(qnorm(wilcox.padj/2)) * sign(avgB - avgA)
  ori <- as.data.frame(cbind(gene, pctA, pctB, avgA, avgB, zstat, zstatadj, wilcox.p, wilcox.padj))
  ori$ntime <- ntime
  ori <- ori[order(ori$wilcox.padj),]
  message(paste(dim(ori), collapse = '\t'))
  return(ori)
}

diff_calc <- function(idx, n, ntime, cells_hip_all, cells_all, mat_ds){
  res <- foreach(i = seq(ntime), .combine = rbind) %dopar% {
    grpA <- get_n_hip_cells_per_control(cells_hip_all, n)
    grpB <- rownames(cells_all)[cells_all$idx == idx]
    if(length(grpB) < 500){
     grpB <- rep(grpB, ceiling(500/length(grpB))) 
    }
    grpB <- sample(grpB, n)
    cal_exp_diff(grpA, grpB, mat_ds, i)
  }
  return(res)
}

diff_stat <- function(diff_res){
  genes_freq <- table(diff_res$gene)
  genes <- names(genes_freq)
  message(paste("final gene counts: ", length(genes)))
  res <- foreach(i = genes, .combine = rbind) %dopar% {
    now <- diff_res[diff_res$gene == i,]
    pctA.mean <- mean(as.numeric(as.vector(now$pctA)))
    pctB.mean <- mean(as.numeric(as.vector(now$pctB)))
    avgA.mean <- mean(as.numeric(as.vector(now$avgA)))
    avgB.mean <- mean(as.numeric(as.vector(now$avgB)))
    zstat.mean <- mean(as.numeric(as.vector(now$zstat)))
    zstatadj.mean <- mean(as.numeric(as.vector(now$zstatadj)))
    wilcox.p.mean <- -mean(log10(as.numeric(as.vector(now$wilcox.p))))
    wilcox.padj.mean <- -mean(log10(as.numeric(as.vector(now$wilcox.padj))))
    pctA.sd <- sd(as.numeric(as.vector(now$pctA)))
    pctB.sd <- sd(as.numeric(as.vector(now$pctB)))
    avgA.sd <- sd(as.numeric(as.vector(now$avgA)))
    avgB.sd <- sd(as.numeric(as.vector(now$avgB)))
    zstat.sd <- sd(as.numeric(as.vector(now$zstat)))
    zstatadj.sd <- sd(as.numeric(as.vector(now$zstatadj)))
    wilcox.p.sd <- sd(log10(as.numeric(as.vector(now$wilcox.p))))
    wilcox.padj.sd <- sd(log10(as.numeric(as.vector(now$wilcox.padj))))
    res <- data.frame(pctA.mean, pctB.mean, avgA.mean, avgB.mean, zstat.mean, zstatadj.mean, wilcox.p.mean, wilcox.padj.mean,
                      pctA.sd, pctB.sd, avgA.sd, avgB.sd, zstat.sd, zstatadj.sd, wilcox.p.sd, wilcox.padj.sd)
    rownames(res) <- i
    res
  }
  res <- res[order(res$wilcox.p.mean, decreasing = T),]
  return(res)
}
  run_analysis <- function(METHOD, mat_ds, cells_all_original) {
  message(paste("\n=== Starting analysis for method:", METHOD, "===\n"))
  
  # base_path <- '/home/projects/amit/annaku/repos/MM_2024_AK'
  # data_path <- file.path(base_path, 'data/zscore_outputs')
  
  cells_all <- cells_all_original[cells_all_original$Method == METHOD,]
  message(paste("Cells after METHOD filter:", nrow(cells_all)))
  
  # split healthy and non-healthy
  cells_hip_all <- cells_all[cells_all$Tissue == 'BM' & cells_all$Populations == 'Normal_PC',]
  # cells_all <- cells_all[cells_all$Tissue != 'Blood' & cells_all$Populations != 'Healthy',]
  cells_all <- cells_all[cells_all$Tissue != 'Blood' & cells_all$Populations == 'Malignant',]
  
  # filter healthy cells
  healthy_counts <- table(cells_hip_all$idx)
  pool_healthy <- sort(names(healthy_counts[healthy_counts > 50]))
  cells_hip_all <- cells_hip_all[cells_hip_all$idx %in% pool_healthy,]
  message(paste("Healthy cells:", nrow(cells_hip_all)))
  message(paste("Non-healthy cells:", nrow(cells_all)))
  
  # get pool
  clone_counts <- table(cells_all$idx)
  message(paste("Unique samples before threshold:", length(unique(cells_all$idx))))
  message(paste("Samples with >50 cells:", sum(clone_counts > 50)))
  pool <- sort(names(clone_counts[clone_counts > 50]))
  pool <- grep('^hip', pool, v = T, invert = T)
  message(paste("Final pool size:", length(pool)))
  
  saveRDS(pool, file = file.path(data_path, paste0('idx_pool_', METHOD, '_v_', VERSION, '_samplelevel.Rds')))
  write.table(pool, file = file.path(data_path, paste0('idx_pool_', METHOD, '_v_', VERSION, '_samplelevel.txt')))
  
  # run de
  n <- 100 # number of cells
  ntime <- 100 # number of sampling
  
  for(idx in pool) {
    output_file <- file.path(data_path, paste0('all_v_', VERSION, '_diff_', idx, '_full.txt'))
    
    if (file.exists(output_file)) {
      message(paste('File already exists for idx:', idx, '. Skipping.'))
      next
    }
    
    message(paste('Processing idx:', idx))
    ind_diff <- diff_calc(idx, n, ntime, cells_hip_all, cells_all, mat_ds)
    ind_diff_stat <- diff_stat(ind_diff)
    write.table(ind_diff_stat, file = output_file, sep = '\t', quote = F, col.names = NA)
  }
  
  message(paste("\n=== Completed analysis for method:", METHOD, "===\n"))
}

########## execution ##########

cleanup_memory <- function() {
    for(i in 1:3) {
        gc(full=TRUE, verbose=TRUE)
    }
}

main <- function() {
    setwd(repo_path)
    cleanup_memory()
    
    total_cores <- detectCores()
    cores_to_use <- min(total_cores, 4)  
    message(paste("Using cores:", cores_to_use))
    registerDoParallel(cores = cores_to_use)
    
    message("Reading input data...")
    cells_all <- read.csv(paste0(data_path, 'cells_all_v_', VERSION, '.csv'),
                         stringsAsFactors = FALSE, check.names = FALSE, row.names = 1)
    message(paste("Total cells:", nrow(cells_all)))
    
    message("Loading matrix...")
    mat_ds <- readMM(paste0(data_path,"matrix_forz_down_300_v_", VERSION, ".mtx"))
    gene_names <- readLines(paste0(data_path,'genes_forz_down_300_v_', VERSION, '.mtx'))
    cell_names <- readLines(paste0(data_path,'cells_forz_down_300_v_', VERSION, '.mtx'))
    colnames(mat_ds) <- cell_names
    rownames(mat_ds) <- gene_names
    rm(gene_names, cell_names)  
    cleanup_memory()
    
    cells <- intersect(colnames(mat_ds), rownames(cells_all))
    message(paste("Cells after intersection:", length(cells)))
    cells_all <- cells_all[cells,]
    rm(cells) 
    cleanup_memory()
    
    cells_all$pID <- cells_all$'Code'
    cells_all$idx <- paste(cells_all$Method, cells_all$Populations, cells_all$pID, sep = '_')
    
    # each method
    for(method in c('MARS', 'SPID')) {
        message(sprintf("\nProcessing method: %s", method))
        run_analysis(method, mat_ds, cells_all)
        cleanup_memory()
    }
    
    rm(mat_ds, cells_all)
    cleanup_memory()
    message("\nall analyses completed")
}

main()

##############

# 3. merging results

pool_MARS <- readRDS(paste0(data_path,'idx_pool_MARS_v_', VERSION, '_samplelevel.Rds')) 
pool_SPID <- readRDS(paste0(data_path,'idx_pool_SPID_v_', VERSION, '_samplelevel.Rds'))
length(pool_MARS)
length(pool_SPID)
pool <- c(pool_MARS, pool_SPID)

g <- list()
## for case if some files are missed (possible because of new version)
missed_files <- character(0)
for(idx in pool){
  fpath <- paste0(data_path,'all_v_', VERSION, '_diff_', idx, '_full.txt') # delete selected!
  
  tryCatch({
    now <- read.table(fpath, header = T)
    res <- data.frame(rownames(now), now$zstat.mean)
    colnames(res) <- c('gene', paste0('z.', idx))
    g[[idx]] <- res
    message(paste("Processed:", fpath))
  }, error = function(e) {
    missed_files <- c(missed_files, fpath)
    message(paste("Error processing:", fpath))
  })
}

if (length(missed_files) > 0) {
  cat("The following files were missed:\n")
  cat(paste(missed_files, collapse = "\n"))
  cat("\n")
}
g_merged_mean <- join_all(g, type = 'full')
g_merged_mean[is.na(g_merged_mean)] <- 0
dim(g_merged_mean)
dat_mean <- g_merged_mean[,2:ncol(g_merged_mean)]
rownames(dat_mean) <- g_merged_mean$gene
dat_mean <- as.matrix(dat_mean)
dat_mean <- dat_mean[apply(dat_mean, 1, function(x) max(abs(x)) > 1),]
message(paste(dim(dat_mean)))
write.table(dat_mean, file = paste0(data_path,'zstat_Atlas_v_', VERSION, '_full_samplelevel.txt'), sep = '\t', quote = F, col.names = NA)
message(paste(paste0('saved to ', data_path, 'zstat_Atlas_v_', VERSION, '_full_samplelevel.txt')))

##############

# 4. preproc zscore table

dat_mean <- read.table(paste0(data_path,'zstat_Atlas_v_', VERSION, '_full_samplelevel.txt'), sep = '\t', header = T, row.names = 1, check.names = F)

dat_mean[is.na(dat_mean)] <- 0

# new addition
ribo_genes <- grep("^RPL|^RPS", rownames(dat_mean), value = TRUE)
cat("Excluding ribosomal genes:", length(ribo_genes), "\n")

dat_mean <- dat_mean[!rownames(dat_mean) %in% ribo_genes, ]
cat("Dimensions after ribo exclusion:", dim(dat_mean), "\n")

fc <- 2 # cutoff for the zscore for genes selecting
dat <- dat_mean
res_mean_pre <- dat[apply(dat, 1, function(x) (sort(x)[5] < -fc | sort(x, decreasing = T)[5] > fc)),]
print(dim(res_mean_pre))

# exclude genes apriori

blk <- c('BX679664.3', 'GPC5-AS1', "KIAA0040", "KIAA1549L", 'RNU6-333P', "Z74021.1")
res_mean_pre <- res_mean_pre[setdiff(rownames(res_mean_pre), blk),]

dat <- res_mean_pre
print(dim(dat))

dat[dat < 0] <- 0 #all negative values to zero

dat <- dat[apply(dat, 1, function(x) sd(x) > 0.5),] # only keep genes with sd > 0.5
print(dim(dat))

# clipping fold change
fc_cut <- 2
dat[dat > fc_cut] <- fc_cut
dat[dat < -fc_cut] <- -fc_cut
print(dim(dat))

# excl hip if any
dat_mean <- dat_mean[, !grepl("hip", colnames(dat_mean), ignore.case = TRUE)]

write.table(dat, paste0(data_path, '/zstat_Atlas_v_', VERSION, '_full_samplelevel_preproc.txt'),
sep = '\t', quote = F, col.names = NA)

print(paste0('saved to ', data_path, 'zstat_Atlas_v_', VERSION, '_full_samplelevel_preproc.txt'))


library(infercnv)
#MARS
infercnv_obj = CreateInfercnvObject(raw_counts_matrix="/home/projects/amit/annaku/repos/Blueprint/data/processed/infercnv_r_input/exp_arch_prolif_MARS.tsv",
                                    annotations_file="/home/projects/amit/annaku/repos/Blueprint/data/processed/infercnv_r_input/ann_arch_prolif_MARS.txt",
                                    delim="\t",
                                    gene_order_file="/home/projects/amit/annaku/repos/Blueprint/data/processed/infercnv_r_input/ann_gene.txt",
                                    ref_group_names=c("Normal_PC"))

infercnv_obj = infercnv::run(infercnv_obj,
                             cutoff=0.1,  
                             out_dir="/home/projects/amit/annaku/repos/Blueprint/data/processed/infercnv_r_output/output_infercnv_r_arch_prolif_MARS",  
                             cluster_by_groups=TRUE,   
                             denoise=TRUE,
                             HMM=TRUE,
                             output_format="pdf"
                             )
#MARS v2
infercnv_obj = CreateInfercnvObject(raw_counts_matrix="/home/projects/amit/annaku/repos/Blueprint/data/processed/infercnv_r_input/exp_arch_prolif_SPID.tsv",
                                    annotations_file="/home/projects/amit/annaku/repos/Blueprint/data/processed/infercnv_r_input/ann_arch_prolif_SPID.txt",
                                    delim="\t",
                                    gene_order_file="/home/projects/amit/annaku/repos/Blueprint/data/processed/infercnv_r_input/ann_gene.txt",
                                    ref_group_names=c("Normal_PC"))

infercnv_obj = infercnv::run(infercnv_obj,
                             cutoff=0.1,  
                             out_dir="/home/projects/amit/annaku/repos/Blueprint/data/processed/infercnv_r_output/output_infercnv_r_arch_prolif_SPID",  
                             cluster_by_groups=TRUE,   
                             denoise=TRUE,
                             HMM=TRUE,
                             output_format="pdf"
                             )
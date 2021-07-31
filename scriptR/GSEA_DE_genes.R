require(tidyverse)
require(plyranges)
require(BSgenome.Hsapiens.UCSC.hg19)
require(clusterProfiler)
DE_DIVA <- PhDfunc::GetDE_DIvA()

DE_DIVA <- DE_DIVA %>% 
  mutate(Type = case_when(
    logFC < 0 & FILTER.FC == 1 & FILTER.P == 1 ~ "Downregulated",
    # logFC < 0 & FILTER.FC == 1 ~ "Downregulated",
    # logFC < -0.15 ~ "Downregulated",
    logFC > 0 & FILTER.FC == 1 & FILTER.P == 1~ "Upregulated",
    # logFC > 0 & FILTER.FC == 1~ "Upregulated",
    # logFC > 0.15 ~ "Upregulated",
    TRUE ~ "None"
  )) 



formula_res_2 <- compareCluster(rowname~Type, data=filter(DE_DIVA,Type !="None"), fun="enrichGO",
                                OrgDb         = org.Hs.eg.db::org.Hs.eg.db,
                                keyType       = 'ENSEMBL',
                                ont           = "BP",
                                pAdjustMethod = "BH",
                                pvalueCutoff  = 0.1,
                                qvalueCutoff  = 0.05,
                                readable=TRUE)




bp2 <- simplify(gofilter(formula_res_2, level=10), cutoff=0.7, by="p.adjust", select_fun=min)
# bp2 <- simplify(formula_res_2, cutoff=0.7, by="p.adjust", select_fun=min)
p2 <- bp2 %>% dotplot(showCategory=30)

ggsave(glue::glue("/home/rochevin/Documents/PROJET_THESE/PAPIER_COLINE_AUDE/results/EnrichGO_on_DE_genes.pdf"),p2,width = 12,height=10)
formula_res_2@compareClusterResult %>% as_tibble() %>% write_tsv(glue::glue("/home/rochevin/Documents/PROJET_THESE/PAPIER_COLINE_AUDE/results/EnrichGO_on_DE_genes.tsv"))

ens.genes <- read.table("/home/rochevin/Documents/PROJET_THESE/DSB_EFFECT_ON_GENE/data/EnsDB.Hsapiens.v75.genes.bed",sep="\t",h=T) %>% GRanges()
ens.genes <- regioneR::filterChromosomes(ens.genes,keep.chr=c(1:22,"X","Y"))
seqlevels(ens.genes) <- paste0("chr",seqlevels(ens.genes))
DE_DIVA.chr <- DE_DIVA %>% filter(rowname %in% filter(ens.genes,seqnames %in% glue::glue("chr{c(1,17,'X')}"))$gene_id)

formula_res_3 <- compareCluster(rowname~Type, data=filter(DE_DIVA.chr,Type !="None"), fun="enrichGO",
                                OrgDb         = org.Hs.eg.db::org.Hs.eg.db,
                                keyType       = 'ENSEMBL',
                                ont           = "BP",
                                pAdjustMethod = "BH",
                                pvalueCutoff  = 0.1,
                                qvalueCutoff  = 0.05,
                                readable=TRUE)




bp3 <- simplify(gofilter(formula_res_3, level=10), cutoff=0.7, by="p.adjust", select_fun=min)
# p3 <- bp3 %>% dotplot(showCategory=30)
p3 <- formula_res_3 %>% dotplot(showCategory=30)

ggsave(glue::glue("/home/rochevin/Documents/PROJET_THESE/PAPIER_COLINE_AUDE/results/EnrichGO_on_DE_genes_chr_1_17_X_only.pdf"),p3,width = 12,height=10)
formula_res_3@compareClusterResult %>% as_tibble() %>% write_tsv(glue::glue("/home/rochevin/Documents/PROJET_THESE/PAPIER_COLINE_AUDE/results/EnrichGO_on_DE_genes_chr_1_17_X_only.tsv"))

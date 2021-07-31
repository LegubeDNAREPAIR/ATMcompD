require(HiTC)
require(tidyverse)
require(plyranges)
require(reshape2)
require(keras)
require(rtracklayer)
require(cowplot)
seqlens <- BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19 %>% seqlengths()
#Functions

loadRData <- function(fileName){
    #loads an RData file, and returns it
    load(fileName)
    get(ls()[ls() != "fileName"])
}



get_index_of_bed <- function(HTC,bed){
    rows.chr <- x_intervals(HTC)
    as(findOverlaps(bed,rows.chr), "List")
}


get_val_from_gene_inside_TADs <- function(file,DSB,genes){
    HTC <- loadRData(file)
    val1 <- HTC %>% intdata %>% sum
    
    res <- lapply(1:length(DSB),function(i){
        one_dsb <- DSB[i] %>% anchor_center() %>% mutate(width=2000000)
        mes_genes <- genes %>% filter_by_overlaps(one_dsb)
        one_dsb <- one_dsb %>% anchor_center() %>% mutate(width=30000)
        mes_genes <- mes_genes %>% filter_by_non_overlaps(one_dsb)
        if(length(mes_genes) > 0){
            
            pos.dsb <- get_index_of_bed(HTC,one_dsb) %>% unlist()
            pos.genes <- get_index_of_bed(HTC,mes_genes)
            
            res <- lapply(1:length(pos.genes),function(i){
                x <- pos.genes[[i]]
                
                intdata(HTC)[x,pos.dsb] %>% sum()
            }) %>% unlist()
            
            mes_genes %>% as_tibble() %>% 
                dplyr::select(seqnames,gene_id,width) %>%
                mutate(score = res) %>%
                mutate(DSB = one_dsb$name)
        }
            
        
            
    }) %>% compact() %>% bind_rows()
    return(
        list(res,val1)
    )
}


get_val_from_gene_inside_TADs_2 <- function(file,DSB,genes){
    HTC <- loadRData(file)
    val1 <- HTC %>% intdata %>% sum
    
    res <- lapply(1:length(DSB),function(i){
        one_dsb <- DSB[i] %>% anchor_center() %>% mutate(width=2000000)
        mes_genes <- genes %>% filter_by_overlaps(one_dsb)
        if(length(mes_genes) > 0){
            
            pos.genes <- get_index_of_bed(HTC,mes_genes)
            
            res <- lapply(1:length(pos.genes),function(i){
                sres <- lapply(1:length(pos.genes),function(j){
                    x <- pos.genes[[i]]
                    y <- pos.genes[[j]]
                    intdata(HTC)[x,y] %>% sum()
                }) %>% compact() %>% unlist()
                mes_genes %>% as_tibble() %>% 
                    dplyr::select(seqnames,width) %>%
                    mutate(score = sres) %>%
                    mutate(Gene1 = mes_genes$gene_id[i]) %>% 
                    mutate(Gene2 = mes_genes$gene_id)
            }) %>% bind_rows()
            return(res)
            
        }
        
        
        
    }) %>% compact() %>% bind_rows()
    return(
        list(res,val1)
    )
}


ens.genes <- read.table("/home/rochevin/Documents/PROJET_THESE/DSB_EFFECT_ON_GENE/data/EnsDB.Hsapiens.v75.genes.bed",sep="\t",h=T) %>% GRanges()
ens.genes <- regioneR::filterChromosomes(ens.genes,keep.chr=c(1:22,"X","Y"))
seqlevels(ens.genes) <- paste0("chr",seqlevels(ens.genes))
ens.genes.ext <- ens.genes %>% anchor_5p() %>% stretch(3000)
genes_RNA_seq <- PhDfunc::GetDE_DIvA()
expressed_genes <- genes_RNA_seq %>% filter(FILTER.FC == 0 & FILTER.P == 0) %>% pull(rowname) 
expressed_genes <- ens.genes.ext[ens.genes.ext$gene_id %in% expressed_genes]
DSB174 <- read_bed("/mnt/NAS/DATA/AsiSI/ASIsites_hg19_174_clived_IQ1.5.bed")
bless80 = read_bed("/mnt/NAS/DATA/AsiSI/BLESS_80best_JunFragPE_Rmdups_pm500bp.bed")

m.w <- "5kb"

HiC.files <- list.files("/home/rochevin/Documents/PROJET_INGE/HiC_Coline/data/HTC",pattern=m.w,full.names=T)

window <- 5000
#Compute interaction between DSB AND genes
res.DSB <- lapply(seqlevels(bless80),function(chrom){
    message(chrom)
    f <- HiC.files[str_detect(HiC.files,str_c("_",chrom,"_"))]
    Replicate <- f %>% map(.%>% basename() %>% str_extract("manip[A-Z]+|rep[0-9]+"))
    diva.chr <- bless80 %>% filter(seqnames ==chrom) %>% 
        anchor_center() %>% mutate(width = window) 
    
    res <- lapply(unique(Replicate),function(j){
        s <- f[grep(j,f)]
        if(length(s)<2){
            NULL
        }else{
            message(s)
            Cond <- s %>% map(.%>% basename() %>% str_extract("DIvA|OHT"))
            
            
            
            c(res1,val1) %<-% get_val_from_gene_inside_TADs(file = s[[1]],DSB = diva.chr,genes = expressed_genes)
            c(res2,val2) %<-% get_val_from_gene_inside_TADs(file = s[[2]],DSB = diva.chr,genes = expressed_genes)
            
            res1 <- res1 %>% mutate(score = (score+1) * (val1/val2)) %>% mutate(Condition = Cond[[1]])
            
            res2 <- res2 %>% mutate(score = (score+1) * (val1/val2)) %>% mutate(Condition = Cond[[2]])
            
            
            full <- list(res1,res2)
            full %>% bind_rows()
        }
        
        
    })
    names(res) <- unique(Replicate)
    res <- res %>% bind_rows(.id = "Replicate")
    
    res
    
}) %>% bind_rows()

data.plot <- res.DSB %>% spread(key = Replicate,value = score) %>%
    rowwise() %>% 
    mutate(meanValue = mean(manipA,manipB,na.rm=TRUE)/width) %>% ungroup()

FCtresh <- -0.20
down_genes <- genes_RNA_seq %>% filter(logFC < FCtresh) %>% pull(rowname)

data.plot <- data.plot %>% mutate(Type = ifelse(gene_id %in% down_genes,"Downregulated","No-Down")) 



n_cat <- data.plot %>% dplyr::count(Type) %>% pull(n)
stats <- str_c(
    c("Downregulated","No-Down"),
    " (n=",
    n_cat,
    ")\n(",
    c("<",">="),
    FCtresh,
    ") "
)

p1 <- data.plot %>% mutate(Type = ifelse(Type == "Downregulated",stats[1],stats[2])) %>% ggplot(aes(x=Type,y = meanValue,fill=Condition)) + geom_boxplot() + theme_classic(base_size = 18) + scale_y_log10() +
    theme(legend.position = "none") +
    scale_fill_manual(values = c("#2980b9","#e74c3c"))


n_cat <- data.plot %>% filter(Condition != "OHT") %>% dplyr::count(Type) %>% pull(n)
stats <- str_c(
    c("Downregulated","No-Down"),
    " (n=",
    n_cat,
    ")\n(",
    c("<",">="),
    FCtresh,
    ") "
)


p1.5 <- data.plot %>% filter(Condition != "OHT") %>% mutate(Type = ifelse(Type == "Downregulated",stats[1],stats[2])) %>% ggplot(aes(x=Type,y = meanValue,fill=Type)) + geom_boxplot() + theme_classic(base_size = 18) + scale_y_log10() +
    theme(legend.position = "none") 


p1.55 <- data.plot  %>% filter(Condition == "OHT")%>% mutate(Type = ifelse(Type == "Downregulated",stats[1],stats[2])) %>% ggplot(aes(x=Type,y = meanValue,fill=Type)) + geom_boxplot() + theme_classic(base_size = 18) + scale_y_log10() +
    theme(legend.position = "none") + facet_wrap(~Condition)

data.plot %>% mutate(Type = ifelse(Type == "Downregulated",stats[1],stats[2])) %>% group_by(Condition) %>% nest() %>% 
    mutate(pval = map_dbl(data,function(x){wilcox.test(meanValue ~ Type,data=x)$p.value})) %>% dplyr:::select(-data)
    

pdf("/media/HDD_ROCHER/PROJET_THESE/PAPIER_COLINE_AUDE/results/Interaction_between_genes_and_DSB_inside_TADs.pdf",height=8,width=6)
print(p1)
print(p1.5)
print(p1.55)
dev.off()
#Compute interaction between genes (two categories, down and no down)
#Compute interaction between DSB AND genes

res.genes <- lapply(seqlevels(bless80),function(chrom){
    message(chrom)
    f <- HiC.files[str_detect(HiC.files,str_c("_",chrom,"_"))]
    Replicate <- f %>% map(.%>% basename() %>% str_extract("manip[A-Z]+|rep[0-9]+"))
    diva.chr <- bless80 %>% filter(seqnames ==chrom) %>% 
        anchor_center() %>% mutate(width = window) 
    
    res <- lapply(unique(Replicate),function(j){
        s <- f[grep(j,f)]
        if(length(s)<2){
            NULL
        }else{
            message(s)
            Cond <- s %>% map(.%>% basename() %>% str_extract("DIvA|OHT"))
            
            
            
            c(res1,val1) %<-% get_val_from_gene_inside_TADs_2(file = s[[1]],DSB = diva.chr,genes = expressed_genes)
            c(res2,val2) %<-% get_val_from_gene_inside_TADs_2(file = s[[2]],DSB = diva.chr,genes = expressed_genes)
            
            res1 <- res1 %>% mutate(score = (score+1) * (val1/val2)) %>% mutate(Condition = Cond[[1]])
            
            res2 <- res2 %>% mutate(score = (score+1) * (val1/val2)) %>% mutate(Condition = Cond[[2]])
            
            
            full <- list(res1,res2)
            full %>% bind_rows()
        }
        
        
    })
    names(res) <- unique(Replicate)
    res <- res %>% bind_rows(.id = "Replicate")
    
    res
    
}) %>% bind_rows()


data.plot.genes <- res.genes %>% 
    filter(Gene1 != Gene2) %>%
    distinct() %>% 
    spread(key = Replicate,value = score) %>%
    rowwise() %>% 
    mutate(meanValue = mean(manipA,manipB,na.rm=TRUE)/width) %>% ungroup()

data.plot.genes <- list(
    "Downregulated"=data.plot.genes %>% filter(Gene1 %in% down_genes) %>% filter(Gene2 %in% down_genes),
    "No-Down"=data.plot.genes %>% filter(!Gene1 %in% down_genes) %>% filter(!Gene2 %in% down_genes)
) %>% bind_rows(.id = "Type") 


p2 <- data.plot.genes %>% ggplot(aes(x=Type,y = meanValue,fill=Condition)) + geom_boxplot() +
    theme_classic(base_size = 18) + scale_y_log10() +
    scale_fill_manual(values = c("#2980b9","#e74c3c")) +
    theme(legend.position = "top")
pdf("/mnt/NAS/DOSSIERS_PERSO/VINCENT/HiC_Coline/Interaction_between_genes_inside_TADs.pdf",height=8,width=6)
print(p2)
dev.off()

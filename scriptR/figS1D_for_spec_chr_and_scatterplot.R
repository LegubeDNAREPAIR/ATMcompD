#From bopxplot /mnt/NAS/PAPIERS/EN COURS/PAPIER COLINE/NATURE/Revision/DATA/Figure 4/FIG4D/boxplot_damaged_notdamaged_verysmall_80best.pdf
require(tidyverse)
require(plyranges)
# bless80 = read_bed("/mnt/NAS/DATA/AsiSI/BLESS_80best_JunFragPE_Rmdups_pm500bp.bed")
bless80 = read_bed("/mnt/NAS/DATA/AsiSI/ASIsites_hg19_174_clived_IQ1.5.bed")
files <- c(
    str_c(list.files("/home/rochevin/Documents/PROJET_INGE/APA/Gaelle/manipA/174clived/loops_verysmall_DIvA_OHT_notdamaged_manipA_DIvA/10000",pattern="[A-Z0-9]v[A-Z0-9]",full.names = T),"/enhancement.txt") %>% setNames(str_c("DIvA",str_extract(.,"/[A-Z0-9]+v[A-Z0-9]+"))),
    str_c(list.files("/home/rochevin/Documents/PROJET_INGE/APA/Gaelle/manipA/174clived/loops_verysmall_DIvA_OHT_notdamaged_manipA_OHT/10000",pattern="[A-Z0-9]v[A-Z0-9]",full.names = T),"/enhancement.txt") %>% setNames(str_c("OHT",str_extract(.,"/[A-Z0-9]+v[A-Z0-9]+")))
) %>%
    map(read_delim," ",col_names = F) %>%
    map(gather) %>%
    bind_rows(.id = "Type") %>%
    separate(Type,into=c("Condition","seqnames")) %>% 
    mutate(seqnames = str_c("chr",str_extract(seqnames,"[A-Z0-9]+"))) 
DSB_n_chr <- str_c("chr",c(1:22,"X")) %>% enframe(value = "seqnames") %>% dplyr::select(-name) %>% left_join(bless80 %>% as_tibble() %>% dplyr::count(seqnames)) %>% mutate(n = replace_na(n,0))

toplot <- files %>% left_join(DSB_n_chr)  %>% na.omit() %>% group_by(Condition,seqnames) %>% summarise(medval = median(value),n=mean(n)) %>% spread(key= Condition,value=medval) %>%
    mutate(ratio = log2(OHT/DIvA))
p <- toplot %>% gather(key = Condition,value = value,-seqnames,-n) %>% ggplot(aes(x=n,y=value,col=Condition)) + geom_point() + theme_classic() +facet_wrap(~Condition,scales="free_y",ncol=1)
pdf("scatterplot_loop_str_by_chr.pdf",width=12,height=12)
print(p)
dev.off()

pp <- files %>% left_join(DSB_n_chr)  %>% na.omit()%>% spread(key= Condition,value=value) %>%
    mutate(ratio = log2(OHT/DIvA)) %>% mutate(n = as.factor(n)) %>% ggplot(aes(x=n,y=ratio)) + geom_boxplot() + theme_classic()
pdf("boxplot_loop_str_by_DSBn_in_chr_DSB174.pdf",width=12,height=6)
print(pp)
dev.off()

p2 <- files %>% filter(seqnames %in% c("chr16","chr1")) %>% spread(key= Condition,value=value) %>%
    mutate(ratio = log2(OHT/DIvA)) %>% ggplot(aes(x=seqnames,y=ratio)) +geom_boxplot() + geom_hline(yintercept = 0,col="red",linetype="dashed")
pdf("boxplot_notdamaged_174best_chr1_chr16_manipA.pdf",width=8,height=8)
print(p2)
dev.off()

#Refaire plot 4D avec manip A et B (et plus tard si CTRL)
files <- c(
    str_c(list.files("/home/rochevin/Documents/PROJET_INGE/APA/Gaelle/manipA/80best/loops_verysmall_DIvA_OHT_notdamaged_manipA_DIvA/10000",pattern="[A-Z0-9]v[A-Z0-9]",full.names = T),"/enhancement.txt") %>% setNames(str_c("manipA/notdamaged/","DIvA",str_extract(.,"/[A-Z0-9]+v[A-Z0-9]+"))),
    str_c(list.files("/home/rochevin/Documents/PROJET_INGE/APA/Gaelle/manipA/80best/loops_verysmall_DIvA_OHT_notdamaged_manipA_OHT/10000",pattern="[A-Z0-9]v[A-Z0-9]",full.names = T),"/enhancement.txt") %>% setNames(str_c("manipA/notdamaged/","OHT",str_extract(.,"/[A-Z0-9]+v[A-Z0-9]+"))),
    str_c(list.files("/home/rochevin/Documents/PROJET_INGE/APA/Gaelle/manipA/80best/loops_verysmall_DIvA_OHT_damaged_manipA_DIvA/10000",pattern="[A-Z0-9]v[A-Z0-9]",full.names = T),"/enhancement.txt") %>% setNames(str_c("manipA/damaged/","DIvA",str_extract(.,"/[A-Z0-9]+v[A-Z0-9]+"))),
    str_c(list.files("/home/rochevin/Documents/PROJET_INGE/APA/Gaelle/manipA/80best/loops_verysmall_DIvA_OHT_damaged_manipA_OHT/10000",pattern="[A-Z0-9]v[A-Z0-9]",full.names = T),"/enhancement.txt") %>% setNames(str_c("manipA/damaged/","OHT",str_extract(.,"/[A-Z0-9]+v[A-Z0-9]+"))),
    str_c(list.files("/home/rochevin/Documents/PROJET_INGE/APA/Gaelle/manipB/80best/loops_verysmall_DIvA_OHT_notdamaged_manipB_DIvA/10000",pattern="[A-Z0-9]v[A-Z0-9]",full.names = T),"/enhancement.txt") %>% setNames(str_c("manipB/notdamaged/","DIvA",str_extract(.,"/[A-Z0-9]+v[A-Z0-9]+"))),
    str_c(list.files("/home/rochevin/Documents/PROJET_INGE/APA/Gaelle/manipB/80best/loops_verysmall_DIvA_OHT_notdamaged_manipB_OHT/10000",pattern="[A-Z0-9]v[A-Z0-9]",full.names = T),"/enhancement.txt") %>% setNames(str_c("manipB/notdamaged/","OHT",str_extract(.,"/[A-Z0-9]+v[A-Z0-9]+"))),
    str_c(list.files("/home/rochevin/Documents/PROJET_INGE/APA/Gaelle/manipB/80best/loops_verysmall_DIvA_OHT_damaged_manipB_DIvA/10000",pattern="[A-Z0-9]v[A-Z0-9]",full.names = T),"/enhancement.txt") %>% setNames(str_c("manipB/damaged/","DIvA",str_extract(.,"/[A-Z0-9]+v[A-Z0-9]+"))),
    str_c(list.files("/home/rochevin/Documents/PROJET_INGE/APA/Gaelle/manipB/80best/loops_verysmall_DIvA_OHT_damaged_manipB_OHT/10000",pattern="[A-Z0-9]v[A-Z0-9]",full.names = T),"/enhancement.txt") %>% setNames(str_c("manipB/damaged/","OHT",str_extract(.,"/[A-Z0-9]+v[A-Z0-9]+")))
) %>%
    map(read_delim," ",col_names = F) %>%
    map(gather) %>%
    bind_rows(.id = "Type") %>%
    separate(Type,into=c("Manip","Type","Condition","seqnames")) %>% 
    mutate(seqnames = str_c("chr",str_extract(seqnames,"[A-Z0-9]+"))) 
files_ratio <- files %>% spread(key= Condition,value=value) %>%
    mutate(ratio = log2(OHT/DIvA)) 
p1 <- files_ratio %>% filter(Manip =="manipA")  %>% mutate(Type = fct_relevel(Type,"notdamaged","damaged")) %>% ggplot(aes(x=Type,y=ratio)) +geom_boxplot() + geom_hline(yintercept = 0,col="red",linetype="dashed")+ facet_wrap(~Manip,ncol=1)
p2 <- files_ratio %>% filter(Manip =="manipB")  %>% mutate(Type = fct_relevel(Type,"notdamaged","damaged"))%>% ggplot(aes(x=Type,y=ratio)) +geom_boxplot() + geom_hline(yintercept = 0,col="red",linetype="dashed")+ facet_wrap(~Manip,ncol=1)

pdf("boxplot_damaged_notdamaged_80best_manipB_manipA_verysmall.pdf",width=8,height=8)
print(p1)
print(p2)
dev.off()
#compute stat
## pvsmOHT
pval1 <- files_ratio %>% group_by(Manip) %>% nest() %>%  mutate(pval = map_dbl(data,. %>% wilcox.test(ratio~Type,data=.) %>%  .$p.value)) %>% dplyr::select(-data)
pval2 <- files_ratio %>% group_by(Manip,Type) %>% nest()
pval2$pval <- map_dbl(pval2$data,function(data){
    wilcox.test(data$ratio)$p.value
})
pval2 <- pval2 %>% dplyr::select(-data)

pval1
pval2

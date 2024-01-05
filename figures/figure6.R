## CNV burden across cell types and LOY for all modalities
#############################
# figure proportion loy by age group
# note that only a minority of kpmp male samples have an age in the meta data (n=11)
library(here)
library(openxlsx)
library(dplyr)
library(ggplot2)

clinical_meta <- openxlsx::read.xlsx(here("ckd","clinical_meta.xlsx"))

kpmp <- readRDS(here("ckd","rna_aggr_prep","kpmp","kpmp_prep.rds"))
rna_meta <- kpmp@meta.data %>%
  dplyr::filter(rkd_genotype %in% c("XY","LOY"), nCount_RNA > 1000) %>%
  dplyr::select(barcode, rkd_genotype, age, library_id) %>%
  dplyr::rename(genotype = rkd_genotype, age_group = age)
rna_meta$Modality <- "RNA"

atac_meta <- readRDS(here("ckd","atac_aggr_prep","step6_meta.rds"))
atac_meta <- atac_meta %>%
  dplyr::filter(akd_genotype %in% c("XY","LOY")) %>%
  dplyr::select(barcode, akd_genotype, library_id) %>%
  dplyr::rename(genotype = akd_genotype) %>% 
  left_join(clinical_meta, by = "library_id")
atac_meta$Modality <- "ATAC"

multi_meta <- readRDS(here("ckd","multi_aggr_prep","step3_meta_multi_loy.rds"))
multi_meta <- multi_meta %>%
  dplyr::filter(gmm_genotype %in% c("XY","LOY")) %>%
  dplyr::select(barcode, gmm_genotype, library_id) %>%
  dplyr::rename(genotype = gmm_genotype) %>% 
  left_join(clinical_meta, by = "library_id")
multi_meta$Modality <- "MULTI"

all_meta <- bind_rows(list(rna_meta, atac_meta, multi_meta))

toplot <- all_meta %>%
  dplyr::select(barcode, library_id, age_group, genotype, Modality) %>%
  na.omit() %>%
  group_by(library_id) %>%
  summarize(genotype, age_group, Modality) %>%
  table() %>%
  as.data.frame() %>%
  group_by(library_id) %>%
  dplyr::mutate(total_cells = sum(Freq)) %>%
  dplyr::mutate(prop_loy = Freq / total_cells) %>%
  dplyr::filter(genotype == "LOY", Freq > 0) 

levels(toplot$Modality) <- c("RNA","ATAC","MULTI")
library(grid)
age <- textGrob("Age(y)", gp=gpar(fontsize=10, fontface="bold"))

plotA <- toplot %>%
  ggplot(aes(age_group, prop_loy, fill=age_group)) +
  geom_boxplot() +
  geom_point(position=position_dodge(width=0.75), aes(group=age_group, shape=Modality)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  ylab("Proportion LOY") +
  xlab("") +
  ggtitle("A)") + 
  theme(plot.title = element_text(face="bold"),
        axis.text = element_text(color="black"),
        legend.title=element_blank(),
        legend.pos = c(0.3,0.7),
        legend.direction="vertical",
        panel.grid = element_blank()) +
  guides(fill = "none") +
  annotation_custom(age,xmin=-0.1,xmax=-0.1,ymin=-0.05,ymax=-0.05) +
  coord_cartesian(ylim=c(0,0.35), clip = "off") 

######################
# figure multiome cnv burden by ATAC and RNA
library(here)
library(ggplot2)
library(here)
library(data.table)
library(dplyr)
library(tidyr)
library(openxlsx)
library(parameters)

# create a plot dir
plots <- here("ckd","figures")
dir.create(here(plots))

# read in clinical meta data (ie age)
clinical_meta <- read.xlsx(here("ckd","clinical_meta.xlsx")) %>%
  dplyr::select(library_id, age, age_group)
clinical_meta$age_group <- factor(clinical_meta$age_group, levels = c("0-9","10-19","20-29","30-39","40-49","50-59","60-69","70-79","80-89","90-99"))

# read in loy anno for multiomes
multi_meta <- readRDS(here("ckd","multi_aggr_prep","step3_meta_multi_loy.rds"))

# read multiome atac cnv burden
multi_atac_cnv <- fread(here("ckd","multi_aggr_prep","step5_atac_cnv_burden.tsv"))

# format multiome meta data
multi_meta <- multi_meta %>%
  dplyr::filter(gmm_genotype %in% c("XY","LOY")) %>%
  dplyr::select(barcode, gmm_genotype, celltype, library_id) %>%
  dplyr::rename(genotype = gmm_genotype) %>% 
  left_join(clinical_meta, by = "library_id")

# join multiome meta data with cnv 
multi_cnv <- multi_atac_cnv %>% left_join(multi_meta, by = "barcode") %>% na.omit()
multi_cnv$Modality <- "MULTI"

# read in snATAC meta data
atac_meta <- readRDS(here("ckd","atac_aggr_prep","step6_meta.rds"))
atac_meta <- atac_meta %>%
  dplyr::filter(akd_genotype %in% c("XY","LOY")) %>%
  dplyr::select(barcode, akd_genotype, celltype, library_id) %>%
  dplyr::rename(genotype = akd_genotype) %>% 
  left_join(clinical_meta, by = "library_id")

# read snatac cnv burden and join to meta data
atac_cnv <- fread(here("ckd","atac_aggr_prep","step7_atac_cnv_burden.tsv"))
atac_cnv <- atac_cnv %>% left_join(atac_meta, by = "barcode") %>% na.omit()
atac_cnv$Modality <- "ATAC"

# merge the snATAC and multiome cnv and meta data
cnv <- bind_rows(atac_cnv, multi_cnv)

# harmonize cell types between experiments and relevel
cnv <- cnv %>% 
  dplyr::mutate(orig_celltype = celltype) %>%
  dplyr::mutate(celltype = case_when(celltype %in% c("TCELL","BCELL","MONO") ~ "LEUK",
                                     celltype %in% c("PCT","PST") ~ "PT",
                                     celltype %in% c("PT_PROM1") ~ "PT_VCAM1",
                                     celltype %in% c("DCT1","DCT2") ~ "DCT",
                                     celltype %in% c("ATL","TAL1","TAL2","MD") ~ "LOH",
                                     celltype %in% c("ICA","ICB") ~ "IC",
                                     TRUE ~ as.character(celltype)))
cnv$celltype <- factor(cnv$celltype, levels = c("PT","PT_VCAM1","PEC","LOH","DCT","PC","IC","PODO","ENDO","FIB_VSMC_MC","LEUK"))

# binarize loy
cnv <- cnv %>% dplyr::mutate(loy = ifelse(genotype == "LOY",1,0))

# make sure there are no duplicate rows
cnv <- cnv %>% distinct()

# compute wilcoxon pval with bonferroni adjust
wilcox <- cnv %>%
  distinct(barcode, celltype, atac_cnv_total_burden, library_id, genotype) %>%
  group_by(celltype) %>%
  summarize(w=wilcox.test(atac_cnv_total_burden[genotype=="XY"],
                          atac_cnv_total_burden[genotype=="LOY"],
                          paired=FALSE)$p.value) %>%
  ungroup()
wilcox <- wilcox %>%
  dplyr::mutate(padj = w * nrow(wilcox))

cnv <- cnv %>% left_join(wilcox, by = "celltype")

# total atac cnv burden by male cell types with and without loy 
library(ggpubr)
library(RColorBrewer)
palette <- brewer.pal(n = 2, name = "Set1")

plotB <- cnv %>%
  distinct(barcode, celltype, atac_cnv_total_burden, library_id, genotype, padj) %>%
  ggplot(aes(celltype, log10(atac_cnv_total_burden + 1), fill=genotype)) +
  geom_boxplot(outlier.shape = NA) +
  scale_fill_manual(values = palette) +
  theme_bw() +
  xlab("") +
  ylab("Log ATAC CNV Burden") +
  ggtitle("B)") +
  theme(panel.grid = element_blank(),
        legend.title = element_blank(),
        legend.pos = c(0.50,-0.35), legend.direction="horizontal",
        plot.title = element_text(face="bold"), 
        axis.text = element_text(color="black"),
        axis.text.x = element_text(angle=90, vjust=1, hjust=1)) +
  guides(alpha = "none") +
  ggpubr::stat_compare_means(method = 'wilcox.test',
                             label = 'p.signif',
                             show.legend = F,
                             label.y.npc = 0.80)

# fold-change atac cnv burden by male cell types with and without loy
# one estimate for each cell type and library_id
plot_fold_change <- cnv %>%
  distinct(barcode, celltype, atac_cnv_total_burden, library_id, genotype) %>%
  group_by(genotype, celltype, library_id) %>%
  summarize(mean_burden = mean(atac_cnv_total_burden)) %>%
  pivot_wider(names_from = genotype, values_from = mean_burden) %>%
  dplyr::mutate(fold_change = LOY / XY) %>%
  ggplot(aes(celltype, fold_change, fill=celltype)) +
  geom_boxplot() +
  geom_point(position=position_dodge(width=0.75), aes(group=celltype))

#######################
# joint snATAC and multiome models
# figure forest plot of odds ratios and predicted probabilities for glmm of loy ~ celltype + cnv + age (1|library_id)
library(lme4)
library(sjPlot)
library(sjlabelled)
library(sjmisc)
library(ggplot2)
library(stringr)
library(effects)
library(ggeffects)

df <- cnv %>%
  dplyr::select(barcode, library_id, loy, celltype, atac_cnv_total_burden, age, age_group, Modality)

# scale and center cnv burden 
df <- df %>% dplyr::mutate(scaled_cnv = as.numeric(scale(atac_cnv_total_burden, center=0)))
df <- distinct(df)

# glm for effect of celltype and cnv on loy
jmod1 <- glm(loy ~ celltype + scaled_cnv, family = binomial(link = "logit"), data = df)
summary(jmod1)
# plot_jmod1 <- plot_model(jmod1)
# plot marginal effects
plot_jmarg1 <- plot_model(jmod1, type = "emm", terms = c("scaled_cnv","celltype"))

# glm for effect of celltype and cnv with age on loy
jmod2 <- glm(loy ~ celltype + scaled_cnv + age_group, family = binomial(link = "logit"), data = df)
summary(jmod2)
# plot_jmod2 <- plot_model(jmod2)
plot_jmarg2 <- plot_model(jmod2, type = "emm", terms = c("scaled_cnv","celltype"))
plot_jmarg2b <- plot_model(jmod2, type = "emm", terms = c("scaled_cnv","age_group"))
plot_jmarg2c <- plot_model(jmod2, type = "emm", terms = c("celltype","age_group"))

# compare first two models
anova(jmod1, jmod2, test = "Chisq")

# glm for effect of celltype, age and cnv with interaction on loy
jmod3 <- glm(loy ~ celltype + scaled_cnv * age_group, family = binomial(link = "logit"), data = df)
summary(jmod3)
# plot_jmod2 <- plot_model(jmod2)
plot_jmarg3 <- plot_model(jmod3, type = "emm", terms = c("scaled_cnv","celltype"))
plot_jmarg3b <- plot_model(jmod3, type = "emm", terms = c("scaled_cnv","age_group"))
plot_jmarg3c <- plot_model(jmod3, type = "emm", terms = c("celltype","age_group"))
plot_jint3 <- plot_model(jmod3, type = "int")

# glmm with mixed effect for sample with celltype and cnv on loy
jmod4 <- lme4::glmer(loy ~ celltype + scaled_cnv + age_group + (1|library_id),
                     family = binomial(link = "logit"),
                     data = df)
summary(jmod4)
saveRDS(jmod4, here("ckd","figures","jmod4.rds"))


# run again to compute OR using age as continuous variable
jmod5 <- lme4::glmer(loy ~ celltype + scaled_cnv + age + (1|library_id),
                     family = binomial(link = "logit"),
                     data = df)

# plot_jmod4 <- plot_model(jmod4)
plot_jmarg4 <- plot_model(jmod4, type = "emm", terms = c("scaled_cnv","celltype"))
plot_jmarg4b <- plot_model(jmod4, type = "emm", terms = c("scaled_cnv","age_group"))
plot_jmarg4c <- plot_model(jmod4, type = "emm", terms = c("celltype","age_group"))

# subset factor levels in plots
plotC <- plot_model(jmod4, type = "emm", terms = c("scaled_cnv","celltype [PT,PT_VCAM1,PEC,DCT,PODO]"))
plotC <- plotC + theme_bw() + ggtitle("C)") +  
  theme(legend.title=element_blank(),
        legend.pos=c(0.2,0.75),
        legend.key.size = unit(1,"line"),
        plot.title = element_text(face="bold"),
        axis.text = element_text(color="black"),
        panel.grid = element_blank()) +
  ylab("Predicted Probability for LOY") +
  xlab("Scaled ATAC CNV Burden") 

plotD <- plot_model(jmod4, type = "emm", terms = c("scaled_cnv","age_group"))
plotD <- plotD + theme_bw() + ggtitle("D)") +  
  theme(legend.title=element_blank(),
        legend.pos=c(0.15,0.75),
        legend.key.size = unit(1,"line"),
        plot.title = element_text(face="bold"),
        axis.text = element_text(color="black"),
        panel.grid = element_blank()) +
  ylab("Predicted Probability for LOY") +
  xlab("Scaled ATAC CNV Burden") +
  guides(color = guide_legend(reverse=TRUE)) 
  
#########################################
# arrange
# library(cowplot)
# pl <- list(plotA, plotB, plotC, plotD)

# pdf(here(plots,"figure5.pdf"),width = 7.5, height = 10)
# plot_grid(plotlist = pl, align="v")
# dev.off()

# pdf(here(plots,"figure5.pdf"))
# plot_grid(plotlist = pl, align="v")
# dev.off()

###################################
# atac density for chr7 and chr10 XY vs LOY vs XX
atac_frags <- readRDS(here("ckd","atac_aggr_prep","atac_frags.rds")) %>% as.data.frame()
atac_frags <- atac_frags %>% dplyr::select(barcode, seqnames, scaled_atac_frags, malesex, library_id)
meta <- readRDS(here("ckd","atac_aggr_prep","step6_meta.rds"))
meta <- dplyr::select(meta, barcode, celltype, library_id, akd_genotype) %>%
  dplyr::rename(Genotype = akd_genotype)
atac_frags <- left_join(atac_frags, meta, by = c("barcode","library_id"))

multi_frags <- readRDS(here("ckd","multi_aggr_prep","step3_multi_counts.rds")) %>%
  dplyr::select(barcode, seqnames, scaled_atac_frags, malesex, library_id, gmm_genotype) %>%
  dplyr::rename(Genotype = gmm_genotype)
multi_meta <- readRDS(here("ckd","multi_aggr_prep","step3_meta_multi_loy.rds")) %>%
  dplyr::select(barcode, celltype, library_id)
multi_frags <- left_join(multi_frags, multi_meta, by = c("barcode","library_id"))
all_frags <- bind_rows(atac_frags, multi_frags)

# harmonize cell types between experiments and relevel
all_frags <- all_frags %>% 
  dplyr::mutate(orig_celltype = celltype) %>%
  dplyr::mutate(celltype = case_when(celltype %in% c("TCELL","BCELL","MONO") ~ "LEUK",
                                     celltype %in% c("PCT","PST") ~ "PT",
                                     celltype %in% c("PT_PROM1") ~ "PT_VCAM1",
                                     celltype %in% c("DCT1","DCT2") ~ "DCT",
                                     celltype %in% c("ATL","TAL1","TAL2","MD") ~ "LOH",
                                     celltype %in% c("ICA","ICB") ~ "IC",
                                     TRUE ~ as.character(celltype)))
all_frags$celltype <- factor(all_frags$celltype, levels = c("PT","PT_VCAM1","PEC","LOH","DCT","PC","IC","PODO","ENDO","FIB_VSMC_MC","LEUK"))

toplot <- all_frags %>%
  dplyr::filter(Genotype %in% c("XY","LOY")) %>%
  dplyr::mutate(Genotype = factor(Genotype, levels = c("LOY","XY"))) %>%
  dplyr::filter(celltype %in% c("PT_VCAM1","PT")) %>%
  dplyr::filter(seqnames %in% c("chr1","chr7","chr10","chrX"))

library(RColorBrewer)
palette <- brewer.pal(n = 2, name = "Set1")
plotE <- toplot %>%
  ggplot(aes(scaled_atac_frags, alpha=0.5, fill=Genotype)) +
  geom_density(aes(y = after_stat(scaled)), adjust=2) +
  guides(alpha = guide_legend(element_blank())) +
  facet_wrap(~celltype + seqnames, nrow = 2, 
             labeller = label_wrap_gen(multi_line=FALSE)) +
  xlim(c(0,2)) +
  scale_fill_manual(values = palette) +
  theme_bw() +
  theme(plot.title = element_text(face="bold"),
        axis.text = element_text(color="black")) +
  guides(alpha = "none") +
  guides(fill = guide_legend(title = "Genotype")) +
  ylab("Density") +
  xlab("Scaled ATAC Coverage") +
  ggtitle("E)")

library(gridExtra)
library(cowplot)
row1 <- cowplot::plot_grid(plotA, plotB, rel_widths = c(1,2), align = "hv")
row2 <- grid.arrange(plotC, plotD, widths = c(1,1))
row3 <- grid.arrange(plotE, widths = 2)

pdf(here(plots,"figure6.pdf"),width = 7.5, height=10)
grid.arrange(row1, row2, row3)
dev.off()

# # scale the cnv burden
# meta <- meta %>%
#   group_by(library_id, chr) %>%
#   mutate(scaled_rna_chrom_cnv = scale(rna_cnv_chrom_burden, center = 0)) %>%
#   mutate(scaled_atac_chrom_cnv = scale(atac_cnv_chrom_burden, center = 0)) %>%
#   group_by(library_id) %>%
#   mutate(scaled_rna_total_cnv = scale(rna_cnv_total_burden, center = 0)) %>%
#   mutate(scaled_atac_total_cnv = scale(atac_cnv_total_burden, center = 0)) %>%
#   ungroup()

# palette <- viridis(9)
# palette[1] <- "white"

# chrom_sel <- "chr7"
# int <- meta %>%
#   dplyr::filter(chr == chrom_sel) %>%
#   dplyr::summarize(median_scaled_atac_chrom_cnv = median(scaled_atac_chrom_cnv),
#                    median_scaled_rna_chrom_cnv = median(scaled_rna_chrom_cnv))
# yint <- int$median_scaled_atac_chrom_cnv
# xint <- int$median_scaled_rna_chrom_cnv

# res <- meta %>%
#   dplyr::group_by(gmm_genotype) %>%
#   dplyr::filter(chr == chrom_sel) %>%
#   summarize(atac = scaled_atac_chrom_cnv > yint, rna = scaled_rna_chrom_cnv > xint) %>%
#   table() %>%
#   as.data.frame() %>% 
#   group_by(gmm_genotype) %>%
#   dplyr::mutate(total_cells = sum(Freq)) %>%
#   dplyr::filter(atac == TRUE, rna == TRUE) %>%
#   dplyr::mutate(prop_multi_cnv = Freq / total_cells)

# # note the addition of a small amount of noise with jitter() function
# # bc o/w geom_density cannot compute band estimates
# p1 <- meta %>%
#   dplyr::filter(chr == chrom_sel) %>%
#   distinct(barcode, celltype, scaled_atac_chrom_cnv, scaled_rna_chrom_cnv, library_id, gmm_genotype) %>%
#   dplyr::filter(gmm_genotype %in% c("XY")) %>%
#   dplyr::mutate(scaled_atac_chrom_cnv = jitter(scaled_atac_chrom_cnv),
#                 scaled_rna_chrom_cnv = jitter(scaled_rna_chrom_cnv)) %>%
#   na.omit() %>%
#   ggplot(aes(scaled_atac_chrom_cnv, scaled_rna_chrom_cnv)) +
#   geom_density2d_filled(aes(fill = after_stat(level)), bins=9) +
#   scale_fill_manual(values = palette) +
#   xlim(c(0,2)) +
#   ylim(c(0,2)) +
#   coord_fixed() +
#   geom_hline(yintercept = yint) +
#   geom_vline(xintercept = xint) +
#   annotate(geom="text",
#            x = 1,
#            y = 1,
#            label = paste0("Prop = ", 100*round(res$prop_multi_cnv[res$gmm_genotype == "XY"], digits=2),"%"))

# p2 <- meta %>%
#   dplyr::filter(chr == chrom_sel) %>%
#   distinct(barcode, celltype, scaled_atac_chrom_cnv, scaled_rna_chrom_cnv, library_id, gmm_genotype) %>%
#   dplyr::filter(gmm_genotype %in% c("LOY")) %>%
#   dplyr::mutate(scaled_atac_chrom_cnv = jitter(scaled_atac_chrom_cnv),
#                 scaled_rna_chrom_cnv = jitter(scaled_rna_chrom_cnv)) %>%
#   na.omit() %>%
#   ggplot(aes(scaled_atac_chrom_cnv, scaled_rna_chrom_cnv)) +
#   geom_density2d_filled(aes(fill = after_stat(level)), bins=9) +
#   scale_fill_manual(values = palette) +
#   xlim(c(0,2)) +
#   ylim(c(0,2)) +
#   coord_fixed() +
#   geom_hline(yintercept = yint) +
#   geom_vline(xintercept = xint) +
#   annotate(geom="text",
#            x = 1,
#            y = 1,
#            label = paste0("Prop = ", 100*round(res$prop_multi_cnv[res$gmm_genotype == "LOY"], digits=2),"%"))

# grid.arrange(p1,p2, ncol=2)
  

# # facet wrap chrom for atac or rna chrom density
# meta <- meta %>%
#   group_by(library_id, chr) %>%
#   mutate(scaled_atac_cnv = scale(atac_cnv_chrom_burden, center=0),
#          scaled_rna_cnv = scale(rna_cnv_chrom_burden, center=0))

# meta %>%
#   na.omit() %>%
#   dplyr::filter(celltype %in% c("PCT","PST")) %>%
#   distinct(barcode, celltype, scaled_atac_cnv, scaled_rna_cnv, library_id, gmm_genotype) %>%
#   dplyr::filter(gmm_genotype %in% c("XY","LOY")) %>%
#   ggplot(aes(scaled_atac_cnv, fill=gmm_genotype, alpha=0.5)) +
#   geom_density(aes(y = ..scaled..)) +
#   facet_wrap(~chr) +
#   xlim(c(0,5))

# meta %>%
#   na.omit() %>%
#   dplyr::filter(celltype %in% c("PT_VCAM1")) %>%
#   distinct(barcode, celltype, scaled_atac_cnv, scaled_rna_cnv, library_id, gmm_genotype) %>%
#   dplyr::filter(gmm_genotype %in% c("XY","LOY")) %>%
#   ggplot(aes(scaled_atac_cnv, fill=gmm_genotype, alpha=0.5)) +
#   geom_density(aes(y = ..scaled..)) +
#   facet_wrap(~chr) +
#   xlim(c(0,5))

# # pairwise wilcox test by cell type for cells with and without loy
# meta$gmm_genotype <- as.factor(meta$gmm_genotype)
# meta %>%
#   dplyr::filter(gmm_genotype %in% c("XY","LOY")) %>%
#   group_by(celltype) %>%
#   do(w = wilcox.test(propcnv~gmm_genotype, data=., paired=FALSE)) %>% 
#   summarise(celltype, Wilcox = w$p.value)

# # wilcox test for all cells with and without loy
# meta %>%
#   dplyr::filter(gmm_genotype %in% c("XY","LOY")) %>%
#   na.omit() %>%
#   ungroup() %>%
#   summarise(w=wilcox.test(propcnv~gmm_genotype, data=., paired=FALSE)$p.value)

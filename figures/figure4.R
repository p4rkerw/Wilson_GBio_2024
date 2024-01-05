# dpcr validation and analysis

library(ggplot2)
library(dplyr)
library(ggplot2)
library(dplyr)
library(here)

df <- read.csv(here("ckd","digital_pcr","Validation run RPTEC & HEK 293T (238,110)-Quantification.csv")) %>%
  dplyr::filter(Well != "Average")
clinical_meta <- read.xlsx(here("ckd","clinical_meta.xlsx"))

sample_order <- c(100,90,80,50,0)
df$group <- factor(df$Name, levels = sample_order)
df <- df %>% dplyr::mutate(expected_copy_number = Name / (Name + 2*(100-Name)))
df$expected_copy_number <- as.character(round(df$expected_copy_number, 2))
df$expected_copy_number <- factor(df$expected_copy_number, levels = c(1, 0.82, 0.67, 0.33, 0))

df <- df %>%
  dplyr::mutate(copy_number_aby_as_ref = Conc..cp.uL / Conc..cp.uL.2) 

means <- aggregate(copy_number_aby_as_ref ~  expected_copy_number, df, mean)
means$copy_number_aby_as_ref <- round(means$copy_number_aby_as_ref, 2)

p1 <- df %>%
  ggplot(aes(expected_copy_number, copy_number_aby_as_ref, label=copy_number_aby_as_ref, fill=expected_copy_number)) +
  geom_boxplot(aes(color = expected_copy_number)) +
  geom_jitter(position=position_jitter(width=0.1)) +
  geom_text(data = means, nudge_y=0.05) +
  theme_bw() +
  theme(text = element_text(size = 14),
        plot.title = element_text(face="bold"),
        axis.text = element_text(color="black")) +
  scale_fill_discrete(labels = c("100:0","90:10","80:20","50:50","0:100")) +
  xlab("Expected Copy Ratio") +
  ylab("Measured Copy Ratio") +
  guides(fill = guide_legend(title = "Input M:F"), color="none") +
  ggtitle("A)")

df <- read.csv(here("ckd","digital_pcr","digital_pcr_loy_single_cell.csv"))

sample_order <- unique(df$Name)
sample_order <- sample_order[sample_order != "CONTROL"]
sample_order <- c("CONTROL", sample_order)
sample_order <- sample_order[1:10]

df$group <- factor(df$Name, levels = sample_order)

df <- df %>%
  dplyr::mutate(copy_number_aby_as_ref = Conc..cp.uL / Conc..cp.uL.2) 

palette <- brewer.pal(6, "Set1")
palette <- c("lightgray", palette)

p2 <- df %>%
  ggplot(aes(group, copy_number_aby_as_ref)) +
  geom_boxplot(aes(fill=group)) +
  scale_fill_manual(values = palette) +
  geom_jitter(position=position_jitter(width=0.1), aes(shape=Run)) +
  theme_bw() +
  theme(text = element_text(size = 14),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        plot.title = element_text(face="bold"),
        axis.text = element_text(color="black")) +
  guides(shape = guide_legend(title = "dPCR Batch", order = 2),
         fill=guide_legend(title="Sample", order = 1)) +
  xlab("") +
  ylab("Measured Copy Ratio") +
  coord_cartesian() +
  ggtitle("B)")

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

# read in loy anno for multiomes
multi_meta <- readRDS(here("ckd","multi_aggr_prep","step3_meta_multi_loy.rds"))

# format multiome meta data
multi_meta <- multi_meta %>%
  dplyr::filter(gmm_genotype %in% c("XY","LOY")) %>%
  dplyr::select(barcode, gmm_genotype, celltype, library_id) %>%
  dplyr::rename(genotype = gmm_genotype) %>% 
  left_join(clinical_meta, by = "library_id")

# join multiome meta data with cnv 
multi_meta$Modality <- "MULTI"

# read in snATAC meta data
atac_meta <- readRDS(here("ckd","atac_aggr_prep","step6_meta.rds"))
atac_meta <- atac_meta %>%
  dplyr::filter(akd_genotype %in% c("XY","LOY")) %>%
  dplyr::select(barcode, akd_genotype, celltype, library_id) %>%
  dplyr::rename(genotype = akd_genotype) %>% 
  left_join(clinical_meta, by = "library_id")
atac_meta$Modality <- "ATAC"


# merge the snATAC and multiome cnv and meta data
cnv <- bind_rows(atac_meta, multi_meta)

# binarize loy
cnv <- cnv %>% dplyr::mutate(loy = ifelse(genotype == "LOY",1,0))

# make sure there are no duplicate rows
cnv <- cnv %>% distinct()

#######################
# joint snATAC and multiome models
library(stringr)

df2 <- cnv %>%
  dplyr::select(barcode, library_id, loy, Modality) %>%
  distinct()

df2 <- table(df2$library_id, df2$loy) %>% 
  as.data.frame() %>% 
  group_by(Var1) %>% 
  mutate(total_cells = sum(Freq)) %>% 
  mutate(prop_loy = 1- Freq / total_cells) %>% 
  dplyr::filter(Var2 == 0)

df2 <- df2 %>% dplyr::mutate(library_id = Var1) 
df <- df %>% dplyr::mutate(library_id = Name)
df <- df %>% left_join(df2, by = "library_id")
df <- df %>% na.omit()

res <- lm(1 - copy_number_aby_as_ref ~ prop_loy, data = df)
summary(res)$r.squared

#define function to extract overall p-value of model
overall_p <- function(my_model) {
  f <- summary(my_model)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}

#extract overall p-value of model
pval <- round(overall_p(res), 8)
palette <- brewer.pal(6, "Set1")

p3 <- df %>%
  ggplot(aes(prop_loy, 1- copy_number_aby_as_ref, color=group)) +
  scale_color_manual(values = palette) +
  geom_point(aes(shape = Run), size = 2) +
  geom_smooth(method = lm, col='darkgray', size=1) +
  theme_bw() + 
  ylab("dPCR LOY") +
  xlab("Single-cell LOY") +
  theme(text = element_text(size = 14),
        plot.title = element_text(face="bold"),
        axis.text = element_text(color="black")) +
  guides(shape = guide_legend(title = "dPCR Batch", order = 2),
         color=guide_legend(title="Sample", order = 1)) +
  annotate(geom="text",
           x = 0.15,
           y = 0.1,
           label = paste0("r^2=",round(summary(res)$r.squared,2),
                          "\np=",pval),
           hjust=0) +
  xlim(c(0,0.3)) +
  ylim(c(0,0.3)) +
  coord_fixed() +
  ggtitle("C)")


# layout figure
lay <- rbind(c(1,2),
             c(3,NA))
gs <- list(p1,p2,p3)


pdf(here(plots,"figure4.pdf"),width = 7.5, height = 7.5)
grid.arrange(grobs = gs, layout_matrix = lay, heights = c(1,1))
dev.off()

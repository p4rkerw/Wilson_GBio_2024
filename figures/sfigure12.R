# odds ratios for loy in full model

# read in model from figure5.R
mod <- readRDS(here("ckd","figures","jmod4.rds"))

plotA <- plot_model(mod, p.adjust=TRUE)

plotA <- plotA + 
  ggtitle("A)") +
  geom_hline(yintercept = 1, linetype = "dotted") + 
  theme_bw() +
  theme(plot.title = element_text(face="bold"), axis.text = element_text(color="black")) + 
  ylab("Odds Ratios for LOY relative to PT")

pdf(here(plots,"sfigure12.pdf"), width = 7.5, height = 5)
  print(plotA)
dev.off()

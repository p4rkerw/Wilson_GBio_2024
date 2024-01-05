# show marginal effects for all cell types
jmod4 <- readRDS(here("ckd","figures","jmod4.rds"))
  
# subset factor levels in plots
plotA <- plot_model(jmod4, type = "emm", terms = c("scaled_cnv","celltype"))
plotA <- plotA + theme_bw() + ggtitle("A)") +  
  theme(legend.title=element_blank(),
        legend.pos=c(0.2,0.75),
        legend.key.size = unit(1,"line"),
        plot.title = element_text(face="bold"),
        axis.text = element_text(color="black")) +
  ylab("Predicted Probability for LOY") +
  xlab("Scaled ATAC CNV Burden") 
  
pdf(here(plots,"sfigure13.pdf"), width = 7.5, height = 6)
 print(plotA)
dev.off()

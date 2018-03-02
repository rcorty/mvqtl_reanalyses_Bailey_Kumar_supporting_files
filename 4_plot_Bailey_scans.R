library(vqtl)
library(gridExtra)
library(grid)
library(ggplot2)

bailey_replication <- readRDS('Bailey_scans_1000_perms_smaller_files.RDS')

sos <- bailey_replication$so
sops <- bailey_replication$sop
sovs <- bailey_replication$sov

plotting_names <- c('Total Distance', 'Ambulatory Episodes', 'Percent Time Resting', 'Percent Time in Center', 'Average Velocity', 'Rearing')


# replication scans
orig_phen_names <- names(sos)[c(2, 4, 6, 1, 5, 3)]
names(plotting_names) <- orig_phen_names

pdf(file = 'bailey_plots/bailey_replication_scanones.pdf',
    width = 8,
    height = 9)
par(mfrow = c(6, 1), mar = c(2, 3, 3, 2))
for (phen_name in orig_phen_names) {

  so <- sos[[phen_name]]
  sop <- sops[[phen_name]]

  plot(x = so, bandcol = 'gray', ylim = c(0, 4.5))
  add.threshold(out = so,
                perms = sop,
                alpha = 0.05)
  abline(h = 3.3, lty = 2)
  mtext(text = 'LOD', side = 2, line = 1.8, cex = 0.8)
  mtext(text = plotting_names[phen_name], side = 3, line = 0.5)

  text(x = 360, y = 3.3, labels = '3.3', adj = c(0.5, 1.2))
  text(x = 360,
       y = quantile(x = sop, probs = 0.95),
       labels = broman::myround(x = quantile(x = sop, probs = 0.95), digits = 2),
       adj = c(0.5, -0.2))
}
dev.off()


# scanones and scanonevars on Box-Cox transformed traits
plot_sov_and_so <- function(n) {
  plot(x = sovs[[n]],
       y = sos[[n]],
       plot.title = plotting_names[n], ymax = 2.5) +
    theme(legend.position = 'none')
}

bc_phen_names <- c('TOTDIST', 'AMBEPIS', 'bc_PCTREST', 'bc_PCTT10', 'bc_AVGVELO', 'bc_TOTREAR')
names(plotting_names) <- bc_phen_names

gs <- lapply(X = bc_phen_names, FUN = plot_sov_and_so)
g <- arrangeGrob(grobs = gs, ncol = 1)
ggsave(plot = g, filename = 'bailey_plots/bailey_scans_bc.pdf', height = 11, width = 8)



# scanones and scanonevars on RINT transformed traits
rint_phen_names <- paste0('rint_', c('TOTDIST', 'AMBEPIS', 'PCTREST', 'PCTT10', 'AVGVELO', 'TOTREAR'))
names(plotting_names) <- rint_phen_names

gs <- lapply(X = rint_phen_names, FUN = plot_sov_and_so)
g <- arrangeGrob(grobs = gs, ncol = 1)
ggsave(plot = g, filename = 'bailey_plots/bailey_scans_rint.pdf', height = 11, width = 8)


# scanonevar for other transforms of TOTREAR
plot(x = sovs$TOTREAR)
ggsave(filename = 'bailey_plots/TOTREAR_no_transform.pdf', height = 1.5, width = 8)

plot(x = sovs$log_TOTREAR)
ggsave(filename = 'bailey_plots/TOTREAR_log.pdf', height = 1.5, width = 8)

plot(x = sovs$sqrt_TOTREAR)
ggsave(filename = 'bailey_plots/TOTREAR_sqrt.pdf', height = 1.5, width = 8)

plot(x = sovs$TOTREAR_poisson)
ggsave(filename = 'bailey_plots/TOTREAR_poisson.pdf', height = 1.5, width = 8)




# phenotype at top marker plot
library(ggplot2)
c$pheno$sex <- factor(x = c$pheno$SEX, labels = c('female', 'male'))
phenotype_at_marker_plot(cross = c,
                         phenotype_name = 'bc_TOTREAR',
                         marker_name = '2.65.484',
                         color_by = 'sex',
                         genotype_labels = c('B6', 'Het', 'C58')) +
  xlab('Chr 2, 65Mb marker') +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5)) +
  ggtitle(label = 'Rearing Behavior (counts)')
ggsave(filename = 'bailey_plots/TOTREAR_bc_phen_plot.pdf', height = 4, width = 5)


# mean_var_plot
c$pheno$chr2_marker <- factor(c$geno$`2`$data[,'2.65.484'])

totrear_lm <- lm(formula = bc_TOTREAR ~ SEX + chr2_marker,
                 data = c$pheno,
                 na.action = na.exclude)

summary(totrear_lm)

lm_pred <- data.frame(sex = c$pheno$SEX, chr2_marker = c$pheno$chr2_marker)
lm_pred$lm_pred <- predict(object = totrear_lm)
lm_pred$sigma <- sd(residuals(totrear_lm), na.rm = TRUE)



mean_var_plot_model_based(cross = c,
                          phenotype.name = 'bc_TOTREAR',
                          focal.groups = c('sex', 'chr2_marker'),
                          point_size = 3,
                          title = 'Mean and Variance Effect Estimates\nfor Rearing Events',
                          genotype.names = c('B6', 'Het', 'C58')) +
  geom_point(data = unique(lm_pred),
             mapping = aes(x = lm_pred, y = sigma), color = 'black', shape = 4, size = 3) +
  scale_shape(name = 'chr 2,\n65Mb marker') +
  scale_color_discrete(name = 'Sex') +
  theme(legend.background = element_rect(fill = 'white'))
ggsave(filename = 'bailey_plots/TOTREAR_bc_mean_var_plot.pdf', height = 4, width = 5)

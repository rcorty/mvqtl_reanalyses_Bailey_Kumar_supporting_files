library(vqtl)
library(ggplot2)
library(dplyr)


MPD_URL <- 'https://phenomedoc.jax.org/QTL_Archive/kumar_2014/Kumar2014_Data.csv'
c <- read.cross(format = 'csv',  file = url(description = MPD_URL))


# read the already-done scans if you don't want to do the perms
sos <- readRDS(file = 'Kumar_scans_1000_perms.RDS')[['so']]
sops <- readRDS(file = 'Kumar_scans_1000_perms.RDS')[['sop']]
sovs <- readRDS(file = 'Kumar_scans_1000_perms.RDS')[['sov']]


# cocaine response scans
plot(x = sovs$X130avg, y = sos$X130avg, plot.title = '30 Minute Cocaine Response', alpha_size = 2)
ggsave(filename = 'Kumar_30min_cocaine_scan.pdf', width = 10, height = 2)

plot(x = sovs$X1sum.60min, y = sos$X1sum.60min, plot.title = '60 Minute Cocaine Response', alpha_size = 2)
ggsave(filename = 'Kumar_60min_cocaine_scan.pdf', width = 10, height = 2)



# average counts scan
plot(x = sovs$Avg_Counts, y = sos$Avg_Counts, plot.title = 'Circadian Wheel Running Activity (revolutions/minute)', alpha_size = 2)
ggsave(filename = 'Kumar_avg_counts_scan.pdf', width = 10, height = 2)


# data for chr 6 QTL
sov$Avg_Counts$result %>% filter(mQTL.empir.p < 0.01)

# zoom to chromosome 6
plot(x = sovs$Avg_Counts, y = sos$Avg_Counts, chrs = 6, plot.title = 'Circadian Wheel Running Activity (revolutions/minute)', legend_pos = 'none', alpha_pos = 'right')
ggsave(filename = 'Kumar_avg_counts_chr6.pdf', width = 3, height = 4)



# Figure 2a
c$pheno$sex <- factor(c$pheno$sex, labels = c('female', 'male'))

actogram_mice <- c('WT010G3NCIB6F20702', 'WT010G3NCIB6F20629', 'WT010G3NCIB6F20654', 'WT010G3NCIB6F20618',
                   'WT010G3NCIB6F20514', 'WT010G3NCIB6F20556', 'WT010G3NCIB6F20810', 'WT010G3NCIB6F20617',
                   'WT010G3NCIB6F20802', 'WT010G3NCIB6F20835', 'WT010G3NCIB6F20441', 'WT010G3NCIB6F20427')
c$pheno$in_actogram <- c$pheno$NA. %in% actogram_mice
set.seed(2)   # so we can get the same jitter again
phenotype_at_marker_plot(cross = c,
                         phenotype_name = 'Avg_Counts',
                         marker_name = 'rs30314218',
                         color_by = 'sex',
                         # shape_by = 'in_actogram',
                         point_size = 'in_actogram',
                         point_alpha = 0.7,
                         genotype_labels = c('C57BL/6J', 'Het', 'C57BL/6N')) +
  xlab('chr6, rs30314218') +
  ylab('average wheel speed') +
  ggtitle('Circadian Wheel Running Activity (revolutions/minute)') +
  scale_shape_manual(values = c(16, 8), guide = FALSE) +
  scale_size_discrete(guide = FALSE) +
  theme(plot.title = element_text(hjust = 0.5))
ggsave(filename = 'Kumar_avg_counts_phenotype_plot.pdf', height = 4, width = 5)


# Figure 2b
c$pheno$chr6_marker <- factor(c$geno$`6`$data[,'rs30314218'])

ac_lm <- lm(formula = Avg_Counts ~ factor(sex) + chr6_marker,
            data = c$pheno,
            na.action = na.exclude)

ac_dglm <- dglm::dglm(formula = Avg_Counts ~ factor(sex) + chr6_marker,
                      dformula = ~ factor(sex) + chr6_marker,
                      data = c$pheno)

summary(ac_lm)
summary(ac_dglm)

lm_pred <- data.frame(sex = c$pheno$sex, chr6_marker = c$pheno$chr6_marker)
lm_pred$lm_pred <- predict(object = ac_lm)
lm_pred$sigma <- sd(residuals(ac_lm), na.rm = TRUE)


mean_var_plot_model_based(cross = c,
                          phenotype.name = 'Avg_Counts',
                          focal.groups = c('sex', 'rs30314218'),
                          point_size = 4,
                          title = 'Mean and Variance Effect Estimates for\n Circadian Wheel Running Activity (revolutions/minute)',
                          genotype.names = c('C57BL/6J', 'Het', 'C57BL/6N')) +
  geom_point(data = unique(lm_pred),
             mapping = aes(x = lm_pred, y = sigma), color = 'black', shape = 4, size = 3)
ggsave(filename = 'Kumar_avg_counts_mean_var_plot.pdf', height = 4, width = 5)



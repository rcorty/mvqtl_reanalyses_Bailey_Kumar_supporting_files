library(qtl)
library(vqtl)
library(evd)

PGEV <- function(q, gev, ...) {

  if ((!'estimate' %in% names(gev))) {
    stop("argument 'gev' must have an element named 'estimate' (all gev objects do)")
  }
  if (length(gev$estimate) != 3) {
    stop("gev$estimate must have three elements")
  }

  return(pgev(q = q,
              loc = gev$estimate[1],
              scale = gev$estimate[2],
              shape = gev$estimate[3], ...))
}

MPD_URL <- 'https://phenomedoc.jax.org/QTL_Archive/kumar_2014/Kumar2014_Data.csv'
c <- read.cross(format = 'csv',  file = url(description = MPD_URL))




c <- calc.genoprob(cross = c, step = 2)
c$pheno %>% glimpse()

# LL at lambda=1 is close enough to maxLL that it's reasonable not to transform
MASS::boxcox(object = c$pheno$Avg_Counts ~ c$pheno$sex)

hist(x = c$pheno$Avg_Counts, breaks = 30)


# scanone
num_perms <- 1000
num_cores <- 4

so <- scanone(cross = c, pheno.col = 'Avg_Counts', addcovar = c$pheno$sex)
sop <- scanone(cross = c, pheno.col = 'Avg_Counts', addcovar = c$pheno$sex, n.perm = num_perms, n.cluster = num_cores, verbose = FALSE)
evd <- fgev(x = sop)
so$empir.p <- PGEV(q = so$lod, gev = evd, lower.tail = FALSE)

# plot LODs
plot(x = so, bandcol = 'gray')

# plot FWER-controlling p-value
plot(x = so, bandcol = 'gray', lodcolumn = 2)

# takes ~30 seconds
sov <- scanonevar(cross = c,
                  mean.formula = Avg_Counts ~ sex + mean.QTL.add + mean.QTL.dom,
                  var.formula = ~ sex + var.QTL.add + var.QTL.dom)

# plot LODs
plot(x = sov, y = so)

# calculate FWER-controlling p-values
# careful -- this will take 0.5 mins * 1000 perms / 4 cores = 2 hours on a typical computer
# recommend to use a many-core computer to make it faster...40 cores -> 12 minutes
# sov <- scanonevar.perm(sov = sov, n.perms = num_perms, n.cores = num_cores)


saveRDS(object = list(so = so, sov = sov),
        file = 'Kumar2014_avg_counts_scan.RDS')

# read the already-done scans if you don't want to do the perms
so <- readRDS(file = 'Kumar2014_avg_counts_scan.RDS')[['so']]
sov <- readRDS(file = 'Kumar2014_avg_counts_scan.RDS')[['sov']]


# Figure 1
plot(x = sov, y = so, plot.title = 'Circadian Wheel Running Activity (revolutions/minute)', alpha_size = 2)
ggsave(filename = 'Kumar_avg_counts_scan.pdf', width = 10, height = 2)

# data for chr 6 QTL
sov$result %>% filter(mQTL.empir.p < 0.01)

# zoom to chromosome 6
plot(x = sov, y = so, chrs = 6, plot.title = 'Circadian Wheel Running Activity (revolutions/minute)', legend_pos = 'none', alpha_pos = 'right')
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


c$pheno$sex <- factor(c$pheno$sex, labels = c('female', 'male'))
mean_var_plot_model_based(cross = c,
                          phenotype.name = 'Avg_Counts',
                          focal.groups = c('sex', 'rs30314218'),
                          point_size = 4,
                          title = 'Mean and Variance Effect Estimates for\n Circadian Wheel Running Activity (revolutions/minute)',
                          genotype.names = c('C57BL/6J', 'Het', 'C57BL/6N')) +
  geom_point(data = unique(lm_pred),
             mapping = aes(x = lm_pred, y = sigma), color = 'black', shape = 4, size = 3)
ggsave(filename = 'Kumar_avg_counts_mean_var_plot.pdf', height = 4, width = 5)



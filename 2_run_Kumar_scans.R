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


# LL at lambda=1 is close enough to maxLL that it's reasonable not to transform
MASS::boxcox(object = c$pheno$Avg_Counts ~ c$pheno$sex)


# cocaine response phenotypes published in Kumar2013
phen_names <- c('X130avg', 'X1sum.60min', 'Avg_Counts')


num_cores <- 40
num_perms <- 25*num_cores

so <- sov <- sop <- list()
for (phen_name in phen_names) {

  message('Starting ', phen_name)

  so[[phen_name]] <- scanone(cross = c, pheno.col = phen_name, addcovar = as.numeric(c$pheno$sex))

  sov[[phen_name]] <- scanonevar(cross = c,
                                 mean.formula = formula(paste0(phen_name, '~ sex + mean.QTL.add + mean.QTL.dom')),
                                 var.formula = ~ sex+ var.QTL.add + var.QTL.dom)

  if (num_perms) {

    message('Starting perms')

    sop[[phen_name]] <- scanone(cross = c, pheno.col = phen_name, addcovar = as.numeric(c$pheno$sex), n.perm = num_perms, n.cluster = num_cores, verbose = FALSE)
    evd <- fgev(x = sop[[phen_name]])
    so[[phen_name]]$empir.p <- PGEV(q = so[[phen_name]]$lod, gev = evd, lower.tail = FALSE)

    sov[[phen_name]] <- scanonevar.perm(sov = sov[[phen_name]], n.perms = num_perms, n.cores = num_cores, random.seed = 27599)
  }

}


saveRDS(object = list(so = so, sop = sop, sov = sov),
        file = paste0('Kumar_scans_', num_perms, '_perms.RDS'))


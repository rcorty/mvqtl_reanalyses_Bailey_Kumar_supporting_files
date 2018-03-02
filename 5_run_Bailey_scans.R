library(qtl)
library(vqtl)
library(evd)
library(MASS)

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
RINT <- function(v, c = 3/8, ...) {
  qnorm(p = (rank(v) - c)/(length(v) - 2*c + 1), ...)
}

MPD_URL <- 'https://phenomedoc.jax.org/QTL_Archive/bailey_2008/Bailey2008_B6xC58J_B37_Data.csv'
c <- read.cross(format = 'csv',  file = url(description = MPD_URL),
                genotypes = c('B', 'H', 'C'),
                na.strings = c('O', '-'))
c <- calc.genoprob(cross = jittermap(object = c), step = 2)

# trim extreme outliers
c$pheno$PCTT10[c$pheno$PCTT10 > 60] <- NA
c$pheno$TOTDIST[c$pheno$TOTDIST > 6000] <- NA
c$pheno$TOTREAR[c$pheno$TOTREAR < 10] <- NA
c$pheno$ANXFACT[c$pheno$ANXFACT > 3.3] <- NA


boxcox(object = c$pheno$PCTT10 ~ c$pheno$SEX)   # 0.75
boxcox(object = c$pheno$TOTDIST ~ c$pheno$SEX)  # 1
boxcox(object = c$pheno$TOTREAR + 1 ~ c$pheno$SEX) # 0.25
boxcox(object = c$pheno$AMBEPIS ~ c$pheno$SEX) # 1
boxcox(object = c$pheno$AVGVELO ~ c$pheno$SEX) # 0
boxcox(object = c$pheno$PCTREST ~ c$pheno$SEX) # 0

c$pheno$bc_PCTT10 <- c$pheno$PCTT10^(0.75)
c$pheno$bc_TOTREAR <- c$pheno$TOTREAR^(0.25)
c$pheno$bc_AVGVELO <- log(c$pheno$AVGVELO)
c$pheno$bc_PCTREST <- log(c$pheno$PCTREST)

c$pheno$rint_PCTT10  <- RINT(c$pheno$PCTT10)
c$pheno$rint_TOTDIST <- RINT(c$pheno$TOTDIST)
c$pheno$rint_TOTREAR <- RINT(c$pheno$TOTREAR)
c$pheno$rint_AMBEPIS <- RINT(c$pheno$AMBEPIS)
c$pheno$rint_AVGVELO <- RINT(c$pheno$AVGVELO)
c$pheno$rint_PCTREST <- RINT(c$pheno$PCTREST)

c$pheno$sqrt_TOTREAR <- sqrt(c$pheno$TOTREAR)
c$pheno$log_TOTREAR <- log(c$pheno$TOTREAR)


orig_phen_names <- names(c$pheno)[4:9]
bc_phen_names <- grep(pattern = 'bc', x = names(c$pheno), value = TRUE)
rint_phen_names <- grep(pattern = 'rint', x = names(c$pheno), value = TRUE)
totrear_additional_transform_names <- c('log_TOTREAR', 'sqrt_TOTREAR')

phen_names <- c(orig_phen_names, bc_phen_names, rint_phen_names, totrear_additional_transform_names)


num_cores <- 40
num_perms <- 25*num_cores

so <- sov <- sop <- list()
for (phen_name in phen_names) {

  message('Starting ', phen_name)

  so[[phen_name]] <- scanone(cross = c, pheno.col = phen_name, addcovar = as.numeric(c$pheno$SEX))

  sov[[phen_name]] <- scanonevar(cross = c,
                                 mean.formula = formula(paste0(phen_name, '~ SEX + mean.QTL.add + mean.QTL.dom')),
                                 var.formula = ~ SEX  + var.QTL.add + var.QTL.dom)

  if (num_perms) {

    message('Starting perms')

    sop[[phen_name]] <- scanone(cross = c, pheno.col = phen_name, addcovar = as.numeric(c$pheno$SEX), n.perm = num_perms, n.cluster = num_cores, verbose = FALSE)
    evd <- fgev(x = sop[[phen_name]])
    so[[phen_name]]$empir.p <- PGEV(q = so[[phen_name]]$lod, gev = evd, lower.tail = FALSE)

    sov[[phen_name]] <- scanonevar.perm(sov = sov[[phen_name]], n.perms = num_perms, n.cores = num_cores, random.seed = 27599)
  }

}

# TOTREAR Poisson model
sov[['TOTREAR_poisson']] <- scanonevar(cross = c,
                                       mean.formula = TOTREAR ~ SEX + mean.QTL.add + mean.QTL.dom,
                                       var.formula = ~ SEX  + var.QTL.add + var.QTL.dom,
                                       glm_family = 'poisson')

sov[['TOTREAR_poisson']] <- scanonevar.perm(sov = sov[['TOTREAR_poisson']],
                                            n.perms = num_perms,
                                            n.cores = num_cores,
                                            random.seed = 27599)

saveRDS(object = list(so = so,
                      sop = sop,
                      sov = sov),
        file = paste0('Bailey_scans_', num_perms, '_permsc.RDS'))


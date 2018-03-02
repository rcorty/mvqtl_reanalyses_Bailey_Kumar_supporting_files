# need to make files smaller to put on github

library(vqtl)

bailey_replication <- readRDS('Bailey_scans_1000_perms.RDS')

sovs <- bailey_replication$sov

remove_cross_and_genoprobs_from_meta <- function(sov) {
  tiny_cross <- qtl::sim.cross(map = qtl::sim.map(len = 0))
  sov$meta$cross <- tiny_cross
  sov$meta$genoprob.df <- NULL
  return(sov)
}

smaller_sovs <- lapply(X = sovs, FUN = remove_cross_and_genoprobs_from_meta)

plot(smaller_sovs$PCTT10)

saveRDS(object = list(so = bailey_replication$so,
                      sop = bailey_replication$sop,
                      sov = smaller_sovs),
        file = 'Bailey_scans_1000_perms_smaller_files.RDS')

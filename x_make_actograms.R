library(tidyverse)
library(actogrammr)



x <- read_clock_lab_files(file_names = list.files(path = 'actogram_data', full.names = TRUE))

x %>% glimpse()

mice <- unique(x$file_name)
short_mice <- sapply(X = str_split(string = mice, pattern = '/'), FUN = function(v) v[2], simplify = TRUE)
shorter_mice <- sapply(X = str_split(string = short_mice, pattern = '-'), FUN = function(v) v[1], simplify = TRUE)
shoterterer_mice <- gsub(pattern = '(M|F)$', replacement = '', x = shorter_mice)

# figure out which mice are which
library(qtl)

MPD_URL <- 'https://phenomedoc.jax.org/QTL_Archive/kumar_2014/Kumar2014_Data.csv'
c <- read.cross(format = 'csv',  file = url(description = MPD_URL))

c$pheno$qtl_geno <- c$geno$`6`$data[, 'rs30314218']
sum(c$pheno$NA. %in% shoterterer_mice)

y <- c$pheno[c$pheno$NA. %in% shoterterer_mice, c(8, 25, 26, 27)] %>% arrange(sex, qtl_geno, Avg_Counts)
match(x = y$NA.[7:12], table = shoterterer_mice)



# the males are
# GT1: m9, m7
# GT2: m10, m12
# GT3: m3, m4

day_to_bin <- function(day, mins_per_bin = 6) {
  (day-1)*24*60/mins_per_bin+1
}

m9 <- x %>%
  filter(file_name == mice[9]) %>%
  bin_data(minutes_per_bin = 6)

m9 %>%
  actotile() +
  scale_fill_continuous(low = 'white', high = 'black', limits = c(0, 700)) +
  labs(fill = 'rotations')

m7 <- x %>%
  filter(file_name == mice[7]) %>%
  bin_data(minutes_per_bin = 6)

m7 %>%
  actotile() +
  # actotile(start_date = m9$date[day_to_bin(3)], end_date = m9$date[day_to_bin(34)]) +
  scale_fill_continuous(low = 'white', high = 'black', limits = c(0, 700)) +
  labs(fill = 'rotations')




m10 <- x %>%
  filter(file_name == mice[10]) %>%
  bin_data(minutes_per_bin = 6)

m10 %>%
  actotile() +
  # actotile(start_date = m9$date[day_to_bin(3)], end_date = m9$date[day_to_bin(34)]) +
  scale_fill_continuous(low = 'white', high = 'black', limits = c(0, 700)) +
  labs(fill = 'rotations')

m12 <- x %>%
  filter(file_name == mice[12]) %>%
  bin_data(minutes_per_bin = 6)

m12 %>%
  actotile() +
  # actotile(start_date = m9$date[day_to_bin(3)], end_date = m9$date[day_to_bin(34)]) +
  scale_fill_continuous(low = 'white', high = 'black', limits = c(0, 700)) +
  labs(fill = 'rotations')





m3 <- x %>%
  filter(file_name == mice[3]) %>%
  bin_data(minutes_per_bin = 6)

m3 %>%
  actotile() +
  # actotile(start_date = m9$date[day_to_bin(3)], end_date = m9$date[day_to_bin(34)]) +
  scale_fill_continuous(low = 'white', high = 'black', limits = c(0, 700)) +
  labs(fill = 'rotations')

m4 <- x %>%
  filter(file_name == mice[4]) %>%
  bin_data(minutes_per_bin = 6)

m4 %>%
  actotile() +
  # actotile(start_date = m9$date[day_to_bin(3)], end_date = m9$date[day_to_bin(34)]) +
  scale_fill_continuous(low = 'white', high = 'black', limits = c(0, 700)) +
  labs(fill = 'rotations')


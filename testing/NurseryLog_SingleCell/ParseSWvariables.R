#' ---
#' title: Parse soil water variables from the Log file
#' author: Haga Chihiro
#' date: 2020.1.14
#' ---

# Load libraries ===============================================================
library(tidyverse)
library(lubridate)


# Define functions =============================================================
ReadLogFile <- function(fname) {
  f <- file(fname, "r")
  while (TRUE) {
    line <- readLines(f, 1)
    if (length(line) == 0) break
    # pase dataset
    data_line <- unlist(str_split(line, pattern = ' - '))
    if (length(data_line) > 1) {
      if (str_detect(data_line[2], pattern = 'currentNurseryLog')) {
        buf <- unlist(str_split(data_line[2], pattern = ':'))[2]
        dat <- data.frame(t(as.numeric(unlist(str_split(buf, pattern = ',')))))
      } else {
        next
      }
      # append lists
      if (exists('dat_df')) {
        dat_df <- bind_rows(dat_df, dat)
      } else {
        dat_df <- dat
      }
    } 
  }
  close(f)
  colnames(dat_df) <- c('year', 'month', paste0('nl', 1:(ncol(dat_df) - 2)))
  return(dat_df)
}



# Main =========================================================================
# Read and parse log file -----------------------------------
root_path <- 'C:/Working/LANDIS-II_repos/Extension-NECN-H-Succession/testing/NurseryLog_SingleCell/' # Modify HERE !!!!!!!!
log_fullpath <- file.path(root_path, 'Landis-Log.txt')
remove(dat_df) # remove old object
dat_df <- ReadLogFile(log_fullpath) # Read data from .txt file


# plot -------------------------------------------------------
plt <- dat_df %>% 
  filter(month == 5) %>%
  tidyr::gather(key = key, val = nl, -year, -month) %>% 
  ggplot(aes(x = year, y = nl, color = key)) +
  geom_line() +
  ylab('Dead wood biomass (gC m-2)') +
  viridis::scale_color_viridis(discrete = T) +
  theme_classic()
plot(plt)

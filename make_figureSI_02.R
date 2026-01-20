rm(list=ls())

library(ggplot2)
library(gridExtra)
library(epitools)
library(lubridate)
library(dplyr)
library(tidyr)
library(ggpubr)
library(broom)

##############################################
# proposed model names
##############################################

#
# 1. linear-scaling value[nn] = value[nn] * 1./initial_reporting_fraction[nn]
# 2. prospective-scaling value[nn] = (1+ratio)/2 where ratio = max_reporting_fraction/initial_reporting_fraction[nn]
# this is what we did during the season
# 3. retrospective-scaling - do a linear fit of current/initial vs initial reporting fraction and use that to 
# create a retrospective time series
# 4. null-model no scaling
#

##############################################
# path to 2024-25 data

current_data_path = "data/PROF_state_hosp.rds"

source("hosp_data_funs.R")

source('utils.R')

# load the current data - this will have the most up-to-date data

current_hosp = readRDS(current_data_path)

# remove some locations and HHS Regions from data

current_hosp = delete_locations(current_hosp, locs = c('AS', 'VI', 'GU','MP', 'USA'))

current_hosp = delete_locations(current_hosp, locs = paste0('Region ', 1:10))

##############################################

current_hosp_data <- current_hosp$conf_inf_admin

current_hosp_dates <- current_hosp$dates

current_perc_hospinf_report <- current_hosp$perc_hospinf_report

# subset to the current season
date_start = as.Date('2024-11-01')
ind <- which(current_hosp_dates >= date_start)

current_hosp_data = current_hosp_data[ind,]
current_hosp_dates = current_hosp_dates[ind]

current_perc_hospinf_report <- current_perc_hospinf_report[ind,]

df_current_hosp = as.data.frame(current_hosp_data)
df_current_hosp$Date = current_hosp_dates

df_current_report = as.data.frame(current_perc_hospinf_report)
df_current_report$Date = current_hosp_dates

tmp1 <- pivot_longer(df_current_hosp, cols = -Date, names_to ='state', values_to = 'flu_adm')
tmp2 <- pivot_longer(df_current_report, cols = -Date, names_to ='state', values_to = 'perc_hosp_rep')

current_data <- tmp1 %>% left_join(tmp2, by = c('Date', 'state'))

current_data <- current_data %>% rename(weekendingdate = Date)

current_data_update_date = as.Date(file.info(current_data_path)$mtime)

current_data$report_date = current_data_update_date

current_data$issue = 'current'

locs = unique(current_data$state)
nloc = length(locs)

##############################################
# find peak admission date for each location
##############################################

loc_max_date <- as.Date(rep(NA, nloc))
for (ii in 1:nloc) {
  tmp = subset(current_data, state == locs[ii])
  ind = which.max(tmp$flu_adm)
  loc_max_date[ii] = tmp$weekendingdate[ind]
}

names(loc_max_date) <- locs

df_vline <- data.frame(
  state = names(loc_max_date),
  vline_date = as.Date(loc_max_date)
)


# List all archived csv files
##############################################

# submission date is Wednesday and so is data update date we use

archive_dir = "data/nhsn_archive/"

archive_files <- list.files(path = archive_dir, pattern = "HHS_weekly-hosp_state__", full.names = TRUE)

# extract dates
archive_dates <- as.Date(sub(".*__(\\d{6})\\d+\\.csv", "\\1", archive_files), format = "%y%m%d")

process_dates <- archive_dates 
process_files <- archive_files 

# locations to remove 
remove_locs = c('AS','GU', 'MP','USA','VI')
# these are the column names we need 
col_names <- c('weekendingdate',  'jurisdiction','totalconfflunewadmperchosprep', 'totalconfflunewadm')
new_col_names <- c('weekendingdate', 'state', 'perc_hosp_rep', 'flu_adm', 'report_date', 'issue','scaling')

for (ifile in 1:length(process_files)) {
  archive_data = read.csv(process_files[ifile])
  archive_data$weekendingdate = as.Date(archive_data$weekendingdate)
  archive_data = subset(archive_data, weekendingdate >= date_start)
  archive_data = subset(archive_data, !jurisdiction %in% remove_locs)
  iorder = order(archive_data$jurisdiction)
  archive_data = archive_data[iorder,]
  last_date = archive_data$weekendingdate[nrow(archive_data)]
  data_small = subset(archive_data, weekendingdate == last_date)
  data_small = data_small[,col_names]
  data_small$report_date = process_dates[ifile]
  data_small$issue = 'initial'
  data_small$scaling = 1.
  colnames(data_small) <- new_col_names
  
  
  current_small = subset(current_data, weekendingdate == last_date)
  current_small$scaling = current_small$flu_adm/data_small$flu_adm
  
  
  new_rows = rbind(data_small, current_small)
  
  
  if (ifile == 1) {
    df = new_rows
  } else {
    df = rbind(df, new_rows)
  }
  
}


df_complete = df

# subset to three data frame

df_cur = subset(df, issue == 'current') 
df_ini  = subset(df, issue == 'initial')

# remove dates/locations where current/initial > scaling_limit and where current % reporting is < initial % reporting
# and where initial % of reporting is less than 30% 
upper_scaling_limit = 3.
lower_scaling_limit = 0.5
lower_initial_perc_limit = 0.3

ind_rmv1 = which(df_cur$scaling > upper_scaling_limit | df_cur$scaling < lower_scaling_limit)

ind_rmv2 = which(df_cur$perc_hosp_rep < df_ini$perc_hosp_rep)

ind_rmv3 = which(df_ini$perc_hosp_rep < lower_initial_perc_limit)

ind_rmv = c(ind_rmv1, ind_rmv2, ind_rmv3)

percent_rmv = round(length(ind_rmv)/nrow(df_ini) * 100.)

cat("Removing ", length(ind_rmv),' data points out of ', nrow(df_ini),' data points ~', percent_rmv,' percent\n\n')

df = rbind(df_cur[-ind_rmv,], df_ini[-ind_rmv, ])

p2 <- ggplot(df, aes(x = weekendingdate, y = perc_hosp_rep, color = issue)) +
  geom_line() + 
  coord_cartesian(ylim = c(NA, 1)) +  # <- replaces ylim()
  labs(y = "Fraction Hospitals Reporting", x = NULL) +
  facet_wrap(~ state, ncol = 5, scales = "free_y") +
  scale_x_date(
    date_breaks = "1 month" ,
    date_labels = "%y-%m-%d"
  ) + 
  theme_minimal(base_size = 12) +
  guides(color = guide_legend(title = NULL)) + 
  theme(
    strip.text = element_text(size = 10),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  )

plot(p2)

ggsave("fig/figure_SI02.pdf", plot = p2, width = 15, height = 11, units = "in")



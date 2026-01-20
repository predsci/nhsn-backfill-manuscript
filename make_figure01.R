rm(list=ls())

library(ggplot2)
library(gridExtra)
library(epitools)
library(lubridate)
library(dplyr)
library(tidyr)
library(ggpubr)
library(broom)
library(grid)

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
new_col_names <- c('weekendingdate', 'state', 'perc_hosp_rep', 'flu_adm', 'report_date', 'issue')

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

  colnames(data_small) <- new_col_names
  
  
  # this is setting the stage for the retrospective linear fit.  We need to know current/initial and we need to
  # know the initial reporting 
  
  current_small = subset(current_data, weekendingdate == last_date)

  new_rows = rbind(data_small, current_small)
  
  # new_rows$scaling[new_rows$scaling == 0] <- 1.
  
  if (ifile == 1) {
    df = new_rows
  } else {
    df = rbind(df, new_rows)
  }
  
}

df_complete = df

# subset to three data frame

df_current = subset(df, issue == 'current') 
df_initial  = subset(df, issue == 'initial')

# remove dates/locations where current/initial > scaling_limit and where current % reporting is < initial % reporting
# and where initial % of reporting is less than 30% 
upper_scaling_limit = 3.
lower_scaling_limit = 0.5
lower_initial_perc_limit = 0.3

ind_rmv1 = which(df_current$scaling > upper_scaling_limit | df_current$scaling < lower_scaling_limit)

ind_rmv2 = which(df_current$perc_hosp_rep < df_initial$perc_hosp_rep)

ind_rmv3 = which(df_initial$perc_hosp_rep < lower_initial_perc_limit)

ind_rmv = c(ind_rmv1, ind_rmv2, ind_rmv3)

percent_rmv = round(length(ind_rmv)/nrow(df_initial) * 100.)

cat("Removing ", length(ind_rmv),' data points out of ', nrow(df_initial),' data points ~', percent_rmv,' percent\n\n')
df_complete = df


df_current = df_current[-ind_rmv,]
df_initial = df_initial[-ind_rmv,]
# re-assemble df

df = rbind(df_current, df_initial)


p0 <- ggplot(df, aes(x = weekendingdate, y = flu_adm, color = issue)) +
  geom_line() +
  facet_wrap(~ state, ncol = 5, scales = "free_y") +
  labs(y = "Hosp. Admin", x = NULL) +
  theme_minimal(base_size = 12) +
  theme(
    strip.text = element_text(size = 10),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "top"
  )

# p0


# subset to the state we want to plot

selected_state = 'AZ'

date_cutoff = c(as.Date("2024-12-22"), as.Date(max(df$weekendingdate)))

df_state = subset(df, state == selected_state)

legend = c('(A)', '(B)')

scale_x_date_list = list()
scale_x_date_list[[1]] =     scale_x_date(
  date_breaks = "1 week" ,
  date_labels = "%y-%m-%d"
)
scale_x_date_list[[2]] =     scale_x_date(
  date_breaks = "1 month" ,
  date_labels = "%y-%m-%d"
)
pl = list()

for (icount in 1:length(date_cutoff)) {

  df_plot = subset(df_state, weekendingdate <= date_cutoff[icount])
  
  # Find the max flu_adm
  max_flu_adm <- max(df_state$flu_adm, na.rm = TRUE)
  
  # Rescale flu_adm
  df_plot <- df_plot %>%
    mutate(flu_adm_scaled = flu_adm / max_flu_adm)
  
  pl[[icount]] <- ggplot(df_plot, aes(x = weekendingdate)) +
    geom_line(aes(y = perc_hosp_rep, color = issue), linetype = 'dashed', linewidth = 1.) +
    geom_line(aes(y = flu_adm_scaled, color = issue), linewidth = 1.) +
    scale_y_continuous(
      name = "Fraction Hospitals Reporting",
      sec.axis = sec_axis(~ . * max_flu_adm, name = "Weekly New Hospital Admission")
    ) +
    scale_x_date_list[[icount]] +
    labs(
      title = legend[icount],
      x = "Week Ending Date",
      color = NULL
    ) +
    theme(legend.position = 'none',
          axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
          axis.text.y = element_text(size = 12),
          axis.title.x= element_text(size = 12),
          axis.title.y= element_text(size = 12)) +
    coord_cartesian(xlim = c(min(df_state$weekendingdate), date_cutoff[icount])) +
    coord_cartesian(ylim = c(0, 1)) 
  
}


grid.arrange(pl[[1]], pl[[2]], ncol = 2)
p_combined <- arrangeGrob(pl[[1]], nullGrob(), pl[[2]], ncol = 3, widths = c(1,0.1,1))
ggsave("fig/figure_01.pdf", plot = p_combined, width = 11, height = 7, units = "in")

## some stats

# Step 1: Merge the two data frames
df_merged <- df_current[,c('weekendingdate', 'state', 'perc_hosp_rep', 'flu_adm')] %>%
  rename(flu_adm_current = flu_adm,
         perc_rep_current = perc_hosp_rep) %>%
  inner_join(
    df_initial[,c('weekendingdate', 'state', 'perc_hosp_rep', 'flu_adm')] %>%
      rename(flu_adm_initial = flu_adm,
             perc_rep_initial = perc_hosp_rep),
    by = c("state", "weekendingdate")
  )

# Step 2: Create comparison flags
df_comparison <- df_merged %>%
  mutate(
    adm_decreased = flu_adm_current < flu_adm_initial,
    rep_change = case_when(
      perc_rep_current > perc_rep_initial ~ "increased",
      perc_rep_current < perc_rep_initial ~ "decreased",
      TRUE ~ "no_change"
    )
  )

# Step 3: Summarize results by state
summary_by_state <- df_comparison %>%
  group_by(state) %>%
  summarise(
    n_total = n(),
    n_adm_decreased = sum(adm_decreased),
    frac_adm_decreased = n_adm_decreased / n_total,
    n_with_rep_increase = sum(adm_decreased & rep_change == "increased"),
    n_with_rep_decrease = sum(adm_decreased & rep_change == "decreased"),
    .groups = "drop"
  )

iorder = order(summary_by_state$n_adm_decreased, decreasing = TRUE)

print(summary_by_state[iorder,], n=52)
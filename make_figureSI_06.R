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
  data_small$issue = 'null-model'
  data_small$scaling = 1.
  colnames(data_small) <- new_col_names
  
 
  max_report_df <- archive_data %>%
    filter(jurisdiction %in% locs) %>%
    group_by(jurisdiction) %>%
    summarize(max_report = max(totalconfflunewadmperchosprep, na.rm = TRUE)) %>%
    arrange(jurisdiction)
  
  # Match order to data_small
  max_report <- max_report_df$max_report[match(data_small$state, max_report_df$jurisdiction)]
  
  data_psi <- data_small
  data_psi$issue = 'prospective-scaling'
  ratio = max_report/data_small$perc_hosp_rep
  data_psi$flu_adm = data_small$flu_adm * (1.+ratio)*0.5
  data_psi$scaling = (1.+ratio)*0.5
  
  # this is setting the stage for the retrospective linear fit.  We need to know current/initial and we need to
  # know the initial reporting 
  
  current_small = subset(current_data, weekendingdate == last_date)
  current_small$scaling = current_small$flu_adm/data_small$flu_adm

  # linear-scaling model if report fraction is x than flu_adm_null = flu_adm * 1/perc_hosp_rep
  
  data_linear <- data_small
  data_linear$issue = 'linear-scaling'
  data_linear$flu_adm = data_small$flu_adm * 1./data_small$perc_hosp_rep
  data_linear$scaling = 1./data_small$perc_hosp_rep
  
  new_rows = rbind(data_small, current_small, data_psi, data_linear)
  
  new_rows$scaling[is.na(new_rows$scaling)] <- 1.
  # new_rows$scaling[new_rows$scaling == 0] <- 1.
  
  if (ifile == 1) {
    df = new_rows
  } else {
    df = rbind(df, new_rows)
  }
  
}

df$scaling[is.infinite(df$scaling)] <- 1. #NA

df_complete = df

# subset to three data frame

df_linear = subset(df, issue == 'linear-scaling') # scaling is 1/initial_perc_hosp_rep
df_current = subset(df, issue == 'current') # scaling is current/initial, it is setting the stage for the linear fit 
df_psi  = subset(df, issue == 'prospective-scaling') # scaling is our model (1+ratio)/2
df_null    = subset(df, issue == 'null-model') # no scaling, scaling == 1

# df_initial - needed at the very end.  It is basically df_psi. For good measure set scaling to 1.  It is never going to be used

df_initial = df_psi
df_initial$scaling = 1.
df_initial$issue = 'initial'



# remove dates/locations where current/initial > scaling_limit and where current % reporting is < initial % reporting
# and where initial % of reporting is less than 30% 
upper_scaling_limit = 3.
lower_scaling_limit = 0.5
lower_initial_perc_limit = 0.3

ind_rmv1 = which(df_current$scaling > upper_scaling_limit | df_current$scaling < lower_scaling_limit)

ind_rmv2 = which(df_current$perc_hosp_rep < df_linear$perc_hosp_rep)

ind_rmv3 = which(df_initial$perc_hosp_rep < lower_initial_perc_limit)

ind_rmv = c(ind_rmv1, ind_rmv2, ind_rmv3)

percent_rmv = round(length(ind_rmv)/nrow(df_linear) * 100.)

cat("Removing ", length(ind_rmv),' data points out of ', nrow(df_linear),' data points ~', percent_rmv,' percent\n\n')
df_complete = df

df_linear = df_linear[-ind_rmv,]
df_current = df_current[-ind_rmv,]
df_psi = df_psi[-ind_rmv,]
df_null   = df_null[-ind_rmv,]
df_initial = df_initial[-ind_rmv,]
# re-assemble df

df = rbind(df_linear, df_current, df_psi, df_null)

# plot the weekly data scaling factor vs initial reporting fraction and do a linear fit

df_lm_fit = df_current[,c('weekendingdate', 'state', 'flu_adm', 'scaling')]
df_lm_fit <- df_lm_fit %>% left_join(df_null[,c('weekendingdate', 'state','perc_hosp_rep')], by = c('weekendingdate', 'state'))


lm_coeffs <- df_lm_fit %>%
  group_by(state) %>%
  do(tidy(lm(scaling ~ perc_hosp_rep, data = .))) %>%
  filter(term %in% c("perc_hosp_rep", "(Intercept)") ) %>%
  select(state, term, estimate) %>%
  pivot_wider(names_from = term, values_from = estimate)

# Rename columns for clarity
colnames(lm_coeffs) <- c("state", "intercept", "slope")

cor_scaling_hosp_rep <- df_lm_fit %>%
  group_by(state) %>%
  summarise(cor = cor(scaling, perc_hosp_rep), .groups = "drop")
df_cor = as.data.frame(cor_scaling_hosp_rep)

df_cor$label <- sprintf("r = %.2f", df_cor$cor)



# df_fit: will hold the results of doing the linear fit 

df_fit <- df_null %>%
  left_join(lm_coeffs, by = "state") %>%
  mutate(
    scaled_flu_adm = (perc_hosp_rep * slope + intercept) * flu_adm
  )

cols_to_remove <- c("flu_adm","intercept","slope")
df_fit <- df_fit %>% select(-all_of(cols_to_remove))
df_fit$issue <- 'retrospective-scaling'
df_fit <- df_fit %>% rename(flu_adm = scaled_flu_adm)

# and now we will add it to df 

df_w_fit <- rbind(df, df_fit)

# start calculating the errors and relative errors of the different procedures

df_err_null <- df_current %>%
  select(weekendingdate, state, flu_adm) %>%
  left_join(df_null %>% select(weekendingdate, state, flu_adm), 
            by = c("weekendingdate", "state"), 
            suffix = c("_current", "_null")) %>%
  mutate(err = (flu_adm_current - flu_adm_null),
         type = "null-model") %>%
  mutate(rel_err = (flu_adm_current - flu_adm_null)/flu_adm_current,
         type = "null-model") %>%
  select(weekendingdate, state, err, rel_err, type)

df_err_psi_model <- df_current %>%
  select(weekendingdate, state, flu_adm) %>%
  left_join(df_psi %>% select(weekendingdate, state, flu_adm), 
            by = c("weekendingdate", "state"), 
            suffix = c("_current", "_psi_model")) %>%
  mutate(err = (flu_adm_current - flu_adm_psi_model),
         type = "prospective-scaling") %>%
  mutate(rel_err = (flu_adm_current - flu_adm_psi_model)/flu_adm_current,
         type = "prospective-scaling") %>%
  select(weekendingdate, state, err, rel_err, type)

df_err_retrospective_model <- df_current %>%
  select(weekendingdate, state, flu_adm) %>%
  left_join(df_fit %>% select(weekendingdate, state, flu_adm), 
            by = c("weekendingdate", "state"), 
            suffix = c("_current", "_retrospective_model")) %>%
  mutate(err = (flu_adm_current - flu_adm_retrospective_model),
         type = "retrospective-scaling") %>%
  mutate(rel_err = (flu_adm_current - flu_adm_retrospective_model)/flu_adm_current,
         type = "retrospective-scaling") %>%
  select(weekendingdate, state, err, rel_err, type)

df_err_linear_model <- df_current %>%
  select(weekendingdate, state, flu_adm) %>%
  left_join(df_linear %>% select(weekendingdate, state, flu_adm), 
            by = c("weekendingdate", "state"), 
            suffix = c("_current", "_linear_model")) %>%
  mutate(err = (flu_adm_current - flu_adm_linear_model),
         type = "linear-scaling") %>%
  mutate(rel_err = (flu_adm_current - flu_adm_linear_model)/flu_adm_current,
         type = "linear-scaling") %>% 
  select(weekendingdate, state, err, rel_err, type)

# Combine
df_error <- bind_rows(df_err_null, df_err_psi_model, df_err_retrospective_model, df_err_linear_model)

# clean

big =100

df_error$err[is.nan(df_error$err)] <- 0.
df_error$rel_err[is.nan(df_error$rel_err)] <- 0.

df_error$err[is.infinite(df_error$err)] <- big
df_error$rel_err[is.infinite(df_error$rel_err)] <- big

# plot relative error 

p6 <- ggplot(df_error, aes(x = weekendingdate, y = rel_err, fill = type)) +
  geom_col(position='dodge') +
  geom_vline(data = df_vline, aes(xintercept = vline_date), 
             linetype = "dashed", color = "black", linewidth = 0.5) +
  scale_x_date(
    date_breaks = "1 month" ,
    date_labels = "%y-%m-%d"
  ) +
  scale_fill_manual(values = c(
    "null-model" = "red",
    "linear-scaling" = "orange",
    "prospective-scaling" = "#00BFC4",
    "retrospective-scaling" = "#C77CFF"
  )) +  
  guides(color = guide_legend(title = NULL)) +
  facet_wrap(~ state, ncol = 5, scales = 'free_y') +
  labs(title = "",
       x = "Week Ending Date",
       y = "Relative Error: (Final Reported Value-Model Estimate)/Final Reported Value") +
  theme_minimal(base_size = 12) +
  theme(
    strip.text = element_text(size = 10),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "top"
  )

plot(p6)

ggsave("fig/figure_SI06.pdf", plot = p6, width = 15, height = 11, units = "in")


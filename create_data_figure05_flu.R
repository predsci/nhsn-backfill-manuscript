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
# path to current, 2024-25 data

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

# linear fit
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
         type = "null-model") 

df_err_psi_model <- df_current %>%
  select(weekendingdate, state, flu_adm) %>%
  left_join(df_psi %>% select(weekendingdate, state, flu_adm), 
            by = c("weekendingdate", "state"), 
            suffix = c("_current", "_psi_model")) %>%
  mutate(err = (flu_adm_current - flu_adm_psi_model),
         type = "prospective-scaling") 

df_err_retrospective_model <- df_current %>%
  select(weekendingdate, state, flu_adm) %>%
  left_join(df_fit %>% select(weekendingdate, state, flu_adm), 
            by = c("weekendingdate", "state"), 
            suffix = c("_current", "_retrospective_model")) %>%
  mutate(err = (flu_adm_current - flu_adm_retrospective_model),
         type = "retrospective-scaling") 

df_err_linear_model <- df_current %>%
  select(weekendingdate, state, flu_adm) %>%
  left_join(df_linear %>% select(weekendingdate, state, flu_adm), 
            by = c("weekendingdate", "state"), 
            suffix = c("_current", "_linear_model")) %>%
  mutate(err = (flu_adm_current - flu_adm_linear_model),
         type = "linear-scaling")

# Combine
df_error <- bind_rows(df_err_null, df_err_psi_model, df_err_retrospective_model, df_err_linear_model)

# clean

big =100

df_error$err[is.nan(df_error$err)] <- 0.

df_error$err[is.infinite(df_error$err)] <- big

# now relative to the null error: measure the error relative to the NULL-Error

df_err_psi_model_over_null <- df_err_psi_model %>%
  select(weekendingdate, state, err) %>%
  left_join(df_err_null %>% select(weekendingdate, state, err), 
            by = c("weekendingdate", "state"), 
            suffix = c("_psi_scaling", "_null")) %>%
  mutate(err = (err_psi_scaling/err_null),
         type = "prospective-scaling") %>%
  select(weekendingdate, state, err, type)

df_err_linear_model_over_null <- df_err_linear_model %>%
  select(weekendingdate, state, err) %>%
  left_join(df_err_null %>% select(weekendingdate, state, err), 
            by = c("weekendingdate", "state"), 
            suffix = c("_linear_model", "_null")) %>%
  mutate(err = (err_linear_model/err_null),
         type = "linear-scaling") %>%
  select(weekendingdate, state, err, type)

df_err_retrospective_model_over_null <- df_err_retrospective_model %>%
  select(weekendingdate, state, err) %>%
  left_join(df_err_null %>% select(weekendingdate, state, err), 
            by = c("weekendingdate", "state"), 
            suffix = c("_retrospective_model", "_null")) %>%
  mutate(err = (err_retrospective_model/err_null),
         type = "retrospective-scaling") %>%
  select(weekendingdate, state, err, type)

# bind them together

df_error_over_null  <- bind_rows(df_err_psi_model_over_null, df_err_linear_model_over_null,df_err_retrospective_model_over_null)

df_error_over_null$err[is.nan(df_error_over_null$err)] <- 0.
df_error_over_null$err[is.infinite(df_error_over_null$err)] <- big

# for each model fraction of weeks for each location that it does better than the NULL model

df_psi_model_frac <- df_err_psi_model_over_null %>%
  group_by(state) %>%
  summarise(frac = length(which(abs(err) < 1))/length(err), .groups = "drop") %>%
  mutate(type = 'prospective-scaling')

df_linear_model_frac <- df_err_linear_model_over_null %>%
  group_by(state) %>%
  summarise(frac = length(which(abs(err) < 1))/length(err), .groups = "drop") %>%
  mutate(type = 'linear-scaling')

df_retrospective_model_frac <-df_err_retrospective_model_over_null %>%
  group_by(state) %>%
  summarise(frac = length(which(abs(err) < 1))/length(err), .groups = "drop") %>%
  mutate(type = 'retrospective-scaling')

df_frac <- rbind(df_psi_model_frac, df_linear_model_frac, df_retrospective_model_frac)
df_frac$state <- factor(df_frac$state)
df_frac$type <- factor(df_frac$type)

p7 <- ggplot(df_frac, aes(x = type, y = state, fill = frac)) +
  geom_tile() +
  geom_text(aes(label = round(frac, 2)), color = "black", size = 3) +  # value on each tile
  scale_fill_gradient2(high ="#6baed6", low = "#d6936b", midpoint = 0.5, na.value = "grey50", 
                       name = "Fraction", 
                       breaks = c(0, 0.5, 1), 
                       labels = c(0, 0.5, 1)) + 
  
  labs(
    title = "(A) Influenza",
    x = "Model Type",
    y = "Location",
    fill = "Fraction"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "top",
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank()
  )

p7

# save to an RDS file

df_flu = df_frac

filename = 'figure_05_flu_data.rds'

saveRDS(file = filename, df_flu)

# 
# ## Book keeping for saving all the results
# # ----------------------------------------------------------
# 
# # Get all objects in the global environment
# all_objects <- ls(envir = .GlobalEnv)
# 
# # Retrieve actual objects
# obj_list <- mget(all_objects, envir = .GlobalEnv)
# 
# # Helper to get the main class of each object
# get_main_class <- function(x) class(x)[1]
# 
# # Categorize by class
# by_type <- split(names(obj_list), sapply(obj_list, get_main_class))
# 
# # Example: access specific types
# dataframes <- by_type[["data.frame"]]     # Includes tibbles
# lists       <- by_type[["list"]]
# ggplots     <- by_type[["gg"]]            # ggplot objects usually have class "gg"
# matrices    <- by_type[["matrix"]]
# 
# # View available types
# names(by_type)
# 
# # save by types - not everything needs to be saved
# 
# ggplot_objects <- mget(by_type[["gg"]], envir = .GlobalEnv)
# data.frame_objects <- mget(by_type[["data.frame"]], envir = .GlobalEnv)
# tbl_df_objects <- mget(by_type[["tbl_df"]], envir = .GlobalEnv)
# grouped_df_objects <- mget(by_type[["grouped_df"]], envir = .GlobalEnv)
# list_objects <- mget(by_type[["list"]], envir = .GlobalEnv)
# 
# list_by_types <- list(data.frame_objects = data.frame_objects, tbl_df_objects = tbl_df_objects,
#                       grouped_df_objects = grouped_df_objects, list_objects = list_objects)
# 
# filename = paste0('flu_eval_scaling-', Sys.Date(),'.rds')
# saveRDS(file = filename, list_by_types)


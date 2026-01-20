# util.R holds functions that will be used for the 2024-25 Influenza Forecasting Season


# accept daily hospitalization data, remove locations, remove data prior to December 2021,
# convert to weekly cadence and return

make_weekly_hosp_state <- function(state_hosp = state_hosp, locs_rmv = c('AS', 'VI')) {

  state_hosp = delete_locations(state_hosp, locs=locs_rmv)

  # get a list of state abbvs
  all_states = state_hosp$states$state

  # I kept all 4 influenza-specific columns, but the one we want for the
  # FluSight challenge is 'conf_inf_admin' confirmed influenza admissions
  #str(state_hosp$conf_inf_admin)

  # to be consistent with the DICE format, this is a data.frame of incidence
  # with each column corresponding to a location. The dates for this incidence
  # are a separate vector state_hosp$dates

  # set start date for the season

  # remove data before December 2021


  #start_date = as.Date("2020-12-04")

  start_EY = 2020
  start_EW = 49
  start_date = mmwr_week_to_date(start_EY, start_EW, day= 7)

  min_date_ind = which(state_hosp$dates==start_date)[1]
  keep_ind = min_date_ind:length(state_hosp$dates)
  state_hosp$dates = state_hosp$dates[keep_ind]
  state_hosp$conf_inf_admin = state_hosp$conf_inf_admin[keep_ind, ]
  state_hosp$conf_inf_admin_cvrg = state_hosp$conf_inf_admin_cvrg[keep_ind, ]
  state_hosp$inf_death = state_hosp$inf_death[keep_ind, ]
  state_hosp$inf_death_cvrg = state_hosp$inf_death_cvrg[keep_ind, ]

  dates = state_hosp$dates
  ntimes = length(dates)

  # this is needed only for the statistical model - we want the week-to-week changes from 2021-12-04

  state_hosp_hist = daily_to_weekly(state_hosp)

  return(state_hosp_hist )

}

# read ILI+ data

read_ilip <- function(filename="~/Dropbox/CSJT01/smh_inf_rd5/data/ILI-plus_timeseries.rds") {

  ili_plus = readRDS(file = filename)

  # First column is nation and we will replace is with a column of NAs for Puerto Rico
  # it is not in the ILI+ data set but we need it for the hospitalization

  ili_plus[,1] = rep(NA, nrow(ili_plus))
  colnames(ili_plus)[1] <- 'Puerto Rico'

  # map locs name to locs abbreviation

  state_names = colnames(ili_plus)

  state_abb <- state.abb[match(state_names, state.name)]
  state_abb[1] <- 'PR'

  state_abb[is.na(state_abb)] <- 'DC'

  # reorder

  iorder = order(state_abb)

  ili_plus = ili_plus[,iorder]

  state_ili = state_abb[iorder]

  colnames(ili_plus) <- state_ili
  return(ili_plus)
}

# replace any identically zero values with interpolated values assuming that
# an identical zero means no reporting
# this is done only for dates in 2024

clean_ilip <- function(ili_plus = NULL) {

  nloc = ncol(ili_plus)
  ndates = nrow(ili_plus)
#
#   ili_EYEW = as.numeric(rownames(ili_plus))
#   ili_EY <-as.numeric(substr(ili_EYEW, 1, 4))
#   ili_EW <- as.numeric(substr(ili_EYEW, 5, 6))
#   ind_2024 <- which(ili_EY == 2024)
  for (iloc in 1:nloc) {
    my_ilip = as.numeric(ili_plus[,iloc])
    ind = which(my_ilip == 0.0 ) #& (1:ndates %in% ind_2024))
    my_ilip[ind] <- NA
    my_ilip <- na.approx(my_ilip, na.rm = FALSE)
    ili_plus[,iloc] <- my_ilip
    # if (my_ilip[ndates] == 0.0) my_ilip[,iloc] =
  }

  return(ili_plus)
}

replace_last_na_ilip <- function (ili_plus) {
  nloc = ncol(ili_plus)
  for (iloc in 1:nloc) {

    ts_data = as.numeric(ili_plus[,iloc])
    if(all(is.na(ts_data))) next

    # Indices of non-NA values
    non_na_indices <- which(!is.na(ts_data))

    last_three_indices <- tail(non_na_indices, 6)

    # Corresponding values
    non_na_values <- ts_data[non_na_indices]


    last_three_values <- tail(non_na_values, 6)
    # Fit linear model
    model <- lm(last_three_values ~ last_three_indices)

    # Index of the missing value
    missing_index <- which(is.na(ts_data))

    # Predict the missing value

    # Index of the missing value
    missing_index <- which(is.na(ts_data))

    # Predict the missing value
    predicted_value <- predict(model, newdata = data.frame(last_three_indices = missing_index))

    # Replace NA with predicted value
    ts_data[missing_index] <- predicted_value


    ili_plus[,iloc] <- ts_data
  }

    return(ili_plus)
}
# function to impute hospitalization data using ILI+ data

impute_using_ilip <- function(ili_plus = NULL, state_hosp = NULL, exception_list = NULL, zero_intercept = NULL) {

  # get a list of state abbvs in hosp data
  all_states = state_hosp$states$state
  locs = all_states
  nloc = length(locs)

  # set start date for data
  # start_date = as.Date("2021-12-04")

  start_EY = 2021
  start_EW = 48

  start_date = mmwr_week_to_date(start_EY, start_EW, day= 7)

  min_date_ind = which(state_hosp$dates==start_date)[1]
  keep_ind = min_date_ind:length(state_hosp$dates)
  state_hosp$dates = state_hosp$dates[keep_ind]
  state_hosp$conf_inf_admin = state_hosp$conf_inf_admin[keep_ind, ]
  state_hosp$conf_inf_admin_cvrg = state_hosp$conf_inf_admin_cvrg[keep_ind, ]
  state_hosp$inf_death = state_hosp$inf_death[keep_ind, ]
  state_hosp$inf_death_cvrg = state_hosp$inf_death_cvrg[keep_ind, ]

  hosp_dates = state_hosp$dates

  ntimes = length(hosp_dates)

  # the ILI+ data: get EY and EW and use that to get date

  ili_EYEW = as.numeric(rownames(ili_plus))
  ili_EY <-as.numeric(substr(ili_EYEW, 1, 4))
  ili_EW <- as.numeric(substr(ili_EYEW, 5, 6))

  # day 7 is the saturday - end of the EW
  ili_dates = mmwr_week_to_date(year = ili_EY, week = ili_EW, day = 7)

  min_date_ind = which(ili_dates >= start_date)[1]

  keep_ind = min_date_ind:length(ili_EYEW)

  ili_plus_long = ili_plus

  ntimes_iliplus_long = nrow(ili_plus_long)

  ili_copy_EYEW = as.numeric(rownames(ili_plus_long))
  ili_copy_EY <-as.numeric(substr(ili_copy_EYEW, 1, 4))
  ili_copy_EW <- as.numeric(substr(ili_copy_EYEW, 5, 6))

  ili_plus = ili_plus[keep_ind,]
  ili_dates = ili_dates[keep_ind]
  ili_EYEW = ili_EYEW[keep_ind]
  ili_EY = ili_EY[keep_ind]
  ili_EW = ili_EW[keep_ind]


  # indices of ILI+ dates that overlap hospitalization dates

  indices <- which(ili_dates %in% hosp_dates)

  proxy_hosp_inc = array(0, c(length(ili_dates), nloc))
  proxy_hosp_inc_long = ili_plus_long

  colnames(proxy_hosp_inc) = locs

  df = data.frame(iloc = 1:nloc, hosp_loc = locs, jloc = 1:nloc, ilip_loc = locs, r_sq = rep(-1, nloc))

  hosp_data = state_hosp$conf_inf_admin

  # loop over all locations and find the scaling between ILI+ and hospitalization
  for (iloc in 1:nloc) {
    state_abb = locs[iloc]

    ili_inc = ili_plus[indices,iloc]
    hosp_inc = hosp_data[,state_abb]
    x = hosp_inc

    # if (any(is.na(ili_inc))) {
      if (any(is.na(ili_plus[,iloc])) || state_abb %in% exception_list) {
      adj_r_squared_model_nn = -1.
      iloc_best = NA
      for (jloc in 1:nloc) {
        if (jloc == iloc) next
        if (locs[jloc] %in% exception_list) next
        y = ili_plus[indices,jloc]
        if (any(is.na(y))) next
        if (zero_intercept) {
          model <- lm(x ~ 0 + y)
        } else {
          model <- lm(x ~ y)
        }
        adj_r_squared_model <- summary(model)$adj.r.squared
        if (adj_r_squared_model > adj_r_squared_model_nn) {
          adj_r_squared_model_nn = adj_r_squared_model
          iloc_best = jloc
        }
      }
      df$jloc[iloc] = iloc_best
      df$ilip_loc[iloc] = locs[iloc_best]
      df$r_sq[iloc] = adj_r_squared_model_nn
      ili_inc = ili_plus[indices,iloc_best]

      if (zero_intercept) {
        model <- lm(hosp_inc ~ 0 + ili_inc)
        intercept <- 0 #as.numeric(coef(model)[1])
        slope <- as.numeric(coef(model)[1])
      } else {
        model <- lm(hosp_inc ~ ili_inc)
        intercept <- as.numeric(coef(model)[1])
        slope <- as.numeric(coef(model)[2])
      }

      proxy_hosp_inc[,iloc] <- as.numeric(ili_plus[,iloc_best]) * slope + intercept

    } else {
      if (zero_intercept) {
        model <- lm(hosp_inc ~ 0 + ili_inc)
        intercept <- 0
        slope <- as.numeric(coef(model)[1])
      } else {
        model <- lm(hosp_inc ~ ili_inc)
        intercept <- as.numeric(coef(model)[1])
        slope <- as.numeric(coef(model)[2])
      }

      df$r_sq[iloc] = summary(model)$adj.r.squared

      proxy_hosp_inc[,iloc] <- as.numeric(ili_plus[,iloc]) * slope + intercept
    }
  }

  proxy_hosp_inc_long[,iloc] <- as.numeric(proxy_hosp_inc_long[,iloc]) * slope + intercept


  df_proxy = as.data.frame(proxy_hosp_inc)
  df_proxy$Date = ili_dates
  # Melt the data into long format
  df_long_proxy <- reshape2::melt(df_proxy, id.vars = "Date", variable.name = "State", value.name = "Value")

  df_long_proxy$Dataset = 'ILI+'

  df_hosp = as.data.frame(hosp_data)
  df_hosp$Date = hosp_dates

  # Melt the data into long format
  df_long_hosp <- reshape2::melt(df_hosp, id.vars = "Date", variable.name = "State", value.name = "Value")
  df_long_hosp$Dataset = 'Hosp.'

  df_combined = bind_rows(df_long_proxy, df_long_hosp)

  p <- ggplot(df_combined, aes(x = Date, y = Value, color=Dataset)) +
    geom_line() +
    facet_wrap(~ State, scales = "free_y") +  # Facet by state, scales can be "free_y" or "fixed"
    theme_minimal() +
    labs(title = "Time Series by State", x = "Date", y = "Value")

  return(list(dates = ili_dates, results = proxy_hosp_inc, df = df, p = p, df_combined = df_combined, proxy_hosp_inc_long = proxy_hosp_inc_long))
}



# interpolate any possible gaps between the imputed hospitalization data and the current, 2024-25
# weekly hospitalization data


interp_imputed_data <- function(impute_date = NULL, impute_data = NULL, current_ts = NULL) {

  all_states = colnames(impute_data)
  nloc = length(all_states)


  for (iloc in 1:nloc) {
    state_abb = all_states[iloc]


    my_impute_data = impute_data[,state_abb]

    impute_date_last = impute_date[length(impute_date)]
    interp_start = impute_date_last + 7
    interp_end = interp_start + 7# this might be something else
    # current date - for now make up to be 7 days ahead
    interp_date = seq(from = interp_start, to = interp_end, by = 7)
    n_inter= length(interp_date)
    date_all = c(impute_date, interp_date)

    #hosp_date and inc what is the fhospitalization date and incidence

    ts1 <- data.frame(date = impute_date, value = my_impute_data)

    # This will be replaced by dates and data from the 2024-25 season

    start_ts2 <- as.Date("2024-10-05") + n * 7 # this will be replaced by the hosp start date
    dates_ts2 <- seq(start_ts2, by = "week", length.out = 10)  # this will be the hospitalization dates
    values_ts2 <- runif(length(dates_ts2))  # Random values for ts2 - to be replaced by hospital incidence
    ts2 <- data.frame(date = dates_ts2, value = values_ts2)

    # Step 1: Create a sequence of dates that includes the gap between ts1 and ts2
    all_dates <- seq(min(ts1$date), max(ts2$date), by = "week")

    # Step 2: Create a combined time series with missing values (NA) in the gap
    combined_ts <- merge(ts1, ts2, by = "date", all = TRUE)
    combined_ts$value <- ifelse(is.na(combined_ts$value.x), combined_ts$value.y, combined_ts$value.x)
    combined_ts <- combined_ts[, c("date", "value")]

    # Step 3: Interpolate the missing values in the gap using the zoo package
    combined_ts$value <- na.approx(combined_ts$value, x = combined_ts$date, na.rm = FALSE)

    # step 4  round the incidence
    combined_ts$value  = round(combined_ts$value)

    if (iloc == 1) {
      interp_ts = array(0, c(length(combined_ts$value), nloc))
      colnames(interp_ts) = all_states
    }
    interp_ts[,state_abb] = combined_ts$value


  }

  return(list(dates = all_dates, interp_ts = interp_ts))

}

# process the vaccination data

process_vax_data <- function(scen_vax = NULL) {

  vax_states = unique(scen_vax$Geography)
  vax_states = vax_states[vax_states!="United States"]

  nstates = length(vax_states)

  vax_abbv = state.abb[match(vax_states, state.name)]

  vax_abbv[nstates] = 'DC'

  # The order of the states in the vaccination time series is NOT identical to that
  # in the HHS data and two are missing! We will use the HHS order, we will see what
  # to do about the missing locations

  # for each location

  age_names = c('0-17','18-49','50-64','65+')

  age_label = age_names

  vax_rate_list = list()
  for (iloc in 1:length(vax_states)) {

    # we can use one state and one scenario to get the array of vax dates since it is the same for all
    state = vax_states[iloc]
    state_abbv = vax_abbv[iloc]
    scenario = "C"

    # mbn: we will use this or some other source
    out_list = round5_data(scen_vax=scen_vax, scenario=scenario, state=state)

    vax_pc_by_age = out_list$vax_pc[,age_label]

    state_pop_ages = out_list$state_pop_ages
    state_pop = state_pop_ages$Population
    state_pop_pc = state_pop / sum(state_pop)
    week_dates = out_list$vax_pc$wk_end
    my_vax_pc = rep(0, length(week_dates))
    for (itime in 1:length(week_dates)) {
      my_vax_pc[itime] = sum(vax_pc_by_age[itime,] * state_pop_pc)
    }
    vax_rate = my_vax_pc/100.0

    vax_df = data.frame(date = week_dates, value = vax_rate)
    # if we want to subset
    #vax_df = subset(vax_df, date %in% state_hosp$dates)
    vax_rate_list[[state_abbv]] = vax_df

  }

  # Now create a vax time-series for PR using the Florida one

  pr_vax_rate = vax_rate_list[['FL']]
  pr_vax_rate[,'value'] = pr_vax_rate[,'value'] * 0.6

  vax_rate_list[['PR']] = pr_vax_rate
  vax_abbv[(nstates + 1)] = 'PR'

  nstates = length(vax_rate_list)

  return(list(vax_rate_list = vax_rate_list, vax_abbv = vax_abbv))

}

# plot 2024-25 data

trim_2024_data <- function(state_hosp = NULL) {

  start_EY = 2021
  start_EW = 48

  start_date = mmwr_week_to_date(start_EY, start_EW, day= 7)
  start_date = mmwr_week_to_date(start_EY, start_EW, day= 7)

  min_date_ind = which(state_hosp$dates==start_date)[1]
  keep_ind = min_date_ind:length(state_hosp$dates)
  state_hosp$dates = state_hosp$dates[keep_ind]
  state_hosp$conf_inf_admin = state_hosp$conf_inf_admin[keep_ind, ]
  state_hosp$perc_hospdays_report = state_hosp$perc_hospdays_report[keep_ind, ]
  state_hosp$perc_hospinf_report = state_hosp$perc_hospinf_report[keep_ind, ]

  # this is needed only for the statistical model - we want the week-to-week changes from 2021-12-04

  return(state_hosp)

}

# find off-season dates
off_season_dates <- function(start_EY=2023, end_EY = 2023,start_EW = 21, end_EW = 39) {

  off_season_start = mmwr_week_to_date(start_EY, start_EW, day= 7)
  off_season_end = mmwr_week_to_date(end_EY, end_EW, day= 7)

  return(list(off_season_start = off_season_start, off_season_end = off_season_end))
}

calc_off_season_baseline <- function(data = NULL, off_season = NULL) {

  start_date = off_season$off_season_start
  end_date = off_season$off_season_end

  hosp_data = data$conf_inf_admin
  hosp_dates = data$dates

  keep_ind <- which(hosp_dates >= start_date & hosp_dates <= end_date)

  off_season_data = hosp_data[keep_ind,]
  locs = data$states$state
  nloc = length(locs)

  baseline = rep(0, nloc)
  for (iloc in 1:nloc) {
    baseline[iloc] <- mean(off_season_data[,iloc], na.rm = TRUE)
  }
  names(baseline) <- locs
  return(baseline)
}

calc_impute_off_season_baseline <- function(data = NULL, date = NULL, off_season = NULL) {

  start_date = off_season$off_season_start
  end_date = off_season$off_season_end

  hosp_data = data
  hosp_dates = date

  keep_ind <- which(hosp_dates >= start_date & hosp_dates <= end_date)

  off_season_data = hosp_data[keep_ind,]
  locs = colnames(data)
  nloc = length(locs)

  baseline = rep(0, nloc)
  for (iloc in 1:nloc) {
    baseline[iloc] <- mean(off_season_data[,iloc], na.rm = TRUE)
  }
  names(baseline) <- locs
  return(baseline)
}



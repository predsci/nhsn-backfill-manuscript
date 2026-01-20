
library(lubridate)

# retrieve state and territory population data (requires DICE)
get_pop_data <- function() {
  library(DICE)
  myDB = OpenCon()
  flu_lut = dbReadTable(myDB, name="flu_lut")
  zika_lut = dbReadTable(myDB, name="zika_lut")
  
  state_ind = flu_lut$ABBV_2=="US" & flu_lut$level==4
  
  state_abbvs = flu_lut$ABBV_4[state_ind]
  state_keys = flu_lut$master_key[state_ind]
  state_sedac_pop = flu_lut$sedac_pop[state_ind]
  
  
  state_pop = data.frame(master_key=state_keys, state=state_abbvs, sedac_pop=state_sedac_pop, stringsAsFactors=F)
  # get Virgin Island data from zika_lut
  # PRI_ind = zika_lut$ABBV_2=="PR" & zika_lut$level==2
  # state_pop = rbind(state_pop, data.frame(master_key=zika_lut$master_key[PRI_ind], state="PRI", sedac_pop=zika_lut$sedac_pop[PRI_ind]))
  VIR_ind = zika_lut$ABBV_2=="VI" & zika_lut$level==2
  state_pop = rbind(state_pop, data.frame(master_key=zika_lut$master_key[VIR_ind], state="VI", sedac_pop=zika_lut$sedac_pop[VIR_ind]))
  # get American Samoa data from flu_lut
  AS_ind = flu_lut$ABBV_2=="AS" & flu_lut$level==2
  state_pop = rbind(state_pop, data.frame(master_key=flu_lut$master_key[AS_ind], state="AS", sedac_pop=flu_lut$sedac_pop[AS_ind]))
  
  # query census data
  census_pop = get.census.data(master_keys=state_pop$master_key, myDB=myDB)
  # record most recent census value available
  state_pop$census_pop = NA
  census_keys = unique(census_pop$master_key)
  for (ii in 1:length(census_keys)) {
    key = census_keys[ii]
    max_date_ind = max(which(census_pop$master_key == key))
    key_ind = state_pop$master_key == key
    state_pop$census_pop[key_ind] = census_pop$pop[max_date_ind]
  }
  # where census data is unavailable, use sedac pop
  state_pop$use_pop = state_pop$census_pop
  census_na = is.na(state_pop$census_pop)
  state_pop$use_pop[census_na] = state_pop$sedac_pop[census_na]
  
  dbDisconnect(myDB)
  
  return(state_pop)
}


delete_locations <- function(dice_data, locs) {
  # assumes dice_data is a list with states and dates entries.
  # locs should be a vector that matches the dice_data$states$state column
  list_names = names(dice_data)
  data_names = list_names[!(list_names %in% c('dates', 'states'))]
  loc_index = dice_data$states$state %in% locs
  print(paste0("Removing the following locations from data: ", paste(dice_data$states$state[loc_index], collapse=", ")))
  
  # remove from state list
  dice_data$states = dice_data$states[!loc_index, ]
  # remove from all data metrics
  for (col_name in data_names) {
    dice_data[[col_name]] = dice_data[[col_name]][, !loc_index]
  }
  
  return(dice_data)
}


trim_timeseries <- function(dice_data, min_keep_date) {
  # assumes dice_data is a list with states and dates entries.
  # min_keep_date should be a scalar of type Date
  list_names = names(dice_data)
  data_names = list_names[!(list_names %in% c('dates', 'states'))]
  keep_index = dice_data$dates >= min_keep_date
  
  # remove from date list
  dice_data$dates = dice_data$dates[keep_index]
  # remove from all data metrics
  for (col_name in data_names) {
    dice_data[[col_name]] = dice_data[[col_name]][keep_index, ]
  }
  
  return(dice_data)
}


daily_to_weekly <- function(dice_data) {
  sat_index = which(wday(dice_data$dates)==7)
  # remove first from list
  sat_index = sat_index[-1]
  # get a list of time-series names
  list_names = names(dice_data)
  ts_names = list_names[3:length(list_names)]
  # initialize new data list
  new_dice = list(dates=dice_data$dates[sat_index],
                  states=dice_data$states)
  # loop through time-series entries
  for (list_name in ts_names) {
    # initialize the data frame
    new_dice[[list_name]] = as.data.frame(matrix(NA, nrow=length(new_dice$dates), ncol=nrow(new_dice$states)))
    names(new_dice[[list_name]]) = new_dice$states$state
    # for each week
    for (ii in 1:length(sat_index)) {
      cur_sat = sat_index[ii]
      prev_sat = cur_sat - 7
      # sum to weekly
      new_dice[[list_name]][ii, ] = apply(
        X=dice_data[[list_name]][(prev_sat+1):cur_sat, ], 
        MARGIN=2, FUN=sum, na.rm=TRUE
      )
    }
  }
  return(new_dice)
}


inf_hosp_2_plot_list <- function(cur_inc_data, data_cols=data.frame(
  data_name='conf_inf_admin', sub_name='inc flu hosp'), days_back=180,
  state_fips=NULL) {
  
  # get all locations that appear in data
  locs = cur_inc_data$states$state
  # remove locations that do not appear in state_fips
  locs = locs[locs %in% state_fips$state]
  # match abbreviations with fips codes
  state_match = match(x=locs, table=state_fips$state)
  
  # generate dataframe that maps submission location column with plot name
  locs_df = data.frame(loc_col=state_fips$state_code[state_match], 
                       plot_lab=state_fips$state_name[state_match], 
                       data_col=locs)
  # initialize plot_list
  plot_list = list(locs=locs_df, time_series=list())
  
  # index dates to keep
  max_date = max(cur_inc_data$dates)
  keep_ind = cur_inc_data$dates > (max_date-days_back)
  keep_dates = cur_inc_data$dates[keep_ind]
  
  # initialize national totals
  national_df = data.frame(date=keep_dates)
  for (jj in 1:nrow(data_cols)) {
    national_df[[data_cols$sub_name[jj]]] = 0
  }
  
  for (ii in 1:nrow(locs_df)) {
    data_abbv = locs_df$data_col[ii]
    sub_col = locs_df$loc_col[ii]
    
    timeseries_df = data.frame(date=keep_dates)
    for (jj in 1:nrow(data_cols)) {
      data_vec = cur_inc_data[[data_cols$data_name[jj]]][[data_abbv]][keep_ind]
      timeseries_df[[data_cols$sub_name[jj]]] = data_vec
      na_index = is.na(data_vec)
      national_df[[data_cols$sub_name[jj]]][!na_index] = 
        national_df[[data_cols$sub_name[jj]]][!na_index] + data_vec[!na_index]
    }
    plot_list$time_series[[sub_col]] = timeseries_df
  }
  
  # add national totals
  plot_list$time_series[["US"]] = national_df
  plot_list$locs = rbind(plot_list$locs, data.frame(loc_col="US", 
                                                    plot_lab="United States", 
                                                    data_col="US"))
  
  return(plot_list)
}



rsvnet_2_plot_list <- function(plot_inc_data, data_cols=data.frame(
  data_name='value', sub_name='RSV Hosp'), weeks_back=25,
  state_fips=NULL) {
  
  # get all locations that appear in data
  locs = unique(plot_inc_data$location)
  # match abbreviations with fips codes
  state_match = match(x=locs, table=state_fips$state_code)
  # generate dataframe that maps submission location column with plot name
  locs_df = data.frame(loc_col=state_fips$state_code[state_match], 
                       plot_lab=state_fips$state_name[state_match], 
                       data_col=locs)
  # initialize plot_list
  plot_list = list(locs=locs_df, time_series=list())
  
  # ensure that the data are in ascending-date order
  plot_inc_data = plot_inc_data[order(plot_inc_data$dates), ]
  
  # index dates to keep
  max_date = max(plot_inc_data$dates)
  keep_ind = plot_inc_data$dates > (max_date - weeks_back*7)
  keep_dates = unique(plot_inc_data$dates[keep_ind])
  
  for (ii in 1:nrow(locs_df)) {
    data_abbv = locs_df$data_col[ii]
    sub_col = locs_df$loc_col[ii]
    
    # index location
    loc_index = plot_inc_data$location==data_abbv
    
    timeseries_df = data.frame(date=keep_dates)
    for (jj in 1:nrow(data_cols)) {
      data_vec = plot_inc_data[[data_cols$data_name[jj]]][loc_index & keep_ind]
      timeseries_df[[data_cols$sub_name[jj]]] = data_vec
      na_index = is.na(data_vec)
    }
    plot_list$time_series[[sub_col]] = timeseries_df
  }
  
  return(plot_list)
}











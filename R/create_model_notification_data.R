create_model_notification_data <- function(
        infections_timeseries,
        full_infection_dates,
        observed_infection_dates,
        timevarying_delay_dist_ext,
        timevarying_proportion,
        observed_data,
        valid_mat = 1,
        model_likelihood = 'negative_binomial') {

    n_days_infection <- nrow(infections_timeseries)
    n_jurisdictions <- ncol(infections_timeseries)

    # get indices for subset of infection days that end up in observed/forecast data
    obs_infection_idx <- which(full_infection_dates %in% observed_infection_dates)

    # subset infection timeseries to these indices
    infection_obs <- infections_timeseries[obs_infection_idx,]
    #note this is the number of infection days that get observed in the observed
    #data, it is not the same length as the actual number of observed data,
    #because we can observe (in the data) infections that had taken place before
    #data collection
    n_days_infection_obs <- length(obs_infection_idx)

    convolution_matrices <- lapply(1:n_jurisdictions, function(x)
        get_convolution_matrix(timevarying_delay_dist_ext[, x],
                               n_days_infection_obs))


    # compute expected cases of the same length
    # note not all of these dates would have been observed
    expected_cases <- do.call(
        cbind,
        lapply(1:n_jurisdictions, function(x) {
            convolution_matrices[[x]] %*% infection_obs[, x] *
                timevarying_proportion[, x]
        }))

    # define likelihood
    if (model_likelihood == 'negative_binomial') {

        #observed days in the data itself
        #basically this backs out extra left and right days
        obs_data_idx <- which(observed_infection_dates %in% as.Date(rownames(observed_data)))

        #truncated infection timeseries matching data dates, just for plotting
        infection_match_data <- infection_obs[obs_data_idx,]

        #get idx for forecast days, basicaly just backing out extra left days
        forecast_data_idx <- which(observed_infection_dates > max(as.Date(rownames(observed_data)))
        )
        # #number of forecast days
        # n_days_forecast <- length(forecast_data_idx)

        # #use validity matrix - default is all valid but can use this to nullify
        # #expected mean cases for dates when reporting is stopped
        # #extend the valid mat to include forecast days by repeating the last row
        # #this ensures for states where data collection has stopped, the forecast is also 0

        if (!is.matrix(valid_mat)) {

            valid_mat <- observed_data
            valid_mat[] <- 1
        }

        valid_idx <- which(as.logical(valid_mat))


        # negative binomial parameters - need to change from mean and variance
        # specification to size and prob
        sqrt_inv_size <- normal(0, 0.5,
                                truncation = c(0, Inf),
                                dim = n_jurisdictions)
        sqrt_inv_size <- sweep(greta::zeros(n_days_infection_obs,
                                            n_jurisdictions),
                               2, sqrt_inv_size,
                               FUN = "+")

        size <- 1 / sqrt(sqrt_inv_size)
        prob <- 1 / (1 + expected_cases / size)
        size_obs <- size[obs_data_idx,]
        prob_obs <- prob[obs_data_idx,]

        #proposal neg binom params for a forecast bit
        size_forecast <- size[forecast_data_idx,]
        prob_forecast <- prob[forecast_data_idx,]
        #take some random draws with these for the forecast bit?
        #array-ify the data
        observed_data_array <- greta::as_data(as.numeric(observed_data)[valid_idx])

        #stick in the forecast bit
        # timeseries_data_array_with_forecast <- rbind(timeseries_data_array,
        #                                              forecast_cases)

        #I guess you could make a version of infection time series with forecast
        #bit too but can't imagine that being useful at all
        greta::distribution(observed_data_array) <- greta::negative_binomial(
            size_obs[valid_idx],
            prob_obs[valid_idx])
    }

    greta_arrays <- module(
        infection_obs,
        infection_match_data,
        size,
        prob,
        observed_data_array,
        size_forecast,
        prob_forecast,
        timevarying_proportion
    )
    # infection_model_objects[[length(infection_model_objects) + 1]] <- module(size)
    # infection_model_objects
    return(greta_arrays)
}

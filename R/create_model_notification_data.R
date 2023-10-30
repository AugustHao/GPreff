create_model_notification_data <- function(
        infections_timeseries,
        full_infection_dates,
        observable_infection_dates,
        timevarying_delay_dist_ext,
        timevarying_proportion,
        dow_option = c("none","constant","proportional"),
        observed_data,
        model_likelihood = 'negative_binomial',
        valid_mat = 1,
        dataID,
        case_type_proportion = 1) {

    n_days_infection <- nrow(infections_timeseries)
    n_jurisdictions <- ncol(infections_timeseries)

    # get indices for subset of infection days that end up in observed/forecast data
    observable_infection_idx <- which(full_infection_dates %in% observable_infection_dates)
    # subset infection timeseries to these indices
    infection_observable <- infections_timeseries[observable_infection_idx,]
    # note this is the number of infection days that get observed in the observed
    # data, it is not the same length as the actual number of observed dates in
    # notification series, which is below

    # observed days in the obs data itself
    # basically this backs out extra left and right days
    obs_data_idx <- which(observable_infection_dates %in% as.Date(rownames(observed_data)))

    # update proportion arg with type proportion
    case_type_proportion <- pmax(case_type_proportion,1e-3)
    timevarying_proportion[obs_data_idx,] <- timevarying_proportion[obs_data_idx,] * case_type_proportion

    # implement dow correction if specified
    dow_option <- match.arg(dow_option)

    if (dow_option != "none") {
        # get day of week index
        # favour base approach over lubridate for manual declaration of factor level
        # so that we know which index is which day of the week
        dweek <- weekdays(observable_infection_dates)
        dweek <- as.integer(
            factor(
                dweek,
                levels = c("Monday",
                           "Tuesday",
                           "Wednesday",
                           "Thursday",
                           "Friday",
                           "Saturday",
                           "Sunday")
            )
        )

        # prior for dweek correction

        # national mean
        dow_nat <- greta::normal(0, 1,
                                 dim = c(6,1))
        # state specific
        dow_state <- greta::normal(0, 1,
                                   dim = c(n_jurisdictions,6))
        # add together for mixed effect
        dow <- sweep(dow_state, 2, dow_nat, FUN = "+")

        # if assuming constant dow effect
        if (dow_option == "constant") {
            # inverse multilogit to get dow weights that sum to 1
            dow_weight <- greta::imultilogit(dow)

            # normalise multiplier to average to 1 over a week
            dow_weight <- dow_weight * 7

            # match weight to date by state matrix
            dow_correction <- t(dow_weight)[dweek,]

            # if assuming dow is proportional to infection count
        } else if (dow_option == "proportional") {

            # get unique week index
            week_indx <- lubridate::week(observable_infection_dates) |>
                as.factor() |>
                as.integer()
            week_unqiue <- unique(week_indx)
            n_week_unqiue <- length(week_unqiue)

            # calculate weekly mean infection
            week_mean <- greta::zeros(length(week_unqiue),n_jurisdictions)
            # Error in match.arg(FUN) : 'arg' must be NULL or a character vector
            # prevents using apply, need to look into this further but for now use for
            # loops
            for (s in 1:n_jurisdictions) {

                week_mean[,s] <-  week_mean[,s] + greta::tapply(infection_observable[,s],
                                                                week_indx,
                                                                FUN = "mean")
            }
            # multiply weekly dow with mean infection counts of that week
            dow_by_week <- lapply(week_unqiue,FUN = function(x) {
                sweep(dow,1,t(week_mean[x,]), FUN = "*")
            })
            # inverse multilogit to get dow weights that sum to 1
            dow_weight <- lapply(dow_by_week,FUN = function(x) greta::imultilogit(x))

            # normalise multiplier to average to 1 over a week
            dow_weight <- lapply(dow_weight,FUN = function(x) {x * 7})

            # match weight to date by state matrix

            #first get dweek to be the same structure
            dweek_list <- tapply(dweek,week_indx,FUN = unique)

            #use it to re arrange dow weight
            dow_correction <- greta::zeros(dim = c(1,8))
            for (w in week_unqiue) {

                weekly_dow_correction <-  t(dow_weight[[w]])[dweek_list[[w]],]

                dow_correction <- rbind(dow_correction,
                                        weekly_dow_correction)
            }

            dow_correction <- dow_correction[-1,]
        }
        #in both options we get a dow_correction to multiply proportion with
        timevarying_proportion <- timevarying_proportion * dow_correction
    }


    # build convolution matrix
    n_days_infection_observable <- length(observable_infection_idx)

    convolution_matrices <- lapply(1:n_jurisdictions, function(x)
        get_convolution_matrix(timevarying_delay_dist_ext[, x],
                               n_days_infection_observable))


    # compute expected cases of the same length
    # note not all of these dates would have been observed
    expected_cases <- do.call(
        cbind,
        lapply(1:n_jurisdictions, function(x) {
            convolution_matrices[[x]] %*% infection_observable[, x] *
                timevarying_proportion[, x]
        }))

    # define likelihood
    if (model_likelihood == 'negative_binomial') {

        # #use validity matrix - default is all valid but can use this to nullify
        # #expected mean cases for dates when reporting is stopped
        # #extend the valid mat to include forecast days by repeating the last row
        # #this ensures for states where data collection has stopped, the forecast is also 0

        if (!is.matrix(valid_mat)) {

            valid_mat <- observed_data
            valid_mat[] <- TRUE
            valid_idx <- which(as.logical(valid_mat))
        } else {
            valid_idx <- which(valid_mat,arr.ind = FALSE)
        }

        # negative binomial parameters - need to change from mean and variance
        # specification to size and prob
        sqrt_inv_size <- normal(0, 0.5,
                                truncation = c(0, Inf),
                                dim = n_jurisdictions)
        sqrt_inv_size <- sweep(greta::zeros(n_days_infection_observable,
                                            n_jurisdictions),
                               2, sqrt_inv_size,
                               FUN = "+")

        size <- 1 / sqrt(sqrt_inv_size)
        prob <- 1 / (1 + expected_cases / size)
        size_obs <- size[obs_data_idx,]
        prob_obs <- prob[obs_data_idx,]

        #array-ify the data
        observed_data_array <- greta::as_data(as.numeric(observed_data)[valid_idx])

        greta::distribution(observed_data_array) <- greta::negative_binomial(
            size_obs[valid_idx],
            prob_obs[valid_idx])
    }

    if (dow_option != "none") {
        greta_arrays <- list(
            dow_weight,
            size,
            prob,
            observed_data_array
        )

        names(greta_arrays) <- c(
            paste0(dataID, '_dow_weight'),
            paste0(dataID, '_size'),
            paste0(dataID, '_prob'),
            paste0(dataID, '_observed_data_array')
        )
    } else {
        greta_arrays <- list(
            size,
            prob,
            observed_data_array
        )

        names(greta_arrays) <- c(
            paste0(dataID, '_size'),
            paste0(dataID, '_prob'),
            paste0(dataID, '_observed_data_array')
        )
    }


    return(greta_arrays)
}

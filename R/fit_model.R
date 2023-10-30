fit_model <- function (model,
                       n_chains = 4,
                       max_convergence_tries = 5,
                       warmup = 1000,
                       init_n_samples = 1000,
                       iterations_per_step = 1000) {

    # get stable inits
    init <- generate_valid_inits(model = model,
                                 chains = n_chains,
                                 max_tries = 1000)

    # first pass at model fitting
    draws <- mcmc(
        m,
        sampler = hmc(Lmin = 25, Lmax = 30),
        chains = n_chains,
        warmup = warmup,
        n_samples = init_n_samples,
        initial_values = init,
        one_by_one = TRUE
    )

    # if it did not converge, try extending it a bunch more times
    finished <- check_convergence(draws)
    tries <- 1
    while(!finished & tries < max_convergence_tries) {
        draws <- greta::extra_samples(
            draws,
            iterations_per_step,
            one_by_one = TRUE
        )
        tries <- tries + 1
        finished <- check_convergence(draws)
    }

    # warn if we timed out before converging successfully
    if (tries == max_convergence_tries) {
        warning("sampling did not converge according to benchmarks")
    }

    return(draws)
}

# has the sampler converged to our standards?
check_convergence <- function(draws, max_r_hat = 1.1, min_n_eff = 1000) {
    stats <- get_convergence_stats(draws)
    all(stats$r_hats < max_r_hat) &
        all(stats$n_eff >= min_n_eff)
}

# check convergence
get_convergence_stats <- function(draws) {

    r_hats <- coda::gelman.diag(draws, autoburnin = FALSE,
                                multivariate = FALSE)$psrf[, 1]
    n_eff <- coda::effectiveSize(draws)

    # sometimes n_eff underflows to 0 simply because the values being traced are
    # very small, so remove these (exactly 0 is not possible)
    n_eff <- n_eff[n_eff != 0]

    cat(sprintf("maximum R-hat: %.2f\nminimum n effective: %.2f",
                max(r_hats, na.rm = TRUE),
                min(n_eff)))

    result <- list(r_hats = r_hats,
                   n_eff = n_eff)

    invisible(result)
}


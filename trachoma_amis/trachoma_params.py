params = {
    "N": 2500,
    "av_I_duration": 2,
    "av_ID_duration": 300 / 7,
    "inf_red": 0.45,
    "min_ID": 11,  # Parameters relating to duration of infection period, including ID period
    "av_D_duration": 200 / 7,
    "min_D": 10.1 / 7,  # Parameters relating to duration of disease period
    "dis_red": 0.3,
    "v_1": 1,
    "v_2": 2.6,
    "phi": 1.4,
    "epsilon": 0.5,  # Parameters relating to lambda function- calculating force of infection
    # Parameters relating to MDA
    "MDA_Cov": 0.8,
    "MDA_Eff": 0.85,  # Efficacy of treatment
    "rho": 0.3,
    "nweeks_year": 52,
    "babiesMaxAge": 0.5,  # Note this is years, need to check it converts to weeks later
    "youngChildMaxAge": 9,  # Note this is years, need to check it converts to weeks later
    "olderChildMaxAge": 15,  # Note this is years, need to check it converts to weeks later
    "b1": 1,  # this relates to bacterial load function
    "ep2": 0.114,
    "n_inf_sev": 38,
    "TestSensitivity": 0.96,
    "TestSpecificity": 0.965,
    "SecularTrendIndicator": 0,
    "SecularTrendYearlyBetaDecrease": 0.01,
    "vacc_prob_block_transmission": 0.5,
    "vacc_reduce_bacterial_load": 0,
    "vacc_reduce_duration": 0,
    "vacc_coverage": 0,
    "vacc_waning_length": 52 * 5,
    "importation_rate": 0.000008,
    "importation_reduction_rate": (0.9) ** (1 / 10),
    "infection_risk_shape": 6.4,  # extra parameter needed for infection risk shape. equivalent to k in STH/sch model.
    # Set to a very high number if you want to assume everyone is the same
    "min_importation_rate": 1
    / (
        20 * 52 * 2500
    ),  # some small number for the minimum importation rate. Can be 0 if you want
    "importation_reduction_length": 25,
}  # time in weeks after performing an MDA which we wait before reducing the importation rate


sim_params = {
    "timesim": (52 * 96) + 1,
    "burnin": (52 * 70) - 1,
    "N_MDA": 0,
}  # Set by main script later


demog = {"tau": 0.0004807692, "max_age": 3120, "mean_age": 1040}


# Projection-specific parameters
projection_params = {
    "timesim": (52 * 100) - 1,  # Different from fitting
    "burnin": (52 * 70) - 1,
    "N_MDA": 0,
}

# Projection configuration
projection_config = {
    "species": "trachoma",
    "species_prefix": "Trachoma_",
    "dist_to_use": "Expo",
    "initial_infect_frac": 0.1,
    "num_sims": 200,
    "subset_keys": [
        "IndI",
        "IndD", 
        "No_Inf",
        "T_latent",
        "T_ID",
        "T_D",
        "Ind_latent",
        "Ind_ID_period_base",
        "Ind_D_period_base",
        "bact_load",
        "Age",
        "vaccinated",
        "time_since_vaccinated",
        "treatProbability",
        "MDA_coverage",
        "systematic_non_compliance",
        "ids",
    ]
}

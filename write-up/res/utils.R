summarizing_table <-
  function(lm_simu_rlt,
           p = 1,
           num_methods = 3,
           only_beta = FALSE,
           latex_opt = "HOLD_position") {
    if (only_beta) {
      tbl_names <- c("", "$\\sigma^2$",
                     paste0("$\\", paste("beta", seq(0, p), sep = "_"), "$"))
    } else {
      tbl_names <- c(
        "",
        paste0("$\\", paste("mu", seq_len(p), sep = "_"), "$"),
        paste0("$\\", paste("sigma^2", seq_len(p), sep = "_"), "$"),
        "$\\sigma^2$",
        paste0("$\\", paste("beta", seq(0, p), sep = "_"), "$")
      )
    }
    est_m <- lm_simu_rlt$estimates %>%
      group_by(estimator_type) %>%
      summarise_all(.funs = mean,) %>%
      mutate_if(is.numeric, ~ as.character(round(., digits = 3)))
    est_v <- lm_simu_rlt$estimates %>%
      group_by(estimator_type) %>%
      summarise_all(.funs = var) %>%
      mutate_if(is.numeric, ~ as.character(round(., digits = 3)))
    est_coverage <- lm_simu_rlt$coverage %>%
      group_by(estimator_type) %>%
      summarise_all(.funs = mean) %>%
      mutate_if(is.numeric, ~ as.character(round(., digits = 3)))
    rbind(est_m, est_v, est_coverage) %>%
      kableExtra::kbl(booktabs = TRUE,
                      col.names = tbl_names,
                      escape = FALSE) %>%
      kableExtra::group_rows("Mean", 1, num_methods) %>%
      kableExtra::group_rows("Variance", 1 + num_methods, 2 * num_methods) %>%
      kableExtra::group_rows("Coverage rate", 2 * num_methods + 1, 3 * num_methods) %>%
      kableExtra::kable_styling(
        latex_options = "HOLD_position",
        font_size = 7,
        full_width = TRUE
      )
  }
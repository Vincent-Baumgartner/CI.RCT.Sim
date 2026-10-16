analyse_diabetes_demediation_fast <- function(separate = TRUE, B = 100){
  lf_fit <- utils::getFromNamespace("logistf.fit", "logistf")
  lf_ctrl <- logistf::logistf.control(maxit = 2000, maxstep = 0.5)
  function(condition, dat, fixed_objects = NULL){

    # Global Variables
    visits <- 0:condition$k
    R_cols <- paste0("R", visits) # idx 1 visit 0, idx 2 visit 1, ...
    yc_cols <- paste0("yc", visits[-1]) # idx 1 visit 1, idx 2 visit 2, ...
    ys_cols <- paste0("y", visits) # idx 1 visit 0, idx 2 visit 1, ...
    j_cols <- paste0("j", visits[-length(visits)]) # idx 1 visit 0, ..., idx k visit k - 1
    cols <- c("trt", "age", ys_cols, R_cols, yc_cols)
    M_idx <- function(v) match(v, cols)

    # Convert the logical of receiving rescue at any point in to a longitudinal measurement in wide format
    dat <- dat |> mutate(
      across(
        starts_with("R") & !R0 & !rescue_start,
        ~ ifelse(is.na(.) & rescue_start <= as.numeric(gsub("R", "", cur_column())), 1, .),
        .names = "R{gsub('R', '', .col)}"
      )
    )

    setup_mice <- function(dat, drop){

      pred <- make.predictorMatrix(dat)
      pred[c(R_cols, "m_start"),] <- 0
      pred[, c(drop, R_cols, "m_start")] <- 0

      meth <- make.method(dat)
      meth[R_cols] <- ""
      return(list(pred = pred, meth = meth))
    }

    run_mice <- function(dat, setup){
      mice(dat, m = 5, printFlag = FALSE, predictorMatrix = setup$pred,
           method = setup$meth, ridge = 1e-5, remove.collinear = FALSE)
    }

    if (separate){
      dat1 <- dat |> filter(trt == 1)
      dat0 <- dat |> filter(trt == 0)

      suppressWarnings(dats <- rbind(
        run_mice(dat1, setup_mice(dat1, drop = c("id", "trt"))),
        run_mice(dat0, setup_mice(dat0, drop = c("id", "trt")))
      ))
    }
    else {
      suppressWarnings(dats <- run_mice(dat, setup_mice(dat, drop = "id")))
    }

    calculate_baseline_change <- function(data, yc_cols, ys_cols){
      for(i in seq_along(yc_cols)){
        data[yc_cols[i]] <- data[, ys_cols[i + 1]] - data[, ys_cols[1]]
      }
      return(data)
    }

    completed_data <- lapply(seq_len(dats$m), function(i){
      dat_comp <- calculate_baseline_change(mice::complete(dats, i), yc_cols, ys_cols)
      Rmat <- as.matrix(dat_comp[R_cols])
      rs <- dat_comp$rescue_start

      # NA & visit < rescue_start -> 0
      mask0 <- is.na(Rmat) & outer(rs, visits, `>`)

      # NA & visit >= rescue_start -> 1
      mask1 <- is.na(Rmat) & outer(rs, visits, `<=`)

      Rmat[mask0] <- 0
      Rmat[mask1] <- 1

      dat_comp[R_cols] <- Rmat
      return(as.matrix(dat_comp[cols]))
    })

    analysis <- function(dat_comp, indicator){
      dat_comp <- dat_comp[indicator,]
      n <- nrow(dat_comp)

      J <- matrix(NA_real_, nrow = n, ncol = length(j_cols),
                  dimnames = list(NULL, j_cols))

      J[, j_cols[length(j_cols)]] <- dat_comp[, M_idx(yc_cols[length(yc_cols)])]
      for (k in rev(visits[-c(1, length(visits))])){

        j_new <- j_cols[k]
        j_old <- j_cols[k + 1]
        r_k <- dat_comp[, M_idx(paste0("R", k))]

        if (all(is.na(r_k)) ||
            sum(r_k, na.rm = TRUE) == 0
        ){
          J[, j_new] <- J[, j_old]
          next
        }

        Xps <- cbind(1, dat_comp[, M_idx(c("trt", "age", paste0("y", 0:k)))])
        mod <- tryCatch(
          suppressWarnings(lf_fit(Xps, r_k, control = lf_ctrl)),
          error = function(e) {
            tryCatch(suppressWarnings(glm.fit(Xps, r_k, family = binomial())),
                     error = function(e) NULL)
          }
        )
        if (is.null(mod)) {
          r_pred <- rep(mean(r_k), n)
        } else {
          r_pred <- if (is.null(mod$pi)) mod$fitted.values else mod$pi
        }

        X <- cbind(1, dat_comp[, M_idx(c("trt", "age", "y0",
                                         R_cols[2:(k + 1)],
                                         yc_cols[1:k]))],
                                         r_pred)
        model <- .lm.fit(X, J[, j_old])
        pos <- match(4 + k, model$pivot) # col R_k
        b <- if (pos <= model$rank) model$coefficients[pos] else 0

       J[,j_new] <- J[, j_old] - b * r_k
      }

      final.model <- lm(J[, "j0"] ~ dat_comp[, M_idx("trt")] + dat_comp[, M_idx("y0")])
      cf <- summary(final.model)$coefficients
      c(
        p = cf[2, "Pr(>|t|)"],
        coef = cf[2, "Estimate"],
        se = cf[2, "Std. Error"]
        )
    }

    res  <- lapply(completed_data, analysis)
    test <- lapply(completed_data, boot::boot, statistic = analysis, R = B)

    variance <- function(bootdata){
      B <- bootdata$R
      mu <- mean(bootdata$t[, 2])
      var <- sum((bootdata$t[,2]-mu)^2)/(B-1)
      var
    }

    effs <- unlist(lapply(res, function(x) x["coef"]))
    var  <- unlist(lapply(test, variance))

    end_res <- pool.scalar(effs, var)
    crit <- qt(0.975, df = end_res$df)
    ci <- c(
      end_res$qbar - crit*sqrt(end_res$t),
      end_res$qbar + crit*sqrt(end_res$t)
    )

    t_stat <- end_res$qbar / sqrt(end_res$t)
    p_onesided <- pt(t_stat, df = end_res$df)
    p_pool <- 2 * pt(-abs(t_stat), df = end_res$df)

    list(
      coef = end_res$qbar,
      p = p_onesided,
      se = sqrt(end_res$t),
      ci_lower = ci[1],
      ci_upper = ci[2],
      two_sided = p_pool
    )
  }
}

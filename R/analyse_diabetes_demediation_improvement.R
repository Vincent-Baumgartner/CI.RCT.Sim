analyse_diabetes_demediation_t <- function(separate = TRUE) {
  function(condition, dat, fixed_objects = NULL) {
    # Convert the logical of receiving rescue at any point in to a longitudinal measurement in wide format
    daat <- dat |> mutate(
      across(
        starts_with("R") & !R0 & !rescue_start,
        ~ ifelse(is.na(.) & rescue_start <= as.numeric(gsub("R", "", cur_column())), 1, .),
        .names = "R{gsub('R', '', .col)}"
      )
    )
    Rcols <- grep("^R[0-9]+$", names(daat), value = TRUE)

    if (separate) {
      dat1 <- daat |> filter(trt == 1)
      pred1 <- make.predictorMatrix(dat1)
      pred1[Rcols, ] <- 0
      pred1["m_start", ] <- 0
      pred1[, c("id", "trt", "m_start", Rcols)] <- 0
      meth1 <- make.method(dat1)
      meth1[Rcols] <- ""
      dat0 <- daat |> filter(trt == 0)
      pred0 <- make.predictorMatrix(dat0)
      pred0[Rcols, ] <- 0
      pred0["m_start", ] <- 0
      pred0[, c("id", "trt", "m_start", Rcols)] <- 0
      meth0 <- make.method(dat0)
      meth0[Rcols] <- ""
      suppressWarnings(dats <- rbind(
        mice(dat1, m = 5, printFlag = FALSE, predictorMatrix = pred1, method = meth1, ridge = 1e-5, remove.collinear = FALSE),
        mice(dat0, m = 5, printFlag = FALSE, predictorMatrix = pred0, method = meth0, ridge = 1e-5, remove.collinear = FALSE)
      ))
    } else {
      pred <- make.predictorMatrix(daat)
      pred[Rcols, ] <- 0
      pred["m_start", ] <- 0
      pred[, c("id", "m_start", Rcols)] <- 0
      meth <- make.method(dat)
      meth[Rcols] <- ""
      suppressWarnings(
        dats <- mice(daat, m = 5, printFlag = FALSE, predictorMatrix = pred, method = meth, ridge = 1e-5, remove.collinear = FALSE)
      )
    }

    pre_process_data <- function(data){
      dat_comp <- data |> mutate(
        across(
          matches("^y[0-9]+$") & !y0,
          ~ .x - y0,
          .names = "yc{gsub('y', '', .col)}"
        )
      )

      Rcols <- grep("^R[0-9]+$", names(dat_comp), value = TRUE)
      visits <- as.numeric(sub("R", "", Rcols))

      Rmat <- as.matrix(dat_comp[Rcols])
      rs <- dat_comp$rescue_start

      # NA & visit < rescue_start -> 0
      mask0 <- is.na(Rmat) & outer(rs, visits, `>`)

      # NA & visit >= rescue_start -> 1
      mask1 <- is.na(Rmat) & outer(rs, visits, `<=`)

      Rmat[mask0] <- 0
      Rmat[mask1] <- 1

      dat_comp[Rcols] <- Rmat
      dat_comp[, paste0("j", condition$k - 1)] <- dat_comp[, paste0("yc", condition$k)]
      return(dat_comp)
    }


    analysis <- function(dataa, indicator = seq_len(nrow(dataa))) {

      dat_comp <- dataa[indicator,]

      for (k in 1:(condition$k - 1)) {
        if (all(is.na(dat_comp[[paste0("R", condition$k-k)]])) ||
            sum(dat_comp[[paste0("R", condition$k-k)]], na.rm=TRUE)==0) {
          dat_comp[, paste0("j", condition$k - k - 1)] <- dat_comp[, paste0("j", condition$k - k)]
          next
        }


      rvar <- paste0("R", condition$k-k)

      # Fit a model to predict the probability of receiving rescue medication at visit condition$k - k
      mod <- tryCatch(
        logistf(
          as.formula(paste0("R", condition$k - k, " ~ trt + age + y0", paste0("+ y", 1:(condition$k - k), collapse = " "))),
          data = dat_comp,
          pl = FALSE,
          control = logistf::logistf.control(maxit = 2000, maxstep = 0.5)
        ),
        error=function(e) {
          tryCatch(
            glm(as.formula(paste0("R", condition$k - k, " ~ trt + age + y0", paste0("+ y", 1:(condition$k - k), collapse = " "))),
                data = dat_comp, family=binomial()),
            error=function(e) NULL
          )
        }
      )

      if (is.null(mod)) {
        pred <- rep(
          mean(dat_comp[[rvar]], na.rm=TRUE),
          nrow(dat_comp)
        )
      } else {
        pred <- tryCatch(
          predict(mod, newdata=dat_comp, type="response"),
          error=function(e)
            rep(mean(dat_comp[[rvar]], na.rm=TRUE),
                nrow(dat_comp))
        )
      }

      dat_comp[,paste0("pred_",rvar)] <- pred

      # Subset the data
      daats <- dat_comp |>
        select(
          trt,
          age,
          y0,
          paste0("R", (condition$k - k):1),
          paste0("yc", (condition$k - k):1),
          paste0("pred_R", condition$k - k),
          paste0("j", condition$k - k)
        )

      # Fit a model for the outcome
      model <- lm(
        as.formula(
          paste0("j", condition$k - k, "~ .")
        ),
        data = daats
      )

      dat_comp[, paste0("j", condition$k - k - 1)] <-
        dat_comp[, paste0("j", condition$k - k)] -
        case_when(
          !is.na(coef(model)[paste0("R", condition$k - k)]) ~ coef(model)[paste0("R", condition$k - k)],
          TRUE ~ 0
        ) *
        dat_comp[, paste0("R", condition$k - k)]
      }


      final.model <- lm(j0 ~ trt + y0, data = dat_comp)
      c(
        p = summary(final.model)$coefficients["trt", "Pr(>|t|)"],
        coef = summary(final.model)$coefficients["trt", "Estimate"],
        se = summary(final.model)$coefficients["trt", "Std. Error"]
      )
    }


    completed <- lapply(seq_len(dats$m), function(i) {
      d <- pre_process_data(mice::complete(dats, i))
      d
    })

    res <- lapply(completed, analysis)
    test <- lapply(completed, boot::boot, statistic = analysis, R = 100)
    variance <- function(bootdata){
      B <- bootdata$R
      mu <- mean(bootdata$t[,2])
      var <- sum((bootdata$t[,2]-mu)^2)/(B-1)
      var
    }

    effs <- unlist(lapply(res, function(x) x["coef"]))
    var <- lapply(test,variance)

    end_res <- pool.scalar(effs, unlist(var))
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


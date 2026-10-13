devtools::load_all()

generate_imp <- function(dat){
  daat <- dat |> mutate(
    across(
      starts_with("R") & !R0 & !rescue_start,
      ~ ifelse(is.na(.) & rescue_start <= as.numeric(gsub("R", "", cur_column())), 1, .),
      .names = "R{gsub('R', '', .col)}"
    )
  )
  Rcols <- grep("^R[0-9]+$", names(daat), value = TRUE)
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
}

generate_example_data <- function(){
  Design <- diabetes_scenario()[1, ] |>
    diabetes_scenario_set_truevalues()
  data <- generate_diabetes(Design)
  imp_data <- generate_imp(data)
  return(list(Design = Design, Data = data))
}

analysis_test <- function(dataa, indicator) {
  dat_comp <- dataa[indicator,] |> mutate(
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


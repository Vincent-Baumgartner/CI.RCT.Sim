library(CI.RCT.Sim)
library(profvis)
library(tidyverse)
library(mice)
analyse_diabetes_demediation()(Design, dat)

Design <- diabetes_scenario()[1,] |>
  diabetes_scenario_set_truevalues()

dat <- generate_diabetes(Design)

profvis(analyse_diabetes_demediation()(Design, dat))


analysis <- function(dataa, indicator, condition) {
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

imp <- generate_imp(dat)
test <- complete(imp, 1)
??logistf
profvis::profvis(analysis(test, 1:nrow(test), Design))

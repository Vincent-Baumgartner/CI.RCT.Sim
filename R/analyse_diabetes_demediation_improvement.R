
analyse_diabetes_demediation <- function(separate, B){
  lf_fit <- utils::getFromNamespace("logistf.fit", "logistf")
  lf_ctrl <- logistf::logistf.control(maxit = 2000, maxstep = 0.5)
  function(condition, dat, fixed_objects = NULL){

    # Global Variables
    visits <- 0:condition$k
    R_cols <- paste0("R", visits)
    yc_cols <- paste0("yc", visits[-1])
    ys_cols <- paste0(y, visits)
    j_cols <- paste0("j", visits[-length(visits)])
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

      pred <- make.predictorMatrix(d)
      pred[c(R_cols, "m_start"),] <- 0
      pred[, c(drop, R_cols, "m_start")]

      meth <- make.method(daat)
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
        run_mice(dat0, setup_mice(dat0, drop = c("id", "trt")))
        run_mice(dat1, setup_mice(dat1, drop = c("id", "trt")))
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

    completed_data <- lapply(seq_len(test$imp_data$m), function(i){
      dat_comp <- calculate_baseline_change(mice::complete(test$imp_data, i), yc_cols, ys_cols)
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

      dat_comp <- cbind(dat_comp, j_cols[length(j_cols)] = dat_comp[, "yc12"])
      for (k in rev(visits[-c(1, length(visits))])){

        R_idx <- M_idx(paste0("R", k))
        if (all(is.na(dat_comp[, R_idx])) ||
            sum(dat_comp[, R_idx], na.rm = TRUE) == 0
        ){
          dat_comp <- cbind(dat_comp, j_cols[k - 1] = dat_comp[, j_cols[k]])
          next
        }
      }
    }




  }
}

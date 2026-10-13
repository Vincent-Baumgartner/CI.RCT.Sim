# devtools::install()
# renv::restore()
#library(CI.RCT.Sim)

rm(list=ls())
rpath<-"/home/robin/EMA_Causal_rescue_med/CI.RCT.Sim/R"
Rfiles<-list.files(path = rpath)
for(ff in Rfiles) {
  if(grepl("diabetes",ff)) source(paste(rpath,"/",ff,sep=""))

}
source(paste(rpath,"/","utility_funtions.R",sep=""))

library(future)
library(future.apply)

#library(parallel)

library(survival)
library(mvtnorm)
library(dplyr)
library(mice)
#install.packages("mice")
#library(trtswitch)
#library(SimNPH)




# Define parameter values and derived quantities -------------------------

sim_parameters <- diabetes_scenario() |>
  diabetes_scenario_set_truevalues()

cond<-sim_parameters[9,]
cond$miss[[1]][1]<- -100
data<-generate_diabetes(cond)
condition=cond
dat=data
fixe_objects=NULL
fun<-analyse_diabetes_demediation_test()
system.time({
for(i in 1:10) dm<-analyse_diabetes_demediation_test(data,condition,imp_m=1,imp_method="pmm")
})
dm

#bootstrap
head(data)
B<-10
est<-rep(NA,B)
n<-dim(data)[1]
ind_set<-1:n
ind<-ind_set
i<-1
for(i in 1:(B+1)) {
  if(i>1) ind<-sample(ind_set,n,TRUE)
  data_bs<-data[ind,]
  dm<-analyse_diabetes_demediation_test(data_bs,condition,imp_m=1,imp_method="norm")
  est[i]<-dm$coef
}
SD<-sd(est[-1])
est
mean(est)
data.frame(est=est[1])

str(dm)


#jackknife

sim_parameters <- diabetes_scenario() |>
  diabetes_scenario_set_truevalues()

cond<-sim_parameters[9,]
cond$miss[[1]][1]<- -100


analysis <- function(dataa) {
  dat_comp <- dataa |> mutate(
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

sim_one<-function(j) {
  data<-generate_diabetes(cond)
  est<-rep(NA,n+1)
  dm<-data.frame(p=rep(NA,n+1),coef=NA,SE=NA)
  for(i in 1:(n+1)) {
    if(i>1) data_bs<-data[-i,] else data_bs<-data
    dm[i,]<-analysis(data_bs)
  }
  est<-dm$coef[-1]
  SE_jn<-sqrt((n-1)/n*sum((est-mean(est))^2))
  z<-est[1]/SE
  p_jn<-pnorm(z)
  data.frame(dm[1,],p_jn=p_jn,SE_jn=SE_jn)
}

start<-Sys.time()
sim_one(1)
ende<-Sys.time()
ende-start
10*1000/60/60/30




sim_one(1)
example<-sim_one(1)
example
exampleNA<-example
exampleNA[1,]<-NA
exampleNA
sim_one_try<-function(i) {
  tryCatch(sim_one(i), error = function(e) exampleNA)
}

sim_one_try(1)
plan(multisession, workers = availableCores() - 1)


R<-10

start<-Sys.time()

res <- future_lapply(
  1:R,
  function(i) sim_one_try(i),
  future.seed = TRUE
)
end<-Sys.time()
end-start #ca. 11.2 sec for R=10


n_methods<-dim(res[[1]])[1]
n_param<-dim(res[[1]])[2]

Ar<-array(unlist(res),dim=c(n_methods,n_param,R),dimnames=list(rownames(res[[1]]),names(res[[1]]),1:R))
#T1E
Ar
dim(Ar)
Ar[1,,]
mean(Ar[,,4])
A<-t(Ar[1,,])
A


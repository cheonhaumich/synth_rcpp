# synth_rcpp
There are two major functions: synthetic for point estimation and inference for uncertainty estimation. 
Additionally, ATTplot provides a graph of ATT trajectory surrounded by grey area of +/- 1SD. 

Synthetic 
- input: panel data and a choice of quadratic programming
- output: optimized V, optimized W, ATT trajectory
- functions needed : fn_v.cpp, admm.cpp, R::optimx

Inference
- input: panel data and the number of threads
- output: placebo ATTs across all control units
- functions needed: synthetic.cpp

ATTplot
- input: ATT trajectory from synthetic, placebo ATTs from placebo
- output: a graph of ATT trajectory with grey shadow for +/- 1SD.
- functions needed: R::ggplot2


Sequence is as follows. 
Synthetic (balanced panel data, choice of quadratic programming){
01.	data preparation step 
a.	separating pre-/post- treatment times
b.	separating treatment/control units
c.	handling NAs in observables 
d.	normalizing observables in pre-treatment times

02.	optimize V using R::optimx
a.	R::optimx minimizes scalar output loss_v from fn_V.cpp with BFGS as default
b.	fn_V.cpp itself optimizes solution_v using quadratic programming
c.	execution of quadratic programming optimizer through admm.cpp

03.	optimize W using optimized V from 02.
a.	Re-build input matrix for quadratic programming with solution_v
b.	Execution of quadratic programming through one of admm.cpp, 
c.	solution_w is the output of quadratic programming optimization

04.	calculate ATT trajectory
a.	Build synthetic control outcome trajectory using solution_w and solution_v
b.	Calculate a vector of ATT by subtracting outcome trajectory of treated group from synthetic control outcome trajectory

05.	Return loss_v, loss_w, solution_v, solution_w, ATT
}

Inference (balanced panel){
01.	 Run synthetic across all control units
a.	Make treated unit as control
b.	Choose one control and assign treatment at the same time as the original treated unit was assigned
c.	Run synthetic across all control units
d.	Save vectors of ATT 
02.	Bring back ATTs across different threads
03.	Return placeboATT
}

ATTplot (ATT from synthetic, placeboATT from placebo){
01.	Calculate standard deviation of placeboATT across the number of control units
02.	Load R::ggplot2 and draw ATT plot across time with grey shadow of +/- 1SD
}




Example code

library(Rcpp)
library(RcppArmadillo)
library(devtools)
library(magrittr)
library(dplyr)
library(Synth)
library(optimx)
library(ggplot2)
library(syntheticRcpp)
help(syntheticRcpp)

## how R::synth operates
data(basque)
# dataprep: prepare data for synth
dataprep.out <-
  dataprep(
    foo = basque
    ,predictors= c("school.illit",
                   "school.prim",
                   "school.med",
                   "school.high",
                   "school.post.high"
                   ,"invest"
    )
    ,predictors.op = c("mean")
    ,dependent     = c("gdpcap")
    ,unit.variable = c("regionno")
    ,time.variable = c("year")
    ,special.predictors = list(
      list("gdpcap",1960:1969,c("mean")),                            
      list("sec.agriculture",seq(1961,1969,2),c("mean")),
      list("sec.energy",seq(1961,1969,2),c("mean")),
      list("sec.industry",seq(1961,1969,2),c("mean")),
      list("sec.construction",seq(1961,1969,2),c("mean")),
      list("sec.services.venta",seq(1961,1969,2),c("mean")),
      list("sec.services.nonventa",seq(1961,1969,2),c("mean")),
      list("popdens",1969,c("mean")))
    ,treatment.identifier  = 17
    ,controls.identifier   = c(2:16,18)
    ,time.predictors.prior = c(1964:1969)
    ,time.optimize.ssr     = c(1960:1969)
    ,unit.names.variable   = c("regionname")
    ,time.plot            = c(1955:1997))


R_synth <- synth(data.prep.obj = dataprep.out)
gaps   <- dataprep.out$Y1plot-(dataprep.out$Y0plot%*%R_synth$solution.w)


## make a balanced panel for syntheticRcpp

basque %>%
  mutate(
    treatment=case_when(year < 1975 ~ 0,
                        regionno != 17 ~ 0,
                        regionno == 17 ~ 1) # Basque after 1975 is treated
  ) %>%
  filter(regionno != 1) -> basque_panel

basque_panel <- basque_panel[,-2]
basque_panel_mat <- as.matrix(basque_panel)


m <- basque_panel_mat
u <- 1
t <- 2
y <- 3
d <- 17
c <- c(4:16)

library(ggplot2)
library(optimx)

data1    <- data(m,u,t,y,d,c)
synthcpp <- synthetic( data1$y0, data1$y1, data1$x0_scaled, data1$x1_scaled, data1$z0, data1$z1 )
gaps_cpp <- synthcpp$ATT

# plot comparison (red is Rcpp while black is R)
plot(gaps, type = "l", col = "black", main = "ATT with R:synth and Rcpp:synth", ylab = "ATT", xlab = "time")
lines(gaps_cpp, type = "l", col = "red")

# inference
placebos <- inference( data1$y0, data1$y1, data1$x0_scaled, data1$x1_scaled, data1$z0, data1$z1)

# graph
ATTplot(synthcpp$ATT, m, d, t, placebos)

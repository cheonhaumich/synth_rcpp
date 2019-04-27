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

library(Rcpp)<br />
library(RcppArmadillo)<br />
library(devtools)<br />
library(magrittr)<br />
library(dplyr)<br />
library(Synth)<br />
library(optimx)<br />
library(ggplot2)<br />
library(syntheticRcpp)<br />
help(syntheticRcpp)<br />

## how R::synth operates
data(basque)<br />
# dataprep: prepare data for synth
dataprep.out <-<br />
  dataprep(<br />
    foo = basque<br />
    ,predictors= c("school.illit",<br />
                   "school.prim",<br />
                   "school.med",<br />
                   "school.high",<br />
                   "school.post.high"<br />
                   ,"invest"<br />
    )<br />
    ,predictors.op = c("mean")<br />
    ,dependent     = c("gdpcap")<br />
    ,unit.variable = c("regionno")<br />
    ,time.variable = c("year")<br />
    ,special.predictors = list(<br />
      list("gdpcap",1960:1969,c("mean")),                     <br />       
      list("sec.agriculture",seq(1961,1969,2),c("mean")),<br />
      list("sec.energy",seq(1961,1969,2),c("mean")),<br />
      list("sec.industry",seq(1961,1969,2),c("mean")),<br />
      list("sec.construction",seq(1961,1969,2),c("mean")),<br />
      list("sec.services.venta",seq(1961,1969,2),c("mean")),<br />
      list("sec.services.nonventa",seq(1961,1969,2),c("mean")),<br />
      list("popdens",1969,c("mean")))<br />
    ,treatment.identifier  = 17<br />
    ,controls.identifier   = c(2:16,18)<br />
    ,time.predictors.prior = c(1964:1969)<br />
    ,time.optimize.ssr     = c(1960:1969)<br />
    ,unit.names.variable   = c("regionname")<br />
    ,time.plot            = c(1955:1997))<br />


R_synth <- synth(data.prep.obj = dataprep.out)
gaps   <- dataprep.out$Y1plot-(dataprep.out$Y0plot%*%R_synth$solution.w)


## make a balanced panel for syntheticRcpp

basque %>%<br />
  mutate(<br />
    treatment=case_when(year < 1975 ~ 0,<br />
                        regionno != 17 ~ 0,<br />
                        regionno == 17 ~ 1)<br /> 
  ) %>%<br />
  filter(regionno != 1) -> basque_panel<br />

basque_panel <- basque_panel[,-2]<br />
basque_panel_mat <- as.matrix(basque_panel)<br />


m <- basque_panel_mat <br />
u <- 1<br />
t <- 2<br />
y <- 3<br />
d <- 17<br />
c <- c(4:16)<br />

library(ggplot2)<br />
library(optimx)<br />

data1    <- data(m,u,t,y,d,c)<br />
synthcpp <- synthetic( data1$y0, data1$y1, data1$x0_scaled, data1$x1_scaled, data1$z0, data1$z1 )<br />
gaps_cpp <- synthcpp$ATT<br />

# plot comparison (red is Rcpp while black is R)
plot(gaps, type = "l", col = "black", main = "ATT with R:synth and Rcpp:synth", ylab = "ATT", xlab = "time")<br />
lines(gaps_cpp, type = "l", col = "red")<br />

# inference
placebos <- inference( data1$y0, data1$y1, data1$x0_scaled, data1$x1_scaled, data1$z0, data1$z1)<br />

# graph
ATTplot(synthcpp$ATT, m, d, t, placebos)<br />

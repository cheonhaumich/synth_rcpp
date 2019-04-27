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

#'Estimation of ATT trajectory

#'@param m : a size (N x T) matrix of balanced panel data
#'@param u : (numeric) a column index for unit variable
#'@param t : (numeric) a column index for time variable 
#'@param d : (numeric) a column index for treatment intervention (0 or 1)
#'@param solution_W : a N_ct length of vector for weights
#'@return : a list containing the following variables
#' * ATT : (numeric) a T-dim vector of average treatment effect on treated across time (from the beginning to the end of study time)
#' * T_trajectory : (numeric) a T-dim vector of treatment group outcome trajectory across time
#' * SC_trajectory : (numeric) a T-dim vector of synthetic control outcome trajectory across time
#'@export

ATT_estimate <- function(m, u, t, y, d, solution_w){
  
  g_tr <- aggregate(m[,d] ~ m[,u], FUN = sum)
  id_tr <- g_tr[which(g_tr[,2] !=0),1]
  id_co <- g_tr[which(g_tr[,2] ==0),1]
  
  pre_times  <- m[which((m[,d] == 0)*(m[,u] %in% id_tr)==1),t]
  post_times <- m[which(m[,d] != 0),t]
  when <- min(post_times)

  T_trajectory <- m[which(m[,u]==id_tr),y]
  SC <- m[which(m[,u] != id_tr), y]
  SC <- matrix(SC, ncol = length(unique(id_co)))
  SC_trajectory <- SC%*%solution_w
  
  ATT <- T_trajectory - SC_trajectory
  
  answer <- list(T_trajectory, SC_trajectory, ATT)
  names(answer) <- c("T_trajectory", "SC_trajectory", "ATT")
  
  return(answer)
}
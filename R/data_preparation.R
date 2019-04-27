#'Data preparation for optimization of V,W (.R)

#'@param m : a size (N x T) matrix of balanced panel data
#'@param u : (numeric) a column index for unit variable
#'@param t : (numeric) a column index for time variable 
#'@param d : (numeric) a column index for treatment intervention (0 or 1)
#'@param c : a vector of column indices for time-unit-varying controls 
#'@return : a list containing the following variables
#'y0 = a (T x N_co) matrix of outcomes from control units
#'y1 = a (T-dim) vector of outcomes from treated units
#'x0_scaled = normalized x0
#'x1_scaled = normalized x1
#'z0 = a size (K x 1) matrix of control group outcomes during pre-treatment times
#'z1 = a (T_pre) length vector of treated group outcomes during pre-treatment times
#'@export

data <- function(m, u, t, y, d, c){
  
  g_tr <- aggregate(m[,d] ~ m[,u], FUN = sum)
  id_tr <- g_tr[which(g_tr[,2] !=0),1]
  id_co <- g_tr[which(g_tr[,2] ==0),1]
  
  pre_times  <- m[which((m[,d] == 0)*(m[,u] %in% id_tr)==1),t]
  post_times <- m[which(m[,d] != 0),t]
  when <- min(post_times)
  
  pre <- m[m[,t]%in%pre_times,]
  pre_tr <- pre[pre[,u] %in% id_tr,]
  pre_co <- pre[!(pre[,u] %in% id_tr), ]
  
  # x1 treated group mean controls during pre-treatment times
  # z1 treated group outcomes during pre-treatment times
  z1 <- matrix(pre_tr[,y], ncol=1)
  x1 <- pre_tr[,c]
  x1 <- as.matrix(apply(x1, 2, function(x){mean(x, na.rm = TRUE)}), col = 1)
  
  # x0 control group mean controls during pre-treatment times
  # z0 control group outcomes during pre-treatment times
  z0 <- matrix(pre_co[,y], nrow = length(pre_times))
  x0 <- aggregate(pre_co[,c]~pre_co[,u], FUN = function(x){mean(x, na.rm = TRUE)} )
  x0 <- t(x0[,-c(1)])
  
  # normalize x0, x1
  big <- cbind(x0, x1)
  divisor <- sqrt(apply(big, 1, var))
  divisor <- diag(1/divisor)
  
  big <- t(t(big) %*%divisor)
  
  x0_scaled <- big[, c(1: dim(big)[2]-1)]
  x1_scaled <- matrix(big[, c(dim(big)[2])], ncol=1)
  
  
  y1 <- as.matrix(m[which(m[,u]==id_tr),y], ncol = 1)
  y0 <- matrix(m[which(m[,u]%in%id_co),y], ncol = length(unique(id_co)))
  
  answer <- list (y0, y1, x0_scaled, x1_scaled, z0, z1)
  names(answer) <- c("y0", "y1", "x0_scaled", "x1_scaled", "z0", "z1")
  
  return (answer) 
}
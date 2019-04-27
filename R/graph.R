#'ATT trend graph

#'@param ATT : a (T) length vector of ATTs
#'@param m : a size (N x T) matrix of balanced panel data
#'@param d : (numeric) a column index for treatment intervention (0 or 1)
#'@param t : (numeric) a column index for time variable 
#'@param placebo : a size (T X N_ct) of placebo ATTs
#'@return : a graph plot
#'@export

ATTplot <- function(ATT, m, d, t, placebo){
  
  time <- unique(m[,t])
  post_times <- m[which(m[,d] != 0),t]
  when <- min(post_times)
  
  SE <- apply(placebo, 1, sd)
  
  plot_data <- cbind.data.frame(time, ATT, SE)
  plot_data$ci_lower <- ATT-SE
  plot_data$ci_upper <- ATT+SE

  p <- ggplot(plot_data, aes (x = time, y = ATT)) +
        geom_line() + 
        geom_ribbon ( aes ( x = time, ymin = ci_lower, ymax = ci_upper , fill = "+/- 1SD"), alpha = 0.3) +
    scale_colour_manual("",values="blue")+
    scale_fill_manual("",values="grey12")+
        geom_vline(xintercept = when) +
        geom_hline(yintercept = 0)
   
  print(p)
  
}
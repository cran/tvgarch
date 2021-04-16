tv <- function (speed, location, xtv, opt, order.g)
{
  if (length(location) != order.g) {
    stop("Check the number of location parameters.")
  }
  xtv <- as.numeric(xtv)
  n <- length(xtv)
  if (opt == 0) {
    G <- (1 + exp(- speed * rowProds(xtv - matrix(rep(location, n), byrow = TRUE, ncol = order.g))))^(-1) 
  }
  if (opt == 1) {
    G <- (1 + exp(- exp(speed) / sd(xtv) * rowProds(xtv - matrix(rep(location, n), byrow = TRUE, ncol = order.g))))^(-1) 
  }
  if (opt == 2) {
    G <- (1 + exp(- exp(speed) * rowProds(xtv - matrix(rep(location, n), byrow = TRUE, ncol = order.g))))^(-1) 
  }
  return(G) 
}	

# function to adjust the color range to quantiles
truncate_P = function(P) {
  P_q = quantile(P, probs = c(0.01, 0.99))
  P[P > P_q[2]] = P_q[2]
  P[P < P_q[1]] = P_q[1]
  return(P)
}

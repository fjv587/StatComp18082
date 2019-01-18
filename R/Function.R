#' @title A Monte Carlo estimate of the Beta(3,3) using R
#' @name MCestimate
#' @description A Monte Carlo estimate of the Beta(3,3) using R
#' @param n the replicates
#' @param x the number of between-sample random numbers
#' @return the estimated values returned by the my_MCestimate function and the value computated by the cdf
#' @examples
#' \dontrun{
#' rnR <- MCestimate(1000,c(1,2,3,4,5,6,7,8,9)/10)
#' }
#' @export

MCestimate=function(n,x){
  u=runif(n)  #Generate iid Uniform(0,1) random numbers
  cdf=numeric(length(x))
  for(i in 1:length(x)){
    g=30*x[i]^3*u^2*(1-u*x[i])^2
    cdf[i]=mean(g)  #the estimated value of the cdf
  }
  phi=pbeta(x,3,3)  #the estimated values returned by the pbeta function
  print(round(rbind(x,cdf,phi),3)) #Compare the estimates with the value F(x) computed (numerically) by the pbeta function.
}







#' @title  A function to compute the chi-square test statistic when the input is two numeric vectors
#' @name chisq_test
#' @description  A function to compute the chi-square test statistic when the input is two numeric vectors
#' @param x the first input vector
#' @param y the second input vector
#' @return the chi-square test statistic
#' @examples
#' \dontrun{
#' rnR <-chisq_test(c(3,8,4,9),c(6,9,3,5))
#' }
#' @export



chisq_test=function(x, y) {
  n=sum(x) + sum(y)
  ni._x=sum(x)
  ni._y=sum(y)
  chisq.statistic=0
  for (i in seq_along(x)) {
    n.j=x[i] + y[i]
    expected_x=(ni._x*n.j)/n
    expected_y=(ni._y*n.j)/n
    chisq.statistic=chisq.statistic+((x[i]-expected_x)^2)/expected_x
    chisq.statistic=chisq.statistic+((y[i]-expected_y)^2)/expected_y
  }
  chisq.statistic
}




#' @title  A function to estimate the bias and the standard error of the correlation statistic by jacknife method
#' @name Jacknife.estimate
#' @description   A function to estimate the bias and the standard error of the correlation statistic by jacknife method
#' @param x the first input variable
#' @param y the second input variable
#' @return the estimate of the bias and the standard error of the correlation statistic
#' @examples
#' \dontrun{
#' rnR <-Jacknife.estimate(c(0.1,0.5,1,1.2,0.2,0.4,0.8,0.6,0.8,1.3),c(4,2.7,5.4,3.6,6.3,3.9,4.7,3.9,2.8,5.6))
#' }
#' @export




Jacknife.estimate=function(x,y){
  set.seed(1)
  n=length(x)
  theta.hat=cor(x,y) #the correlation estimate
  print(theta.hat)

  #compute the jackknife replicates, leave-one-out estimates
  theta.jack=numeric(n)
  for(i in 1:n){
    theta.jack[i]=cor(x[-i],y[-i])
  }
  bias=(n-1)*(mean(theta.jack)-theta.hat)
  print(bias)  #a jackknife estimate of the bias of the correlation statistic
  se=sqrt((n-1) *mean((theta.jack - mean(theta.jack))^2))
  print(se)  #a jackknife estimate of the standard error of the correlation statistic
}

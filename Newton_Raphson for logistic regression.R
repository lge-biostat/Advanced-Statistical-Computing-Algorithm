## Newton-Raphson algorithm for
## logistic regression model

fun0 <- function(theta,Y,X){   # target function
  n = length(Y)
  fn = 0
  for(i in 1:n){
    x_i = matrix(X[i,],ncol = 1)
    y_i = Y[i]
    bx = t(theta)%*%x_i
    fn = fn + y_i*bx - log(1+exp(bx))
  }
  return(fn)
}
grd0 <- function(theta,Y,X){      # gradient
  n = length(Y)
  p = length(theta)
  fn = matrix(rep(0,p),ncol = 1)
  for(i in 1:n){
    x_i = matrix(X[i,],ncol = 1)
    y_i = Y[i]
    exp_bx = as.numeric(exp(t(theta)%*%x_i))
    fn = fn + (y_i - exp_bx/(1+exp_bx))*x_i
  }
  return(fn)
}   
hes0 <- function(theta,Y,X){       # Hessian
  n = length(Y)
  p = length(theta)
  fn = matrix(rep(0,p*p),ncol = p)
  for(i in 1:n){
    x_i = matrix(X[i,],ncol = 1)
    y_i = Y[i]
    exp_bx = as.numeric(exp(t(theta)%*%x_i))
    fn = fn - exp_bx/(1+exp_bx)^2*x_i%*%t(x_i)
  }
  return(fn)
}

# Newton-Raphson Algorithm
Newton_Raphson <- function(Y,X, fun=fun0, grd=grd0, hes=hes0, kmax=1000, 
                            tol1=1e-6){
  diff <- 10
  k <- 0
  theta = matrix(glm(Y~X-1)$coef,ncol = 1)
  #theta = matrix(X[1,],ncol = 1)
  while( sum(diff^2) > tol1 & k <= kmax){
    g_x <- grd(theta,Y,X)
    h_x <- hes(theta,Y,X)     # calculate the second derivative (Hessian)
    diff <- -solve(h_x)%*%g_x  # calculate the difference used by the stop criteria
    theta_new <- theta + diff
    theta = theta_new
    k <- k + 1
    cat("Iter", k, ": theta=", theta,  "\n")
  }
  return(list(iteration = k, theta = t(theta_new)))
}


## testing
N = 1000       # sample size
n_covar = 2    # number of covariates
X = matrix(rnorm(N*n_covar),ncol=n_covar)
theta_true = 0.3
p = exp(X*theta_true)/(1+exp(X*theta_true))
Y = rbinom(N,1,p)

(NR_result = Newton_Raphson(Y,X))
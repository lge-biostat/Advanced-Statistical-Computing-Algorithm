## IRLS (Iteratively reweighted least squares) algorithm
## for Poison Regression 

poisreg = function(y, X, kmax=1000, tol=1e-6){
  beta = rep(1, dim(X)[2])
  cnt = 0
  diff = 1e5
  while(diff > tol & cnt<kmax){
    mu_k = exp(X%*%beta)
    y_k = log(mu_k) + (y-mu_k)/mu_k
    W_k = diag(as.vector(mu_k))
    beta_k = solve((t(X)%*%W_k%*%X))%*%(t(X)%*%W_k%*%y_k)
    diff = sum((beta_k-beta)^2)
    beta = beta_k
    cnt = cnt +1
  }
  var_beta = solve((t(X)%*%W_k%*%X))
  return(list(beta = beta_k,var_beta = var_beta))
}

#############################################################

## Comparing the "poisreg" with "glm" from R

# parameters
n = 100 ## number of observations
p = 3 ## number of covariates
## generate X, the covariates
X = cbind(1, matrix(rnorm(n*p), ncol=p))
beta = c(1, .5, 1, 2)
mu = exp(X %*% beta)
## generate Y, the outcome
Y = rpois(n, mu)


### use R’s glm function to fit
fit = glm(Y~X-1, family=poisson)
coef(fit) ## estimated coefficients
vcov(fit) ## estimated variance/covariance matrix of the estimates

### use poisreg function to fit
poisreg(Y,X)
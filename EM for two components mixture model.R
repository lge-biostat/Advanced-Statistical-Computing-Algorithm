## EM (expectation–maximization) algorithm for
## two components mixture model

# Part I: Two components mixture of Gaussian model

EM_MoG = function(p, mu1, mu2, sd1, sd2, X, maxiter=1000, tol=1e-5)
{
  diff=1
  iter=0
  while (diff>tol & iter<maxiter) {
    ## E-step: compute omega:
    d1=dnorm(X, mean=mu1, sd=sd1)
    d2=dnorm(X, mean=mu2, sd=sd2)
    omega=d1*p/(d1*p+d2*(1-p))
    # compute density in two groups
    ## M-step: update p, mu and sd
    p.new=mean(omega)
    mu1.new=sum(X*omega) / sum(omega)
    mu2.new=sum(X*(1-omega)) / sum(1-omega)
    resid1=X-mu1
    resid2=X-mu2
    sd1.new=sqrt(sum(resid1^2*omega) / sum(omega))
    sd2.new=sqrt(sum(resid2^2*(1-omega)) / sum(1-omega))
    
    ## calculate diff to check convergence
    diff=sqrt(sum((mu1.new-mu1)^2+(mu2.new-mu2)^2
                  +(sd1.new-sd1)^2+(sd2.new-sd2)^2))
    
    p=p.new
    mu1=mu1.new
    mu2=mu2.new
    sd1=sd1.new
    sd2=sd2.new
    iter=iter+1
    
    cat("Iter", iter, ": mu1=", mu1.new, ", mu2=",mu2.new, ", sd1=",sd1.new,
        ", sd2=",sd2.new, ", p=", p.new, ", diff=", diff, "\n")
  }
  
}

## testing
# true values
p0=0.3
n = 10000
X1 = rnorm(n*p0)
X2 = rnorm(n*(1-p0),mean=4)
X = c(X1,X2)
#hist(X,50)

# initial values
p=0.5
mu1=quantile(X,0.1)
mu2=quantile(X,0.9)
sd1=sd2=sd(X)

# model fiting
EM_MoG(p, mu1, mu2, sd1, sd2, X)


# Part II: Two componets mixture of Poison model

EM_MoPoison = function(y,lambda1,lambda2,p,maxiter=1000, tol=1e-5){
  diff=1
  iter=0
  while(diff>tol & iter<maxiter){
    ## E-step: compute omega:
    d1=dpois(y, lambda1)
    d2=dpois(y, lambda2)
    omega=d1*p/(d1*p+d2*(1-p))
    
    ## M-step: update p, mu and sd
    p.new=mean(omega)
    lambda1.new=sum(y*omega) / sum(omega)
    lambda2.new=sum(y*(1-omega)) / sum(1-omega)
    
    ## calculate diff to check convergence
    diff=sqrt(sum((lambda1.new-lambda1)^2+(lambda2.new-lambda2)^2))
    
    p = p.new
    lambda1 = lambda1.new
    lambda2 = lambda2.new
    iter = iter + 1
    
    cat("Iter", iter, ": lambda1=", lambda1.new, ", lambda2=",lambda2.new, 
        ", p=", p.new, ", diff=", diff, "\n")
  }
}

## testing
# true values
p0 = 0.3
n  = 10000
y1 = rpois(n*p0, lambda=5)
y2 = rpois(n*(1-p0), lambda=15)
Y = c(y1,y2)

# initial values
lambda1 = 1
lambda2 = 3
p0 = 0.5

EM_MoPoison(Y,lambda1,lambda2, p0)
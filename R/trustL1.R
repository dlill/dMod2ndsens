
norm <- function(x) sqrt(sum(x^2))

########## REFERENCES ##########
#####
##### Fletcher, R. (1987)
##### Practical Methods of Optimization, second edition.
##### John Wiley, Chichester.
#####
##### Nocedal, J. and Wright, S. J. (1999)
##### Numerical Optimization.
##### Springer-Verlag, New York.
#####
##### See Section 5.1 of Fletcher
##### See Section 4.2 of Nocedal and Wright
#####
################################

########## COMMENT ##########
##### Our method using one eigendecomposition per iteration is not fastest.
##### Both books recommend using multiple Cholesky decompositions instead.
##### But the eigendecomposition method is simpler to program, easier to
##### understand (which is why both books use it for their theoretical
##### explanation), and hopefully more bulletproof.
#####
##### Our idea for this comes from the way mvrnorm in the MASS package also
##### uses eigendecomposition rather than Cholesky -- also because bulletproof
##### is better than fast.
#############################

trustL1 <- function(objfun, parinit, mu = 0*parinit, one.sided=FALSE, lambda = 1, rinit, rmax, parscale,
    iterlim = 100, fterm = sqrt(.Machine$double.eps),
    mterm = sqrt(.Machine$double.eps),
    minimize = TRUE, blather = FALSE, blather2 = FALSE, ...)
{
    if (! is.numeric(parinit))
       stop("parinit not numeric")
    if (! all(is.finite(parinit)))
       stop("parinit not all finite")
    d <- length(parinit)
    if (missing(parscale)) {
        rescale <- FALSE
    } else {
        rescale <- TRUE
        if (length(parscale) != d)
           stop("parscale and parinit not same length")
        if (! all(parscale > 0))
           stop("parscale not all positive")
        if (! all(is.finite(parscale) & is.finite(1 / parscale)))
           stop("parscale or 1 / parscale not all finite")
    }
    if (! is.logical(minimize))
       stop("minimize not logical")

    r <- rinit
    theta <- parinit
    out <- try(objfun(theta, ...))
    grad0 <- out$gradient
    outL1 <- try(constraintL1(theta, mu, lambda))
    gradL1 <- outL1$gradient
    
    if (inherits(out, "try-error")) {
        warning("error in first call to objfun")
        return(list(error = out, argument = theta, converged = FALSE,
            iterations = 0))
    }
    
    ## Fix L1-parameters on prior if they would be drawn back after step
    is.fixed.theta <- match(names(mu), names(theta))[which(theta[names(mu)] == mu & abs(grad0[names(mu)]) <= lambda)]
    if(one.sided){
      is.fixed.theta <- match(names(mu), names(theta))[which(theta[names(mu)] == mu & -(grad0[names(mu)]) <= lambda)]
    }
    if(length(is.fixed.theta) > 0) {
      
      out$gradient <- out$gradient[-is.fixed.theta]
      out$hessian <- out$hessian[-is.fixed.theta,-is.fixed.theta]  
      outL1$gradient <- outL1$gradient[-is.fixed.theta]
      outL1$hessian <- outL1$hessian[-is.fixed.theta, -is.fixed.theta]
      
    }
    
    

    
    out <- out + outL1
    check.objfun.output(out, minimize, d - length(is.fixed.theta))
    if (! is.finite(out$value))
        stop("parinit not feasible")
    accept <- TRUE

    if (blather) {
        theta.blather <- NULL
        theta.try.blather <- NULL
        type.blather <- NULL
        accept.blather <- NULL
        r.blather <- NULL
        stepnorm.blather <- NULL
        rho.blather <- NULL
        val.blather <- NULL
        val.try.blather <- NULL
        preddiff.blather <- NULL
    }

    for (iiter in 1:iterlim) {
      if(blather2)
        print(paste(iiter,out$value,accept))

        if (blather) {
            theta.blather <- rbind(theta.blather, theta)
            r.blather <- c(r.blather, r)
            if (accept)
                val.blather <- c(val.blather, out$value)
            else
                val.blather <- c(val.blather, out.value.save)
        }

        if (accept) {
            B <- out$hessian
            g <- out$gradient
            f <- out$value
            out.value.save <- f
            if (rescale) { 
                B <- B / outer(parscale, parscale)
                g <- g / parscale
            }
            if (! minimize) {
                B <- (- B)
                g <- (- g)
                f <- (- f)
            }
            eout <- eigen(B, symmetric = TRUE)
            gq <- as.numeric(t(eout$vectors) %*% g)
        }

        ########## solve trust region subproblem ##########

        ##### try for Newton #####
        is.newton <- FALSE
        if (all(eout$values > 0)) {
            ptry <- as.numeric(- eout$vectors %*% (gq / eout$values))
            if (norm(ptry) <= r)
                is.newton <- TRUE
        }

        ##### non-Newton #####
        if (! is.newton) {
            lambda.min <- min(eout$values)
            beta <- eout$values - lambda.min
            imin <- beta == 0
            C1 <- sum((gq / beta)[! imin]^2)
            C2 <- sum(gq[imin]^2)
            C3 <- sum(gq^2)
            if (C2 > 0 || C1 > r^2) {
                is.easy <- TRUE
                is.hard <- (C2 == 0)
                ##### easy cases #####
                beta.dn <- sqrt(C2) / r
                beta.up <- sqrt(C3) / r
                fred <- function(beep) {
                    if (beep == 0) {
                        if (C2 > 0)
                            return(- 1 / r)
                        else
                            return(sqrt(1 / C1) - 1 / r)
                    }
                    return(sqrt(1 / sum((gq / (beta + beep))^2)) - 1 / r)
                }
                if (fred(beta.up) <= 0) {
                    uout <- list(root = beta.up)
                } else if (fred(beta.dn) >= 0) {
                    uout <- list(root = beta.dn)
                } else {
                    uout <- uniroot(fred, c(beta.dn, beta.up))
                }
                wtry <- gq / (beta + uout$root)
                ptry <- as.numeric(- eout$vectors %*% wtry)
            } else {
                is.hard <- TRUE
                is.easy <- FALSE
                ##### hard-hard case #####
                wtry <- gq / beta
                wtry[imin] <- 0
                ptry <- as.numeric(- eout$vectors %*% wtry)
                utry <- sqrt(r^2 - sum(ptry^2))
                if (utry > 0) {
                    vtry <- eout$vectors[ , imin, drop = FALSE]
                    vtry <- vtry[ , 1]
                    ptry <- ptry + utry * vtry
                }
            }
        }

        
        ########## predicted versus actual change ##########

        preddiff <- sum(ptry * (g + as.numeric(B %*% ptry) / 2))
        
        
        ## Compute theta.try
        
        ## Fix prior parameters which are on prior (catch-up from above)
        if(length(is.fixed.theta) > 0) {
          ptry.new <- structure(rep(0, length(parinit)), names = names(parinit))
          ptry.new[names(parinit)[-is.fixed.theta]] <- ptry
          ptry <- ptry.new
        }
        
        if (rescale) {
            theta.try <- theta + ptry / parscale
        } else {
            theta.try <- theta + ptry
        }
        
        ## Set on prior value if step-over
        chgsgn <- (theta[names(mu)]-mu)*(theta.try[names(mu)]-mu)
        theta.try[names(mu)][chgsgn < 0] <- mu[chgsgn < 0]
        if(one.sided){  
          theta.try[names(mu)][theta.try[names(mu)] < mu] <- mu[theta.try[names(mu)] < mu]
        }
        ## Scale down step length if step over prior
        #         chgsgn <- (theta[names(mu)]-mu)*(theta.try[names(mu)]-mu)
        #         steplength.red <- (mu - theta[names(mu)])[chgsgn < 0]
        #         steplength.full <- (theta.try[names(mu)] - theta[names(mu)])[chgsgn < 0]
        #         if(length(steplength.red) > 0) {
        #           fact <- abs(steplength.red/steplength.full)
        #           theta.try <- theta + min(fact)*(theta.try-theta)
        #           theta.try[names(mu)][chgsgn < 0][which.min(fact)] <- mu[chgsgn < 0][which.min(fact)]
        #           ptry.red <- theta.try - theta
        #           if(length(is.fixed.theta) > 0) ptry.red <- ptry.red[-is.fixed.theta]
        #           preddiff <- sum(ptry.red * (g + as.numeric(B %*% ptry.red) / 2))
        #           
        #         }
                
        
        out <- try(objfun(theta.try, ...))
        outL1 <- try(constraintL1(theta.try, mu, lambda))
        if (inherits(out, "try-error"))
            break
        
        
        ## Fix L1-parameters on prior if they would be drawn back after step (theta.try)
      
       # is.fixed.theta.try <- which(names(theta.try)%in%names(mu))[which(theta.try[names(mu)] == mu & abs(out$gradient[names(mu)]) <= lambda)]
      is.fixed.theta.try <- match(names(mu), names(theta.try))[which(theta.try[names(mu)] == mu & abs(out$gradient[names(mu)]) <= lambda)]
      
        if(one.sided){
          is.fixed.theta.try <- match(names(mu), names(theta.try))[which(theta.try[names(mu)] == mu & -(out$gradient[names(mu)]) <= lambda)]
        }
      
        if(length(is.fixed.theta.try) > 0) {
          
          out$gradient <- out$gradient[-is.fixed.theta.try]
          out$hessian <- out$hessian[-is.fixed.theta.try,-is.fixed.theta.try]  
          outL1$gradient <- outL1$gradient[-is.fixed.theta.try]
          outL1$hessian <- outL1$hessian[-is.fixed.theta.try, -is.fixed.theta.try]
          
        }
                
        out <- out + outL1
        
        check.objfun.output(out, minimize, d - length(is.fixed.theta.try))
        ftry <- out$value
        if (! minimize)
            ftry <- (- ftry)
        rho <- (ftry - f) / preddiff

        ########## termination test ##########
        if (ftry < Inf) {
            is.terminate <- abs(ftry - f) < fterm || abs(preddiff) < mterm
        } else {
            is.terminate <- FALSE
            rho <- (- Inf)
        }

        ##### adjustments #####
        if (is.terminate) {
            if (ftry < f) {
                accept <- TRUE
                theta <- theta.try
                is.fixed.theta <- is.fixed.theta.try
            }
        } else {
            if (rho < 1 / 4) {
                accept <- FALSE
                r <- r / 4
            } else {
                accept <- TRUE
                theta <- theta.try
                is.fixed.theta <- is.fixed.theta.try
                if (rho > 3 / 4 && (! is.newton))
                    r <- min(2 * r, rmax)
            }
        }

        if (blather) {
            theta.try.blather <- rbind(theta.try.blather, theta.try)
            val.try.blather <- c(val.try.blather, out$value)
            accept.blather <- c(accept.blather, accept)
            preddiff.blather <- c(preddiff.blather, preddiff)
            stepnorm.blather <- c(stepnorm.blather, norm(ptry))
            if (is.newton) {
                mytype <- "Newton"
            } else {
                if (is.hard) {
                    if (is.easy) {
                        mytype <- "hard-easy"
                    } else {
                        mytype <- "hard-hard"
                    }
                } else {
                    mytype <- "easy-easy"
                }
            }
            type.blather <- c(type.blather, mytype)
            rho.blather <- c(rho.blather, rho)
        }

        if (is.terminate)
            break
    }

    if (inherits(out, "try-error")) {
        out <- list(error = out, argument = theta.try, converged = FALSE)
    } else {
        out <- try(objfun(theta, ...))
        outL1 <- try(constraintL1(theta.try, mu, lambda))
        if (inherits(out, "try-error")) {
          out <- list(error = out)
          warning("error in last call to objfun")
        } else {
          out <- out + outL1
          check.objfun.output(out, minimize, d)
        }
        out$argument <- theta
        out$converged <- is.terminate
    }
    out$iterations <- iiter
    if (blather) {
        dimnames(theta.blather) <- NULL
        out$argpath <- theta.blather
        dimnames(theta.try.blather) <- NULL
        out$argtry <- theta.try.blather
        out$steptype <- type.blather
        out$accept <- accept.blather
        out$r <- r.blather
        out$rho <- rho.blather
        out$valpath <- val.blather
        out$valtry <- val.try.blather
        if (! minimize)
            preddiff.blather <- (- preddiff.blather)
        out$preddiff <- preddiff.blather
        out$stepnorm <- stepnorm.blather
    }
    return(out)
}

check.objfun.output <- function(obj, minimize, dimen)
{
    if (! is.list(obj))
        stop("objfun returned object that is not a list")
    foo <- obj$value
    if (is.null(foo))
        stop("objfun returned list that does not have a component 'value'")
    if (! is.numeric(foo))
        stop("objfun returned value that is not numeric")
    if (length(foo) != 1)
        stop("objfun returned value that is not scalar")
    if (is.na(foo) || is.nan(foo))
        stop("objfun returned value that is NA or NaN")
    if (minimize && foo == (-Inf))
        stop("objfun returned -Inf value in minimization")
    if ((! minimize) && foo == Inf)
        stop("objfun returned +Inf value in maximization")
    if (is.finite(foo)) {
        bar <- obj$gradient
        if (is.null(bar))
            stop("objfun returned list without component 'gradient' when value is finite")
        if (! is.numeric(bar))
            stop("objfun returned gradient that is not numeric")
        if (length(bar) != dimen)
            stop(paste("objfun returned gradient that is not vector of length", dimen))
        if (! all(is.finite(bar)))
            stop("objfun returned gradient not having all elements finite")
        baz <- obj$hessian
        if (is.null(baz))
            stop("objfun returned list without component 'hessian' when value is finite")
        if (! is.numeric(baz))
            stop("objfun returned hessian that is not numeric")
        if (! is.matrix(baz))
            stop("objfun returned hessian that is not matrix")
        if (! all(dim(baz) == dimen))
            stop(paste("objfun returned hessian that is not", dimen, "by", dimen, "matrix"))
        if (! all(is.finite(baz)))
            stop("objfun returned hessian not having all elements finite")
    }
    return(TRUE)
}


#' Soft L1 constraint on parameters
#' 
#' @param p Namec numeric, the parameter value
#' @param mu Named numeric, the prior values
#' @param lambda Named numeric of length of mu or numeric of length one.
#' @param fixed Named numeric with fixed parameter values (contribute to the prior value
#' but not to gradient and Hessian)
#' @return List of class \code{obj}, i.e. objective value, gradient and Hessian as list.
#' @details Computes the constraint value 
#' \deqn{\lambda\|p-\mu\|}{lambda*abs(p-mu)}
#' and its derivatives with respect to p.
#' @seealso \link{wrss}, \link{summation}, \link{constraintL2}, \link{constraintExp2}
#' @examples
#' p <- c(A = 1, B = 2, C = 3)
#' mu <- c(A = 0, B = 0)
#' lambda <- c(A = 0.1, B = 1)
#' constraintL1(p, mu, lambda)
constraintL1 <- function(p, mu, lambda = 1, fixed = NULL) {
   
  ## Augment sigma if length = 1
  if(length(lambda) == 1) 
    lambda <- structure(rep(lambda, length(mu)), names = names(mu)) 
  
  ## Extract contribution of fixed pars and delete names for calculation of gr and hs  
  par.fixed <- intersect(names(mu), names(fixed))
  sumOfFixed <- 0
  if(!is.null(par.fixed)) sumOfFixed <- sum(lambda[par.fixed]*abs(fixed[par.fixed] - mu[par.fixed]))
  
  ## Compute constraint value and derivatives
  parameters <- intersect(names(p), names(mu))
    
  value <- sum(lambda[parameters]*abs(p[parameters] - mu[parameters])) + sumOfFixed
  
  gradient <- rep(0, length(p)); names(gradient) <- names(p)
  gradient[parameters][p[parameters] >  mu[parameters]] <-  lambda[parameters][p[parameters] >  mu[parameters]]
  gradient[parameters][p[parameters] <  mu[parameters]] <- -lambda[parameters][p[parameters] <  mu[parameters]]
  
  hessian <- matrix(0, length(p), length(p), dimnames = list(names(p), names(p)))
  diag(hessian)[parameters] <- 0
  
  dP <- attr(p, "deriv") 
  if(!is.null(dP)) {
    gradient <- as.vector(gradient %*% dP)
    names(gradient) <- colnames(dP)
    hessian <- t(dP) %*% hessian %*% dP
    colnames(hessian) <- colnames(dP)
    rownames(hessian) <- colnames(dP)
  }
  
  out <- list(value = value, gradient = gradient, hessian = hessian)
  class(out) <- c("obj", "list")
  
  return(out)
  
  
}



#' Soft L1 prior on parameters
#' 
#' @param p Namec numeric, the parameter value
#' @param mu Named numeric, the prior values
#' @param lambda Named numeric of length of mu or numeric of length one.
#' @param fixed Named numeric with fixed parameter values (contribute to the prior value
#' but not to gradient and Hessian)
#' @return List of class \code{obj}, i.e. objective value, gradient and Hessian as list.
#' @details Computes the constraint value 
#' \deqn{\lambda\|p-\mu\|}{lambda*abs(p-mu)}
#' and its derivatives with respect to p.
#' @seealso \link{wrss}, \link{summation}, \link{constraintL2}, \link{constraintExp2}
#' @examples
#' p <- c(A = 1, B = 2, C = 3)
#' mu <- c(A = 0, B = 0)
#' lambda <- c(A = 0.1, B = 1)
#' constraintL1(p, mu, lambda)
# priorL1 <- function(p, mu, lambda = "lambda", fixed = NULL) {
#    
#   ## Extract contribution of fixed pars and delete names for calculation of gr and hs  
#   par.fixed <- intersect(names(mu), names(fixed))
#   sumOfFixed <- 0
#   if(!is.null(par.fixed)) sumOfFixed <- sum(c(fixed, p)[lambda]*abs(fixed[par.fixed] - mu[par.fixed]))
#   
#   ## Compute constraint value and derivatives
#   par <- intersect(names(p), names(mu))
#   par0 <- setdiff(par, lambda)
#     
#   value <- sum(c(fixed, p)[lambda]*abs(p[par] - mu[par])) + sumOfFixed
#   
# 
#   direction <- rep(0, length(p)); names(direction) <- names(p)
#   direction[par][p[par] >  mu[par]] <-  1
#   direction[par][p[par] <  mu[par]] <- -1
#   gradient <- c(fixed, p)[lambda]*direction
#   if(lambda %in% names(p)) 
#     gradient[lambda] <- value
#   
#   hessian <- matrix(0, length(p), length(p), dimnames = list(names(p), names(p)))
#   diag(hessian)[par] <- 0
#   if(lambda %in% names(p)) {
#     hessian[lambda, lambda] <- 0 
#     hessian[lambda, par0] <- hessian[par0, lambda] <- direction[par0] 
#   }
#   
#   dP <- attr(p, "deriv") 
#   if(!is.null(dP)) {
#     gradient <- as.vector(gradient %*% dP)
#     names(gradient) <- colnames(dP)
#     hessian <- t(dP) %*% hessian %*% dP
#     colnames(hessian) <- colnames(dP)
#     rownames(hessian) <- colnames(dP)
#   }
#   
#   out <- list(value = value, gradient = gradient, hessian = hessian)
#   class(out) <- c("obj", "list")
#   
#   return(out)
#   
#   
# }
#
#
# constraintLeins <- function(p, mu, lambda = 1, tol = 1e-3) {
#   
#   parameters <- intersect(names(p), names(mu))
#   
#     
#   parameters.l1 <- parameters[abs(p[parameters] - mu[parameters]) >  tol]
#   parameters.l2 <- parameters[abs(p[parameters] - mu[parameters]) <= tol]
#   
#   sigma <- sqrt(tol/lambda)
#   offset <- lambda*tol - .5*(tol/sigma)^2 
#   
#   prior1 <- prior2 <- as.obj(p)
#   
#   if(length(parameters.l1) > 0) {
#     prior1 <- constraintL1(p, mu[parameters.l1], lambda)
#     if(prior1$value > 0) prior1$value <- prior1$value - offset
#   }
#   if(length(parameters.l2) > 0) {
#     prior2 <- constraintL2(p, mu[parameters.l2], sigma)
#   }
#     
#   prior1 + prior2
#   
# }


#' Compute the L1 norm of residuals
#' 
#' @param nout data.frame (result of \link{res})
#' @return list with entries value (numeric, the weighted residual sum of squares), 
#' gradient (numeric, gradient) and 
#' hessian (matrix of type numeric).
l1norm <- function(nout) {
  
  sign.res <- sign(nout$weighted.residual)
  obj <- sum(abs(nout$weighted.residual))
  grad <- NULL
  hessian <- NULL
  
  if(!is.null(attr(nout, "deriv"))) {
    nout$sigma[is.na(nout$sigma)] <- 1 #replace by neutral element
    
    sens <- as.matrix(attr(nout, "deriv")[,-(1:2)])
    grad <- apply(-sign.res*sens/nout$sigma, 2, sum)
    names(grad) <- colnames(sens)
    hessian <- t(sens/nout$sigma)%*%(sens/nout$sigma)
    
    
  }
  
  out <- list(value=obj, gradient=grad, hessian=hessian)
  class(out) <- c("obj", "list")
  
  return(out)
  
}


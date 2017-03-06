#' Get the sdot expression from the output of the output of the integration of the approximated geodesic eqn with constant christoffel symbols.
#' this I do because bvptwpC has problems if i add the sdot=sqrt(gmunu zmu znu) term to the equations. 
#' It somehow has problems with accuracy even though it's just a value to be evaluated, it doesn't affect the solution of the other terms.
#' 
#' The expression for sdot is again just an approximation, that replaces the integral by a discrete sum.
#' 
#' @param bestfit bestfit
#' 
#' @param out Argument of returned function, the output of the numerical integration of the geodesic
#' @param line Argument of the returned function, the line to be accessed in "out".
#' 
#' 
#' @export
sdotfun<- function(bestfit) {
  
  myresidual  <- resi(bestfit)
  myderiv     <- data.matrix(attr(myresidual, "deriv")[,-c(1,2)])
  mysderiv    <- data.matrix(attr(myresidual, "sderiv")[,-c(1,2)])
  gmunu       <- t(myderiv)%*%myderiv
  # covariant derivative of gmunu
  gmunu.alpha <- array(t(myderiv)%*%mysderiv, dim=rep(length(bestfit), 3)) #indices: (mu, nu, alpha). this term is t(dr/dmu)%*%(d^2r/dnu dalpha)
  gmunu.alpha <- array(gmunu.alpha + aperm(gmunu.alpha, c(2,1,3)), dim=c(length(bestfit)^2, length(bestfit))) #add with swapped mu<->nu # indices: (munu,alpha)
  ##################################################
  # !!!! substract the christoffel symbol part of the covariant derivative###
  ##################################################
  
  
  sdot <- function(out, line) {
    Theta         <- out[line,-c(1)]
    Thalpha       <- Theta[1:length(bestfit)]
    gmunu.alpha   <- matrix(gmunu.alpha%*%Thalpha, nrow=length(bestfit))
    gmunu         <- gmunu+gmunu.alpha
    
    value         <- sqrt(Theta[(length(bestfit)+1):(2*length(bestfit))]%*%gmunu%*%Theta[(length(bestfit)+1):(2*length(bestfit))])
    names(value)  <- "sdot"
    return(value)
  }
  
  return(sdot)
}




#' Integrate sdot and square the total distance travelled to get an expression for Delta chi^2
#' 
#' It is assumed that the sqrt-expression of sdot is constant for each time point in "out",
#' so the integral becomes the sum of the sqrts times the DeltaTaus provided by "out", the solution of the geodesic bvp
#' 
#' The function sdot <- sdotfun(bestfit) has to be evaluated in advance
#' 
#' Giving the residual explicitly is necessary to avoid conflicting dlls.
#' 
#' @param out the output of the numerical integration of the geodesic
#' @param sigma sigma of the data
#' @param bestfit
#' 
#' @export
ssquared<-function(out, sigma, residual) {
  chi0      <- residual[,7]%*%residual[,7]
  Deltatau  <- diff(out[,1])
  val       <- do.call(rbind, lapply(1:length(Deltatau), function(i) {Deltatau[i]*sdot(out,i)}))
  val       <- diffinv(val)
  val       <- val[length(val)]
  return(val^2/(sigma^2)+as.numeric(chi0))
}



#' Jacobian dp/dTh in approximate metric approach
#' 
#' @param Thend 
#' @param out 
#' @param sigma 
#' @param bestfit 
#' @param residual 
#' @param stepwidth 
#'
#' @export
RNCjacssfinite <- function(Thend, out, sigma, bestfit, residual, stepwidth) {
  
  # stepwidth <- 0.001
  xguess    <- out[,1]
  yguess    <- yguessout(out)
  
  mybd        <- geodesicboundaries(bestfit, Thend)
  yini        <- mybd[1:length(bestfit),2]
  names(yini) <- names(bestfit)
  yend        <- mybd[1:length(bestfit),3]
  names(yend) <- names(bestfit)
  
  x           <- xguess
  vini        <- out[1,-(1:(length(Thend)+1))]
  
  
  myjac <- do.call(rbind, lapply(1:length(Thend), function(n) {
    # Variate the parameters
    newyend      <- yend
    newyend[n]   <- newyend[n] + stepwidth
    newvini      <- bvptwpC(x = x, func = eqns, yini = yini, yend = newyend, parms = pars, nmax=50000, xguess = xguess, yguess = yguess)[1,-(1:(length(Thend)+1))]
    
    return((newvini - vini)/stepwidth)
  }) ) # matrix (dp^mu/dTh^alpha)^mu_alpha
  return(myjac)
}




#' Compute the gradient of s^2 with finite differences
#' 
#' This function computes the gradient of s^2 with finite differences. s^2 is the squared arclength of the geodesic. 
#' It takes really long to evaluate but I know of no other expression for the derivative.
#' 
#' @param Thend
#' @param out
#' @param sigma
#' @param bestfit
#' @param residual
#' @param stepwidth
#' 
#' 
#' @export
gradssfinite <- function(Thend, out, sigma = 0.02, bestfit, residual, stepwidth) {
  
  # stepwidth <- 0.001
  xguess    <- out[,1]
  yguess    <- yguessout(out)
  
  mybd        <- geodesicboundaries(bestfit, Thend)
  yini        <- mybd[1:length(bestfit),2]
  names(yini) <- names(bestfit)
  yend        <- mybd[1:length(bestfit),3]
  names(yend) <- names(bestfit)
  
  if(is.null(out)) {
    x <- seq(0,1, by= 0.1)
    myssquared <- 0
  } else {
    x <- xguess
    myssquared <- ssquared(out, sigma, residual)[length(xguess)]
  }
  
  mygrad <- unlist(lapply(names(Thend), function(n) {
    # Variate the parameters
    newyend       <- yend
    newyend[n]    <- newyend[n] + stepwidth
    newout        <- bvptwpC(x = x, func = eqns, yini = yini, yend = newyend, parms = pars, nmax=50000, xguess = xguess, yguess = yguess)
    newssquared   <- ssquared(newout, sigma, residual)[nrow(newout)]
    
    # Take finite difference
    return((newssquared-myssquared)/stepwidth)
  }) )
  
  names(mygrad) <- names(Thend)
  return(mygrad)
}



#' Hessian of ssquared with finite differences
#' 
#' This function computes the hessian of s^2 with finite differences. s^2 is the squared arclength of the geodesic. 
#' It takes really long to evaluate but I know of no other expression for the derivative.
#' 
#' @param Thend
#' @param out
#' @param sigma
#' @param bestfit
#' @param residual
#' @param stepwidth
#' 
#' @export
hessianssfinite <- function(Thend, out, sigma = 0.02, bestfit, residual, stepwidth) {
  
  # stepwidth <- 0.001
  xguess    <- out[,1]
  yguess    <- yguessout(out)
  
  mybd        <- geodesicboundaries(bestfit, Thend)
  yini        <- mybd[1:length(bestfit),2]
  names(yini) <- names(bestfit)
  yend        <- mybd[1:length(bestfit),3]
  names(yend) <- names(bestfit)
  
  if(is.null(out)) {
    x <- seq(0,1, by= 0.1)
    myssquared <- 0
  } else {
    x <- xguess
    myssquared <- ssquared(out, sigma, residual)[length(xguess)]
  }
  
  #compute hessian
  hessianij<-NULL
  myhessian <- unlist(lapply(1:length(Thend), function(i) {
    lapply(1:i, function(j) {
      # Variate the parameters
      inewyend      <- jnewyend     <- ijnewyend     <- yend
      inewyend[i]   <- inewyend[i] + stepwidth
      jnewyend[j]   <- jnewyend[j] + stepwidth
      ijnewyend[i]  <- yend[i] + stepwidth
      ijnewyend[j]  <- ijnewyend[j] + stepwidth
      
      
      inewout       <-  bvptwpC(
        x = x,
        func = eqns,
        yini = yini,
        yend = inewyend,
        parms = pars,
        nmax = 50000,
        xguess = xguess,
        yguess = yguess
      )
      inewssquared    <- ssquared(inewout, sigma, residual)[nrow(inewout)]
      
      
      ijnewout        <-  bvptwpC(
        x = x,
        func = eqns,
        yini = yini,
        yend = ijnewyend,
        parms = pars,
        nmax = 50000,
        xguess = xguess,
        yguess = yguess
      )
      ijnewssquared <- ssquared(ijnewout, sigma, residual)[nrow(ijnewout)]
      
      
      
      if(i!=j) {
        jnewout       <-  bvptwpC(
          x = x,
          func = eqns,
          yini = yini,
          yend = jnewyend,
          parms = pars,
          nmax = 50000,
          xguess = xguess,
          yguess = yguess
        )
        jnewssquared <- ssquared(jnewout, sigma, residual)[nrow(newout)]
        # Take finite difference
        hessianij <- (ijnewssquared-inewssquared-jnewssquared+myssquared)/(stepwidth^2)
        names(hessianij) <- paste(as.character(i), as.character(j), sep=".")
      } else {
        # Take finite difference
        hessianij <- (ijnewssquared-2*inewssquared+myssquared)/(stepwidth^2)
        names(hessianij) <- paste(as.character(i), as.character(j), sep=".")
      }
      
      
      return(hessianij)
    })
  }) )
  
  
  
  myhs                             <- diag(nrow = length(bestfit))
  myhs[lower.tri(myhs, diag = T)] <- myhessian
  myhs                             <- myhs + t(myhs) - diag(diag(myhs))
  
  dimnames(myhs) <- list(names(bestfit), names(bestfit))
  
  return(myhs)
}


#' Objective function for the bvp-solutions
#' 
#' Objective function for the approach with the approximated metric. The value is the integrated arclength computed by ssquared(), gradient and hessian are computed by finite differences.
#' Takes really long. Residual has to be given to the function because else one has two conflicting DLLs.
#' 
#' 
#' @param Thend
#' @param out
#' @param sigma
#' @param bestfit
#' @param residual
#' @param stepwidth
#' 
#' 
#' @export
ssquaredobj <- function(Thend, fixed=NULL, sigma=0.02, bestfit, residual) {  
  
  
  out       <- geodesicbvp(Thend = Thend, bestfit)
  
  value     <- ssquared(out = out, sigma=sigma, residual = residual)[nrow(out)]
  gradient  <- gradssfinite(Thend = Thend, out = out, sigma = sigma, bestfit = bestfit, residual = residual, stepwidth = 0.001)
  hessian   <- hessianssfinite(Thend = Thend, out = out, sigma = sigma, bestfit = bestfit, residual = residual, stepwidth = 0.001)
  
  out                     <- list(value = value, gradient = gradient, hessian = hessian)
  attr(out, "class")      <- "objlist"
  attr(out, "valueData")  <- value
  attr(out, "valuePrior") <- 0
  
  return(out)
}

#' Get the difference of initial and end velocity to see how much they change.
#' This is to evaluate whether the approximation in RNC works well.
#' Actually this is nonsense, because gmunuvmuvnu=const along the geodesic. Thus RNCssquared is a perfect measure.
#' 
#' @param Thend 
#' @param bestfit 
#' 
#' 
#' 
#' @export
velocitydifference <- function(Thend, bestfit) {
  out <- geodesicbvp(Thend = Thend, bestfit = bestfit)
  vini <- out[1, -(1:(length(Thend)+1))]
  vend <- out[nrow(out), -(1:(length(Thend)+1))]
  return(vend-vini)
}
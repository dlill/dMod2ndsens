#' Translate geodesic equation with constant christoffels into first order differential equation
#' 
#' @param bestfit The point of best fit.
#' @details A function resi() to evaluate the residuals at a point must be loaded. 
#' The function computes automatically the christoffel symbols at the point of best fit and then translates it into a geodesic equation.
#' 
#' @export
firstordergeodesiceqn <- function(bestfit, chris = NULL) {
  
  mychristoffels <- NULL
  if(is.null(chris))  { 
    myresidual      <- resi(bestfit)
    mychristoffels  <- christoffels(myresidual)
  } else {mychristoffels <- chris}
  
  # set z=Thdot.
  # generate the zsymb-vector
  z           <- t(outer("z", 1:length(bestfit), paste, sep=""))
  zdot        <- paste("-1*(",as.character(prodSymb(apply(mychristoffels, 2, function(x) prodSymb(x, z)), z)),")" )
  names(zdot) <- z
  
  Thdot         <- z
  names(Thdot)  <- names(bestfit)
  
  equations     <- as.eqnvec(c(Thdot, zdot))
  return(equations)
}


#' Boundaries for geodesic equation
#' @description  Define boundaries for the geodesic equation to be solved with bvptwpC. These boundaries are to be supplied to funC.
#' @param bestfit The starting point
#' @param Thend The end-point
#' 
#' 
#' @export
geodesicboundaries <- function(bestfit, Thend) {
  return(data.frame(name=c(names(equations) ),#  , "dummy1"),            # insert the dummy that needs to be introduced because of the sensitivities
                    yini = c(bestfit, rep(NA, length(bestfit)) ),# ,0), # set dummy-value to 0 
                    yend = c(Thend, rep(NA, length(bestfit))   ))# ,0))   # set dummy-value to 0
         )
}


#' Get the yguesses for neighboured points
#' @param out The output of the numerical solution of the geodesic bvp
#' 
#' @export
yguessout <- function(out) {
  if(is.null(out)) return(NULL)
  else return(t(out[,-1]))
}



#' Solve the geodesic bvp by only specifying Thend and bestfit
#' 
#' @param Thend
#' @param bestfit
#' @details needs "eqns", where eqns is the compiled odemodel, to be evaluated
#' 
#' @export
geodesicbvp <- function(Thend, bestfit, lobatto = FALSE) {
  
  mybd        <- geodesicboundaries(bestfit, Thend)
  yini        <- mybd[1:length(bestfit),2]
  names(yini) <- names(bestfit)
  yend        <- mybd[1:length(bestfit),3]
  names(yend) <- names(bestfit)
  
  # eqns<- funC(equations, jacobian = "full", boundary = mybd)
  
  pars        <- NULL
  x           <- seq(0, 1)
  
  out <- bvptwpC(yini=yini, x = x, func = eqns, yend=yend, parms = pars, lobatto = lobatto)
  return(out)
}





#' Title
#'
#' @param out
#' @param Thend
#' @param bestfit
#' @param sigma
#' @param residual
#'
#' @return
#' @export
#'
#' @examples
RNCdataspacecorrection <- function(out, Thend, bestfit, sigma, residual) {

  # Get the residual vector r
  r     <- residual$residual/sigma

  # get derivatives
  myderiv  <- data.matrix(attr(residual, "deriv")[,-c(1,2)])

  mysderiv <- data.matrix(attr(residual, "sderiv")[,-c(1,2)])
  mysderiv <- array(mysderiv, dim = c(nrow(myderiv), ncol(myderiv), ncol(myderiv))) #array with indices (r, mu, nu)

  # construct the projection operator PN
  S   <- svd(myderiv)
  PN  <- diag(nrow = nrow(S$u)) - S$u%*%t(S$u)

  # get the velocities, ie the coordinates p.
  vini      <- out[1, -(1:(length(Thend)+1))]

  # multiply the second derivatives with the initivial velocities according to formula 26
  mysderiv <- apply(mysderiv, 2, "%*%" , vini) # multiply 2nd index
  mysderiv <- mysderiv %*% vini

  # apply the projection operator onto r and sderiv*v*v. correction =  r%*%PN%*%(sderiv*v*v)
  r <- r%*%PN%*%mysderiv

  return(r)
}




#' Get the chi^2 value in RNC but with vini(Thend) obtained with the approximate metric.
#' 
#' chi^2 = chi0+ pmu pnu(gmunu0+r_m %*% PN_mn %*% d^2_munu r)
#' 
#' This formular corresponds to the formulat 26 of the Transtrum paper and accounts also for the extrinsic curvature. 
#' In this formula Tau=1 and the factor of 1/2 is killed by my other definition of the cost, which doesn't include the 1/2.
#' 
#' Since the extrinsic curvature part is not quadratic, this formula is only an approximation up to second order in the geodesic parameter Tau.
#' See the sheet "corrections of DeltaChi^2 as predicted by s^2" from the 1st of June, 2016.
#' 
#' @param out
#' @param Thend
#' @param bestfit
#' @param sigma
#' @param residual
#' 
#' 
#' @export
RNCssquared <- function(out, Thend, bestfit, sigma = 1, residual, corr = FALSE)  {
  
  # compute the quadratic contribution of the geodesic.
  mysigmas <- residual$sigma
  
  mymetric  <- data.matrix(attr(residual, "deriv")[,-c(1,2)])/mysigmas
  mymetric  <- t(mymetric)%*%mymetric
  

  
  vini      <- out[1, -(1:(length(Thend)+1))]
  
  quadratic      <- as.vector(t(vini)%*%mymetric%*%vini)
  
  # compute the correction that accounts for extrinsic curvature
  correction <- 0
  if(corr)  correction <- RNCdataspacecorrection(out, Thend, bestfit, sigma, residual)

  # constant part, the offset
  chi0      <- (residual$residual/mysigmas)%*%(residual$residual/mysigmas)
  
  return(quadratic + correction + chi0) #gmunu.p.p
}





#' Jacobian via sensitivities
#' 
#' Compute the Jacobian dTh/dp via sensitivity equations that are then inverted to obtain dp/dTh
#' 
#' 
#' @param Thend 
#' @param out 
#' @param sigma 
#' @param bestfit 
#' @param residual 
#'
#' @export
RNCjacsssens <- function(Thend, out, sigma, bestfit, residual) {
  
  # get initial values to integrate the GE.
  pars      <- c(out[1,-1], dummy1 = 0, dummy2 = 0)
  times     <- 0:1
  

  mypred    <- GEivp(0:1, pars = pars, deriv = 1)

  myderiv   <- attr(mypred, "deriv")
  mysnames  <- dimnames(myderiv)[[2]]
  
  mynames   <- names(bestfit) 
  n         <- length(bestfit)
  z         <- as.character(outer("z", 1:n, paste, sep=""))
  mynames   <- as.character(outer(mynames, z, paste, sep="."))
  
  myderiv   <- myderiv[2,mysnames%in%mynames] # take only the derivatives of "Th_i.z_j" at timepoint t = 1
  myderiv   <- matrix(myderiv, nrow = n)
  dimnames(myderiv) <- list(names(bestfit), z)
  
  mysvd <- svd(myderiv)
  myjac <- mysvd$v%*%diag(1/(mysvd$d))%*%t(mysvd$u) # invert to get dp/dTh
  

  
  return(myjac)
}






#' Grad in approximate metric approach
#' RNCgradssf
#' 
#' @export
RNCgradss <- function(Thend, out, sigma, bestfit, residual, stepwidth, jacpTh) {  
  
  vini      <- out[1,-(1:(length(Thend)+1))]
  mysigmas <- residual$sigma
  
  mymetric  <- data.matrix(attr(residual, "deriv")[,-c(1,2)])/mysigmas
  mymetric  <- t(mymetric)%*%mymetric
  
  mygrad <- 2*as.vector(t(vini)%*%mymetric%*%jacpTh)
  
  
  
  names(mygrad) <- names(Thend)
  return(mygrad)
}

#' Hessian dobj(p)/dTh^2 in approximate metric approach
#' 
#' @export
RNChessianss <- function(Thend, out, sigma = 1, bestfit, residual, stepwidth, jacpTh) {
  
  vini      <- out[1,-(1:(length(Thend)+1))]
  
  mysigmas <- residual$sigma
  mymetric  <- data.matrix(attr(residual, "deriv")[,-c(1,2)])/  mysigmas 
  mymetric  <- t(mymetric)%*%mymetric
  
  myhessian <- 2*(t(jacpTh)%*%mymetric%*%jacpTh)# +p*d^2p/dTh^2
  dimnames(myhessian) <- list(names(Thend), names(Thend))
  return(myhessian)
}



#' Title
#'
#' @param Thend 
#'
#' @return
#' @export
#'
#' @examples
RNCssquaredobj <- function(Thend, fixed = NULL, with_prior = FALSE){#}, sigma = 1){ #}, fixed=NULL, sigma=0.02, bestfit = bestfit, residual = residual) {  
  # solve bvp for vini
  loadDLL(eqns) 
  Thend <- c(Thend, fixed) # Achtung, Thend und fixed dürfen nicht überlappen
  Thend <- Thend[names(bestfit)]
  out       <- geodesicbvp(Thend = Thend, bestfit = bestfit)
  # solve ivp for sensitivities
  loadDLL(GEodemodel$extended)
  jacpTh    <- RNCjacsssens(Thend = Thend , out = out, sigma = 1, bestfit = bestfit, residual = residual)
  
  value     <- RNCssquared(out = out, Thend = Thend, sigma = 1, residual = residual, corr = FALSE)
  # value     <- ssquared(out = out, sigma = sigma, residual = residual)
  gradient  <- RNCgradss(Thend = Thend, out = out, sigma = 1, bestfit = bestfit, residual = residual, stepwidth = 0.001, jacpTh = jacpTh)
  hessian   <- RNChessianss(Thend = Thend, out = out, sigma = 1, bestfit = bestfit, residual = residual, stepwidth = 0.001, jacpTh = jacpTh)
  
  myobjlist                     <- list(value = value , gradient = gradient, hessian = hessian)
  attr(myobjlist, "class")      <- c("list", "objlist")
  if(with_prior) { 
    priorsig <- 10
    out.prior <-  constraintL2(Thend, prior, sigma = priorsig)
    myobjlist <- myobjlist  + out.prior
    attr(myobjlist, "valuePrior") <- out.prior$value
    
  } else {
    attr(myobjlist, "valuePrior") <- 0  
  } 

    
    if (!is.null(fixed)) {
      myobjlist$gradient <- myobjlist$gradient[!names(myobjlist$gradient) %in% names(fixed)]
      myobjlist$hessian <- myobjlist$hessian[!rownames(myobjlist$hessian) %in% names(fixed), !colnames(myobjlist$hessian) %in% names(fixed)]
    }
    
    attr(myobjlist, "valueData")  <- value
    # attr(myobjlist, "vini")       <- out[1,-(1:(length(Thend)+1))]


  return(myobjlist)
}
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
# for models with analytically given christoffel symbols. basically everything that changes, needs a "true_" before
# true_firstordergeodesiceqn
# true_eqns
# true_geodesicbvp
# true_RNCjacsssens
# true_GEivp
# true_RNCssquaredobj
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################

#' Create the geodesic eqn odemodel for analytically known christoffel symbols
#'
#' @param ana.chris The analytically knonw christoffel symbols
#'
#' @return
#' @export
true_firstordergeodesiceqn <- function(bestfit, ana.chris) {
  
  mychristoffels  <- ana.chris
  
  # set z=Thdot.
  # generate the zsymb-vector
  z           <- t(outer("z", 1:length(bestfit), paste, sep=""))
  zdot        <- paste("-1*(",as.character(prodSymb(apply(mychristoffels, 2, function(x) prodSymb(x, z)), z)),")" )
  names(zdot) <- z
  
  Thdot         <- z
  names(Thdot)  <- names(bestfit)
  
  equations     <- as.eqnvec(c(Thdot, zdot))
  return(equations)
}

#' Solve the geodesic bvp by only specifying Thend and bestfit
#' 
#' @param Thend
#' @param bestfit
#' @details needs "eqns", where eqns is the compiled odemodel, to be evaluated
#' 
#' @export
true_geodesicbvp <- function(Thend, bestfit, out_guess=NULL, ...) {
  
  mybd        <- geodesicboundaries(bestfit, Thend)
  yini        <- mybd[1:length(bestfit),2]
  names(yini) <- names(bestfit)
  yend        <- mybd[1:length(bestfit),3]
  names(yend) <- names(bestfit)
  
  # eqns<- funC(equations, jacobian = "full", boundary = mybd)
  
  pars        <- NULL
  x           <- seq(0, 1, by = 0.1)
  
  out <- bvptwpC(yini=yini, x = x, func = true_eqns, yend=yend, parms = pars, xguess = out_guess[,1] , yguess = yguessout(out_guess), ...)
  return(out)
}


#' Title
#'
#' @param Thend 
#' @param out 
#' @param sigma 
#' @param bestfit 
#' @param residual 
#'
#' @return
#' @export
true_RNCjacsssens <- function(Thend, out, sigma, bestfit, residual) {
  
  # get initial values to integrate the GE.
  pars      <- c(out[1,-1], dummy1 = 0, dummy2 = 0)
  times     <- 0:1
  
  
  mypred    <- true_GEivp(0:1, pars = pars, deriv = 1)
  
  myderiv   <- attr(mypred, "deriv")
  mysnames  <- dimnames(myderiv)[[2]]
  
  mynames   <- names(bestfit) 
  n         <- length(bestfit)
  z         <- as.character(outer("z", 1:n, paste, sep=""))
  mynames   <- as.character(outer(mynames, z, paste, sep="."))
  
  myderiv   <- myderiv[2,mysnames%in%mynames] # take only the derivatives of "Th_i.z_j" at timepoint t = 1
  myderiv   <- matrix(myderiv, nrow = n)
  dimnames(myderiv) <- list(names(bestfit), z)
  
  mysvd <- svd(myderiv)
  myjac <- mysvd$v%*%diag(1/(mysvd$d))%*%t(mysvd$u) # invert to get dp/dTh
  
  
  
  return(myjac)
}


#' Title
#'
#' @param Thend 
#'
#' @return
#' @export
true_RNCssquaredobj <- function(Thend, fixed = NULL, out_guess = NULL){#}, sigma = 1){ #}, fixed=NULL, sigma=0.02, bestfit = bestfit, residual = residual) {  
  # solve bvp for vini
  loadDLL(true_eqns) 
  Thend <- c(Thend, fixed) # Achtung, Thend und fixed dürfen nicht überlappen
  Thend <- Thend[names(bestfit)]
  out       <- true_geodesicbvp(Thend = Thend, bestfit = bestfit, out = out_guess)
  # solve ivp for sensitivities
  loadDLL(true_GEodemodel$extended)
  jacpTh    <- true_RNCjacsssens(Thend = Thend , out = out, sigma = sigma, bestfit = bestfit, residual = residual)
  
  value     <- RNCssquared(out = out, Thend = Thend, sigma = sigma, residual = residual, corr = FALSE)
  # value     <- ssquared(out = out, sigma = sigma, residual = residual)
  gradient  <- RNCgradss(Thend = Thend, out = out, sigma = sigma, bestfit = bestfit, residual = residual, stepwidth = 0.001, jacpTh = jacpTh)
  hessian   <- RNChessianss(Thend = Thend, out = out, sigma = sigma, bestfit = bestfit, residual = residual, stepwidth = 0.001, jacpTh = jacpTh)
  
  myobjlist                     <- list(value = value , gradient = gradient, hessian = hessian)
  attr(myobjlist, "class")      <- "objlist"
  if(with_prior) { 
    out.prior <-  constraintL2(Thend, prior, sigma = 10)
    myobjlist <- myobjlist  + out.prior
    
    if (!is.null(fixed)) {
      myobjlist$gradient <- myobjlist$gradient[!names(myobjlist$gradient) %in% names(fixed)]
      myobjlist$hessian <- myobjlist$hessian[!rownames(myobjlist$hessian) %in% names(fixed), !colnames(myobjlist$hessian) %in% names(fixed)]
    }
    
    attr(myobjlist, "valueData")  <- value
    attr(myobjlist, "valuePrior") <- out.prior$value
    # attr(myobjlist, "vini")       <- out[1,-(1:(length(Thend)+1))]
  } else {
    out.prior <-  list(value = 0, gradient = 0, hessian = 0)    
    class(out.prior)<- c("objlist", "list")
    
    myobjlist <- myobjlist  #+ out.prior
    
    if (!is.null(fixed)) {
      myobjlist$gradient <- myobjlist$gradient[!names(myobjlist$gradient) %in% names(fixed)]
      myobjlist$hessian <- myobjlist$hessian[!rownames(myobjlist$hessian) %in% names(fixed), !colnames(myobjlist$hessian) %in% names(fixed)]
    }
    
    attr(myobjlist, "valueData")  <- value
    attr(myobjlist, "valuePrior") <- out.prior$value
    # attr(myobjlist, "vini")       <- out[1,-(1:(length(Thend)+1))]
  }
  
  return(myobjlist)
}

###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################

#' Comparison of the results of obj() and of the approximation
#' 
#' This compares the values of the real objective function to the integration of the arclength of one geodesic that is the solution of the geodesic bvp
#' 
#' @param bestfit
#' @param Thend
#' @param sigma
#' 
#' 
#' @export
compareobj<- function(bestfit, Thend, sigma=0.02) {
  
  residual <- resi(bestfit)
  equations   <- firstordergeodesiceqn(bestfit)
  mybd        <- geodesicboundaries(bestfit, Thend)
  eqns        <- funC(equations, jacobian = "full", boundary = mybd)
  
  times       <- seq(0, 1, by = 1)
  pars        <- NULL
  
  out         <- bvptwpC(x = times, func = eqns, parms = pars)#, xguess = xguess, yguess = yguess, nmax=1000)
  
  # The real objective function
  myobj       <- do.call(c, lapply(1:nrow(out), function(i) obj(out[i,1:length(bestfit)+1])$value))
  
  # The arclength of the geodesic path
  sdot    <-sdotfun(bestfit)
  ssquare <-ssquared(out, sigma = sigma, residual = residual)
  
  RNCssquare <- RNCssquared(out, Thend, bestfit, sigma=sigma, residual)
  return(cbind(out, obj=myobj, ssquared = ssquare, RNCssquared = RNCssquare))
}



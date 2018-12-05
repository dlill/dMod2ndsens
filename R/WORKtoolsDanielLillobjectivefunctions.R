#' Objective function Theta(p) in RNC - Riemann Normal Coordinates
#' 
#' @description
#' This is a function returning the objective function that approximates the real objective function via transforming Theta to p and then squaring the absolute value of the p-vector.
#' A residual function resi() for the residuals has to be defined before using it.
#' 
#' @param sigma The estimated error of the data, to normalize the residuals.
#' @param bestfit The point of best fit around which the RNC (ie. p-coordinates) are constructed
#' @param rootsolve Boolean. Determines if the Coordinate trafo p(Theta) should be done via the RNC or by inverting the 2nd order geodesic path.
#' 
#' @param Theta Argument of the returned function.
#' 
#' @export
geodesicobjfun <- function(sigma = 1 , bestfit, rootsolve = FALSE) {
  myresi <- resi(bestfit)
  mysigmas <- residual$sigma
  myderiv <-  data.matrix(attr(myresi, "deriv")[,-c(1,2)]) /mysigmas
  mymetric <- 2*t(myderiv)%*%myderiv
  
 
  
  Thetaobj <- NULL
  
  if(!rootsolve) {
    
    # p<-pfun(bestfit, residual,christoffelsymbols = christoffels(myresi))
    Thetaobj <- function(Theta){
      
      mynames <- names(Theta)
      # Theta <- matrix(Theta, nrow = 1)
      myp <- p(Theta)
      # myJTp <- JTp(myp)
      myJpT <- JpT(Theta)

      # angle between the proposed step and the radial vector

      
      
      value <- 1/(2)*as.numeric(myp%*%mymetric%*%myp)+as.vector(t(myresi[,6])%*%myresi[,6])
      
      gradient <- as.vector(myp%*%mymetric%*%myJpT)
      names(gradient) <- mynames
      
      hessian <- t(myJpT)%*%mymetric%*%myJpT  +apply(christoffels(myresi), c(2,3), function(i) 2*myp%*%mymetric%*%i)
      
      out <- list(value=value, gradient=gradient, hessian=hessian)
      attr(out, "class")<- c("objlist", "list")
      attr(out, "prior") <- 0
      return(out)
    }
    
  } else {
    
    # pnum<-pfunnum(bestfit, residual,christoffelsymbols = christoffels(myresi))
    Thetaobj <- function(Theta){
      
      mynames <- names(Theta)
      # Theta <- matrix(Theta, nrow = 1)
      myp <- pnum(Theta)
      myJTp <- JTp(myp)
      # myJpT <- JpT(Theta)
      myJpT <- svd(myJTp)
      myJpT <- myJpT$v%*%diag(1/myJpT$d)%*%t(myJpT$u)
      # bfobj <- obj(bestfit, deriv = 0)$value
      # value <- as.numeric(myp%*%t(t(mymetric%*%myJTp)%*%myJTp)%*%t(myp)/(sigma*sigma)) #+ bfobj
      
      
      value <- 1/(2)*as.numeric(myp%*%mymetric%*%myp)+as.vector(t(myresi[,6])%*%myresi[,6])
      
      gradient <- as.vector(myp%*%mymetric%*%myJpT)
      names(gradient) <- mynames
      
      # hessian <- 2*apply(christoffels(myresi), c(2,3), sum)/(sigma*sigma)
      hessian <- t(myJpT)%*%mymetric%*%myJpT # +höhere ableitungen von p nach theta. schwierig aus der inversen koordinatentrafo zu bekommen.
      out <- list(value=value, gradient=gradient, hessian=hessian)
      attr(out, "class")<- c("objlist", "list")
      attr(out, "prior") <- 0
      return(out)
    }
    
  }
  return(Thetaobj)
}

#' (True) Objective function in Theta-Space for analytical functions
#' 
#' @description The objective function for functions that are analytically known without data points, like the Rosenbrock function
#' 
#' @param Theta The point at which the objective function shall be evaluated.
#' 
#' 
#' @export
Thetaobj <- function(Theta, fixed = NULL){
  if(!is.null(fixed)) Theta <- c(Theta, fixed)
  mynames <- names(Theta)
  residual <- resi(Theta)
  value <- as.numeric(t(residual$residual)%*%residual$residual) #+bfobj=0 in diesem Falle.
  gradient <- as.vector(2*t(data.matrix(attr(residual, "deriv")[,-c(1,2)]))%*%residual$residual)
  names(gradient) <- mynames
  hessian <- 2*t(data.matrix(attr(residual, "deriv")[,-c(1,2)]))%*%data.matrix(attr(residual, "deriv")[,-c(1,2)])#+Terme mit 2.Ableitungen
  out <- list(value=value, gradient=gradient, hessian=hessian)
  attr(out, "class")<- c("objlist", "list")
  attr(out, "prior") <- 0
  return(out)
}


#' Quadratic approximation to the Log Likelihood
#' 
#' 
#' @export
hessianobjfun <- function(sigma = 1 , bestfit, residual = NULL) {
  
  if (!is.null(residual))
    residual <- combineresiduals(resi(bestfit))
  
  myresidual <- as.vector(residual$residual)
  mysigmas <- residual$sigma
  myderiv <- data.matrix(attr(residual, "deriv")[,-c(1,2)]) / mysigmas
  mysderiv <- data.matrix(attr(residual, "sderiv")[,-c(1,2)]) /mysigmas
  
  gmunu   <- t(myderiv)%*%myderiv
  rsderiv <- matrix((myresidual/mysigmas)%*%mysderiv, nrow=length(bestfit))
  
  Thetaobj<-function(Theta, fixed = NULL) {
    Theta <- c(Theta, fixed)
    DeltaTheta <- Theta-bestfit
    
    value     <- ((myresidual/mysigmas)%*%(myresidual/mysigmas) + 2*(myresidual/mysigmas)%*%myderiv%*%DeltaTheta + DeltaTheta%*%gmunu%*%DeltaTheta + DeltaTheta%*%rsderiv%*%DeltaTheta)
    gradient  <- as.numeric(2*(as.vector(gmunu%*%DeltaTheta) + (myresidual/mysigmas)%*%myderiv + as.vector(rsderiv%*%DeltaTheta)))
    names(gradient) <- names(bestfit)
    hessian   <- 2*(gmunu + rsderiv)
    dimnames(hessian) <- dimnames(gmunu)
    
    out <- list(value=value, gradient=gradient, hessian=hessian)
    attr(out, "class")<- c("objlist", "list")
    attr(out, "prior") <- 0
    return(out)
  }
  return(Thetaobj)
}



#' Objective function with gmunu as an approximation to Hmunu
#'
#' @param sigma 
#' @param bestfit 
#'
#' @return
#' @export
#'
#' @examples
gmunuobjfun <- function(sigma = 1 , bestfit, residual = NULL) {
  
  if (is.null(residual))
    residual <- combineresiduals(resi(bestfit))
  
  myresidual <- as.vector(residual$residual)
  mysigmas <- residual$sigma
  myderiv <- data.matrix(attr(residual, "deriv")[,-c(1,2)])/mysigmas
  mysderiv <- data.matrix(attr(residual, "sderiv")[,-c(1,2)])
  
  gmunu   <- t(myderiv)%*%myderiv
  
  gobj<-function(Theta, fixed = NULL) {
    Theta <- c(Theta, fixed)
    # rsderiv <- matrix(myresidual%*%mysderiv, nrow=length(Theta))
    DeltaTheta <- Theta-bestfit
    
    value     <- 1*((myresidual/mysigmas)%*%(myresidual/mysigmas) + 2*(myresidual/mysigmas)%*%myderiv%*%DeltaTheta + DeltaTheta%*%gmunu%*%DeltaTheta) # + DeltaTheta%*%rsderiv%*%DeltaTheta)
    gradient  <- 2*((myresidual/mysigmas)%*%myderiv + as.vector(gmunu%*%DeltaTheta)) # + as.vector(rsderiv%*%DeltaTheta))
    hessian   <- 2*(gmunu) # + rsderiv)
    dimnames(hessian) <- dimnames(gmunu)
    
    out <- list(value=value, gradient=gradient, hessian=hessian)
    attr(out, "class")<- c("objlist", "list")
    attr(out, "prior") <- 0
    return(out)
  }
  return(gobj)
}


#' The normal objective function of dMod description
#' 
#' Evaulate like obj<-objfun(priorsigma)
#' 
#' @return obj(pouter, fixed=NULL, deriv=2)
#' @export
#'
#' @examples
objfun <- function(priorsigma){
  obj <- function(pouter, fixed=NULL, deriv=2) {
    
    prediction <- x(timesD, pouter, fixed = fixed, deriv = deriv) #pred(pinner), deriv= pred(pinner).pouter, sderiv= pred(pinner).pouter.pouter,
    
    # Apply res() and wrss() to compute residuals and the weighted residual sum of squares
    out.data <- lapply(names(data), function(cn) wrss(res(data[[cn]], prediction[[cn]]))) #res(obs(pinner)) => out.data(pinner)
    out.data <- Reduce("+", out.data)
    
    # -------------------------------------------------------------------------------------------------------------------------------
    if(deriv == 0) {
      out.data$gradient <- 0
      out.data$hessian <- 0
    }
    # -------------------------------------------------------------------------------------------------------------------------------
    
    # Working with weak prior (helps avoiding runaway solutions of the optimization problem)
    out.prior <-  constraintL2(pouter, prior, sigma = priorsigma)
    
    if(deriv == 0) {
      out.data$gradient <- out.prior$gradient <- 0
      out.data$hessian <- out.prior$hessian <- 0
    }
    out <- out.prior + out.data
    
    # Comment in if you want to see the contribution of prior and data
    # e.g. in profile() and plotProfile()
    attr(out, "valueData") <- out.data$value
    attr(out, "valuePrior") <- out.prior$value
    
    return(out)
    
  }
  return(obj)
}








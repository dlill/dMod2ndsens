#' Refill symmetric second sensitivities
#' 
#' @description To save time, derivatives that behave symmetrical in the derivative indices are only computed once.
#' Thus one has to replace the ones that were set to zero by their actual value.
#' @param myssensitivities The second sensitivities as given in P2X of the prediction function Xs
#' @param dv number of variables
#' @param ds number of states
#' @param dp number of parameters
#' @export
srefill <- function(myssensitivities = NULL, variables = 0, ssvariables.2=0, ssparameters.2=0) {
  
  dv <- length(variables)
  ds <- length(ssvariables.2)
  dp <- length(ssparameters.2)
  
  myssens<- submatrix(myssensitivities, cols = -1) #entfernt "times"
  mysnames<-colnames(myssens)
  mysnamesarray <- array(mysnames, dim = c(dv,dp+ds,dp+ds))
  
  indexmat <- matrix(1:((ds+dp)^2), ncol=(ds+dp)) #erstelle matrix, die die indizes beinhaltet
  indexmat[upper.tri(indexmat)] <- t(indexmat)[upper.tri(indexmat)] #mache die matrix so symmetrisch, wie es für srefill benötigt ist, also upper.tri <- lower.tri
  indexnames <- t(apply(mysnamesarray, 1, function(x) #wende die indexmat auf die names an und transponiere diesen array so, dass die originalindexreihenfolge wiederhergestellt ist.
    x[as.vector(indexmat)]
  ))
  
  myssens2 <- myssens[,indexnames] #nehme die Spalten, die vorer rausgefunden wurden
  myssensitivities <- cbind(matrix(myssensitivities[,1], ncol = 1), matrix(myssens2, ncol = (ncol(myssensitivities)-1)))
  colnames(myssensitivities) <- c("time", colnames(myssens))
  
  return(myssensitivities)
}


#' Reduce the output of the christoffel symbols by its symmetries
#'
#' @param chris An array of christoffel symbols
#'
#' @return
#' @export
#'
#' @examples
#' chris <- christoffels(resi(pouter))
reducechristoffels <- function(chris) {
  
  n         <- dim(chris)[1]
  mylower   <- lower.tri(matrix(NA, nrow= n, ncol = n), diag = TRUE)
  namesmat  <- dimnames(chris)[[2]]
  namesmat  <- outer(namesmat, namesmat, paste, sep=".")
  namesmat  <- as.character(namesmat[mylower])
  
  chris <- apply(chris, 1, function(m) {
    as.numeric(m[mylower])
  })
  
  rownames(chris) <- namesmat
  return(chris)
}


#' combine a list of residuals of different conditions
#' 
#' @param residual A list of residuals where entries of the list correspond to an experimental condition
#'
#' @export
combineresiduals <- function(residual){
  if(is.null(residual$time)) {
    residualvalues <- do.call(rbind, residual)
    
    deriv.data <- lapply(1:length(residual), function(cn) {deriv <- attr(residual[[cn]], "deriv")}) # hier noch, falls manche spalten durch paratrafos rausfallen, sollten, diese spalten wieder einfügen
    deriv.data <- do.call(rbind, deriv.data)
    
    sderiv.data <- lapply(1:length(residual), function(cn) {sderiv <- attr(residual[[cn]], "sderiv")})# hier noch, falls manche spalten durch paratrafos rausfallen, sollten, diese spalten wieder einfügen
    sderiv.data <- do.call(rbind, sderiv.data)  
    
    # Prior the jacobians
    residual<-objframe(residualvalues, deriv = deriv.data, sderiv = sderiv.data)
  }
  return(residual)
}



#' Function for residual at point pouter
#' 
#' @details Returns a function of pouter, to compute the residual at that point.
#' @param data To get the names of the conditions.
#' @param timesD The times at which the data points are taken
#' @param fixed Fixed (whatever that means :D)
#' @param deriv (normally 2)
#' @param pouter PARAMETER OF THE FUNCTION RETURNED. The point at which the residual shall be evaluated.
#' @export
resfun <- function(data, timesD, bestfit, fixed = fixed, deriv = 2) {
  return(function(pouter){
    residual <- lapply(names(data), function(cn) res(data[[cn]], x(timesD, pouter, fixed = fixed, deriv = 2)[[cn]])) #create list of residuals of all conditions
    combineresiduals(residual)
  })
}


#' Test Cost function for non-quadratic behaviour
#' 
#' @description If the cost function shows quadratic behaviour around a point, a coordinate transformation towards geodesic coordinates can be found.
#'  This tool is supposed to analyse the cost function and tell the user if a coordinate trafo is reasonable to execute.
#' @param residual The residuals-objframe returned by dMod::res, including first and second derivatives
#' @param vmu,vnu Direction vectors in which the nonquadratic parts shall be investigated. The length has to be the same as the number of parameters
#' 
#' @export
nonquadparttest <- function(residual){
  
  # get residuals
  myresidual <- residual[,"residual"]
  # mytimes <- residual[,1]
  # mynames <- residual[,2]
 
  # get 1st and 2nd derivatives
  myderiv <- data.matrix(attr(residual, "deriv")[,-c(1,2)])
  mysderiv <- data.matrix(attr(residual, "sderiv")[,-c(1,2)])

  
  # Construct the projection operator PN, which projects data (residual) vectors to the space
  # normal to the tangent space of the model manifold. PN= 1 - PT
  # consider that the matrix U is expressed in the basis of the svd:
  # the svd returns something like svd(myderiv) = list(u,d,v) with myderiv = u %*% d %*% t(v)
  # where v is a unitary matrix with the columns containing the eigenvectors of J^T %*% J
  # If everything is regular, the projection operator looks like PN = diag(0,...0,1,...,1)
  # where the length of the 0's is length(pars), and the length of the 1's is m-length(pars)
  
  # Execute the svd
  S <- svd(myderiv)
  d <- S$d
  u <- S$u
  v <- S$v
  
  # # Delete dimensions with too small eigenvalues.
  # tol <-10^(-8)
  # index <- which(abs(d) > max(dim(myderiv))*max(d)*tol) #take only dimensions with singular values bigger than a given tolerance tol, this is copied from statistics.R, pseudoinverse()
  # index <- diag((1:length(d))%in%index) 
  # u <- u%*%index # this is not correct, because v transforms the vectors into eigendimensions of u.delete the dimensions with small singular values in u respectively. Doesn't change output of the svd, nor the idempotency of u*t(u)
  
  # Construct the projection operator PN and then test if parts of the manifold have high extrinsic curvature
  PN<-diag(nrow = nrow(u)) - u%*%t(u)
  if(!any(eigen(PN)$values>0.9)) { #arbitrary border, since all eigenvalues are either very close to 1 or 0. 
    # If there are no eigenvalues that are 1, the normal space to the manifold is zero, therefore no extrinsic curvature exists.
    out <- "non quadratic terms are negligible, PN == 0"
  } else {
    # compute the term d²C/dtau² which is responsible for the non-quadratic behaviour
    mysderivarray <- array(mysderiv, dim=c(nrow(mysderiv), dim(myderiv)[2], dim(myderiv)[2])) #first index switches between the rm, the second and third go through the parameter indices
    mysderivvec <- apply(mysderivarray, 1, sum) #executes d_mu d_nu r_m*v^mu v^nu, as v^mu v^nu is just a matrix full of 1s
    mysderivvec <- PN%*%mysderivvec #Apply the projection operator "wich part of the manifold will be bent out of the tangent plane?"
    nonquadmagnitude <- sum(myresidual*as.vector(mysderivvec)) # compute the scalar product with r_m, as in the formula
    
    # now the nonquadmagnitude must be compared to gmunu vmu vnu, aka "quadmagnitude"
    quadmagnitude <- sum(t(myderiv)%*%myderiv)
    
    
    out <-as.character(signif(abs(nonquadmagnitude/quadmagnitude), digits = 2))
    out <- c("The ratio of non-quadratic terms/quadratic terms at the point, where the residual has been computed, is:", out)
  }
  
  
  return(out)
}

#' Analyse the Jacobian for singular values and their corresponding directions in parameter space
#' @param SVDJac The SVD of the jacobian (attr(residual, "deriv")) 
#' @param tolerance Tolerance level for the singular value
#' @export
analyseJ <- function(residual, tolerance=10^-8) {
  
  myderiv <- data.matrix(attr(residual, "deriv")[,-c(1,2)])
  S <- svd(myderiv)

  
  singularindex <- which(abs(S$d) < max(dim(S$u))*max(S$d)*tolerance) #take only dimensions with singular values smaller than a given tolerance tol, this is copied from statistics.R, pseudoinverse()
  singulardirections <- submatrix(S$v, cols = singularindex)
  singularvalue <- d[singularindex]
  
  return(list(singularindex=singularindex, singularvalue = singularvalue, singulardirections = singulardirections))
}




#' Compute the christoffel symbols
#'
#' @description Compute the christoffel symbols of the solutions of the ode at the point given by the parameter estimates given by the best fit.
#' @param residual The combined residuals of all conditions, regularized if necessary.
#' @param tolerance Tolerance level for the singular value
#' @return Christoffel symbols. A 3-dimensional array with dimension names. The upper index of the christoffel smybols marked with a dash.
#' @export
christoffels <- function(residual) {
  
  # 1st deriv
  myderiv <- data.matrix(attr(residual, "deriv")[,-c(1,2)])
  S <- svd(myderiv)
  myinverse<-S$v%*%diag(1/(S$d))%*%t(S$u)
  
  # 2nd deriv
  mysderiv <- data.matrix(attr(residual, "sderiv")[,-c(1,2)])
  
  # christoffel symbols
  mychristoffel <- myinverse %*% mysderiv
  mychristoffelarray <- array(mychristoffel, dim = rep(dim(myderiv)[2], 3))
  #symmetrize
  mychristoffelarray <- array(apply(mychristoffelarray, 1, function(i) (i+t(i))/2), dim=dim(mychristoffelarray))
  mychristoffelarray <- aperm(mychristoffelarray, c(3,2,1))
  
  
  #naming
  mypars <- colnames(myderiv)
  mydashedpars <- outer(mypars, "'", paste, sep="")
  dimnames(mychristoffelarray) <- list(mydashedpars, mypars, mypars)
  
  return(mychristoffelarray)
}


#' Regularize Jacobian 
#' @description Regularize Jacobian matrix with weak prior if necessary. 
#' Why not regularize the Hessian matrix directly? Because the svd simplifies the formula of the Christoffel symbols so that one just 
#' needs the "inverse" of the Jacobian.
#' @param residual The residual at a given point
#' @param pouter The parameters at which the residual is computed
#' @param prior A vector of length(pouter) with mean values around which the prior is centered
#' @param priorthreshold The threshold of min(singularvalues^2) for which a prior is to be applied. If any singular value is smaller than it, apply the prior.
#' @param sigma The sigma of the prior
#'@export
regularizejacobian <- function(residual, pouter, prior, priorthreshold= 10^(-9), sigma = 10){

  # format residuals, make svd
  myresidual <- combineresiduals(residual)
  S <- svd(attr(myresidual, "deriv")[,-c(1,2)])

  # if at least one direction has a small singular value, regularize
  if(any((S$d)^2 < priorthreshold)) {#priorthreshold is a border which has to be found by experience
  
    out.prior <- constraintL2(pouter, prior, sigma = sigma) #use L2-constraint, in general try to make sigma as large as possible minimize prior bias
    
    priorderiv <- data.frame(time = rep(0, length(out.prior$gradient)), name = rep("prior", length(out.prior$gradient)), diag(out.prior$gradient))
    colnames(priorderiv) <- c("time", "name", names(out.prior$gradient))
    deriv.data <- rbind(attr(myresidual, "deriv"), priorderiv)
    
    priorsderiv <- data.frame(matrix(0, nrow = length(out.prior$gradient), ncol = ncol(attr(myresidual, "sderiv"))))
    priorsderiv[,2] <- "prior"
    colnames(priorsderiv)<- colnames(attr(myresidual, "sderiv"))
    sderiv.data <- rbind(attr(myresidual, "sderiv"), priorsderiv)
    
    myresidual<-objframe(myresidual, deriv = deriv.data, sderiv = sderiv.data)
    
  }
  
  return(myresidual)
}


#' Riemanntensor
#' 
#' @description Compute the Riemann tensor out of the second derivatives, to gain an even better approximation of the geodesics.
#' The underlying formula is found in the papaer of Transtrum, section 7A.
#' @param residual The residual at a given point
#' @param raise Logical. Raise first index (TRUE) or not (FALSE)?
#' 
#'@export
riemanntensor <- function(residual, raise=FALSE) {
  
  myderiv <- data.matrix(attr(residual, "deriv")[,-c(1,2)])
  
  mysderiv <- data.matrix(attr(residual, "sderiv")[,-c(1,2)])
  mysderiv <- array(mysderiv, dim = c(nrow(myderiv), ncol(myderiv), ncol(myderiv)))
  
  S <- svd(myderiv)
  d <- S$d
  u <- S$u
  v <- S$v
  
  PN<-diag(nrow = nrow(u)) - u%*%t(u)
  
  namesarray <- outer(dimnames(mysderiv)[[1]], dimnames(mysderiv)[[2]], paste, sep=".")
  namesarray <- array(namesarray, c(nrow(myderiv), ncol(myderiv), ncol(myderiv)))

  PNsderiv <- apply(mysderiv, c(2,3), function(i) PN%*%i)
  riemann1 <- apply(mysderiv, c(2,3), function(mysd) {
    apply(PNsderiv, c(2,3), function(PNsd) {sum(mysd*PNsd)})
  })
  riemann1 <- array(riemann1, dim= rep(ncol(myderiv), 4))
  riemann2 <- aperm(riemann1, c(1,4,3,2))
  
  riemann <- riemann1 - riemann2
  
  # raise fist index, ignore singular directions
  if(raise) {
    ginverse <- v%*%diag(1/S$d^2)%*%t(v)
    riemann <- apply(riemann, c(2,3,4), function(x) ginverse%*%x)
  }
  return(riemann)
}

#' Christoffel derivative
#' Compute derivatives of the Christoffel symbols wrt the parameters via the Riemann tensor
#' It's a 1,3 tensor, given out like [mu alpha beta nu], it means the christoffel symbol ^mu_alpa_beta derived wrt _nu
#' @export
christoffelderiv <- function(residual) {
  riemann <- riemanntensor(residual)
  chrisderiv <- -1/3*(riemann + aperm(riemann, c(1,3,2,4)))
  return(chrisderiv)
}


#' compute the relative absolute value of christoffelderiv/chris
#'  This function sums over delta: 
#'  relativechrisderiv = sum_delta (sqrt(christoffelderiv^mu_alpha_beta_delta ^2 )/ abs(chris ^mu _alpha_beta))
#' @param residual 
#'
#' @return
#' @export
#'
#' @examples
relativechrisderiv <- function(residual) {
  mychris <- christoffels(residual)
  mychrisderiv <- christoffelderiv(residual)
  
  abschrisderiv <- apply(mychrisderiv, c(1,2,3), function(i) {
    i %>% raise_to_power(2) %>% sum() %>% sqrt()
  })
  
  return(abschrisderiv/abs(mychris))
}





#' Compute the second order derivative via finite differences
#' 
#' @param residual The residual evaluated at pouter
#' @param pouter The parameter vector at which the residual is evaluated
#' @param stepwidth Just for debugging and to experiment which stepsize is good.
#' 
#' @export
reswithfinitesderiv <- function(pouter, stepwidth = 0.001) {
  
  prediction <- function(pouter){
    x(timesD, pouter, fixed = fixed, deriv = 1) #pred(pinner), deriv= pred(pinner).pouter, sderiv= pred(pinner).pouter.pouter,
  }
  residual <- lapply(names(data), function(cn) res(data[[cn]], prediction(pouter)[[cn]])) #create list of residuals of all conditions
  
  residual0 <- combineresiduals(residual)
  deriv0 <- data.matrix(attr(residual0, "deriv")[,-c(1,2)])
  
  directions <- diag(1, nrow= length(pouter), ncol= length(pouter))
  directions <- lapply(1:length(pouter), function(i) directions[,i])  
  
  sderiv <- lapply(directions, function(direction) {
    newpouter<- pouter+stepwidth*direction
    newresidual <- combineresiduals(lapply(names(data), function(cn) res(data[[cn]], prediction(newpouter)[[cn]])))
    newderiv <- data.matrix(attr(newresidual, "deriv")[,-c(1,2)])
    mysderiv <- (newderiv - deriv0)/stepwidth
  })
  sderiv <- array(unlist(sderiv), dim = c(nrow(deriv0), ncol(deriv0), ncol(deriv0)))
  dimnames(sderiv) <- list(as.character(1:nrow(deriv0)), colnames(deriv0), colnames(deriv0))
  sderiv<-data.frame(residual0[,c(1,2)], sderiv)
  attr(residual0, "sderiv") <- sderiv
  return(residual0)
}


#' Project the path travelled onto the radial vector "e_r" in data-space
#' actually an unnecessary function, since one can just directly evaluate the chi^2 and subtract...
#' 
#' @param out the output of the numerical integration of the geodesic
#' 
#' 
#' @export
erprojection<-function(out) {
  
  myrvecs<- do.call(cbind,lapply(1:nrow(out), function(i) {
    myresidual<-resi(out[i,2:5])[7]
  }))
  
  Deltachisquare<-do.call(rbind, lapply(1:(ncol(myrvecs)-1), function(i) {
    delchi<-myrvecs[,i+1]%*%myrvecs[,i+1]+myrvecs[,i+1]%*%myrvecs[,i+1]-2*myrvecs[,i]%*%myrvecs[,i+1]
  }))
  
  chisquare<-diffinv(Deltachisquare, xi=myrvecs[,1]%*%myrvecs[,1])
}



#' Print out the christoffels arranged so that the lower indices define the matrix.
#' could be made nicer with a class called christoffels and then to write a print method for this class. But on the other hand this is " mit Kanonen auf Spatzen"...
#' @param arr The array with the christoffel symbols
#'
#' @return
#' @examples
#'
#' @export
printarray <- function(arr) {
  n<-dim(arr)[1]
  wup <- lapply(1:n, function(i) arr[i,,])
  return(wup)
}



#' Commutator of two christoffel matrices
#' Compute the commutator of two christoffel matrices with one lower index kept constant:
#' commutator = [expm(chris[,a,]*Tha),expm(chris[,b,]*Thb)]
#'
#' @param chris 
#' @param Th1 
#' @param Th2 
#'
#' @return
#' @export
#'
#' @examples
chriscommutator <- function(chris, Th1, Th2) {
  mychris <- lapply(1:(dim(chris)[1]), function(i) t(chris[,i,]))
  commutator <- lapply(1:(dim(chris)[1]), function(i) {
    lapply(i:(dim(chris)[1]), function(j) {
      if(i==j) return(NULL)
      else 
        ei <- expm(mychris[[i]]*Th1)
        ej <- expm(mychris[[j]]*Th2)
        comm <- ei%*%ej-ej%*%ei
      return(comm)
    })
  })
  return(commutator)
}



#' Create a 2d grid with colnames specified by outer_pars
#'
#' @param bestfit 
#' @param grid_pars 
#' @param seqx 
#' @param seqy 
#'
#' @return
#' @export
#'
#' @examples
create_outer_pars_grid <- function(bestfit, grid_pars, seqx, seqy) {
  Thgrid <- rep(bestfit, length(seqx)*length(seqy)) %>% matrix(nrow = length(bestfit)) %>% t() %>% as.data.frame() %>% structure(names = names(bestfit))
  mygrid <- expand.grid(seqx, seqy)
  Thgrid[,grid_pars] = mygrid
  return(Thgrid)
}


#' Compute the metric of a given point
#'
#' @param bestfit 
#'
#' @return gmunu, the metric
#' @export
#'
#' @examples
gmunu <- function(bestfit) {
  residual <- combineresiduals(resi(bestfit))
  myresidual <- as.vector(residual$residual)
  myderiv <- data.matrix(attr(residual, "deriv")[,-c(1,2)])
  mysderiv <- data.matrix(attr(residual, "sderiv")[,-c(1,2)])
  
  return(t(myderiv)%*%myderiv)
}


#' Compute the true hessian matrix of a given point
#'
#' @param bestfit 
#'
#' @return hmunu the hessian
#' @export
#'
#' @examples
hmunu <- function(bestfit) {
  residual <- combineresiduals(resi(bestfit))
  myresidual <- as.vector(residual$residual)
  myderiv <- data.matrix(attr(residual, "deriv")[,-c(1,2)])
  mysderiv <- data.matrix(attr(residual, "sderiv")[,-c(1,2)])
  
  gmunu   <- t(myderiv)%*%myderiv
  rsderiv <- matrix(myresidual%*%mysderiv, nrow=length(bestfit))
  return(gmunu + rsderiv)
}

#' For a given identifier, set the object bestfit and profile_objfun to the respective objects of the identifier
#'
#' @param identifier The identifying phrase, ie "_1"
#'
#' @return
#' @export
#'
#' @examples
set_bestfit_and_profile_objfun <- function(identifier) {
  bf_string <- paste0("bestfit", identifier)
  obj_string <- paste0("obj", identifier)
  bestfit <<- eval(as.symbol(bf_string))
  profile_objfun <<- obj_string # what about RNCssquaredobj und hobj?
}

#' Define a global object with an identifiable name out of an object with given objname
#'
#' @param objtypes a vector of strings naming an object like "obj", "bestfit" or "profiles"
#' @param identifier The identifying phrase, ie "_1"
#'
#' @details if obj exists, but is a special obj that is supposed to be stored with a descrptive name, like obj_1 (for obj(dataset 1)), 
#'  one can use define_global_object("obj","_1")
#' @return
#' @export
#'
#' @examples
define_global_objects <- function(objtypes, identifiers) {
  objname <- outer(objtypes, identifiers, paste0) %>% as.character()
  objtypes <- rep(objtypes, length(identifiers))
  lapply(1:length(objname), function(i) assign(objname[i], eval(as.name(objtypes[i])), envir = .GlobalEnv))
}

#' Return object type and symbol as name
#'
#' @param objtype 
#' @param identifier 
#'
#' @return
#' @export
#'
#' @examples
eval_object <- function(objtype, identifier) {
  objname <- paste0(objtype, identifier)
  eval(as.name(objname))
}

#' Title
#'
#' @param ... 
#'
#' @return
#' @export
#'
#' @examples
pasteDash <- function(...) {
  paste(..., sep = "_")
}




#' Some basic definitions
#'
#' @param withtrue 
#'
#' @return
#' @export
#'
#' @examples
dlIinit <- function(withtrue = F) {
  data_set <- ""
  all_objtypes <<- c("obj", "bestfit", "profiles")
  all_identifiers <- paste0("_" ,c("data", "RNC", "h", "g", "_p", "_p_rootsolve"))
  all_objfuns <<- paste0("obj", all_identifiers)
  namestring <<- "_data"
  with_prior <<- F
  assign("sigma" , 0.02, envir = .GlobalEnv)
  tol <<- 1e5
}




#' Weighted residual, defaulting to "residual", if weighted does not exist...
#'
#' @param pouter 
#'
#' @return
#' @export
#'
#' @examples
weighted_resi <- function(pouter) {
  myresi <- resi(pouter)
  myweighted <- myresi$weighted.residual
  if(is.null(myresi$weighted.residual)) myweighted <- myresi$residual
  return(myweighted)
}




#' ggtheme for my Master thesis.
#'
#' @param base_size 
#' @param base_family 
#'
#' @description This is a modification of Daniel Kaschek's theme_dMod, which can be found in the dMod-Package in CRAN
#'
#' @return
#' @export
#'
#' @examples
theme_eladin <- function (base_size = 12, base_family = "") {
  colors <- list(
    medium = c(gray = '#737373', red = '#F15A60', green = '#7AC36A', blue = '#5A9BD4', orange = '#FAA75B', purple = '#9E67AB', maroon = '#CE7058', magenta = '#D77FB4'),
    dark = c(black = '#010202', red = '#EE2E2F', green = '#008C48', blue = '#185AA9', orange = '#F47D23', purple = '#662C91', maroon = '#A21D21', magenta = '#B43894'),
    light = c(gray = '#CCCCCC', red = '#F2AFAD', green = '#D9E4AA', blue = '#B8D2EC', orange = '#F3D1B0', purple = '#D5B2D4', maroon = '#DDB9A9', magenta = '#EBC0DA')
  )
  gray <- colors$medium["gray"]
  black <- colors$dark["black"]
  
  theme_bw(base_size = base_size, base_family = base_family) + 
    theme(line = element_line(colour = gray), 
          rect = element_rect(fill = "white", colour = NA), 
          text = element_text(colour = black), 
          axis.ticks = element_line(colour = gray), 
          legend.key = element_rect(colour = NA), 
          panel.border = element_rect(colour = gray), 
          panel.grid = element_line(colour = gray, size = 0.2), 
          strip.background = element_rect(fill = "white", colour = NA))
}


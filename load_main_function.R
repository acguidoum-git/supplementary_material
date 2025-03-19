###########################################
### Simulation model Eq.(1): R function ###
###########################################

####################
## A: alpha
## B: beta
## G: gamma
## g_t: g(t)
####################
## ncpu: number of CPU cores
####################

Sim_mod <- function(x0,gt,k=NULL,Alpha=1,Beta=1,Gamma=1,H=0.75,
                    T=1,N=1000,M=100,ncpu=6)
           {
   if (H <= 0) stop("'H' must be in '[0 . 1]'")
   if (H >= 1) stop("'H' must be in '[0 . 1]'")
   if (Gamma <= 0) stop("'Gamma' must be numeric > 0")
   if (!is.expression(gt)) stop(" 'g(t)' must be expressions in 't'")
   samp <- yuima::setSampling(Terminal=T, n=N)
   cl <- parallel::makeCluster(ncpu)
   parallel::clusterEvalQ(cl, library(yuima))
   parallel::clusterExport(cl, varlist = c("H"))
   W_h <- parallel::parSapply(cl,1:M,function(i) 
             yuima::simulate(yuima::setModel(drift = 0, diffusion = 1 ,solve.var = "x"),
             hurst = H, sampling = yuima::setSampling(Terminal=T, n=N))@data@original.data)
   parallel::stopCluster(cl)
   temps <- samp@grid[[1]]  
   if (is.null(k)){ Gt <- gt}else{
   if (ttutils::isInteger(k)==FALSE) stop("k must be non-negative integer.")
   Gt  <- Deriv::Deriv(gt,"t",nderiv = k)}
   G_t  <- function(t) eval(Gt)
   X <- (G_t(temps))^(Alpha)*(x0*G_t(0)^(-Alpha)+Beta*temps+Gamma*W_h)
   name <- "X"
   name <- if(M > 1) paste("X",1:M,sep="")
   X <- ts(X, start = 0, deltat = samp@delta, names=name)
   return(invisible(X))
}


###########################################
###       Empirical moments             ###
###########################################

moment <- function(x, order = 1,center = TRUE,...)
{
  if (!is.numeric(x)) 
    stop("argument 'x' must be numeric")
  x = x[!is.na(x)]
  if (center) x <- x - mean(x)
  return(sum(x^order)/length(x))
}

###################################################
###      Exact moments Eq.(5): R function       ###
###################################################

HO_mom <- function(x0, gt, k=NULL, Alpha=1, Beta=1, Gamma=1, H=0.75, p=1, t) {
  if (!inherits(gt, "expression")) stop("'gt' must be an expression in 't'")
  if (H <= 0 | H > 1) stop("'H' must be in '(0,1)'")
  if (p <= 0) stop("'p' must be numeric > 0'")
  if (Gamma <= 0) stop("'Gamma' must be numeric > 0'")
  t <- as.numeric(t)
  is_single_t <- length(t) == 1  
  if (all(t == 0)) {
    return(x0^p)  
  }
  if (is.null(k)) {
    Gt <- gt
  } else {
    if (!ttutils::isInteger(k)) stop("k must be a non-negative integer.")
    Gt <- Deriv(gt, "t", nderiv = k)
  }
  G_t <- function(t) sapply(t, function(ti) eval(Gt, envir = list(t = ti)))
  G_t_0 <- G_t(0)
  if (length(G_t_0) == 0) stop("Evaluation of g(t) at t=0 failed.")
  A   <- function(t) ifelse(t == 0, NA, 2^(-0.5 * p) * ((Gamma * G_t(t)^Alpha) / 1i)^p * t^(H * p))
  B   <- function(t) ifelse(t == 0, NA, 0.5 * sqrt(2) * 1i * ((G_t_0^(-Alpha) * x0 + Beta * t) / Gamma) * t^(-H))
  compute_moment <- function(t) {
    if (t == 0) {
      return(x0^p)  
    }
    if (ttutils::isInteger(p)) {
      Herm <- hermite(x = B(t), n = p, prob = FALSE)
    } else {
      M1 <- gsl::hyperg_1F1(a = -0.5 * p, b = 0.5, x = B(t)^2)
      M2 <- gsl::hyperg_1F1(a = 0.5 - 0.5 * p, b = 1.5, x = B(t)^2)
      Herm <- 2^p * sqrt(pi) * ((M1 / gamma(0.5 - 0.5 * p)) - 2 * B(t) * (M2 / gamma(-0.5 * p)))
    }
    return(Re(A(t) * Herm))
  }
  result <- sapply(t, compute_moment)
  return(if (is_single_t) result else invisible(result))
}

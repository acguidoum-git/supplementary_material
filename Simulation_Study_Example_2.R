########################################################################
########################################################################
########################################################################
### Example 2: case g(t) is nonlinear time-dependent function       ####
########################################################################
########################################################################
########################################################################
library(future)
future::plan(cluster, workers = 12)

set.seed(12345789)
x0=20;T=10;
g_t <- expression(cos(t+pi)) ## g(t) == cos(t+pi)
A=2;B=2;G=2;k=4 
H=0.75 # H=0.55

## Simulation 
mod2 <- Sim_mod(x0=x0,gt=g_t,Alpha=A,Beta=B,Gamma=G,H=H,
                    T=T,N=1000,M=1000,ncpu=12)

## Monte-Carlo moments approximation

p=c(1/6,22/5,67/7,78/5,93/4)
MC_mom_p <- do.call("cbind",furrr::future_map(1:length(p),function(i)apply(mod2,1,moment,order=p[i],center=FALSE)))
MC_mom_scaled_p <- sapply(1:length(p),function(i) MC_mom_p[,i]^(1/p[i]))


## Exact moments

Ex_mom_p <- sapply(1:length(p),function(i)
  HO_mom(x0=x0,gt=g_t,Alpha=A,Beta=B,Gamma=G,H=H,p=p[i],t=as.vector(time(mod2))))
Ex_mom_scaled_p <- sapply(1:length(p),function(i) Ex_mom_p[,i]^(1/p[i]))

## Plots

Color <- ggsci::pal_jco("default",alpha=0.7)(5)
Expr <- expression(M[frac(1,6)](t),M[frac(22,5)](t),M[frac(67,7)](t),
                   M[frac(78,5)](t),M[frac(93,4)](t))
par(mar = c(4, 3, 1, 1), xaxs = "i")
matplot(as.vector(time(mod2)),MC_mom_scaled_p,type="l",lty=1,las=1,col=Color,lwd=2,
        ylab="Moment",xlab="Time")
grid(col = "gray", lty = "dotted")
matplot(as.vector(time(mod2)),Ex_mom_scaled_p,type="l",lty=2,add=TRUE,col=1,lwd=2)
legend("topleft",Expr,inset = .01,fill=Color,lty=NA,
        lwd=2,cex=1.4)
legend("topleft",expression(m[p](t),E(X[t]^{p})),inset = c(0.35,0),
       col=1,lty=c(1,2),lwd=2,cex=1.4,horiz = FALSE, merge = TRUE,bty="n")
box()

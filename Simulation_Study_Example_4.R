##################################################################
##################################################################
##################################################################
### Example 4: fractional Ornstein-Uhlenbeck process          ####
##################################################################
##################################################################
##################################################################
library(future)
future::plan(cluster, workers = 12)


set.seed(12345)
x0=10;T=5;
g_t <- expression(t+1) ## g(t) == t+1
A=-2;B=0;G=1;k=0
H=0.6 

## Simulation 
mod4 <- Sim_mod(x0=x0,gt=g_t,k=k,Alpha=A,Beta=B,Gamma=G,H=H,
                    T=T,N=1000,M=10000,ncpu=12)

## Monte-Carlo moments approximation

p=seq(1/2,9/2,by=1)
MC_mom_p <- do.call("cbind",furrr::future_map(1:length(p),function(i)apply(mod4,1,moment,order=p[i],center=FALSE)))
MC_mom_scaled_p <- sapply(1:length(p),function(i) MC_mom_p[,i]^(1/p[i]))

## Exact moments

Ex_mom_p <- sapply(1:length(p),function(i)
  HO_mom(x0=x0,gt=g_t,k=k,Alpha=A,Beta=B,Gamma=G,H=H,p=p[i],t=as.vector(time(mod4))))
Ex_mom_scaled_p <- sapply(1:length(p),function(i) Ex_mom_p[,i]^(1/p[i]))

## Plots

Color <- ggsci::pal_jco("default",alpha=0.7)(5)
Expr <- expression(M[frac(1,2)](t),M[frac(3,2)](t),M[frac(5,2)](t),
                    M[frac(7,2)](t),M[frac(9,2)](t))
par(mar = c(4, 3, 1, 1)+0.5)
matplot(as.vector(time(mod4)),MC_mom_p,type="l",lty=1,las=1,col=Color,lwd=2,
        ylab="",xlab="Time",log="y")
grid(col = "gray", lty = "dotted")
matplot(as.vector(time(mod4)),Ex_mom_p,type="l",lty=2,add=TRUE,col=1,lwd=2,log="y")
legend("topright",Expr,inset = .01,fill=Color,lty=NA,
        lwd=2,cex=1.5)
legend("topleft",expression(m[p](t),E(X[t]^{p})),inset = c(0.3,0),
       col=1,lty=c(1,2),lwd=2,cex=1.45,horiz = FALSE, merge = TRUE,bty="n")
box()



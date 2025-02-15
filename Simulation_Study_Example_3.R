##################################################################
##################################################################
##################################################################
### Example 3: fractional Brownian bridge process             ####
##################################################################
##################################################################
##################################################################
library(future)
future::plan(cluster, workers = 12)

set.seed(123456)

x0=10;T=20;
g_t <- expression(log(t-T)) ## g(t) == ln(t-T)
A=-1;B=0;G=1;k=1
H=0.75 # H=0.55

## Simulation 
mod3 <- Sim_mod(x0=x0,gt=g_t,k=k,Alpha=A,Beta=B,Gamma=G,H=H,
                    T=T,N=1000,M=100000,ncpu=12)


## Monte-Carlo moments approximation

p=seq(1,4,by=1)
MC_mom_p <- do.call("cbind",furrr::future_map(1:length(p),function(i)apply(mod3,1,moment,order=p[i],center=FALSE)))
MC_mom_scaled_p <- sapply(1:length(p),function(i) MC_mom_p[,i]^(1/p[i]))

## Exact moments

Ex_mom_p <- sapply(1:length(p),function(i)
  HO_mom(x0=x0,gt=g_t,k=k,Alpha=A,Beta=B,Gamma=G,H=H,p=p[i],t=as.vector(time(mod3))))
Ex_mom_scaled_p <- sapply(1:length(p),function(i) Ex_mom_p[,i]^(1/p[i]))

## Plots

Color <- ggsci::pal_jco("default",alpha=0.7)(4)
Expr <- expression(M[1](t),M[2](t),M[3](t),M[4](t))
par(mar = c(4, 3, 1, 1), xaxs = "i")
matplot(as.vector(time(mod3)),MC_mom_scaled_p,type="l",lty=1,las=1,col=Color,lwd=2,
        ylab="Moment",xlab="Time")
grid(col = "gray", lty = "dotted")
matplot(as.vector(time(mod3)),Ex_mom_scaled_p,type="l",lty=2,add=TRUE,col=1,lwd=2)
legend("topright",Expr,inset = .01,fill=Color,lty=NA,
        lwd=2,cex=1.35)
legend("topleft",expression(m[p](t),E(X[t]^{p})),inset = c(0.75,0.25),
       col=1,lty=c(1,2),lwd=2,cex=1.35,horiz = FALSE, merge = TRUE,bty="n")
box()

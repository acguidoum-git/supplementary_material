###############################################################
###############################################################
###############################################################
### Example 1: case g(t) is linear time-dependent function ####
###############################################################
###############################################################
###############################################################

set.seed(12345)
x0=30;T=5
g_t <- expression(2*t+1) ## g(t) == 2t+1
A=-0.5;B=14;G=6;k=0 
H=0.55 # H=0.75

## Simulation 
mod1 <- Sim_mod(x0=x0,gt=g_t,k=k,Alpha=A,Beta=B,Gamma=G,H=H,
                    T=T,N=1000,M=100000,ncpu=12)

## Monte-Carlo moments approximation

p=seq(1/2,5,by=1/2)
MC_mom_p <- sapply(1:length(p),function(i)apply(mod1,1,moment,order=p[i],center=FALSE))
MC_mom_scaled_p <- sapply(1:length(p),function(i) MC_mom_p[,i]^(1/p[i]))

## Exact moments

Ex_mom_p <- sapply(1:length(p),function(i)
  HO_mom(x0=x0,gt=g_t,k=k,Alpha=A,Beta=B,Gamma=G,H=H,p=p[i],t=as.vector(time(mod1))))
Ex_mom_scaled_p <- sapply(1:length(p),function(i) Ex_mom_p[,i]^(1/p[i]))

## Plots

Color <- ggsci::pal_jco("default", alpha = 0.7)(10)
Expr <- expression(M[frac(1,2)](t), M[1](t), M[frac(3,2)](t), M[2](t),
                   M[frac(5,2)](t), M[3](t), M[frac(7,2)](t), M[4](t),
                   M[frac(9,2)](t), M[5](t))
par(mar = c(4, 3, 1, 1), xaxs = "i")
matplot(as.vector(time(mod1)), MC_mom_scaled_p, type = "l", lty = 1, las = 1,
        col = Color, lwd = 2, ylab = "Moment", xlab = "Time")
grid(col = "gray", lty = "dotted")
matplot(as.vector(time(mod1)), Ex_mom_scaled_p, type = "l", lty = 2,
        add = TRUE, col = 1, lwd = 2)
legend("topleft", legend = Expr, inset = c(0.09, 0.01),
       fill = Color, lty = NA, lwd = 2, cex = 1.45)
legend("topleft", legend = expression(m[p](t), E(X[t]^{p})),
       inset = c(0.45, 0), col = 1, lty = c(1, 2), lwd = 2,
       cex = 1.45, horiz = FALSE, merge = TRUE, bty = "n")
box()

########################################################################
########################################################################
########################################################################
###            case g(t)  nonlinear time-dependent function       ####
########################################################################
########################################################################
########################################################################
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

## plots with ggplot2
time_vec <- as.vector(time(mod2))
p_vals <- c(1/6,22/5,67/7,78/5,93/4)

## data for h = 0.55
# Rename MC_mom_scaled_p to MC_mom_scaled_p_55 for H=0.55
df_MC_55 <- as.data.frame(MC_mom_scaled_p_55)
colnames(df_MC_55) <- paste0("p_", p_vals)
df_MC_55$time <- time_vec
df_MC_55$h <- "0.55"
df_MC_long_55 <- pivot_longer(df_MC_55, 
                              cols = starts_with("p_"), 
                              names_to = "p", 
                              values_to = "value") %>%
  mutate(p = as.numeric(sub("p_", "", p)),
         type = "MC")
# Rename Ex_mom_scaled_p to Ex_mom_scaled_p_55 for H=0.55
df_Ex_55 <- as.data.frame(Ex_mom_scaled_p_55)  
colnames(df_Ex_55) <- paste0("p_", p_vals)
df_Ex_55$time <- time_vec
df_Ex_55$h <- "0.55"
df_Ex_long_55 <- pivot_longer(df_Ex_55, 
                              cols = starts_with("p_"), 
                              names_to = "p", 
                              values_to = "value") %>%
  mutate(p = as.numeric(sub("p_", "", p)),
         type = "Exact")

## data for h = 0.75
# Rename MC_mom_scaled_p to MC_mom_scaled_p_75 for H=0.75
df_MC_75 <- as.data.frame(MC_mom_scaled_p_75)
colnames(df_MC_75) <- paste0("p_", p_vals)
df_MC_75$time <- time_vec
df_MC_75$h <- "0.75"
df_MC_long_75 <- pivot_longer(df_MC_75, 
                              cols = starts_with("p_"), 
                              names_to = "p", 
                              values_to = "value") %>%
  mutate(p = as.numeric(sub("p_", "", p)),
         type = "MC")
# Rename Ex_mom_scaled_p to Ex_mom_scaled_p_75 for H=0.75
df_Ex_75 <- as.data.frame(Ex_mom_scaled_p_75)
colnames(df_Ex_75) <- paste0("p_", p_vals)
df_Ex_75$time <- time_vec
df_Ex_75$h <- "0.75"
df_Ex_long_75 <- pivot_longer(df_Ex_75, 
                              cols = starts_with("p_"), 
                              names_to = "p", 
                              values_to = "value") %>%
  mutate(p = as.numeric(sub("p_", "", p)),
         type = "Exact")

## Combinez data
df_all <- bind_rows(df_MC_long_55, df_Ex_long_55, df_MC_long_75, df_Ex_long_75)
## labels expression
num <- c(1, 22, 67, 78, 93)
den <- c(6, 5, 7, 5, 4)
labels_expr <- sapply(seq_along(num), function(i) {
  deparse(bquote(M[frac(.(num[i]),.(den[i]))](t)))
})
labels_expr <- parse(text = labels_expr)
Color <- ggsci::pal_jco("default", alpha = 0.65)(length(p_vals))
names(Color) <- sort(unique(df_MC_long_75$p))

## plot with ggplot

plot_2 <- ggplot() +
  geom_line(data = df_all %>% filter(type == "MC"),
            aes(x = time, y = value, group = factor(p), 
                color = factor(p), linetype = "MC"),
            linewidth = 1) +
  geom_line(data = df_all %>% filter(type == "Exact"),
            aes(x = time, y = value, group = factor(p), linetype = "Exact"),
            color = "black", linewidth = 0.6) +
  scale_linetype_manual(name="",values = c("MC" = "solid", "Exact" = "dashed"),
                        labels = c("MC" = expression(m[p](t)), "Exact" = expression(E(X[t]^p)))) +
  scale_color_manual(name = "",
                     values = Color,
                     labels = labels_expr) +
  scale_x_continuous(breaks = pretty_breaks(n = 12),
                     labels = number_format(accuracy = 0.1)) +
  scale_y_continuous(breaks = pretty_breaks(n = 10),
                     labels = number_format(accuracy = 0.1)) +
  labs(x = "Time", y = "Moment") +
  facet_wrap(~ h, labeller = label_bquote(~ H == .(h)), nrow = 1) +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        panel.grid.major = element_line(color = "grey80", linewidth = 0.5),
        panel.grid.minor = element_line(color = "grey90", linewidth = 0.25),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.75),
        panel.spacing = unit(0.01, "lines"),
        legend.position = "top",
        legend.title = element_text(size = 18, face = "bold", hjust = 1.5, margin = margin(b = 0)), 
        legend.text = element_text(size = 18, face = "bold"),
        legend.box.margin = margin(t = -11), 
        legend.box.spacing = unit(0.01, "cm"), 
        legend.spacing.y = unit(0.01, "cm"),
        axis.title.x = element_text(face = "bold"),
        axis.title.y = element_text(face = "bold"),
        plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm"),
        strip.text = element_text(size = 15, face = "bold.italic",margin = margin(1, 1, 1, 1)),
        strip.background = element_rect(fill = "#f0f0f0", color = "black", linewidth = 0.75))
plot_2
                          

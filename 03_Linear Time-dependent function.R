###############################################################
###############################################################
###############################################################
###         Case g(t) linear time-dependent function      ####
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

## Plots with ggplot

time_vec <- as.vector(time(mod1))
p_vals <- seq(1/2, 5, by = 1/2)

## data for h = 0.55
# rename MC_mom_scaled_p by MC_mom_scaled_p_55 for H=0.55
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
# rename Ex_mom_scaled_p by Ex_mom_scaled_p_55 for H=0.55
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
# rename MC_mom_scaled_p by MC_mom_scaled_p_75 for H=0.75
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
# rename Ex_mom_scaled_p by Ex_mom_scaled_p_75 for H=0.75
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

# Combinez data
df_all <- bind_rows(df_MC_long_55, df_Ex_long_55, df_MC_long_75, df_Ex_long_75)
# labels expression
labels_expr <- sapply(sort(unique(df_MC_long$p)), function(x) {
  if(x %% 1 == 0){
    deparse(bquote(M[.(x)](t)))
  } else {
        deparse(bquote(M[frac(.(as.integer(2*x)),2)](t)))
  }
})
labels_expr <- parse(text = labels_expr)
Color <- ggsci::pal_jco("default", alpha = 0.65)(length(p_vals))
names(Color) <- sort(unique(df_MC_long_75$p))

## plot

plot_1 <- ggplot() +
  geom_line(data = df_all %>% filter(type == "MC"),
            aes(x = time, y = value, group = factor(p), 
                color = factor(p), linetype = "MC"),
            linewidth = 1) +
  geom_line(data = df_all %>% filter(type == "Exact"),
            aes(x = time, y = value, group = factor(p), linetype = "Exact"),
            color = "black", linewidth = 0.6) +
  scale_color_manual(name = "",
                     values = Color,
                     labels = labels_expr) +
  scale_linetype_manual(name="",values = c("MC" = "solid", "Exact" = "dashed"),
                        labels = c("MC" = expression(m[p](t)), "Exact" = expression(E(X[t]^p)))) +
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
plot_1

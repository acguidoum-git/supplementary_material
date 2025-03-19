##################################################################
##################################################################
##################################################################
###            Fractional Ornstein-Uhlenbeck process          ####
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
matplot(as.vector(time(mod4)),MC_mom_scaled_p,type="l",lty=1,las=1,col=Color,lwd=2,
        ylab="",xlab="Time")
grid(col = "gray", lty = "dotted")
matplot(as.vector(time(mod4)),Ex_mom_scaled_p,type="l",lty=2,add=TRUE,col=1,lwd=2,log="y")
legend("topright",Expr,inset = .01,fill=Color,lty=NA,
        lwd=2,cex=1.5)
legend("topleft",expression(m[p](t),E(X[t]^{p})),inset = c(0.3,0),
       col=1,lty=c(1,2),lwd=2,cex=1.45,horiz = FALSE, merge = TRUE,bty="n")
box()

## Plots whit ggplot2

time_vec <- as.vector(time(mod4))
p_vals=seq(1/2,9/2,by=1)

## data for h = 0.55

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

## Combinez data

df_all <- bind_rows(df_MC_long_55, df_Ex_long_55)

## labels expression

labels_expr <- sapply(sort(unique(df_MC_long_55$p)), function(x) {
  if(x %% 1 == 0){
    deparse(bquote(M[.(x)](t)))
  } else {
    deparse(bquote(M[frac(.(as.integer(2*x)),2)](t)))
  }
})
labels_expr <- parse(text = labels_expr)
Color <- ggsci::pal_jco("default", alpha = 0.65)(length(p_vals))
names(Color) <- sort(unique(df_MC_long_55$p))

plot_4 <- ggplot() +
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
  scale_y_log10(breaks = scales::log_breaks(n = 10),
                 labels = scales::label_number(accuracy = 0.1)) +
  labs(x = "Time", y = "Moment") +
  #coord_cartesian(xlim = c(3, 5),ylim=c(0.27,0.7)) +  # Adjust these values as needed
  #facet_wrap(~ h, labeller = label_bquote(~ H == .(h)~ ", k" == 2), nrow = 1) +
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

plot_5 <- ggplot() +
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
  scale_y_log10(breaks = scales::log_breaks(n = 6),
                 labels = scales::label_number(accuracy = 0.05)) +
  labs(x = "", y = "") +
  coord_cartesian(xlim = c(3.5, 5),ylim=c(0.28,0.55)) +  # Adjust these values as needed
  #facet_wrap(~ h, labeller = label_bquote(~ H == .(h)~ ", k" == 2), nrow = 1) +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        #panel.grid.major = element_line(color = "grey80", linewidth = 0.5),
        #panel.grid.minor = element_line(color = "grey90", linewidth = 0.25),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.75),
        panel.spacing = unit(0.01, "lines"),
        legend.position = "nono",
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
final_plot <- ggdraw() +
  draw_plot(plot_4) +
  draw_plot(plot_5, x = 0.3, y = 0.466, width = 0.7, height = 0.485) 
final_plot

############
## Table 1 #
############
                          
g_t <- expression(t+1)
p_vals <- c(1/4,seq(1/2, 9/2, by=1))
H_vals <- seq(0.5, 1, by=0.1)
t_vals <- seq(0, 5, by=1)
results <- expand.grid(p=p_vals, H=H_vals, t=t_vals) %>%
  arrange(p, H)
results$Moment <- round(mapply(HO_mom, x0=10, gt=list(g_t),k=0, Alpha=-2, Beta=0, Gamma=1, 
                         H=results$H, p=results$p, t=results$t)^(1/results$p),4)
results_wide <- results %>%
  pivot_wider(names_from = "t", values_from = "Moment") %>%
  arrange(p) %>%
  mutate(across(where(is.numeric), ~ round(.x, 8)))
print(results_wide,n=30)

##################
## plot: Figure 5
###################
results_filtered <- results_wide %>%
  filter(p == 0.25, H %in% c(0.5, 0.6, 0.7, 0.8, 0.9, 1))

results_long <- results_wide %>%
  pivot_longer(cols = -c(p, H), names_to = "t", values_to = "Moment") %>%
  mutate(t = as.numeric(t))  # Convertir t en numérique

unique_H <- unique(results_long$H)  # Extraire les valeurs uniques de H
Color <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#666666")
names(Color) <- as.character(unique_H)

format_p_as_fraction <- function(p) {
  fractions <- c("1/4", "1/2", "3/2", "5/2", "7/2", "9/2")
  p_values <- c(1/4, 1/2, 3/2, 5/2, 7/2, 9/2)
  formatted_p <- setNames(fractions, p_values)[as.character(p)]
  return(paste0("p = ", formatted_p))  # Ajoute "p = " devant chaque fraction
}

format_H <- function(H) {
  return(sprintf("%.1f", as.numeric(H)))  # Convertit H en nombre avec 1 chiffre après la virgule
}

ggplot(results_long, aes(x = t, y = Moment, color = factor(H))) +
  geom_line(linewidth = 0.6) +
  scale_color_manual(values = Color, labels = format_H) +  
  scale_x_continuous(breaks = pretty_breaks(n = 12),
                     labels = number_format(accuracy = 0.1)) +
  scale_y_log10(breaks = scales::log_breaks(n = 10),
                labels = scales::label_number(accuracy = 0.05)) +
  coord_cartesian(xlim = c(3, 5),ylim=c(0.25,0.72)) +
  labs(x = "Time", y = "Moment", color = "H") +
  facet_wrap(~p, labeller = labeller(p=format_p_as_fraction), nrow = 2) +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        panel.grid.major = element_line(color = "grey80", linewidth = 0.5),
        panel.grid.minor = element_line(color = "grey90", linewidth = 0.25),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.75),
        panel.spacing = unit(0.01, "lines"),
        legend.position = "top",
        legend.title = element_text(size = 18, face = "bold", hjust = 1, margin = margin(t = 1)), 
        legend.text = element_text(size = 18, face = "bold"),
        legend.box.margin = margin(t = -11), 
        legend.box.spacing = unit(0.01, "cm"), 
        legend.spacing.y = unit(0.01, "cm"),
        axis.title.x = element_text(face = "bold"),
        axis.title.y = element_text(face = "bold"),
        plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm"),
        strip.text = element_text(size = 15, face = "bold.italic", margin = margin(1, 1, 1, 1)),
        strip.background = element_rect(fill = "#f0f0f0", color = "black", linewidth = 0.75))

###################
## plot: Figure 6 #
###################

g_t <- expression(t+1)
p_vals <- c(1/4,1/3,1/2,3/4,1,3/2,5/2,7/2,9/2)
H_vals <- seq(0.5, 1, by=0.1)
t_vals <- seq(0, 5, by=1)
results <- expand.grid(p=p_vals, H=H_vals, t=t_vals) %>%
  arrange(p, H)
results$Moment <- round(mapply(HO_mom, x0=10, gt=list(g_t),k=0, Alpha=-2, Beta=0, Gamma=1, 
                               H=results$H, p=results$p, t=results$t)^(1/results$p),4)

results_t5 <- results %>% filter(t == 5)
format_p_as_fraction <- function(p) {
  fractions <- c("1/4","1/3", "1/2","3/4", "1.0", "3/2", "5/2", "7/2", "9/2")
  p_values <- c(1/4,1/3, 1/2,3/4,1.0, 3/2, 5/2, 7/2, 9/2)
  formatted_p <- setNames(fractions, p_values)[as.character(p)]
  return(paste0(formatted_p)) 
}

ggplot(results_t5, aes(x = H, y = Moment, color = factor(p), group = factor(p))) +
  geom_line(linewidth = 0.8) +  
  geom_point(size = 1.2) +  
  scale_color_manual(
    values = RColorBrewer::brewer.pal(n = length(unique(results_t5$p)), name = "Set1"),
    name = expression(p),
    labels = format_p_as_fraction(unique(results_t5$p))
  ) +
  scale_y_continuous(breaks = pretty_breaks(n = 6),
                     labels = number_format(accuracy = 0.02)) +
  labs(x = "H (Hurst Parameter)",y = "Moment") +
  theme_minimal()+
  theme(panel.grid = element_blank(),
        panel.grid.major = element_line(color = "grey80", linewidth = 0.5),
        panel.grid.minor = element_line(color = "grey90", linewidth = 0.25),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.75),
        panel.spacing = unit(0.01, "lines"),
        legend.position = "top",
        legend.title = element_text(size = 18, face = "bold", hjust = 1, margin = margin(t = 1)), 
        legend.text = element_text(size = 18, face = "bold"),
        legend.box.margin = margin(t = -11), 
        legend.box.spacing = unit(0.01, "cm"), 
        legend.spacing.y = unit(0.01, "cm"),
        axis.title.x = element_text(face = "bold"),
        axis.title.y = element_text(face = "bold"),
        plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm"),
        strip.text = element_text(size = 15, face = "bold.italic", margin = margin(1, 1, 1, 1)),
        strip.background = element_rect(fill = "#f0f0f0", color = "black", linewidth = 0.75))                          

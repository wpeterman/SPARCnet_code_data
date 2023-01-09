## Code for projecting and plotting salamander growth curves
## Written by Bill Peterman
## Updated 9 January 2023
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# ** Growth Curves -----------------------------------------------------------
## Import growth model results
growth_df <- readRDS("fitted_model_df.rds")

n_ind <- 25000 # Number of individuals
int <- 0.05
yrs <- seq(1,10 + 1, by = int)

svl <- array(NA, c(n_ind, length(yrs)))
out_list <- vector('list', 4)
names(out_list) <- c('slope', 'ridge',
                     'slope_m', 'ridge_m')

for(i in 1:4){
  if(i == 1){ # Slope female
    sex <- 2
    pos_sex <- 7
    pos <- 3
    male <- 0
  } else if(i == 2)  { # Ridge female
    sex <- 2
    pos_sex <- 9
    pos <- 4
    male <- 0
  } else if(i == 3){ # Slope Male
    sex <- 1
    pos_sex <- 6 
    pos <- 3
    male <- 1
  } else {       # Ridge Male
    sex <- 1
    pos_sex <- 8
    pos <- 4
    male <- 1
  }
  cnt <- 0
  for(y in yrs){
    cnt <- cnt + 1
    samp <- sort(sample(1:nrow(growth_df), n_ind, replace = F))
    
    if(cnt == 1){
      svl[,cnt] <- rnorm(n_ind, 13.5, 1/growth_df$sigma.SVL[samp]^2)
      
    } else {
      ## Apply growth model
      K_ <- exp(growth_df[samp,pos] + growth_df[samp,5] * male)
      
      
      ## Data SD
      svl[,cnt] <- rnorm(n = n_ind, svl[,cnt-1], 1/growth_df$sigma.SVL[samp]^2) + (growth_df[samp,sex] - svl[,cnt-1]) * (1 - exp(-K_ * int))
      
      
    } 
  } # End year loop
  out_list[[i]] <- list(svl = svl)
} # End location loop

# Growth Plot -------------------------------------------------------------

# *** Data prep -----------------------------------------------------------

slope_growth_m <- out_list$slope_m$svl
slope_growth_f <- out_list$slope$svl
ridge_growth_m <- out_list$ridge_m$svl
ridge_growth_f <- out_list$ridge$svl

apply(slope_growth_m, 2, mean) # 13 * 0.05 = 0.65
apply(slope_growth_f, 2, mean) # 28 * 0.05 = 1.4
apply(ridge_growth_m, 2, mean) # 9 *0.05 = 0.45
apply(ridge_growth_f, 2, mean) # 20 * 0.05 = 1

slope_m <- data.frame(sex = 'male',
                      loc = 'Mature forest',
                      sex_loc = 'male_slope',
                      year = yrs,
                      svl_mn = apply(slope_growth_m, 2, mean),
                      lci = apply(slope_growth_m, 2, function(x) quantile(x, 0.025)),
                      uci = apply(slope_growth_m, 2, function(x) quantile(x, 0.975)))

slope_f <- data.frame(sex = 'female',
                      loc = 'Mature forest',
                      sex_loc = 'female_slope',
                      year = yrs,
                      svl_mn = apply(slope_growth_f, 2, mean),
                      lci = apply(slope_growth_f, 2, function(x) quantile(x, 0.025)),
                      uci = apply(slope_growth_f, 2, function(x) quantile(x, 0.975)))

ridge_m <- data.frame(sex = 'male',
                      loc = 'Successional forest',
                      sex_loc = 'male_ridge',
                      year = yrs,
                      svl_mn = apply(ridge_growth_m, 2, mean),
                      lci = apply(ridge_growth_m, 2, function(x) quantile(x, 0.025)),
                      uci = apply(ridge_growth_m, 2, function(x) quantile(x, 0.975)))

ridge_f <- data.frame(sex = 'female',
                      loc = 'Successional forest',
                      sex_loc = 'female_ridge',
                      year = yrs,
                      svl_mn = apply(ridge_growth_f, 2, mean),
                      lci = apply(ridge_growth_f, 2, function(x) quantile(x, 0.025)),
                      uci = apply(ridge_growth_f, 2, function(x) quantile(x, 0.975)))

growth_plot <- rbind(slope_m, slope_f,
                     ridge_m, ridge_f)

library(ggplot2)
library(cowplot)


# *** Male Plot -----------------------------------------------------------

(p_m <- ggplot(data = growth_plot[growth_plot$sex == 'male' & growth_plot$year <= 8,],
               aes(x = year, y = svl_mn, fill = loc, linetype = loc)) +
   geom_line(size = 1.) +
   scale_color_manual(values=c("#666666", "#bb0000")) +
   geom_ribbon(aes(ymin=lci, ymax=uci),
               alpha = 0.5) +
   scale_fill_manual(values = c("#666666", "#bb0000")) +
   geom_hline(yintercept = 34, linetype="dotted", 
              color = "black", size = 0.75) +
   scale_x_continuous(n.breaks = 8) +
   scale_y_continuous(n.breaks = 7) +
   labs(x = 'Years', y = 'SVL (mm)') +
   theme_cowplot() +
   theme(legend.position = 'none')
)

(male_growth.plot <- p_m + geom_segment(aes(x = 2.2, y = 13, xend = 2.2, yend = 34),
                                        linetype = "dashed", color = "black", size = 0.5) +
    geom_segment(aes(x = 2.75, y = 13, xend = 2.75, yend = 34),
                 linetype = "dashed", color = "black", size = 0.5) +
    annotate('text', 
             x = 1.65,
             y = 14,
             label= '2.25 years',
             size = 4.5) +
    annotate('text', 
             x = 3.4,
             y = 14,
             label= '2.75 years',
             size = 4.5) +
    annotate('text', 
             x = 5,
             y = 35,
             label= 'Minimum size at maturity',
             size = 4.5))



# *** Female plot ---------------------------------------------------------

(p_f <- ggplot(data = growth_plot[growth_plot$sex == 'female' & growth_plot$year <= 8,],
               aes(x = year, y = svl_mn, fill = loc, linetype = loc)) +
   geom_line(size = 1.) +
   scale_color_manual(values=c("#666666", "#bb0000")) +
   geom_ribbon(aes(ymin=lci, ymax=uci),
               alpha = 0.5) +
   scale_fill_manual(values = c("#666666", "#bb0000")) +
   geom_hline(yintercept = 34, linetype="dotted", 
              color = "black", size = 0.75) +
   scale_x_continuous(n.breaks = 8) +
   scale_y_continuous(n.breaks = 8) +
   labs(x = 'Years', y = 'SVL (mm)') +
   theme_cowplot() +
   theme(legend.position = c(0.01, 0.95),
         legend.title = element_blank(),
         axis.title.x = element_blank(),
         axis.text.x = element_blank())
)

(female_growth.plot <- p_f + geom_segment(aes(x = 3.3, y = 13, xend = 3.3, yend = 34),
                                          linetype = "dashed", color = "black", size = 0.5) +
    geom_segment(aes(x = 4.3, y = 13, xend = 4.3, yend = 34),
                 linetype = "dashed", color = "black", size = 0.5) +
    annotate('text', 
             x = 2.75,
             y = 14,
             label= '3.30 years',
             size = 4.5) +
    annotate('text', 
             x = 4.9,
             y = 14,
             label= '4.30 years',
             size = 4.5) +
    annotate('text', 
             x = 1.95,
             y = 36.5,
             label= 'Minimum size at maturity',
             size = 4.5))


## Combine plots
(mf_plot <- plot_grid(female_growth.plot, male_growth.plot,
                      labels = c('A', 'B'),
                      nrow = 2,
                      label_size = 14))



# Ridgeline Plot ----------------------------------------------------------

# >> Growth Rate ----------------------------------------------------------

summary(growth_df)

ridge_df <- data.frame('Mature Male' = exp(growth_df$K_ridge_m),
                       'Mature Female' = exp(growth_df$K_ridge_f),
                       'Successional Male' = exp(growth_df$K_slope_m),
                       'Successional Female' = exp(growth_df$K_slope_f))
ridge_long <- reshape::melt(ridge_df[1:10000,])
colnames(ridge_long) <- c('sex_loc', 'K')

ridge_long <- tidyr::separate(ridge_long,
                              col = sex_loc,
                              into = c("loc", "sex"))

library(ggridges)
library(ggplot2)
library(cowplot)

(ridge_plot_sex.loc <- ggplot(ridge_long, aes(y = loc, x = K,
                                              fill = sex)) +
    geom_density_ridges(alpha = 0.5, scale = 3) + 
    theme_cowplot() +
    theme(legend.position = 'none',
          panel.grid.major.x = element_line(colour = 'lightgray')) + 
    scale_fill_manual(values = c("#666666", "#bb0000")) +
    labs(x = expression(paste(bold("Growth coefficient ("), bolditalic("K"),bold(")"))),
         y = 'Location',
         fill = "Sex") +
    theme(axis.title = element_text(size = 13, face = "bold"), 
          axis.text = element_text(size = 12)) 
) 

ridge_plot_sex.loc + theme(legend.position = c(0.88, 0.9)) +labs(fill = "Sex")

# >> Asymptotic size ----------------------------------------------------------

ridge_size_df <- data.frame('Male' = (growth_df$size_male),
                            'Female' = (growth_df$size_female))
# ridge_size_df <- data.frame('Male' = (growth_df$size_slope),
#                             'Female' = (growth_df$size_ridge))
ridge_size_long_df <- reshape::melt(ridge_size_df[1:10000,])
colnames(ridge_size_long_df) <- c('Sex', 'Size')

library(ggridges)
library(ggplot2)
library(cowplot)

(ridge_plot_size <- ggplot(ridge_size_long_df, aes(y = Sex, x = Size, fill = Sex)) +
    geom_density_ridges(alpha = 0.5, scale = 3) + 
    theme_cowplot() +
    scale_fill_manual(values = c("#bb0000", "#666666")) +
    theme(legend.position = c(0.78, 0.85),,
          panel.grid.major.x = element_line(colour = 'lightgray')) + 
    theme(axis.title = element_text(size = 13, face = "bold"), 
          axis.text = element_text(size = 12)) +
    labs(x = "Asymptotic SVL (mm)",
         y = 'Sex') 
) 

## Combine plots
(ridge_plots <- plot_grid(ridge_plot_size, ridge_plot_sex.loc,
                          labels = c('A', 'B'),
                          nrow = 2,
                          label_size = 14))
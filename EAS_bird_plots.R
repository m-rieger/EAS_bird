########################-
#### EAS bird plots ####
########################-

## plot functions for posterior predictions

## 1) plot abundance

## plot abundance trend
plot.ab <- function(df.raw = NULL, newdat = NULL, group = NULL) {
  
  g <- ggplot() +
    
    # add horizontal line at y = 0
    geom_hline(yintercept = 0, lty = "dashed", linewidth = 0.5) +
    
    # add raw data (colored dots)
    geom_beeswarm(data = df.raw, aes(x = year, y = ab2/2, color = nat.region),
                  size = 0.5, cex = 0.3, alpha = 0.5) +
    
    # add trend
    geom_ribbon(data = tmp,
                aes(x = year, ymin = lwr, ymax = upr, 
                    color = nat.region, fill = nat.region),
                alpha = 0.3) +
    geom_line(data = tmp,
              aes(x = year, y = fit, color = nat.region),
              lwd = 1) +
    
    # add sign. slope
    geom_line(aes(year, fit.slope, group = nat.region),
              data = tmp,
              col = "black", lty = 1, linewidth = 1, na.rm = TRUE) +
    
    # add labs
    ylim(0, max(c(df.raw$ab2/2, tmp$upr))) +
    labs(x = "year", y = "abundance per km²",
         title = paste0(sp, " (", mod.fam, ", ", zi.type, ")"),
         subtitle = "black = period of robust increase/decrease") +
    facet_wrap(~nat.region) +
    theme_classic() +
    theme(legend.position = "none")

  ## add colors
  if(group == "nat.region") {
    g <- g +
      scale_color_viridis_d("nat. region", end = 0.8, option = "plasma") +
      scale_fill_viridis_d("nat. region", end = 0.8, option = "plasma")
  }

  if(group == "region") {
    g <- g +
      scale_color_manual("region", values = c("atl" = "blue", "kon" = "orange", 
                                              "overall" = "grey40")) +
      scale_fill_manual("region",  values = c("atl" = "blue", "kon" = "orange", 
                                              "overall" = "grey40"))
  }
  
  return(print(g))
  
}


## plot index

plot.in <- function(newdat = NULL, group = NULL) {
  
  g <- ggplot(newdat, aes(x = year, y = fit.i, ymin = lwr.i, ymax = upr.i)) +
    
    # add horizontal lines at y = 0 and y = 1
    geom_hline(yintercept = 0, lty = "dashed", lwd = 0.5) +
    geom_hline(yintercept = 1, lty = "dashed", lwd = 0.5) +
    
    # add trend
    geom_ribbon(aes(color = nat.region, fill = nat.region), alpha = 0.2) +
    geom_line(aes(color = nat.region), lwd = 1) +
    
    # add sign. slope
    geom_line(aes(year, fit.slope.i, group = nat.region),
              col = "black", lty = 1, linewidth = 1, na.rm = TRUE) +
    
    # add labs
    ylim(0, max(tmp$upr.i)) +
    labs(y = "Index", x = "survey year",
         title = paste0(sp, " (", mod.fam, ", ", zi.type, ")")) + 
    facet_wrap(~nat.region) +
    theme_classic() +
    theme(legend.position = "none")
  
  ## add colors
  if(group == "nat.region") {
    g <- g +
      scale_color_viridis_d("nat. region", end = 0.8, option = "plasma") +
      scale_fill_viridis_d("nat. region", end = 0.8, option = "plasma")
  }
  
  if(group == "region") {
    g <- g +
      scale_color_manual("region", values = c("atl" = "blue", "kon" = "orange", 
                                              "overall" = "grey40")) +
      scale_fill_manual("region",  values = c("atl" = "blue", "kon" = "orange", 
                                              "overall" = "grey40"))
  }
  
  return(print(g))
}


## plot longterm trend

plot.lt <- function(newdat = NULL, group = NULL) {
  
  g <- ggplot(newdat) +
    
    # add horizontal line at 0
    geom_hline(yintercept = 0, lty = "dashed", linewidth = 0.5) +
    
    # add longterm trend
    geom_pointrange(aes(x = nat.region, y = fit, ymin = lwr, ymax = upr, color = nat.region),
                    linewidth = 1, fatten = 8) +
    geom_pointrange(aes(x = nat.region, y = fit, ymin = lwr25, ymax = upr75, color = nat.region),
                    linewidth = 2, fatten = 8) +
    
    # add labs
    labs(x = group,
         y = "mean annual change in abundance per km²",  
         title = paste0(sp, " (", mod.fam, ", ", zi.type, ")"),
         subtitle = paste0("Longterm trend (", period, " years, ", yearP - period, "-", yearP, ") with 50% (thick) and 95% (thin) CrI")) +
    theme_classic() +
    theme(legend.position = "none")

  ## add colors
  if(group == "nat.region") {
    g <- g +
      scale_color_viridis_d("nat. region", end = 0.8, option = "plasma")
  }
  
  if(group == "region") {
    g <- g +
      scale_color_manual("region", values = c("atl" = "blue", "kon" = "orange", 
                                              "overall" = "grey40"))
  }
  
  return(print(g))
}

library(tidyverse)
options(stringsAsFactors = F)
library(plotly)
source("plot_theme.R")

landScape <- function(dfTraj, nam, pThresh = -3)
{
    kd <- with(dfTraj, MASS::kde2d(EMTScore,Immune,n = 200))
    names(kd) <- c("EMTScore", "Immune", "Probability")
    kd$Probability <- log10(kd$Probability)
    kd$Probability[kd$Probability < pThresh] <- NA
    colorAxis <- kd$Probability
    colorAxis[kd$EMTScore < -0.25, ] <- -1
    colorAxis[kd$EMTScore > 0.5, ] <- 1
    colorAxis[kd$EMTScore <=0.5& kd$EMTScore >=-0.25] <- 0
    colorAxis <- apply(colorAxis, 2, function(x){
        as.character(x)
    })
    
    colorAxis2 <- kd$Probability
    colorAxis2[kd$Immune < 0, ] <- -1
    colorAxis2[kd$Immune >= 0, ] <- 1
    colorAxis2 <- apply(colorAxis2, 2, function(x){
        as.character(x)
    })
    tryCatch({
        p <- plot_ly(x = kd$EMTScore, y = kd$Immune, z = -kd$Probability,
                     colors = c("grey", "orange", "red")) %>% 
            add_surface(surfacecolor = colorAxis, 
                        contours = list(
                            z = list(
                                show=TRUE,
                                usecolormap=TRUE,
                                highlightcolor="#ff0000",
                                project=list(z=TRUE)
                            )
                        )) %>%
            layout(
                scene = list(
                    camera=list(
                        eye = list(x=1.1, y=-1, z=0.64)
                    )
                )
            )
        
        orca(p, paste0(nam, "_EMTCol.png"), width = 1280, height = 720)
    }, error = function(e){conditionMessage(e)})
}

# Generated trajectory plots for multiple initial conditions
trajectoriesnIC <- function(dfTraj, nam)
{
    dfInit <- dfTraj %>% filter(Time == max(dfTraj$Time))
    dfInit$State <- "Epithelial"
    dfInit$State[dfInit$EMTScore > 0.5] <- "Mesenchymal"
    dfInit$State[dfInit$EMTScore >= -0.25 & dfInit$EMTScore <= 0.5] <- "Hybrid"
    dfTraj <- merge(dfTraj, 
                    dfInit %>% select(Id, State), by = "Id",all = T)
    dfNew <- dfTraj %>% select(Id, Time, EMTScore, Immune, State) %>% 
        gather(key = "Metric", value = "Score", - Time, -Id, -State)
    dfNew$State <- factor(dfNew$State, levels = c("Epithelial", "Hybrid", "Mesenchymal"))
    p <- ggplot()
    for (i in dfNew$Id %>% unique)
    {
        p <- p + geom_line(data = dfNew %>% filter(Id == i), mapping = aes(x = Time, y = Score, color = State))
    }
    p <- p + facet_wrap(~Metric) + 
        theme_Publication() + ylim(-2,1) + 
        scale_color_manual(values = c('#7F7F7F', "#FFD17F", '#FF7F7F')) +
        theme(legend.position = "top", panel.grid = element_blank())
            
    ggsave(plot = p, filename = paste0(nam, "_trajecs.png"))
}

trajPlots <- function(topoFile, workers)
{
    net <- topoFile %>% str_remove(".topo")
    setwd(net)
    filz <- list.files(".",".csv")
    plan(multiprocess, workers = workers)
    dummy <- future_lapply(filz, function(x){
        nam <- x %>% str_remove("_simDat.csv")
        dfTraj <- read.csv(x)
        trajectoriesnIC(dfTraj, nam)
    })
    future:::ClusterRegistry("stop")
    setwd("..")
    
}

#Generates trajectory plots for switching.
trajectoriesSwitching <- function(dfTraj, id, tMin = 0, tMax = NULL)
{
    d <- dfTraj %>% filter(Id == id) %>% select(Time, EMTScore, Immune) %>% 
        gather(key = "Metric", value = "Score", -Time)
    p <- ggplot(d, aes(x = Time, y = Score, color = Metric)) + geom_line() + theme_Publication() + 
        theme(legend.position = "top")
    ggsave(paste0(id,"Normal.png"))
    print(p)
    if (!is.null(tMax))
    {
        p <- ggplot(d %>% filter(Time >=tMin, Time<=tMax), aes(x = Time, y = Score, color = Metric)) + 
            geom_line() + theme_Publication() + 
            theme(legend.position = "top")
        ggsave(paste0(id, "Zoom.png"))
    }
}






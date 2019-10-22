### Plotting some Magnetization
args <- commandArgs(trailingOnly=TRUE)
subdir <- args[1]
setwd(subdir)

plot_M <- function(title, df) {
  t <- seq(0, nrow(df)-1)
  plot(x = t, y = df$Mps, type = "l", col="green", main = title, ylim=c(-1,1),
       ylab = "Instantaneous quantity", xlab = "time steps")
  #legend(x="topleft", y=0.95, legend=c("M", "E"), col=c("green", "red"), pch=c("-","-"))
}

files <- list.files(pattern=".csv")

for (f in files) {
  df <- read.csv(f)
  png(paste0(substr(f, 1, nchar(f)-4), ".png"))
  plot_M(paste0("Plot of Magnetisation for ", f), df)
  dev.off()
}


### old: M & E
# get_ylim <- function(M, E) {c(min(min(M),min(E)), max(max(M),max(E)))}
# 
# plot_ME <- function(title, df) {
#   t <- seq(0, nrow(df)-1)
#   plot(x = t, y = df$M, type = "l", col="green", 
#        ylim = get_ylim(df$M, df$E), main = title, ylab = "Instantaneous quantity", xlab = "time steps")
#   lines(x = t, y = df$E, type = "l", col="red")
#   legend(x="topleft", y=0.95, legend=c("M", "E"), col=c("green", "red"), pch=c("-","-"))
# }
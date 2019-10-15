### Plotting some E and M
get_ylim <- function(M, E) {c(min(min(M),min(E)), max(max(M),max(E)))}

plot_ME <- function(title, df) {
  t <- seq(0, nrow(df)-1)
  plot(x = t, y = df$M, type = "l", col="green", 
       ylim = get_ylim(df$M, df$E), main = title, ylab = "Instantaneous quantity", xlab = "time steps")
  lines(x = t, y = df$E, type = "l", col="red")
  legend(x="topleft", y=0.95, legend=c("M", "E"), col=c("green", "red"), pch=c("-","-"))
}

files <- list.files(pattern=".csv")

for (f in files) {
  df <- read.csv(f)
  png(paste0(substr(f, 1, nchar(f)-4), ".png"))
  plot_ME(paste0("Plot of Magnetisation and Energy for ", f), df)
  dev.off()
}
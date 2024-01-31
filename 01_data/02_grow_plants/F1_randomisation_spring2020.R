set.seed(8)

lines <- read.csv("03_randomisation/F1_lines_to_grow.csv")

ix <- lapply(1:6, function(x) sample(1:16, replace = F))
ix <- do.call('c', ix)

d <- data.frame(
  tray = rep(1:4, each=24),
  pos  = rep(1:24, 4),
  lines[ix,]
)

write.csv(d, file="data_processed/F1_randomisation.csv")
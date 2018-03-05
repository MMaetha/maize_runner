data <- read.csv("~/GitHub/maize_runner/TE_search07_chr1.csv", colClasses = "numeric", header=FALSE)
hist(data$v1)
dt2 <- table(data)

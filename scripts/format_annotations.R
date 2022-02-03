full_annotations <- read.csv("./annotdata.csv")

full_annotations$package_name <- as.factor(full_annotations$package_name)
full_annotations$version <- as.factor(full_annotations$version)

rownames(full_annotations) <- full_annotations$ID
full_annotations$ID <- NULL

save(full_annotations, file = "./full_annotations.RData")
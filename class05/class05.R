#' ---
#' title: "Class 5: R Graphics"
#' author: "Mary Garcia"
#' date: "April 18. 2019"
#' ---

# Class 5 R graphics

# 2A. Line plot
weight <- read.table("bimm143_05_rstats/weight_chart.txt", header = TRUE)
# above to use console: command+return
# default is false for header -> need to change to TRUE if 1st line is header

weight$Age
weight$Weight

plot(weight, pch=15, cex=1.5, lwd=2, ylim=c(2,10), xlab="Age (months)", ylab="Weight (kg)", main="Baby Weight by Month")

plot(weight$Age, weight$Weight, type = "b", pch=15, cex=1.5, lwd=2, ylim=c(2,10), xlab="Age (months)", ylab="Weight (kg)", main="Baby Weight by Month")
# Above is more concise -> favorable if there are more variables

#2B. Line plot
feat <- read.table("bimm143_05_rstats/feature_counts.txt", 
                   header= TRUE, sep = "\t")
# sep="\t" separates info by tab

barplot(feat$Count)
# I need to argue w this plot to make it nicer

old.par <- par(mar = c(0, 0, 0, 0))
par(mar=c(4,11,1,1))
barplot(feat$Count, horiz = TRUE, xlab="", names.arg = feat$Feature,
        main="Number of features in the mouse GRCm38 genome", las=1, xlim = c(0,80000))
par(mar=old.par)

#?? use par (the parameter) in the mouse plot

# Section 3
counts <- read.table("bimm143_05_rstats/male_female_counts.txt",
                     sep="\t", header=TRUE)

# use read.delim() to get same result as above
counts <-read.delim("bimm143_05_rstats/male_female_counts.txt")

barplot(counts$Count, names.arg = counts$Sample, las=2,
        col=rainbow(10))

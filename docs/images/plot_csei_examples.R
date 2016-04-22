#! /usr/bin/Rscript

library(ggplot2)
library(grid)

#https://stackoverflow.com/questions/25106508/ggplot2-is-there-an-easy-way-to-wrap-annotation-text
wrapper <- function(x, ...) paste(strwrap(x, ...), collapse = "\n")

#Setup plot skeleton
p <- ggplot() + ylim(50, 110)  + theme_classic() + theme(axis.line=element_blank(), axis.text = element_blank()) + xlab("") + ylab("") + theme(axis.ticks = element_blank()) + ggtitle("")

#Topmost transcript
p <- p + geom_rect(fill = "#3363CC", aes(xmin = 5, ymin = 92.5, xmax = 10, ymax = 97.5)) + geom_rect(fill = "#3363CC", aes(xmin = 10, ymin = 90, xmax = 30, ymax = 100)) + geom_rect(fill = "#3363CC", aes(xmin = 40, ymin = 90, xmax = 60, ymax = 100)) + geom_rect(fill = "#3363CC", aes(xmin = 70, ymin = 90, xmax = 90, ymax = 100)) + geom_rect(fill = "#3363CC", aes(xmin = 100, ymin = 90, xmax = 120, ymax = 100))    + geom_rect(fill = "#3363CC", aes(xmin = 130, ymin = 90, xmax = 150, ymax = 100)) + geom_rect(fill = "#3363CC", aes(xmin = 150, ymin = 92.5, xmax = 155, ymax = 97.5)) + geom_segment(aes(x = 30, y = 95, xend = 40, yend = 95), arrow = arrow(length = unit(0.5, "cm"))) + geom_segment(aes(x = 60, y = 95, xend = 70, yend = 95), arrow = arrow(length = unit(0.5, "cm")))+ geom_segment(aes(x = 90, y = 95, xend = 100, yend = 95), arrow = arrow(length = unit(0.5, "cm")))+ geom_segment(aes(x = 120, y = 95, xend = 130, yend = 95), arrow = arrow(length = unit(0.5, "cm")))

#Plot variants
p <- p + geom_rect(fill = "#006600", aes(xmin = 28, ymin = 74, xmax = 30, ymax = 79))

#Plot junctions
p <- p + geom_rect(aes(xmin = 20, ymin = 55, xmax = 30, ymax = 59), fill = "red") + geom_rect(aes(xmin = 70, ymin = 55, xmax = 80, ymax = 59), fill = "red") + geom_segment(aes(x = 30, y = 57, xend = 70, yend = 57), arrow = arrow(length = unit(0.5, "cm")))

#Plot dotted vertical lines
p <- p + geom_vline(xintercept = 30, linetype = "longdash") + geom_vline(xintercept = 70, linetype = "longdash")

#Add labels
p <- p + annotate("text", label = "Transcript annotation", x = 120, y = 105, size = 5)
p <- p + annotate("text", label = "Variant called in the genomic data", x = 120, xmin = 10, y = 85, size = 5)
p <- p + annotate("text", label = "Spliced junction identified by regtools", x = 120, y = 65, size = 5)


print(p)
ggsave("csei_examples.png")
ggsave("csei_examples.pdf")

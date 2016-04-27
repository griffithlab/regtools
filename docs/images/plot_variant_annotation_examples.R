#! /usr/bin/Rscript

library(ggplot2)
library(grid)

#https://stackoverflow.com/questions/25106508/ggplot2-is-there-an-easy-way-to-wrap-annotation-text
wrapper <- function(x, ...) paste(strwrap(x, ...), collapse = "\n")

#Setup plot skeleton
p <- ggplot() + ylim(30, 50) + theme_classic() + theme(axis.line=element_blank(), axis.text = element_blank()) + xlab("") + ylab("")+ theme(axis.ticks = element_blank()) + scale_x_continuous(breaks=-10:70) + ggtitle("")

#Topmost transcript
p <- p + geom_rect(fill = "#3363CC", aes(xmin = 10, ymin = 40, xmax = 30, ymax = 45)) + geom_rect(fill = "#3363CC", aes(xmin = 40, ymin = 40, xmax = 60, ymax = 45)) + geom_segment(aes(x = 30, y = 42.5, xend = 40, yend = 42.5), arrow = arrow(length = unit(0.5, "cm")))

#Plot variants
p <- p + geom_rect(fill = "#006600", aes(xmin = 22, ymin = 34, xmax = 23, ymax = 35), fill = "red")
p <- p + geom_rect(fill = "#006600", aes(xmin = 27, ymin = 33, xmax = 28, ymax = 34), fill = "red")
p <- p + geom_rect(fill = "#006600", aes(xmin = 31, ymin = 32, xmax = 32, ymax = 33), fill = "red")
p <- p + geom_rect(fill = "#006600", aes(xmin = 35, ymin = 31, xmax = 36, ymax = 32), fill = "red")

#Add labels
p <- p + annotate("text", label = "Transcript annotation", x = 20, xmin = 10, y = 46.5, size = 5)
p <- p + annotate("text", label = "Variants", x = 20, xmin = 10, y = 36.5, size = 5)
p <- p + annotate("text", label = "Annotations", x = 45, xmin = 44, y = 36.5, size = 5)
p <- p + annotate("text", label = "Distance=6;", x = 44, xmin = 44, y = 34.6, size = 4)
p <- p + annotate("text", label = "Distance=2;splicing_exonic", x = 48.5, xmin = 40, y = 33.6, size = 4)
p <- p + annotate("text", label = "Distance=1;splicing_intronic", x = 48.8, xmin = 40, y = 32.6, size = 4)
p <- p + annotate("text", label = "Distance=4;", x = 44, xmin = 44, y = 31.6, size = 4)

#Add lines to make it clear
p <- p + geom_vline(xintercept = 30, linetype = "longdash") + geom_vline(xintercept = 40, linetype = "longdash")

print(p)
ggsave("variant_annotation_examples.png")
ggsave("variant_annotation_examples.pdf")

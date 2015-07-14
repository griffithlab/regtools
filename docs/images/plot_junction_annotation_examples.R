#! /usr/bin/Rscript

library(ggplot2)
library(grid)

#https://stackoverflow.com/questions/25106508/ggplot2-is-there-an-easy-way-to-wrap-annotation-text
wrapper <- function(x, ...) paste(strwrap(x, ...), collapse = "\n")

#Setup plot skeleton
p <- ggplot() + ylim(10, 110) + xlim(-10, 200) + theme(axis.text = element_blank()) + xlab("") + ylab("")+ theme(axis.ticks = element_blank()) + ggtitle("Junction-annotation examples")

#Topmost transcript
p <- p + geom_rect(aes(xmin = 10, ymin = 90, xmax = 30, ymax = 100)) + geom_rect(aes(xmin = 40, ymin = 90, xmax = 60, ymax = 100)) + geom_rect(aes(xmin = 70, ymin = 90, xmax = 90, ymax = 100)) + geom_segment(aes(x = 30, y = 95, xend = 40, yend = 95), arrow = arrow(length = unit(0.5, "cm"))) + geom_segment(aes(x = 60, y = 95, xend = 70, yend = 95), arrow = arrow(length = unit(0.5, "cm")))
p <- p + geom_rect(aes(xmin = 10, ymin = 70, xmax = 30, ymax = 80)) + geom_rect(aes(xmin = 45, ymin = 70, xmax = 60, ymax = 80)) + geom_rect(aes(xmin = 70, ymin = 70, xmax = 90, ymax = 80)) + geom_segment(aes(x = 30, y = 75, xend = 45, yend = 75), arrow = arrow(length = unit(0.5, "cm"))) + geom_segment(aes(x = 60, y = 75, xend = 70, yend = 75), arrow = arrow(length = unit(0.5, "cm")))
p <- p + geom_rect(aes(xmin = 10, ymin = 50, xmax = 30, ymax = 60)) + geom_rect(aes(xmin = 40, ymin = 50, xmax = 65, ymax = 60)) + geom_rect(aes(xmin = 70, ymin = 50, xmax = 90, ymax = 60)) + geom_segment(aes(x = 30, y = 55, xend = 40, yend = 55), arrow = arrow(length = unit(0.5, "cm"))) + geom_segment(aes(x = 65, y = 55, xend = 70, yend = 55), arrow = arrow(length = unit(0.5, "cm")))

#Plot junctions
p <- p + geom_rect(aes(xmin = 20, ymin = 30, xmax = 30, ymax = 35), fill = "red") + geom_rect(aes(xmin = 70, ymin = 30, xmax = 80, ymax = 35), fill = "red") + geom_segment(aes(x = 30, y = 32.5, xend = 70, yend = 32.5), arrow = arrow(length = unit(0.5, "cm")))

#Add labels
p <- p + annotate("text", label = "Transcripts", x = 7, y = 105, size = 5)
p <- p + annotate("text", label = "Junction", x = 5, y = 42.5, size = 5)
p <- p + annotate("text", label = "Donors skipped = 1", x = 140, y = 32.5, size = 5)
p <- p + annotate("text", label = "Acceptors skipped = 1", x = 140, y = 22.5, size = 5)
p <- p + annotate("text",
                  label = wrapper("Exons skipped = 1 or 3 depending on if exons are merged or not", width = 40),
                  x = 140, y = 12.5, size = 5)


#Add lines to make it clear
p <- p + geom_vline(xintercept = 30, linetype = "longdash") + geom_vline(xintercept = 40, linetype = "longdash") + geom_vline(xintercept = 45, linetype = "longdash")+ geom_vline(xintercept = 65, linetype = "longdash")+ geom_vline(xintercept = 70, linetype = "longdash")

print(p)
ggsave("junction_annotation_examples.png")


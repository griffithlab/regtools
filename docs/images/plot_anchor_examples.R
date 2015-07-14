#! /usr/bin/Rscript

library(ggplot2)
library(grid)

p <- ggplot() + ylim(-10, 110) + xlim(-10, 200) + theme(axis.text = element_blank()) + xlab("") + ylab("")+ theme(axis.ticks = element_blank()) + ggtitle("Anchor-annotation examples")

#Topmost transcript
p <- p + geom_rect(aes(xmin = 10, ymin = 90, xmax = 30, ymax = 100)) + geom_rect(aes(xmin = 40, ymin = 90, xmax = 60, ymax = 100)) + geom_rect(aes(xmin = 70, ymin = 90, xmax = 90, ymax = 100)) + geom_rect(aes(xmin = 100, ymin = 90, xmax = 120, ymax = 100))    + geom_rect(aes(xmin = 130, ymin = 90, xmax = 150, ymax = 100)) + geom_segment(aes(x = 30, y = 95, xend = 40, yend = 95), arrow = arrow(length = unit(0.5, "cm"))) + geom_segment(aes(x = 60, y = 95, xend = 70, yend = 95), arrow = arrow(length = unit(0.5, "cm")))+ geom_segment(aes(x = 90, y = 95, xend = 100, yend = 95), arrow = arrow(length = unit(0.5, "cm")))+ geom_segment(aes(x = 120, y = 95, xend = 130, yend = 95), arrow = arrow(length = unit(0.5, "cm")))
p <- p + geom_rect(aes(xmin = -10, xmax = 10, ymin = 70, ymax = 80)) + geom_rect(aes(xmin = 70, xmax = 90, ymin = 70, ymax = 80)) + geom_rect(aes(xmin = 100, ymin = 70, xmax = 120, ymax = 80)) + geom_rect(aes(xmin = 130, ymin = 70, xmax = 150, ymax = 80)) + geom_segment(aes(x = 10, xend = 70, y = 75, yend = 75), arrow = arrow(length = unit(0.5, "cm"))) + geom_segment(aes(x = 90, y = 75, xend = 100, yend = 75), arrow = arrow(length = unit(0.5, "cm"))) + geom_segment(aes(x = 120, y = 75, xend = 130, yend = 75), arrow = arrow(length = unit(0.5, "cm")))
p <- p + geom_rect(aes(xmin = 10, ymin = 50, xmax = 30, ymax = 60)) + geom_rect(aes(xmin = 100, ymin = 50, xmax = 120, ymax = 60)) + geom_rect(aes(xmin = 130, ymin = 50, xmax = 150, ymax = 60)) + geom_segment(aes(x = 30, y = 55, xend = 100, yend = 55), arrow = arrow(length = unit(0.5, "cm"))) + geom_segment(aes(x = 120, y = 55, xend = 130, yend = 55), arrow = arrow(length = unit(0.5, "cm")))

#Plot junctions
p <- p + geom_rect(aes(xmin = 20, ymin = 30, xmax = 30, ymax = 35), fill = "red") + geom_rect(aes(xmin = 40, ymin = 30, xmax = 50, ymax = 35), fill = "red") + geom_segment(aes(x = 30, y = 32.5, xend = 40, yend = 32.5), arrow = arrow(length = unit(0.5, "cm")))
p <- p + geom_rect(aes(xmin = 20, ymin = 20, xmax = 30, ymax = 25), fill = "red") + geom_rect(aes(xmin = 70, ymin = 20, xmax = 80, ymax = 25), fill = "red") + geom_segment(aes(x = 30, y = 22.5, xend = 70, yend = 22.5), arrow = arrow(length = unit(0.5, "cm")))
p <- p + geom_rect(aes(xmin = 20, ymin = 10, xmax = 30, ymax = 15), fill = "red") + geom_rect(aes(xmin = 45, ymin = 10, xmax = 55, ymax = 15), fill = "red") + geom_segment(aes(x = 30, y = 12.5, xend = 45, yend = 12.5), arrow = arrow(length = unit(0.5, "cm")))
p <- p + geom_rect(aes(xmin = 10, ymin = 0, xmax = 20, ymax = 5), fill = "red") + geom_rect(aes(xmin = 70, ymin = 0, xmax = 80, ymax = 5), fill = "red") + geom_segment(aes(x = 20, y = 2.5, xend = 70, yend = 2.5), arrow = arrow(length = unit(0.5, "cm")))
p <- p + geom_rect(aes(xmin = -5, ymin = -10, xmax = 5, ymax = -5), fill = "red") + geom_rect(aes(xmin = 45, ymin = -10, xmax = 55, ymax = -5), fill = "red") + geom_segment(aes(x = 5, y = -7.5, xend = 45, yend = -7.5), arrow = arrow(length = unit(0.5, "cm")))

#Add labels
p <- p + annotate("text", label = "Transcripts", x = 7, y = 105, size = 5)
p <- p + annotate("text", label = "Junctions", x = 5, y = 42.5, size = 5)
p <- p + annotate("text", label = "DA", x = 90, y = 32.5, size = 5)
p <- p + annotate("text", label = "NDA", x = 90, y = 22.5, size = 5)
p <- p + annotate("text", label = "D", x = 90, y = 12.5, size = 5)
p <- p + annotate("text", label = "A", x = 90, y = 2.5, size = 5)
p <- p + annotate("text", label = "N", x = 90, y = -7.5, size = 5)


#Add lines to make it clear
p <- p + geom_vline(xintercept = 30, linetype = "longdash") + geom_vline(xintercept = 40, linetype = "longdash") + geom_vline(xintercept = 70, linetype = "longdash")

print(p)
ggsave("anchor_examples.png")


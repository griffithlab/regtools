#! /usr/bin/Rscript

library(ggplot2)
library(grid)

p <- ggplot() + ylim(-10, 110) + xlim(-20, 220) + theme_classic() + theme(axis.line=element_blank(), axis.text = element_blank()) + xlab("") + ylab("")+ theme(axis.ticks = element_blank()) + ggtitle("")

#Topmost transcript
p <- p + geom_rect(fill = "#3363CC", aes(xmin = 0, ymin = 92.5, xmax = 10, ymax = 97.5)) + geom_rect(fill = "#3363CC", aes(xmin = 10, ymin = 90, xmax = 30, ymax = 100)) + geom_rect(fill = "#3363CC", aes(xmin = 40, ymin = 90, xmax = 60, ymax = 100)) + geom_rect(fill = "#3363CC", aes(xmin = 70, ymin = 90, xmax = 90, ymax = 100)) + geom_rect(fill = "#3363CC", aes(xmin = 100, ymin = 90, xmax = 120, ymax = 100))    + geom_rect(fill = "#3363CC", aes(xmin = 130, ymin = 90, xmax = 150, ymax = 100)) + geom_rect(fill = "#3363CC", aes(xmin = 150, ymin = 92.5, xmax = 160, ymax = 97.5)) + geom_segment(aes(x = 30, y = 95, xend = 40, yend = 95), arrow = arrow(length = unit(0.5, "cm"))) + geom_segment(aes(x = 60, y = 95, xend = 70, yend = 95), arrow = arrow(length = unit(0.5, "cm")))+ geom_segment(aes(x = 90, y = 95, xend = 100, yend = 95), arrow = arrow(length = unit(0.5, "cm")))+ geom_segment(aes(x = 120, y = 95, xend = 130, yend = 95), arrow = arrow(length = unit(0.5, "cm")))
p <- p + geom_rect(fill = "#3363CC", aes(xmin = -20, ymin = 72.5, xmax = -10, ymax = 77.5)) + geom_rect(fill = "#3363CC", aes(xmin = -10, xmax = 10, ymin = 70, ymax = 80)) + geom_rect(fill = "#3363CC", aes(xmin = 70, xmax = 90, ymin = 70, ymax = 80)) + geom_rect(fill = "#3363CC", aes(xmin = 100, ymin = 70, xmax = 120, ymax = 80)) + geom_rect(fill = "#3363CC", aes(xmin = 130, ymin = 70, xmax = 150, ymax = 80)) + geom_rect(fill = "#3363CC", aes(xmin = 150, ymin = 72.5, xmax = 160, ymax = 77.5)) + geom_segment(aes(x = 10, xend = 70, y = 75, yend = 75), arrow = arrow(length = unit(0.5, "cm"))) + geom_segment(aes(x = 90, y = 75, xend = 100, yend = 75), arrow = arrow(length = unit(0.5, "cm"))) + geom_segment(aes(x = 120, y = 75, xend = 130, yend = 75), arrow = arrow(length = unit(0.5, "cm")))
p <- p + geom_rect(fill = "#3363CC", aes(xmin = 0, ymin = 52.5, xmax = 10, ymax = 57.5)) + geom_rect(fill = "#3363CC", aes(xmin = 10, ymin = 50, xmax = 30, ymax = 60)) + geom_rect(fill = "#3363CC", aes(xmin = 100, ymin = 50, xmax = 120, ymax = 60)) + geom_rect(fill = "#3363CC", aes(xmin = 130, ymin = 50, xmax = 150, ymax = 60)) + geom_rect(fill = "#3363CC", aes(xmin = 150, ymin = 52.5, xmax = 160, ymax = 57.5)) + geom_segment(aes(x = 30, y = 55, xend = 100, yend = 55), arrow = arrow(length = unit(0.5, "cm"))) + geom_segment(aes(x = 120, y = 55, xend = 130, yend = 55), arrow = arrow(length = unit(0.5, "cm")))

#Plot junctions
p <- p + geom_rect(aes(xmin = 20, ymin = 30, xmax = 30, ymax = 35), fill = "red") + geom_rect(aes(xmin = 40, ymin = 30, xmax = 50, ymax = 35), fill = "red") + geom_segment(aes(x = 30, y = 32.5, xend = 40, yend = 32.5), arrow = arrow(length = unit(0.5, "cm")))
p <- p + geom_rect(aes(xmin = 20, ymin = 20, xmax = 30, ymax = 25), fill = "red") + geom_rect(aes(xmin = 70, ymin = 20, xmax = 80, ymax = 25), fill = "red") + geom_segment(aes(x = 30, y = 22.5, xend = 70, yend = 22.5), arrow = arrow(length = unit(0.5, "cm")))
p <- p + geom_rect(aes(xmin = 20, ymin = 10, xmax = 30, ymax = 15), fill = "red") + geom_rect(aes(xmin = 45, ymin = 10, xmax = 55, ymax = 15), fill = "red") + geom_segment(aes(x = 30, y = 12.5, xend = 45, yend = 12.5), arrow = arrow(length = unit(0.5, "cm")))
p <- p + geom_rect(aes(xmin = 10, ymin = 0, xmax = 20, ymax = 5), fill = "red") + geom_rect(aes(xmin = 70, ymin = 0, xmax = 80, ymax = 5), fill = "red") + geom_segment(aes(x = 20, y = 2.5, xend = 70, yend = 2.5), arrow = arrow(length = unit(0.5, "cm")))
p <- p + geom_rect(aes(xmin = -5, ymin = -10, xmax = 5, ymax = -5), fill = "red") + geom_rect(aes(xmin = 45, ymin = -10, xmax = 55, ymax = -5), fill = "red") + geom_segment(aes(x = 5, y = -7.5, xend = 45, yend = -7.5), arrow = arrow(length = unit(0.5, "cm")))

#Add labels
p <- p + annotate("text", label = "Transcripts", x = 7, y = 105, size = 5)
p <- p + annotate("text", label = "Junctions", x = 5, y = 42.5, size = 5)
p <- p + annotate("text", label = "Annotations from regtools", x = 120, y = 42.5, size = 5)
p <- p + annotate("text", label = "'DA': Known donor-acceptor", x = 140, y = 32.5, size = 5)
p <- p + annotate("text", label = "'NDA': Novel donor-acceptor combination", x = 140, y = 22.5, size = 5)
p <- p + annotate("text", label = "'D': Novel donor", x = 90, y = 12.5, size = 5)
p <- p + annotate("text", label = "'A': Novel acceptor", x = 90, y = 2.5, size = 5)
p <- p + annotate("text", label = "'N': Novel acceptor and novel donor.", x = 90, y = -7.5, size = 5)


#Add lines to make it clear
p <- p + geom_vline(xintercept = 30, linetype = "longdash") + geom_vline(xintercept = 40, linetype = "longdash") + geom_vline(xintercept = 70, linetype = "longdash")

print(p)
ggsave("anchor_examples.png")
ggsave("anchor_examples.pdf")

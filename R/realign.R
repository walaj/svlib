library(data.table)
library(ggplot2)

.load.results <- function(x, name, mode) {
  f <- fread(x)
  f[, min_size := pmin(V12, V13)]
  f[, round_size := 10* round(min_size / 10)]
  f[, ratio := sum(V9==1 & V10==1 & V15 == 60 & V11 == 2) / nrow(.SD), by=round_size]
  f[, class := name]
  f[, mode := mode]
  setkey(f,round_size)
  return(unique(f))
}

f <- rbindlist(list(.load.results("/xchip/gistic/Jeremiah/GIT/svlib/tenk.out", "bwa", "e0"),
                    .load.results("/xchip/gistic/Jeremiah/GIT/svlib/tenk.int.out", "bwa -intractg", "e0"),
                    .load.results("/xchip/gistic/Jeremiah/GIT/svlib/tenk.pac.out", "bwa -pacbio", "e0"),
                    .load.results("/xchip/gistic/Jeremiah/GIT/svlib/tenk.e50.out", "bwa", "e50"),
                    .load.results("/xchip/gistic/Jeremiah/GIT/svlib/tenk.e50.int.out", "bwa -intractg", "e50"),
                    .load.results("/xchip/gistic/Jeremiah/GIT/svlib/tenk.e50.pac.out", "bwa -pacbio", "e50")))

g <- ggplot(data=f[!(ratio == 0 & round_size == 500)], aes(x=round_size, y=ratio, color=class)) + geom_point() + geom_line() +
  theme_bw() + xlab("Minimum fragment size") + ylab("% correct multi-part mappings") +
  scale_y_continuous(breaks=seq(0,1,by=0.2), labels=seq(0,100,by=20)) +
  scale_color_manual(name="Aligner", values=c("darkred", "darkgreen", "blue")) + facet_wrap(~ mode) +
  
pdf("~/public_html/realign.pdf", width=8, height=4)
print(g)
dev.off()

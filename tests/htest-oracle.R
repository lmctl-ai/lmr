args <- commandArgs(TRUE)
library(stats)
out <- args[[1]]
save_test <- function(name, h) {
  scalar <- function(x) if(length(x)) sprintf("%.17g", as.numeric(x)[[1L]]) else ""
  vector <- function(x) if(length(x)) paste(sprintf("%.17g", as.numeric(x)), collapse=",") else ""
  fields <- c(scalar(h$statistic), vector(h$parameter), scalar(h$p.value), vector(h$conf.int), vector(h$estimate), h$method, if(is.null(h$alternative)) "" else h$alternative)
  writeLines(fields, file.path(out,paste0(name,".txt")))
}
save_test("t",t.test(c(1,2,3,4)))
save_test("t-paired",t.test(c(1,2,4,6),c(0,1,2,3),paired=TRUE))
save_test("wilcox",wilcox.test(1:8,3:10))
save_test("chisq",chisq.test(c(10,20,30),p=c(.2,.3,.5)))
save_test("chisq-table",chisq.test(matrix(c(10,20,30,40),2,2)))
save_test("shapiro",shapiro.test(c(1,2,3,5,7)))
save_test("ks",ks.test(c(-1,0,1,2),pnorm))
save_test("fisher",fisher.test(matrix(c(1,9,11,3),2,2)))
save_test("binom",binom.test(7,10,p=.5))
save_test("prop",prop.test(c(5,7),c(10,10),p=c(.5,.5)))
save_test("cor",cor.test(c(1,2,3,4,5),c(2,1,4,3,7),method="spearman"))
save_test("var",var.test(c(1,2,3,4,5),c(2,4,5,7,8)))

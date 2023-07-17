
# 提取分析 .out 文件
out <- read.delim("main_test2.out")
library(stringr)
info_qc <- str_extract_all(out, '(?<=q_c_seed )\\d')
info_aa <- str_extract_all(out, '(?<= a )[\\d\\.]*')
info_l1 <- str_extract_all(out, '(?<= l1 )[\\d\\.]*')
info_l2 <- str_extract_all(out, '(?<= l2 )[\\d\\.]*')
info_l3 <- str_extract_all(out, '(?<= l3 )[\\d\\.]*')
info_tau <- str_extract_all(out, '(?<= tau )[\\d\\.e-]*')
info_sc <- str_extract_all(out, '(?<=\\[1\\] )[\\d]*[\\.]+[\\d]*')



outres <- cbind(data.frame(q_c_seed = info_qc[[1]], 
                           aa = info_aa[[1]],
                           l1 = info_l1[[1]],
                           l2 = info_l2[[1]],
                           l3 = info_l3[[1]],
                           tau = info_tau[[1]]), t(matrix(info_sc[[1]], nrow = 2)))
colnames(outres) <- c("q_c_seed", 'aa', "l1", "l2", "l3", "tau", "cdist", "mse")
write.csv(outres, file = "out.csv", row.names = F)


plot(outres$mse, outres$cdist)
par(mfrow = c(2,5))
plot(outres$aa, outres$mse)
plot(outres$l1, outres$mse)
plot(outres$l2, outres$mse)
plot(outres$l3, outres$mse)
plot(outres$tau, outres$mse)

plot(outres$aa, outres$cdist)
plot(outres$l1, outres$cdist)
plot(outres$l2, outres$cdist)
plot(outres$l3, outres$cdist)
plot(outres$tau, outres$cdist)
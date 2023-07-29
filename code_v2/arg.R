library(argparse)

# 创建参数解析对象
parser <- ArgumentParser()

# 设置第一个参数verbose，缩写为v，其作用是告诉脚本是否打印完整的计算过程，其缺省值为TRUE
parser$add_argument("-n", "--num", default=200,
                    help="sample size")
parser$add_argument("-p", default = 8, help = "dim of X")
parser$add_argument("-q", default = 4, help = "dim of Z")

args <- parser$parse_args()
print(paste(args$num, args$p, args$q))
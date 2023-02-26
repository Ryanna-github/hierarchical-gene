# 评估预测准确度
right_class_ratio <- function(ci_true, ci_est){
  return(sum(ci_true == ci_est)/length(ci_true))
}

# 绘图观察样本归属类别随迭代轮次的变化
# input:
#   idx: 样本 id
#   coef_est: 维度*类别数
#   rho_used: 长度为类别数的数组（各类别 rho）
plot_classdis <- function(idx, coef_est, rho_used){
  x_sim <- seq(-20,20,length.out = 1000)
  df_plot <- data.frame(X%*%coef_est) %>% 
    cbind(data.frame(x = x_sim)) %>% 
    rename(c(y1="X1",y2="X2",y3="X3",y4="X4")) %>%
    mutate(den1 = dnorm(x, mean = y1[idx], sd = 1/rho_used[1]),
           den2 = dnorm(x, mean = y2[idx], sd = 1/rho_used[2]),
           den3 = dnorm(x, mean = y3[idx], sd = 1/rho_used[3]),
           den4 = dnorm(x, mean = y4[idx], sd = 1/rho_used[4])) %>%
    select(x, den1, den2, den3, den4)
  
  df_plot <- melt(df_plot, id.vars = c("x")) %>% 
    rename(Class = variable, y = value)
  
  p <- ggplot(df_plot, aes(x = x, ymin = 0, ymax = y, y = y, 
                           col = Class, fill = Class)) +
    geom_ribbon(alpha = 0.5) +
    geom_line() +
    geom_vline(xintercept = y[idx]) +
    theme_minimal() +
    labs(title = idx) +
    annotate("text", x = y[idx], y = 2*max(df_plot$y)/3, label = paste0("y_", idx, ":", round(y[idx],2)))
  return(p)
}
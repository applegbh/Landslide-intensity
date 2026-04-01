# ========================================
# 滑坡流动性 GAM 建模（统一 k 参数设置 - 完整修复版）
# ========================================

library(mgcv)
library(Metrics)
library(sf)
library(hydroGOF)  # 用于NSE计算
# ========================================
# 0. 统一参数设置
# ========================================
num.baseFun  <- 45     # 所有一维平滑项的 k =30, r2=0.56。40==》45#
num.spaceFun <- 180    # 空间二维平滑项的 k。 100==》   200

# ========================================
# 1. 读取数据
# ========================================
WORK <- "/Volumes/SU_Codes/Mobility/LuDmob/Imax_org"
setwd(WORK)

cat("读取 Shapefile 数据...\n")
su_data <- st_read("su_100k_02_clean_with_phy_mobility_cut1.shp")

data.raw <- st_drop_geometry(su_data)
cat("原始数据维度:", dim(data.raw), "\n")

# ========================================
# 2. 数据预处理（只保留流动性和面积 > 0 的单元）
# ========================================
cat("\n数据预处理...\n")
#MaxMOb_LH,I_max_org,I_max_sq,I1_orig,I_star_or
data4fit <- subset(data.raw, I_max_org > 0 & Area_sum > 0)

cat("筛选后样本数:", nrow(data4fit), "\n")

# 分类变量转换
categorical_vars <- c("Landforms", "Lithology", "LULC")
for (var in categorical_vars) {
  if (var %in% names(data4fit)) {
    data4fit[[var]] <- as.factor(data4fit[[var]])
    cat("转换", var, "为因子, 水平数:", nlevels(data4fit[[var]]), "\n")
  }
}

# 创建对数变量
data4fit$log_Mobility <- log1p(data4fit$I_max_org)
data4fit$log_Area     <- log(data4fit$Area_sum)

# 清除关键变量的 NA
data4fit <- data4fit[complete.cases(data4fit[, c("log_Mobility", "I_max_org", "Area_sum")]), ]
cat("最终样本数:", nrow(data4fit), "\n")
data4fit$PGA_total <- sqrt(data4fit$PGA_EW_m^2 + data4fit$PGA_NS_m^2)
# ========================================
# 3. 构建 GAM 公式（使用 bquote 注入 k 值）
# ========================================
cat("\n构建 GAM 公式...\n")

# 含空间项的公式
Formula.mobility.v2 <- bquote(
  log_Mobility ~
    s(East_m,    k = .(num.baseFun)) +
    s(Fault_m,   k = .(num.baseFun)) +
    s(North_m,   k = .(num.baseFun)) +
    s(PGA_total,  k = .(num.baseFun)) +
    s(PGA_UD_s,  k = .(num.baseFun)) +
    s(PLC_s,     k = .(num.baseFun)) +
    s(PRC_m,     k = .(num.baseFun)) +
    #s(PRC_s,     k = .(num.baseFun)) +
    s(River_m,   k = .(num.baseFun)) +
    s(River_s,   k = .(num.baseFun)) +
    s(Road_m,    k = .(num.baseFun)) +
    s(Slope_m, bs = "ad",   k = .(num.baseFun)) +
    s(TWI_s,     k = .(num.baseFun)) +
    #s(ShRelief_s,k = .(num.baseFun)) +
    #s(travel_ang,k = 5) + #
    #s(travel_dis,k = 5) + #
    s(delta_z,   k = .(num.baseFun)) +
    #s(V_Std,     k = .(num.baseFun)) +
    s(V_Mean,     k = .(num.baseFun)) +
    ti(PLC_m,PRC_m, k=15)+
    ti(Fault_m,PGA_UD_s,k=15)+
    ti(Fault_m,PGA_EW_m,k=15)+
    ti(PGA_UD_s,PGA_EW_m,k=10)+
    ti(delta_z,travel_dis,k=10)+
    #ti(PGA_total,V_Std,k=10)+
    ti(PGA_UD_s,V_Std,k=15)+  
    ti(PGA_EW_m,V_Std,k=15)+
    
    ti(delta_z,V_Std,k=8)+ # #k=10
    ti(PGA_total,Slope_m,k=10)+
    s(CoordX, CoordY, k = .(num.spaceFun), bs = "ds") +
    as.factor(Landforms) +
    as.factor(Lithology) +
    as.factor(LULC)
)
cat("使用包含空间项的公式\n")
use_formula <- as.formula(Formula.mobility.v2)

# 显示公式
cat("\n模型公式:\n")
print(use_formula)

# ========================================
# 4. 检查变量是否存在
# ========================================
cat("\n检查变量...\n")

# 提取公式中的变量名
formula_vars <- all.vars(terms(use_formula))
formula_vars <- setdiff(formula_vars, "log_Mobility")  # 移除响应变量

# 检查缺失变量
missing_vars <- setdiff(formula_vars, names(data4fit))

if (length(missing_vars) > 0) {
  cat("警告: 以下变量不存在于数据中:\n")
  print(missing_vars)
  
  # 尝试查找相似变量
  cat("\n尝试查找相似变量:\n")
  for (var in missing_vars) {
    similar <- names(data4fit)[grep(var, names(data4fit), ignore.case = TRUE)]
    if (length(similar) > 0) {
      cat("  ", var, "可能是:", paste(similar, collapse = ", "), "\n")
    }
  }
  
  cat("\n可用的数值型变量:\n")
  numeric_vars <- names(data4fit)[sapply(data4fit, is.numeric)]
  print(numeric_vars)
  
  stop("请修正变量名后重新运行")
}

cat("✓ 所有变量都存在\n")

# ========================================
# 5. 拟合 BAM 模型（使用 bam() 替代 gam()）
# ========================================
cat("\n开始拟合 BAM 模型...\n")

Fit.mobility <- tryCatch({
  mgcv::bam(use_formula, family = gaussian(), data = data4fit, method = "fREML", select = TRUE, nthreads = 4)#, discrete = TRUE
}, error = function(e) {
  cat("模型拟合失败:", e$message, "\n")
  return(NULL)
})

if (is.null(Fit.mobility)) {
  stop("模型拟合失败")
}

cat("✓ 模型拟合成功\n")

# ========================================
# 6. 模型摘要
# ========================================
cat("\n", paste(rep("=", 70), collapse = ""), "\n", sep = "")
cat("模型摘要\n")
cat(paste(rep("=", 70), collapse = ""), "\n", sep = "")
print(summary(Fit.mobility))

# ========================================
# 7. 性能评估
# ========================================
cat("\n", paste(rep("=", 70), collapse = ""), "\n", sep = "")
cat("模型性能指标\n")
cat(paste(rep("=", 70), collapse = ""), "\n", sep = "")

fitted_values   <- Fit.mobility$fitted.values
observed_values <- data4fit$log_Mobility
SSres <- sum((observed_values - fitted_values)^2)
SStot <- sum((observed_values - mean(fitted_values))^2)

r2_value <- 1 - (SSres / SStot)
mae_value      <- mae(fitted_values, observed_values)
rmse_value     <- rmse(fitted_values, observed_values)
#r2_value       <- cor(fitted_values, observed_values, method = "pearson")^2
nse            <- NSE(fitted_values, observed_values)
Pearson_value <- cor(fitted_values, observed_values, method = "pearson")

cat(sprintf("MAE:              %.4f\n", mae_value))
cat(sprintf("RMSE:             %.4f\n", rmse_value))
cat(sprintf("R² (Pearson):     %.4f\n", r2_value))
cat(sprintf("NSE:             %.4f\n", nse))
cat(sprintf("Pearson相关:     %.4f\n", Pearson_value))
cat(sprintf("解释偏差:         %.2f%%\n", summary(Fit.mobility)$dev.expl * 100))

# ========================================
# 8. 可视化
# ========================================
cat("\n生成可视化图表...\n")

# 8.1 观测值 vs 拟合值
png("mobilityImax_org_fitted_vs_observed.png",
    width = 800, height = 800, res = 120)

smoothScatter(observed_values, fitted_values,
              xlab = expression("Observed log(H/L)"),
              ylab = expression("Fitted log(H/L)"),
              main = sprintf("Mobility Model\nR² = %.4f, RMSE = %.4f",
                             r2_value, rmse_value),
              xlim = range(observed_values, na.rm = TRUE),
              ylim = range(fitted_values, na.rm = TRUE),
              family = "serif",
              cex.axis = 1.2,
              cex.lab  = 1.4,
              cex.main = 1.3)

abline(a = 0, b = 1, lwd = 2, lty = 2, col = "red")
abline(lm(fitted_values ~ observed_values), lwd = 2, col = "blue")

legend("topleft",
       legend = c("1:1 Line", "Fitted Line"),
       col    = c("red", "blue"),
       lty    = c(2, 1),
       lwd    = 2,
       cex    = 1.1)

dev.off()
cat("✓ 散点图已保存: mobilityI_max_org_fitted_vs_observed.png\n")

# 8.2 残差诊断
png("mobilityImax_org_residuals.png",
    width = 1200, height = 900, res = 120)

par(mfrow = c(2, 2))
gam.check(Fit.mobility)

dev.off()
cat("✓ 残差诊断图已保存: mobilityI_max_org_residuals.png\n")

# 8.3 偏效应图
png("mobility_partial_effects.png",
    width = 1600, height = 1200, res = 120)

plot(Fit.mobility,
     pages    = 1,
     scheme   = 2,
     shade    = TRUE,
     shade.col= "lightblue",
     cex.main = 1.2,
     cex.lab  = 1.1)

dev.off()
cat("✓ 偏效应图已保存: mobilityImax_org_partial_effects.png\n")

# ========================================
# 9. 保存结果
# ========================================
cat("\n保存结果...\n")

# 保存模型
saveRDS(Fit.mobility, "mobilityImax_org_gam_model.rds")
cat("✓ 模型已保存: mobilityI_max_org_gam_model.rds\n")

# 保存预测结果
predictions_df <- data.frame(
  Observed        = observed_values,
  Fitted          = fitted_values,
  Residual        = observed_values - fitted_values,
  Actual_Mobility = data4fit$H_L_sum
)

write.csv(predictions_df, "mobilityImax_org_predictions.csv", row.names = FALSE)
cat("✓ 预测结果已保存: mobilityI_max_org_predictions.csv\n")

# 保存性能指标
metrics_df <- data.frame(
  Metric = c("MAE", "RMSE", "R2", "NSE", "Pearson", "Deviance_Explained", "N"),
  Value  = c(mae_value, rmse_value, r2_value,nse, Pearson_value,
             summary(Fit.mobility)$dev.expl, nrow(data4fit))
)

write.csv(metrics_df, "mobilityI_max_org_metrics.csv", row.names = FALSE)
cat("✓ 性能指标已保存: mobilityI_max_org_metrics.csv\n")

# 保存模型摘要
sink("mobilityI_max_org_model_summary.txt")
cat(paste(rep("=", 70), collapse = ""), "\n", sep = "")
cat("滑坡流动性 GAM 模型摘要\n")
cat(paste(rep("=", 70), collapse = ""), "\n\n", sep = "")
cat("样本数:", nrow(data4fit), "\n")
cat("筛选条件: H_L_sum > 0 且 Area_sum > 0\n")
cat("响应变量: log(H_L_sum)\n")
cat("基函数数 k:", num.baseFun, "\n")
if("CoordX" %in% names(data4fit)) {
  cat("空间项 k:", num.spaceFun, "\n")
}
cat("\n性能指标:\n")
print(metrics_df)
cat("\n\n模型详细摘要:\n")
print(summary(Fit.mobility))
sink()
cat("✓ 模型摘要已保存: mobility_Imax_org_model_summary.txt\n")

# ========================================
# 10. 最终报告
# ========================================
cat("\n", paste(rep("=", 70), collapse = ""), "\n", sep = "")
cat("分析完成!\n")
cat(paste(rep("=", 70), collapse = ""), "\n", sep = "")
cat("\n输出文件:\n")
cat("  1. mobility_gam_model.rds\n")
cat("  2. mobility_predictions.csv\n")
cat("  3. mobility_metrics.csv\n")
cat("  4. mobility_model_summary.txt\n")
cat("  5. mobility_fitted_vs_observed.png\n")
cat("  6. mobility_residuals.png\n")
cat("  7. mobility_partial_effects.png\n")
cat("\n工作目录:", getwd(), "\n")
cat(paste(rep("=", 70), collapse = ""), "\n", sep = "")

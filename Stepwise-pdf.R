
# ========================================
# 滑坡流动性 GAM 建模 - Stepwise AIC 变量选择（含交互项）
# ========================================

library(mgcv)
library(Metrics)
library(sf)
library(hydroGOF)

# ========================================
# 0. 统一参数设置
# ========================================
num.baseFun  <- 45
num.spaceFun <- 180

# AIC 改善阈值
AIC_threshold_main <- 2    # 主效应选择阈值
AIC_threshold_ti   <- 2    # 交互项选择阈值

# ========================================
# 1. 读取数据
# ========================================
WORK <- "/Volumes/SU_Codes/Mobility/LuDmob/Imax_org_pdf"
setwd(WORK)

cat("读取 Shapefile 数据...\n")
su_data <- st_read("su_100k_02_clean_with_phy_mobility_cut1.shp")
data.raw <- st_drop_geometry(su_data)
cat("原始数据维度:", dim(data.raw), "\n")

# ========================================
# 2. 数据预处理
# ========================================
cat("\n数据预处理...\n")

data4fit <- subset(data.raw, I_max_org > 0 & Area_sum > 0)
cat("筛选后样本数:", nrow(data4fit), "\n")

# 分类变量转换
categorical_vars <- c("Landforms", "Lithology", "LULC")
for (var in categorical_vars) {
  if (var %in% names(data4fit)) {
    data4fit[[var]] <- as.factor(data4fit[[var]])
  }
}

# 创建对数变量
data4fit$log_Mobility <- log1p(data4fit$I_max_org)
data4fit$log_Area     <- log1p(data4fit$Area_sum)

# 派生变量
data4fit$PGA_total <- sqrt(data4fit$PGA_EW_m^2 + data4fit$PGA_NS_m^2)

# 移除 NA
data4fit <- na.omit(data4fit)
cat("最终样本数:", nrow(data4fit), "\n")

# ========================================
# 3. 定义候选变量
# ========================================

# 🔑 候选的连续型变量（一维平滑项）
candidate_smooth_vars <- c(
  "East_m", "Fault_m", "North_m", "PGA_total", "PGA_UD_s",
  "PLC_s", "PRC_m", "PRC_s", "River_m", "River_s",
  "Road_m", "Slope_m", "TWI_s", "ShRelief_s", "delta_z",
  "V_Std", "V_Mean", "PLC_m", "PGA_EW_m", "travel_dis"
)

# 🔑 候选的交互项（二维张量积）
candidate_ti_terms <- list(
  list(vars = c("PLC_m", "PRC_m"),     k = 15),
  list(vars = c("Fault_m", "PGA_UD_s"), k = 15),
  list(vars = c("Fault_m", "PGA_EW_m"), k = 15),
  list(vars = c("PGA_UD_s", "PGA_EW_m"), k = 10),
  list(vars = c("delta_z", "travel_dis"), k = 10),
  list(vars = c("PGA_total", "V_Std"),  k = 10),
  list(vars = c("PGA_UD_s", "V_Std"),   k = 15),
  list(vars = c("PGA_EW_m", "V_Std"),   k = 15),
  list(vars = c("delta_z", "V_Std"),    k = 8),
  list(vars = c("PGA_total", "Slope_m"), k = 10),
  list(vars = c("Slope_m", "ShRelief_s"), k = 10),
  list(vars = c("TWI_s", "Slope_m"),    k = 10),
  list(vars = c("River_m", "Slope_m"),  k = 10),
  list(vars = c("Fault_m", "Slope_m"),  k = 10)
)

# 🔑 候选的分类变量
candidate_factor_vars <- c("Landforms", "Lithology", "LULC")

# 🔑 固定项（空间项）
spatial_term <- sprintf("s(CoordX, CoordY, k = %d, bs = 'ds')", num.spaceFun)

# 检查变量是否存在
available_smooth <- candidate_smooth_vars[candidate_smooth_vars %in% names(data4fit)]
available_factor <- candidate_factor_vars[candidate_factor_vars %in% names(data4fit)]

# 检查交互项中的变量是否存在
available_ti_terms <- list()
for (ti_item in candidate_ti_terms) {
  if (all(ti_item$vars %in% names(data4fit))) {
    available_ti_terms <- c(available_ti_terms, list(ti_item))
  }
}

cat("\n可用的连续变量:", length(available_smooth), "\n")
print(available_smooth)
cat("\n可用的分类变量:", length(available_factor), "\n")
print(available_factor)
cat("\n可用的交互项:", length(available_ti_terms), "\n")
for (i in seq_along(available_ti_terms)) {
  cat(sprintf("  %d. ti(%s, %s, k=%d)\n", i, 
              available_ti_terms[[i]]$vars[1], 
              available_ti_terms[[i]]$vars[2],
              available_ti_terms[[i]]$k))
}

# ========================================
# 4. 辅助函数
# ========================================

# 构建公式（含交互项）
build_formula_full <- function(selected_smooth, selected_factor, selected_ti, 
                               k_smooth, spatial_term) {
  terms <- c()
  
  # 添加平滑项
  for (var in selected_smooth) {
    if (var == "Slope_m") {
      terms <- c(terms, sprintf("s(%s, bs = 'ad', k = %d)", var, k_smooth))
    } else {
      terms <- c(terms, sprintf("s(%s, k = %d)", var, k_smooth))
    }
  }
  
  # 添加交互项
  for (ti_item in selected_ti) {
    terms <- c(terms, sprintf("ti(%s, %s, k = %d)", 
                              ti_item$vars[1], ti_item$vars[2], ti_item$k))
  }
  
  # 添加分类变量
  for (var in selected_factor) {
    terms <- c(terms, sprintf("as.factor(%s)", var))
  }
  
  # 添加空间项（固定）
  terms <- c(terms, spatial_term)
  
  # 构建公式
  formula_str <- paste("log_Mobility ~", paste(terms, collapse = " + "))
  return(as.formula(formula_str))
}

# 安全拟合模型
safe_fit <- function(formula, data) {
  tryCatch({
    model <- mgcv::bam(
      formula,
      family = gaussian(),
      data = data,
      method = "fREML",
      select = TRUE,
      nthreads = 4
    )
    return(model)
  }, error = function(e) {
    return(NULL)
  })
}

# 交互项转字符串（用于显示和记录）
ti_to_string <- function(ti_item) {
  sprintf("ti(%s,%s)", ti_item$vars[1], ti_item$vars[2])
}

# ========================================
# 5. 第一阶段：主效应变量选择
# ========================================

cat("\n", paste(rep("=", 80), collapse = ""), "\n", sep = "")
cat("第一阶段：主效应变量选择（基于 AIC）\n")
cat(paste(rep("=", 80), collapse = ""), "\n", sep = "")

# 初始化
selected_smooth <- c()
selected_factor <- c()
selected_ti     <- list()  # 交互项暂时为空
remaining_smooth <- available_smooth
remaining_factor <- available_factor

# 存储选择历史
selection_history <- data.frame(
  Step = integer(),
  Stage = character(),
  Variable = character(),
  Type = character(),
  AIC = numeric(),
  Delta_AIC = numeric(),
  R2 = numeric(),
  Dev_Expl = numeric(),
  stringsAsFactors = FALSE
)

# Step 0: 仅空间项的基准模型
cat("\n--- Step 0: 基准模型（仅空间项） ---\n")
base_formula <- as.formula(paste("log_Mobility ~", spatial_term))
base_model <- safe_fit(base_formula, data4fit)

if (is.null(base_model)) {
  stop("基准模型拟合失败")
}

current_AIC <- AIC(base_model)
cat(sprintf("基准 AIC: %.2f\n", current_AIC))

selection_history <- rbind(selection_history, data.frame(
  Step = 0,
  Stage = "Baseline",
  Variable = "Spatial_Only",
  Type = "baseline",
  AIC = current_AIC,
  Delta_AIC = NA,
  R2 = summary(base_model)$r.sq,
  Dev_Expl = summary(base_model)$dev.expl * 100
))

current_model <- base_model
step <- 0
continue_selection <- TRUE

# 前向逐步选择主效应
while (continue_selection && (length(remaining_smooth) > 0 || length(remaining_factor) > 0)) {
  step <- step + 1
  cat(sprintf("\n--- Step %d (主效应) ---\n", step))
  cat(sprintf("当前: %d 连续 + %d 分类\n", 
              length(selected_smooth), length(selected_factor)))
  cat(sprintf("剩余: %d 连续 + %d 分类\n", 
              length(remaining_smooth), length(remaining_factor)))
  
  best_AIC <- Inf
  best_var <- NULL
  best_type <- NULL
  best_model <- NULL
  
  # 尝试添加每个剩余的连续变量
  for (var in remaining_smooth) {
    test_smooth <- c(selected_smooth, var)
    test_formula <- build_formula_full(test_smooth, selected_factor, selected_ti,
                                       num.baseFun, spatial_term)
    
    test_model <- safe_fit(test_formula, data4fit)
    
    if (!is.null(test_model)) {
      test_AIC <- AIC(test_model)
      if (test_AIC < best_AIC) {
        best_AIC <- test_AIC
        best_var <- var
        best_type <- "smooth"
        best_model <- test_model
      }
    }
  }
  
  # 尝试添加每个剩余的分类变量
  for (var in remaining_factor) {
    test_factor <- c(selected_factor, var)
    test_formula <- build_formula_full(selected_smooth, test_factor, selected_ti,
                                       num.baseFun, spatial_term)
    
    test_model <- safe_fit(test_formula, data4fit)
    
    if (!is.null(test_model)) {
      test_AIC <- AIC(test_model)
      if (test_AIC < best_AIC) {
        best_AIC <- test_AIC
        best_var <- var
        best_type <- "factor"
        best_model <- test_model
      }
    }
  }
  
  # 检查是否有改善
  if (!is.null(best_var)) {
    delta_AIC <- current_AIC - best_AIC
    
    cat(sprintf("最佳候选: %s (%s)\n", best_var, best_type))
    cat(sprintf("AIC: %.2f -> %.2f (Δ = %.2f)\n", current_AIC, best_AIC, delta_AIC))
    
    if (delta_AIC > AIC_threshold_main) {
      # 接受这个变量
      if (best_type == "smooth") {
        selected_smooth <- c(selected_smooth, best_var)
        remaining_smooth <- setdiff(remaining_smooth, best_var)
      } else {
        selected_factor <- c(selected_factor, best_var)
        remaining_factor <- setdiff(remaining_factor, best_var)
      }
      
      current_AIC <- best_AIC
      current_model <- best_model
      
      selection_history <- rbind(selection_history, data.frame(
        Step = step,
        Stage = "Main",
        Variable = best_var,
        Type = best_type,
        AIC = best_AIC,
        Delta_AIC = delta_AIC,
        R2 = summary(best_model)$r.sq,
        Dev_Expl = summary(best_model)$dev.expl * 100
      ))
      
      cat(sprintf("✓ 接受 %s (R² = %.4f, Dev.Expl = %.2f%%)\n", 
                  best_var, summary(best_model)$r.sq, summary(best_model)$dev.expl * 100))
      
    } else {
      cat(sprintf("✗ AIC 改善不显著 (Δ = %.2f < %.2f)\n", 
                  delta_AIC, AIC_threshold_main))
      continue_selection <- FALSE
    }
  } else {
    cat("没有找到可以改善模型的变量\n")
    continue_selection <- FALSE
  }
}

cat("\n", paste(rep("-", 60), collapse = ""), "\n", sep = "")
cat("第一阶段完成\n")
cat(sprintf("选中连续变量 (%d): %s\n", length(selected_smooth), paste(selected_smooth, collapse = ", ")))
cat(sprintf("选中分类变量 (%d): %s\n", length(selected_factor), paste(selected_factor, collapse = ", ")))
cat(sprintf("当前 AIC: %.2f\n", current_AIC))

# ========================================
# 6. 第二阶段：交互项选择
# ========================================

cat("\n", paste(rep("=", 80), collapse = ""), "\n", sep = "")
cat("第二阶段：交互项选择（基于 AIC）\n")
cat(paste(rep("=", 80), collapse = ""), "\n", sep = "")

# 筛选候选交互项：至少一个变量已被选中（可选：要求两个都被选中）
# 这里采用"至少一个变量在 selected_smooth 中"的策略
remaining_ti <- list()
for (ti_item in available_ti_terms) {
  # 检查交互项的变量是否与已选变量相关
  # 策略1：至少一个变量已选中
  # 策略2：两个变量都已选中（更保守）
  # 这里用策略1
  if (any(ti_item$vars %in% selected_smooth) || 
      any(ti_item$vars %in% available_smooth)) {
    remaining_ti <- c(remaining_ti, list(ti_item))
  }
}

cat(sprintf("\n候选交互项: %d\n", length(remaining_ti)))
for (i in seq_along(remaining_ti)) {
  cat(sprintf("  %d. %s\n", i, ti_to_string(remaining_ti[[i]])))
}

# 记录主效应完成时的步数
step_main_end <- step

continue_ti_selection <- TRUE

while (continue_ti_selection && length(remaining_ti) > 0) {
  step <- step + 1
  cat(sprintf("\n--- Step %d (交互项) ---\n", step))
  cat(sprintf("当前交互项: %d, 剩余候选: %d\n", 
              length(selected_ti), length(remaining_ti)))
  
  best_AIC <- Inf
  best_ti <- NULL
  best_ti_idx <- NULL
  best_model <- NULL
  
  # 尝试添加每个剩余的交互项
  for (i in seq_along(remaining_ti)) {
    ti_item <- remaining_ti[[i]]
    test_ti <- c(selected_ti, list(ti_item))
    
    test_formula <- build_formula_full(selected_smooth, selected_factor, test_ti,
                                       num.baseFun, spatial_term)
    
    test_model <- safe_fit(test_formula, data4fit)
    
    if (!is.null(test_model)) {
      test_AIC <- AIC(test_model)
      cat(sprintf("  测试 %s: AIC = %.2f\n", ti_to_string(ti_item), test_AIC))
      
      if (test_AIC < best_AIC) {
        best_AIC <- test_AIC
        best_ti <- ti_item
        best_ti_idx <- i
        best_model <- test_model
      }
    } else {
      cat(sprintf("  测试 %s: 拟合失败\n", ti_to_string(ti_item)))
    }
  }
  
  # 检查是否有改善
  if (!is.null(best_ti)) {
    delta_AIC <- current_AIC - best_AIC
    
    cat(sprintf("\n最佳候选: %s\n", ti_to_string(best_ti)))
    cat(sprintf("AIC: %.2f -> %.2f (Δ = %.2f)\n", current_AIC, best_AIC, delta_AIC))
    
    if (delta_AIC > AIC_threshold_ti) {
      # 接受这个交互项
      selected_ti <- c(selected_ti, list(best_ti))
      remaining_ti <- remaining_ti[-best_ti_idx]
      
      current_AIC <- best_AIC
      current_model <- best_model
      
      selection_history <- rbind(selection_history, data.frame(
        Step = step,
        Stage = "Interaction",
        Variable = ti_to_string(best_ti),
        Type = "ti",
        AIC = best_AIC,
        Delta_AIC = delta_AIC,
        R2 = summary(best_model)$r.sq,
        Dev_Expl = summary(best_model)$dev.expl * 100
      ))
      
      cat(sprintf("✓ 接受 %s (R² = %.4f, Dev.Expl = %.2f%%)\n", 
                  ti_to_string(best_ti), 
                  summary(best_model)$r.sq, 
                  summary(best_model)$dev.expl * 100))
      
    } else {
      cat(sprintf("✗ AIC 改善不显著 (Δ = %.2f < %.2f)\n", 
                  delta_AIC, AIC_threshold_ti))
      continue_ti_selection <- FALSE
    }
  } else {
    cat("没有找到可以改善模型的交互项\n")
    continue_ti_selection <- FALSE
  }
}

# ========================================
# 7. 输出选择结果
# ========================================

cat("\n", paste(rep("=", 80), collapse = ""), "\n", sep = "")
cat("变量选择完成\n")
cat(paste(rep("=", 80), collapse = ""), "\n", sep = "")

cat("\n📊 完整选择历史:\n")
print(selection_history)

cat("\n✓ 最终选中的连续变量 (", length(selected_smooth), "):\n", sep = "")
print(selected_smooth)

cat("\n✓ 最终选中的分类变量 (", length(selected_factor), "):\n", sep = "")
print(selected_factor)

cat("\n✓ 最终选中的交互项 (", length(selected_ti), "):\n", sep = "")
for (ti_item in selected_ti) {
  cat(sprintf("  - %s (k=%d)\n", ti_to_string(ti_item), ti_item$k))
}

cat(sprintf("\n最终模型 AIC: %.2f\n", current_AIC))
cat(sprintf("最终模型 R²: %.4f\n", summary(current_model)$r.sq))
cat(sprintf("最终模型 Dev.Expl: %.2f%%\n", summary(current_model)$dev.expl * 100))

# 保存选择历史
write.csv(selection_history, "stepwise_selection_history.csv", row.names = FALSE)
cat("\n✓ 选择历史已保存: stepwise_selection_history.csv\n")

# ========================================
# 8. 可视化 AIC 变化
# ========================================

png("stepwise_AIC_progression.png", width = 1200, height = 700, res = 120)

par(mar = c(10, 5, 4, 2))

# 颜色区分阶段
colors <- ifelse(selection_history$Stage == "Baseline", "gray50",
                 ifelse(selection_history$Stage == "Main", "blue", "red"))

plot(1:nrow(selection_history), selection_history$AIC,
     type = "b", pch = 19, col = colors, lwd = 2,
     xlab = "", ylab = "AIC",
     main = "Stepwise GAM Variable Selection (Main Effects + Interactions)",
     xaxt = "n", cex.axis = 1.2, cex.lab = 1.3)

# 添加变量名作为 x 轴标签
axis(1, at = 1:nrow(selection_history), 
     labels = selection_history$Variable, 
     las = 2, cex.axis = 0.8)

mtext("Variables Added", side = 1, line = 8, cex = 1.2)

# 添加阶段分隔线
if (step_main_end > 0 && step_main_end < nrow(selection_history) - 1) {
  abline(v = step_main_end + 1.5, lty = 3, col = "gray50", lwd = 2)
  text(step_main_end + 1.5, max(selection_history$AIC) * 0.95, 
       "← Main | Interaction →", pos = 4, cex = 0.9)
}

# 添加 Delta AIC 标注
for (i in 2:nrow(selection_history)) {
  delta <- selection_history$Delta_AIC[i]
  if (!is.na(delta) && delta > 0) {
    text(i, selection_history$AIC[i],
         sprintf("Δ=%.1f", delta), pos = 3, cex = 0.7, col = "darkgreen")
  }
}

legend("topright", 
       legend = c("Baseline", "Main Effects", "Interactions"),
       col = c("gray50", "blue", "red"), 
       pch = 19, lwd = 2, cex = 1.0)

dev.off()
cat("✓ AIC 变化图已保存: stepwise_AIC_progression.png\n")

# ========================================
# 9. 拟合最终模型
# ========================================

cat("\n", paste(rep("=", 80), collapse = ""), "\n", sep = "")
cat("拟合最终模型\n")
cat(paste(rep("=", 80), collapse = ""), "\n", sep = "")

# 构建最终公式
final_formula <- build_formula_full(selected_smooth, selected_factor, selected_ti,
                                    num.baseFun, spatial_term)

cat("\n最终模型公式:\n")
print(final_formula)

# 拟合最终模型
Fit.mobility <- mgcv::bam(
  final_formula,
  family = gaussian(),
  data = data4fit,
  method = "fREML",
  select = TRUE,
  nthreads = 4
)

cat("\n✓ 最终模型拟合成功\n")

# ========================================
# 10. 模型评估
# ========================================

cat("\n", paste(rep("=", 70), collapse = ""), "\n", sep = "")
cat("最终模型性能指标\n")
cat(paste(rep("=", 70), collapse = ""), "\n", sep = "")

print(summary(Fit.mobility))

fitted_values   <- Fit.mobility$fitted.values
observed_values <- data4fit$log_Mobility
SSres <- sum((observed_values - fitted_values)^2)
SStot <- sum((observed_values - mean(observed_values))^2)

r2_value <- 1 - (SSres / SStot)
mae_value <- mae(fitted_values, observed_values)
rmse_value <- rmse(fitted_values, observed_values)
nse <- NSE(fitted_values, observed_values)
Pearson_value <- cor(fitted_values, observed_values, method = "pearson")

cat(sprintf("\n性能指标汇总:\n"))
cat(sprintf("  MAE:              %.4f\n", mae_value))
cat(sprintf("  RMSE:             %.4f\n", rmse_value))
cat(sprintf("  R²:               %.4f\n", r2_value))
cat(sprintf("  NSE:              %.4f\n", nse))
cat(sprintf("  Pearson:          %.4f\n", Pearson_value))
cat(sprintf("  Dev.Expl:         %.2f%%\n", summary(Fit.mobility)$dev.expl * 100))
cat(sprintf("  AIC:              %.2f\n", AIC(Fit.mobility)))

# ========================================
# 11. 保存结果
# ========================================

cat("\n保存结果...\n")

# 保存模型
saveRDS(Fit.mobility, "stepwise_final_model.rds")
cat("✓ 最终模型已保存: stepwise_final_model.rds\n")

# 保存预测结果
predictions_df <- data.frame(
  Observed = observed_values,
  Fitted   = fitted_values,
  Residual = observed_values - fitted_values
)
write.csv(predictions_df, "stepwise_predictions.csv", row.names = FALSE)
cat("✓ 预测结果已保存: stepwise_predictions.csv\n")

# 保存性能指标
metrics_df <- data.frame(
  Metric = c("MAE", "RMSE", "R2", "NSE", "Pearson", "Deviance_Explained", "AIC", "N",
             "N_smooth", "N_factor", "N_interaction"),
  Value  = c(mae_value, rmse_value, r2_value, nse, Pearson_value,
             summary(Fit.mobility)$dev.expl * 100, AIC(Fit.mobility), nrow(data4fit),
             length(selected_smooth), length(selected_factor), length(selected_ti))
)
write.csv(metrics_df, "stepwise_metrics.csv", row.names = FALSE)
cat("✓ 性能指标已保存: stepwise_metrics.csv\n")

# 保存模型摘要
sink("stepwise_final_model_summary.txt")
cat(paste(rep("=", 80), collapse = ""), "\n")
cat("Stepwise GAM 变量选择 - 最终模型摘要\n")
cat(paste(rep("=", 80), collapse = ""), "\n\n")

cat("选择参数:\n")
cat(sprintf("  主效应 AIC 阈值: %d\n", AIC_threshold_main))
cat(sprintf("  交互项 AIC 阈值: %d\n", AIC_threshold_ti))
cat(sprintf("  基函数 k: %d\n", num.baseFun))
cat(sprintf("  空间项 k: %d\n", num.spaceFun))

cat("\n选中的连续变量:\n")
print(selected_smooth)

cat("\n选中的分类变量:\n")
print(selected_factor)

cat("\n选中的交互项:\n")
for (ti_item in selected_ti) {
  cat(sprintf("  - ti(%s, %s, k=%d)\n", ti_item$vars[1], ti_item$vars[2], ti_item$k))
}

cat("\n最终模型公式:\n")
print(final_formula)

cat("\n\n选择历史:\n")
print(selection_history)

cat("\n\n性能指标:\n")
print(metrics_df)

cat("\n\n完整模型摘要:\n")
print(summary(Fit.mobility))
sink()
cat("✓ 模型摘要已保存: stepwise_final_model_summary.txt\n")

# ========================================
# 12. 可视化（续）
# ========================================

# 12.1 观测值 vs 拟合值
png("stepwise_fitted_vs_observed.png", width = 800, height = 800, res = 120)

smoothScatter(observed_values, fitted_values,
              xlab = "Observed log(Mobility)",
              ylab = "Fitted log(Mobility)",
              main = sprintf("Stepwise GAM Final Model\nR² = %.4f, RMSE = %.4f, AIC = %.0f",
                             r2_value, rmse_value, AIC(Fit.mobility)),
              cex.axis = 1.2, cex.lab = 1.4, cex.main = 1.2)

abline(a = 0, b = 1, lwd = 2, lty = 2, col = "red")
abline(lm(fitted_values ~ observed_values), lwd = 2, col = "blue")

legend("topleft",
       legend = c("1:1 Line", "Fitted Line"),
       col = c("red", "blue"),
       lty = c(2, 1), lwd = 2, cex = 1.1)

dev.off()
cat("✓ 散点图已保存: stepwise_fitted_vs_observed.png\n")

# 12.2 残差诊断图
png("stepwise_residuals.png", width = 1200, height = 900, res = 120)

par(mfrow = c(2, 2))
gam.check(Fit.mobility)

dev.off()
cat("✓ 残差诊断图已保存: stepwise_residuals.png\n")

# 12.3 偏效应图
png("stepwise_partial_effects.png", width = 1600, height = 1200, res = 120)

plot(Fit.mobility,
     pages    = 1,
     scheme   = 2,
     shade    = TRUE,
     shade.col= "lightblue",
     cex.main = 1.2,
     cex.lab  = 1.1)

dev.off()
cat("✓ 偏效应图已保存: stepwise_partial_effects.png\n")

# 12.4 变量重要性图（基于 EDF）
png("stepwise_variable_importance.png", width = 1000, height = 800, res = 120)

# 提取平滑项信息
smooth_summary <- summary(Fit.mobility)$s.table

if (!is.null(smooth_summary) && nrow(smooth_summary) > 0) {
  # 按 EDF 排序
  smooth_df <- as.data.frame(smooth_summary)
  smooth_df$Term <- rownames(smooth_df)
  smooth_df <- smooth_df[order(smooth_df$edf, decreasing = TRUE), ]
  
  par(mar = c(5, 12, 4, 2))
  
  barplot(smooth_df$edf,
          names.arg = smooth_df$Term,
          horiz = TRUE,
          las = 1,
          main = "Variable Importance (Effective Degrees of Freedom)",
          xlab = "EDF",
          col = ifelse(smooth_df$`p-value` < 0.05, "steelblue", "gray70"),
          cex.names = 0.8,
          cex.axis = 1.0,
          cex.lab = 1.2)
  
  legend("bottomright", 
         legend = c("p < 0.05", "p >= 0.05"),
         fill = c("steelblue", "gray70"),
         cex = 1.0)
}

dev.off()
cat("✓ 变量重要性图已保存: stepwise_variable_importance.png\n")

# 12.5 R² 和 Dev.Expl 变化图
png("stepwise_performance_progression.png", width = 1000, height = 600, res = 120)

par(mar = c(10, 5, 4, 5))

# 左轴：R²
plot(1:nrow(selection_history), selection_history$R2,
     type = "b", pch = 19, col = "blue", lwd = 2,
     xlab = "", ylab = expression(R^2),
     main = "Model Performance Progression",
     xaxt = "n", cex.axis = 1.2, cex.lab = 1.3,
     ylim = c(0, max(selection_history$R2, na.rm = TRUE) * 1.1))

axis(1, at = 1:nrow(selection_history), 
     labels = selection_history$Variable, 
     las = 2, cex.axis = 0.8)

mtext("Variables Added", side = 1, line = 8, cex = 1.2)

# 右轴：Dev.Expl
par(new = TRUE)
plot(1:nrow(selection_history), selection_history$Dev_Expl,
     type = "b", pch = 17, col = "red", lwd = 2,
     xlab = "", ylab = "", xaxt = "n", yaxt = "n",
     ylim = c(0, max(selection_history$Dev_Expl, na.rm = TRUE) * 1.1))

axis(4, col = "red", col.axis = "red", cex.axis = 1.2)
mtext("Deviance Explained (%)", side = 4, line = 3, col = "red", cex = 1.2)

legend("bottomright", 
       legend = c(expression(R^2), "Deviance Explained"),
       col = c("blue", "red"), 
       pch = c(19, 17), lwd = 2, cex = 1.0)

dev.off()
cat("✓ 性能变化图已保存: stepwise_performance_progression.png\n")

# ========================================
# Stepwise GAM – Publication-quality Plots
# ========================================
library(grDevices)

# 检查 Ghostscript
#if (gs_path == "") {
#  gs_path <- Sys.which("gswin64c")  # Windows
#}

#if (gs_path == "") {
#  warning("未找到 Ghostscript，将使用普通 PDF（字体可能为轮廓）")
#  use_embed <- FALSE
#} else {
#  cat("Ghostscript 路径:", gs_path, "\n")
#  use_embed <- TRUE
#}
Sys.setenv(R_GSCMD = "/opt/homebrew/bin/gs")

# 定义保存函数
save_pdf <- function(filename, width, height, plot_code) {
  
  if (use_embed) {
    # 使用嵌入字体
    temp_file <- tempfile(fileext = ".pdf")
    pdf(temp_file, width = width, height = height, useDingbats = FALSE)
    eval(parse(text = plot_code))
    dev.off()
    
    # 嵌入字体
    tryCatch({
      embedFonts(temp_file, outfile = filename,
                 options = "-dSubsetFonts=true -dEmbedAllFonts=true -dCompressFonts=true")
      unlink(temp_file)
      cat("✓ 已保存（字体已嵌入）:", filename, "\n")
    }, error = function(e) {
      file.copy(temp_file, filename, overwrite = TRUE)
      unlink(temp_file)
      cat("⚠ 已保存（嵌入失败）:", filename, "\n")
    })
    
  } else {
    # 普通 PDF
    pdf(filename, width = width, height = height, useDingbats = FALSE)
    eval(parse(text = plot_code))
    dev.off()
    cat("✓ 已保存:", filename, "\n")
  }
}

# ========================================
# 8. AIC 变化图
# ========================================

pdf("stepwise_AIC_progression_temp.pdf", width = 10, height = 5.8, useDingbats = FALSE)

par(mar = c(10, 5, 4, 2))

colors <- ifelse(selection_history$Stage == "Baseline", "gray50",
                 ifelse(selection_history$Stage == "Main", "blue", "red"))

plot(1:nrow(selection_history), selection_history$AIC,
     type = "b", pch = 19, col = colors, lwd = 2,
     xlab = "", ylab = "AIC",
     main = "Stepwise GAM Variable Selection (Main Effects + Interactions)",
     xaxt = "n", cex.axis = 1.2, cex.lab = 1.3)

axis(1, at = 1:nrow(selection_history), 
     labels = selection_history$Variable, 
     las = 2, cex.axis = 0.8)

mtext("Variables Added", side = 1, line = 8, cex = 1.2)

if (step_main_end > 0 && step_main_end < nrow(selection_history) - 1) {
  abline(v = step_main_end + 1.5, lty = 3, col = "gray50", lwd = 2)
  text(step_main_end + 1.5, max(selection_history$AIC) * 0.95, 
       "← Main | Interaction →", pos = 4, cex = 0.9)
}

for (i in 2:nrow(selection_history)) {
  delta <- selection_history$Delta_AIC[i]
  if (!is.na(delta) && delta > 0) {
    text(i, selection_history$AIC[i],
         sprintf("Δ=%.1f", delta), pos = 3, cex = 0.7, col = "darkgreen")
  }
}

legend("topright", 
       legend = c("Baseline", "Main Effects", "Interactions"),
       col = c("gray50", "blue", "red"), 
       pch = 19, lwd = 2, cex = 1.0)

dev.off()

# 嵌入字体
embedFonts("stepwise_AIC_progression_temp.pdf", 
           outfile = "stepwise_AIC_progression.pdf",
           options = "-dSubsetFonts=true -dEmbedAllFonts=true")
unlink("stepwise_AIC_progression_temp.pdf")
cat("✓ AIC 变化图已保存: stepwise_AIC_progression.pdf\n")


# ========================================
# 12.1 观测值 vs 拟合值
# ========================================

pdf("stepwise_fitted_vs_observed_temp.pdf", width = 6.67, height = 6.67, useDingbats = FALSE)

smoothScatter(observed_values, fitted_values,
              xlab = "Observed log(Mobility)",
              ylab = "Fitted log(Mobility)",
              main = sprintf("Stepwise GAM Final Model\nR² = %.4f, RMSE = %.4f, AIC = %.0f",
                             r2_value, rmse_value, AIC(Fit.mobility)),
              cex.axis = 1.2, cex.lab = 1.4, cex.main = 1.2)

abline(a = 0, b = 1, lwd = 2, lty = 2, col = "red")
abline(lm(fitted_values ~ observed_values), lwd = 2, col = "blue")

legend("topleft",
       legend = c("1:1 Line", "Fitted Line"),
       col = c("red", "blue"),
       lty = c(2, 1), lwd = 2, cex = 1.1)

dev.off()
embedFonts("stepwise_fitted_vs_observed_temp.pdf", 
           outfile = "stepwise_fitted_vs_observed.pdf",
           options = "-dSubsetFonts=true -dEmbedAllFonts=true")
unlink("stepwise_fitted_vs_observed_temp.pdf")
cat("✓ 散点图已保存: stepwise_fitted_vs_observed.pdf\n")


# ========================================
# 12.2 残差诊断图
# ========================================

residuals_model <- residuals(Fit.mobility)
fitted_model <- fitted(Fit.mobility)
linpred <- predict(Fit.mobility, type = "link")

# Q-Q plot
pdf("temp_QQ.pdf", width = 6, height = 6, useDingbats = FALSE)
qq.gam(Fit.mobility, main = "Normal Q-Q Plot", cex.main = 1.2, cex.lab = 1.1)
dev.off()
embedFonts("temp_QQ.pdf", outfile = "stepwise_residuals_QQ.pdf",
           options = "-dSubsetFonts=true -dEmbedAllFonts=true")
unlink("temp_QQ.pdf")
cat("✓ 残差诊断图(Q-Q)已保存\n")

# Residuals vs Linear Predictor
pdf("temp_LP.pdf", width = 6, height = 6, useDingbats = FALSE)
plot(linpred, residuals_model,
     xlab = "Linear Predictor", ylab = "Residuals",
     main = "Residuals vs Linear Predictor",
     pch = 20, col = rgb(0, 0, 0, 0.3),
     cex.main = 1.2, cex.lab = 1.1, cex.axis = 1.0)
abline(h = 0, lty = 2, col = "red", lwd = 1.5)
dev.off()
# Residuals vs Linear Predictor (续)
embedFonts("temp_LP.pdf", outfile = "stepwise_residuals_vs_linear_pred.pdf",
           options = "-dSubsetFonts=true -dEmbedAllFonts=true")
unlink("temp_LP.pdf")
cat("✓ 残差诊断图(Residuals vs LP)已保存\n")

# Histogram of Residuals
pdf("temp_hist.pdf", width = 6, height = 6, useDingbats = FALSE)
hist(residuals_model, breaks = 50, 
     main = "Histogram of Residuals", xlab = "Residuals",
     col = "lightblue", border = "white",
     cex.main = 1.2, cex.lab = 1.1, cex.axis = 1.0)
dev.off()
embedFonts("temp_hist.pdf", outfile = "stepwise_residuals_histogram.pdf",
           options = "-dSubsetFonts=true -dEmbedAllFonts=true")
unlink("temp_hist.pdf")
cat("✓ 残差诊断图(Histogram)已保存\n")

# Response vs Fitted Values
pdf("temp_resp.pdf", width = 6, height = 6, useDingbats = FALSE)
plot(fitted_model, observed_values,
     xlab = "Fitted Values", ylab = "Response",
     main = "Response vs Fitted Values",
     pch = 20, col = rgb(0, 0, 0, 0.3),
     cex.main = 1.2, cex.lab = 1.1, cex.axis = 1.0)
abline(a = 0, b = 1, lty = 2, col = "red", lwd = 1.5)
dev.off()
embedFonts("temp_resp.pdf", outfile = "stepwise_response_vs_fitted.pdf",
           options = "-dSubsetFonts=true -dEmbedAllFonts=true")
unlink("temp_resp.pdf")
cat("✓ 残差诊断图(Response vs Fitted)已保存\n")


# ========================================
# 12.3 偏效应图 - 每个平滑项单独保存
# ========================================
cat("\n保存偏效应图（单独保存每个项）...\n")

n_smooth_terms <- length(Fit.mobility$smooth)

if (n_smooth_terms > 0) {
  for (i in 1:n_smooth_terms) {
    term_label <- Fit.mobility$smooth[[i]]$label
    safe_label <- gsub("[^a-zA-Z0-9_]", "_", term_label)
    
    file_name <- sprintf("stepwise_partial_effect_%02d_%s.pdf", i, safe_label)
    
    # 判断维度
    term_dim <- Fit.mobility$smooth[[i]]$dim
    
    # 开始绘图
    pdf(file_name, width = 6, height = 5, useDingbats = FALSE)
    
    if (term_dim == 2) {
      # 二维项：使用自定义配色
      plot(Fit.mobility,
           select = i,
           scheme = 2,
           shade = TRUE,
           rug = FALSE,
           cex.main = 1.2,
           cex.lab = 1.1,
           cex.axis = 1.0,
           # 配色方案（可根据需要修改）
           hcolors = colorRampPalette(c("#2C7BB6", "#FFFFBF", "#D7191C"))(50))
      
      # 其他配色选项：
      # hcolors = rev(heat.colors(50))  # 反转热力图
      # hcolors = colorRampPalette(c("blue", "white", "red"))(50)  # 蓝-白-红
      
    } else {
      # 一维项：白色背景
      plot(Fit.mobility,
           select = i,
           scheme = 1,
           shade = TRUE,
           shade.col = "lightblue",
           cex.main = 1.2,
           cex.lab = 1.1,
           cex.axis = 1.0)
    }
    
    dev.off()
    cat(sprintf("  ✓ %s (dim=%d)\n", file_name, term_dim))
  }
}
# ========================================
# 因子变量效应图 - 误差线图版本
# ========================================

cat("\n--- 因子变量（误差线图）---\n")

model_data <- Fit.mobility$model
response_var <- names(model_data)[1]
formula_terms <- attr(terms(Fit.mobility), "term.labels")

# 提取因子变量
factor_vars <- c()
for (term in formula_terms) {
  if (grepl("^as\\.factor\\(", term)) {
    var_name <- gsub("^as\\.factor\\((.*)\\)$", "\\1", term)
    factor_vars <- c(factor_vars, var_name)
  }
}
factor_vars <- unique(factor_vars)

if (length(factor_vars) > 0) {
  cat(sprintf("发现 %d 个因子变量: %s\n", 
              length(factor_vars), 
              paste(factor_vars, collapse = ", ")))
  
  for (j in seq_along(factor_vars)) {
    var_name <- factor_vars[j]
    model_col_name <- paste0("as.factor(", var_name, ")")
    
    safe_name <- gsub("[^a-zA-Z0-9_]", "_", var_name)
    file_name <- sprintf("stepwise_partial_effect_factor_%02d_%s.pdf", j, safe_name)
    
    if (!model_col_name %in% names(model_data)) next
    
    actual_data <- model_data[[model_col_name]]
    factor_levels <- levels(as.factor(actual_data))
    n_levels <- length(factor_levels)
    
    # 从原始数据复制模板
    pred_list <- list()
    for (level in factor_levels) {
      template <- model_data[1, , drop = FALSE]
      template[[model_col_name]] <- factor(level, levels = factor_levels)
      pred_list[[level]] <- template
    }
    pred_data <- do.call(rbind, pred_list)
    rownames(pred_data) <- NULL
    
    tryCatch({
      pred <- predict(Fit.mobility, newdata = pred_data, type = "response", se.fit = TRUE)
      
      pred_df <- data.frame(
        level = factor(factor_levels, levels = factor_levels),
        fit = pred$fit,
        se = pred$se.fit,
        lower = pred$fit - 1.96 * pred$se.fit,
        upper = pred$fit + 1.96 * pred$se.fit
      )
      
      # ========================================
      # 误差线图（点 + 置信区间）
      # ========================================
      pdf(file_name, width = max(6, n_levels * 0.6), height = 5)
      par(mar = c(6, 5, 4, 2))
      
      # 设置 x 坐标
      x_pos <- 1:n_levels
      
      # 绘制空白画布
      plot(x_pos, pred_df$fit,
           type = "n",
           xlim = c(0.5, n_levels + 0.5),
           ylim = c(min(pred_df$lower) * 0.9, max(pred_df$upper) * 1.1),
           xaxt = "n",
           xlab = "",
           ylab = "Predicted Response",
           main = paste("Partial Effect:", var_name),
           cex.main = 1.3,
           cex.lab = 1.2)
      
      # 添加水平参考线（整体均值）
      abline(h = mean(pred_df$fit), lty = 2, col = "gray50", lwd = 1)
      
      # 绘制误差线
      arrows(x_pos, pred_df$lower, x_pos, pred_df$upper,
             angle = 90, code = 3, length = 0.1, lwd = 2, col = "steelblue")
      
      # 绘制点
      points(x_pos, pred_df$fit, pch = 19, cex = 2, col = "steelblue")
      
      # 添加 x 轴标签
      axis(1, at = x_pos, labels = pred_df$level, las = 2, cex.axis = 0.9)
      
      # 添加数值标签
      text(x_pos, pred_df$upper + (max(pred_df$upper) - min(pred_df$lower)) * 0.05,
           labels = round(pred_df$fit, 2), cex = 0.8, col = "gray30")
      
      mtext(var_name, side = 1, line = 4.5, cex = 1.1)
      
      dev.off()
      cat(sprintf("  ✓ [误差线图] %s (levels=%d)\n", file_name, n_levels))
      
    }, error = function(e) {
      cat(sprintf("  ✗ %s 失败: %s\n", var_name, e$message))
      if (dev.cur() > 1) dev.off()
    })
  }
}

cat("✓ 偏效应图已全部保存\n")


# ========================================
# 12.4 变量重要性图（基于卡方统计量）
# ========================================
# ========================================
# 12.4 变量重要性图（基于 EDF）- 修正版
# ========================================

pdf("temp_importance.pdf", width = 12, height = 8, useDingbats = FALSE)

# 提取平滑项信息
smooth_summary <- summary(Fit.mobility)$s.table

if (!is.null(smooth_summary) && nrow(smooth_summary) > 0) {
  
  # 转换为数据框
  smooth_df <- as.data.frame(smooth_summary)
  smooth_df$Term <- rownames(smooth_df)
  
  # 确保 edf 是数值向量
  smooth_df$edf <- as.numeric(smooth_df$edf)
  
  # 按 EDF 排序
  smooth_df <- smooth_df[order(smooth_df$edf, decreasing = TRUE), ]
  
  # 简化变量名（如果太长）
  smooth_df$Term_short <- ifelse(
    nchar(smooth_df$Term) > 30,
    paste0(substr(smooth_df$Term, 1, 27), "..."),
    smooth_df$Term
  )
  
  # 动态调整左边距
  max_name_len <- max(nchar(smooth_df$Term_short))
  left_margin <- min(max_name_len * 0.4, 12)
  
  par(mar = c(5, left_margin, 4, 2))
  
  # 查找 p-value 列（可能是 "p-value" 或 "p.value"）
  pval_col <- NULL
  if ("p-value" %in% colnames(smooth_df)) {
    pval_col <- "p-value"
  } else if ("p.value" %in% colnames(smooth_df)) {
    pval_col <- "p.value"
  }
  
  # 设置颜色
  if (!is.null(pval_col)) {
    pvals <- as.numeric(smooth_df[[pval_col]])
    bar_colors <- ifelse(pvals < 0.05, "steelblue", "gray70")
    show_legend <- TRUE
  } else {
    bar_colors <- "steelblue"
    show_legend <- FALSE
  }
  
  # 绘制条形图
  barplot(smooth_df$edf,
          names.arg = smooth_df$Term_short,
          horiz = TRUE,
          las = 1,
          main = "Variable Importance (Effective Degrees of Freedom)",
          xlab = "EDF",
          col = bar_colors,
          cex.names = 0.55,
          cex.axis = 1.0,
          cex.lab = 1.2)
  
  # 添加图例（如果有 p-value）
  if (show_legend) {
    legend("bottomright", 
           legend = c("p < 0.05", "p ≥ 0.05"),
           fill = c("steelblue", "gray70"),
           cex = 0.9)
  }
}

dev.off()

# 嵌入字体
if (use_embed) {
  tryCatch({
    embedFonts("temp_importance.pdf", 
               outfile = "stepwise_variable_importance.pdf",
               options = "-dSubsetFonts=true -dEmbedAllFonts=true")
    unlink("temp_importance.pdf")
    cat("✓ 已保存（字体已嵌入）: stepwise_variable_importance.pdf\n")
  }, error = function(e) {
    file.rename("temp_importance.pdf", "stepwise_variable_importance.pdf")
    cat("⚠ 嵌入失败，使用原始文件: stepwise_variable_importance.pdf\n")
  })
} else {
  file.rename("temp_importance.pdf", "stepwise_variable_importance.pdf")
  cat("✓ 已保存: stepwise_variable_importance.pdf\n")
}


# ========================================
# 12.5 R² 和 Dev.Expl 变化图
# ========================================

pdf("temp_performance.pdf", width = 8.33, height = 5, useDingbats = FALSE)

par(mar = c(10, 5, 4, 5))

plot(1:nrow(selection_history), selection_history$R2,
     type = "b", pch = 19, col = "blue", lwd = 2,
     xlab = "", ylab = expression(R^2),
     main = "Model Performance Progression",
     xaxt = "n", cex.axis = 1.2, cex.lab = 1.3,
     ylim = c(0, max(selection_history$R2, na.rm = TRUE) * 1.1))

axis(1, at = 1:nrow(selection_history), 
     labels = selection_history$Variable, 
     las = 2, cex.axis = 0.8)

mtext("Variables Added", side = 1, line = 8, cex = 1.2)

par(new = TRUE)
plot(1:nrow(selection_history), selection_history$Dev_Expl,
     type = "b", pch = 17, col = "red", lwd = 2,
     xlab = "", ylab = "", xaxt = "n", yaxt = "n",
     ylim = c(0, max(selection_history$Dev_Expl, na.rm = TRUE) * 1.1))

axis(4, col = "red", col.axis = "red", cex.axis = 1.2)
mtext("Deviance Explained (%)", side = 4, line = 3, col = "red", cex = 1.2)

legend("bottomright", 
       legend = c(expression(R^2), "Deviance Explained"),
       col = c("blue", "red"), 
       pch = c(19, 17), lwd = 2, cex = 1.0)

dev.off()
embedFonts("temp_performance.pdf", outfile = "stepwise_performance_progression.pdf",
           options = "-dSubsetFonts=true -dEmbedAllFonts=true")
unlink("temp_performance.pdf")
cat("✓ 性能变化图已保存: stepwise_performance_progression.pdf\n")


# ========================================
# 确认文件生成
# ========================================
cat("\n========================================\n")
cat("所有图片保存完成！\n")
cat("========================================\n")

pdf_files <- list.files(getwd(), pattern = "\\.pdf$")
cat(sprintf("共 %d 个 PDF 文件:\n", length(pdf_files)))
print(pdf_files)

cat("\n保存位置:", getwd(), "\n")



# ========================================
# 13. 生成最终公式字符串（方便复制）
# ========================================

cat("\n", paste(rep("=", 80), collapse = ""), "\n", sep = "")
cat("最终模型公式（可直接复制使用）\n")
cat(paste(rep("=", 80), collapse = ""), "\n", sep = "")

# 生成可复制的公式代码
formula_code <- "Formula.final <- bquote(\n  log_Mobility ~\n"

# 添加连续变量
for (i in seq_along(selected_smooth)) {
  var <- selected_smooth[i]
  if (var == "Slope_m") {
    formula_code <- paste0(formula_code, 
                           sprintf("    s(%s, bs = 'ad', k = .(num.baseFun))", var))
  } else {
    formula_code <- paste0(formula_code, 
                           sprintf("    s(%s, k = .(num.baseFun))", var))
  }
  formula_code <- paste0(formula_code, " +\n")
}

# 添加交互项
for (i in seq_along(selected_ti)) {
  ti_item <- selected_ti[[i]]
  formula_code <- paste0(formula_code,
                         sprintf("    ti(%s, %s, k = %d)", 
                                 ti_item$vars[1], ti_item$vars[2], ti_item$k))
  formula_code <- paste0(formula_code, " +\n")
}

# 添加分类变量
for (var in selected_factor) {
  formula_code <- paste0(formula_code,
                         sprintf("    as.factor(%s) +\n", var))
}

# 添加空间项
formula_code <- paste0(formula_code,
                       "    s(CoordX, CoordY, k = .(num.spaceFun), bs = 'ds')\n)")

cat("\n", formula_code, "\n", sep = "")

# 保存公式代码
writeLines(formula_code, "stepwise_final_formula.R")
cat("\n✓ 公式代码已保存: stepwise_final_formula.R\n")

# ========================================
# 14. 最终报告
# ========================================

cat("\n", paste(rep("=", 80), collapse = ""), "\n", sep = "")
cat("Stepwise GAM 变量选择完成！\n")
cat(paste(rep("=", 80), collapse = ""), "\n", sep = "")

cat("\n📊 选择结果汇总:\n")
cat(sprintf("  候选变量总数: %d\n", 
            length(available_smooth) + length(available_factor) + length(available_ti_terms)))
cat(sprintf("    - 连续变量: %d\n", length(available_smooth)))
cat(sprintf("    - 分类变量: %d\n", length(available_factor)))
cat(sprintf("    - 交互项:   %d\n", length(available_ti_terms)))

cat(sprintf("\n  选中变量总数: %d\n", 
            length(selected_smooth) + length(selected_factor) + length(selected_ti)))
cat(sprintf("    - 连续变量: %d / %d\n", length(selected_smooth), length(available_smooth)))
cat(sprintf("    - 分类变量: %d / %d\n", length(selected_factor), length(available_factor)))
cat(sprintf("    - 交互项:   %d / %d\n", length(selected_ti), length(available_ti_terms)))

cat("\n📈 最终模型性能:\n")
cat(sprintf("  R²:           %.4f\n", r2_value))
cat(sprintf("  RMSE:         %.4f\n", rmse_value))
cat(sprintf("  MAE:          %.4f\n", mae_value))
cat(sprintf("  NSE:          %.4f\n", nse))
cat(sprintf("  Dev.Expl:     %.2f%%\n", summary(Fit.mobility)$dev.expl * 100))
cat(sprintf("  AIC:          %.2f\n", AIC(Fit.mobility)))

cat("\n📁 输出文件:\n")
cat("  • stepwise_selection_history.csv    (选择历史)\n")
cat("  • stepwise_AIC_progression.png      (AIC 变化图)\n")
cat("  • stepwise_performance_progression.png (R²/DevExpl 变化图)\n")
cat("  • stepwise_final_model.rds          (最终模型对象)\n")
cat("  • stepwise_final_model_summary.txt  (模型摘要)\n")
cat("  • stepwise_final_formula.R          (最终公式代码)\n")
cat("  • stepwise_predictions.csv          (预测结果)\n")
cat("  • stepwise_metrics.csv              (性能指标)\n")
cat("  • stepwise_fitted_vs_observed.png   (拟合散点图)\n")
cat("  • stepwise_residuals.png            (残差诊断图)\n")
cat("  • stepwise_partial_effects.png      (偏效应图)\n")
cat("  • stepwise_variable_importance.png  (变量重要性图)\n")

cat("\n💡 使用说明:\n")
cat("  1. 查看 stepwise_selection_history.csv 了解变量选择过程\n")
cat("  2. 使用 stepwise_final_formula.R 中的公式代码进行后续建模\n")
cat("  3. 加载模型: model <- readRDS('stepwise_final_model.rds')\n")

cat("\n🔧 可调参数:\n")
cat(sprintf("  AIC_threshold_main = %d  (主效应选择阈值)\n", AIC_threshold_main))
cat(sprintf("  AIC_threshold_ti   = %d  (交互项选择阈值)\n", AIC_threshold_ti))
cat("  可根据需要调整阈值重新运行\n")

cat("\n工作目录:", getwd(), "\n")
cat(paste(rep("=", 80), collapse = ""), "\n", sep = "")

# ========================================
# 15. 可选：与原始全变量模型对比
# ========================================

cat("\n", paste(rep("-", 60), collapse = ""), "\n", sep = "")
cat("可选：与原始全变量模型对比\n")
cat(paste(rep("-", 60), collapse = ""), "\n", sep = "")

run_comparison <- FALSE  # 设为 TRUE 启用对比

if (run_comparison) {
  cat("\n拟合全变量模型进行对比...\n")
  
  # 构建全变量公式
  full_formula <- build_formula_full(
    available_smooth, 
    available_factor, 
    available_ti_terms,
    num.baseFun, 
    spatial_term
  )
  
  full_model <- tryCatch({
    mgcv::bam(
      full_formula,
      family = gaussian(),
      data = data4fit,
      method = "fREML",
      select = TRUE,
      nthreads = 4
    )
  }, error = function(e) {
    cat("全变量模型拟合失败:", e$message, "\n")
    return(NULL)
  })
  
  if (!is.null(full_model)) {
    full_AIC <- AIC(full_model)
    full_r2 <- summary(full_model)$r.sq
    full_dev <- summary(full_model)$dev.expl * 100
    
    cat("\n模型对比:\n")
    cat(sprintf("                    Stepwise      Full\n"))
    cat(sprintf("  变量数:           %d            %d\n", 
                length(selected_smooth) + length(selected_factor) + length(selected_ti),
                length(available_smooth) + length(available_factor) + length(available_ti_terms)))
    cat(sprintf("  AIC:              %.2f        %.2f\n", AIC(Fit.mobility), full_AIC))
    cat(sprintf("  R²:               %.4f        %.4f\n", r2_value, full_r2))
    cat(sprintf("  Dev.Expl:         %.2f%%       %.2f%%\n", 
                summary(Fit.mobility)$dev.expl * 100, full_dev))
    cat(sprintf("  AIC 差异:         %.2f (正值表示 Stepwise 更优)\n", 
                full_AIC - AIC(Fit.mobility)))
  }
} else {
  cat("跳过全变量模型对比（设置 run_comparison <- TRUE 启用）\n")
}

cat("\n", paste(rep("=", 80), collapse = ""), "\n", sep = "")
cat("全部完成！\n")
cat(paste(rep("=", 80), collapse = ""), "\n", sep = "")

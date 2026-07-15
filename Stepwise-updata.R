

# ========================================
# 滑坡流动性 GAM 建模
# Stepwise AIC 变量选择稳定性分析
# 每个 SCV 训练折内重复执行变量筛选
# 最终仅输出稳定性分析结果和图表
# ========================================

library(viridis)

theme_set(
  theme_bw(base_size = 12) +
    theme(
      panel.grid = element_blank(),
      plot.title = element_text(face = "bold"),
      axis.text = element_text(color = "black")
    )
)

color_stable <- "#2C7FB8"   # 稳定变量蓝色
color_weak   <- "#F5F5F5"   # 弱变量浅灰
color_alert  <- "#D62728"   # 阈值红线
rm(list = ls())

library(mgcv)
library(sf)
library(ggplot2)

# ========================================
# 0. 参数设置
# ========================================

num.baseFun  <- 45
num.spaceFun <- 180

# AIC 改善阈值
AIC_threshold_main <- 2
AIC_threshold_ti   <- 2

# 交叉验证折数
K_fold <- 10
set.seed(20240625)

# 如果已有空间交叉验证折叠列，请在这里填写列名
# 例如 fold_column <- "SCV_fold"
# 如果为 NULL，代码会自动检测常见列名；若未发现，则基于 CoordX/CoordY 自动生成空间折叠
fold_column <- NULL

# 二次精简参数
# 你前面描述中使用的是“重要性较低且显著性大于 0.5 的因子”
# 因此这里默认 p > 0.50 才删除；如果你想使用传统显著性标准，可改为 0.05
p_remove_threshold <- 0.50

# EDF 低重要性阈值：默认将 EDF 排名最低的 25% 视为低贡献项
# 如果你想删除所有 p > p_remove_threshold 的项，可设为 1.00
edf_low_quantile <- 0.25

# 稳定性阈值：>=80% folds 被选中，认为较稳定
stability_threshold <- 0.80

# 交互项候选规则
# "all_predefined"：使用所有预先定义的物理交互项，与你当前代码实际效果最接近
# "at_least_one_selected_main"：至少一个变量已被主效应选择
# "both_selected_main"：两个变量都已被主效应选择，更保守
ti_candidate_rule <- "all_predefined"

# 输出目录
OUT_DIR <- "stability_selection_outputs"
if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)

# ========================================
# 1. 读取数据
# ========================================

WORK <- "/Volumes/E/SU_Codes/Mobility/LuDmob/I1_org_pdf/I1_org-stablity"
setwd(WORK)

cat("读取 Shapefile 数据...\n")
su_data <- st_read("su_100k_02_clean_with_phy_mobility_cut1.shp")
data.raw <- st_drop_geometry(su_data)

cat("原始数据维度:", dim(data.raw), "\n")

# ========================================
# 2. 数据预处理
# ========================================

cat("\n数据预处理...\n")

data4fit <- subset(data.raw, I1_orig > 0 & Area_sum > 0)####

categorical_vars <- c("Landforms", "Lithology", "LULC")
for (var in categorical_vars) {
  if (var %in% names(data4fit)) {
    data4fit[[var]] <- as.factor(data4fit[[var]])
  }
}

data4fit$log_Mobility <- log1p(data4fit$I1_orig) #####
data4fit$log_Area     <- log1p(data4fit$Area_sum)

data4fit$PGA_total <- sqrt(data4fit$PGA_EW_m^2 + data4fit$PGA_NS_m^2)

data4fit <- na.omit(data4fit)

cat("最终样本数:", nrow(data4fit), "\n")

if (!all(c("CoordX", "CoordY") %in% names(data4fit))) {
  stop("数据中缺少 CoordX 或 CoordY，无法构建空间项或空间折叠。")
}

# ========================================
# 3. 定义候选变量
# ========================================

candidate_smooth_vars <- c(
  "East_m", "Fault_m", "North_m", "PGA_total", "PGA_UD_s",
  "PLC_s", "PRC_m", "PRC_s", "River_m", "River_s",
  "Road_m", "Slope_m", "TWI_s", "ShRelief_s", "delta_z",
  "V_Std", "V_Mean", "PLC_m", "PGA_EW_m", "travel_dis"
)

candidate_ti_terms <- list(
  list(vars = c("PLC_m", "PRC_m"),        k = 15),
  list(vars = c("Fault_m", "PGA_UD_s"),   k = 15),
  list(vars = c("Fault_m", "PGA_EW_m"),   k = 15),
  list(vars = c("PGA_UD_s", "PGA_EW_m"),  k = 10),
  list(vars = c("delta_z", "travel_dis"), k = 10),
  list(vars = c("PGA_total", "V_Std"),    k = 10),
  list(vars = c("PGA_UD_s", "V_Std"),     k = 15),
  list(vars = c("PGA_EW_m", "V_Std"),     k = 15),
  list(vars = c("delta_z", "V_Std"),      k = 8),
  list(vars = c("PGA_total", "Slope_m"),  k = 10),
  list(vars = c("Slope_m", "ShRelief_s"), k = 10),
  list(vars = c("TWI_s", "Slope_m"),      k = 10),
  list(vars = c("River_m", "Slope_m"),    k = 10),
  list(vars = c("Fault_m", "Slope_m"),    k = 10)
)

candidate_factor_vars <- c("Landforms", "Lithology", "LULC")

spatial_term <- sprintf("s(CoordX, CoordY, k = %d, bs = 'ds')", num.spaceFun)

# ========================================
# 4. 工具函数
# ========================================

smooth_to_string <- function(var) {
  sprintf("s(%s)", var)
}

factor_to_string <- function(var) {
  sprintf("factor(%s)", var)
}

ti_to_string <- function(ti_item) {
  sprintf("ti(%s,%s)", ti_item$vars[1], ti_item$vars[2])
}

normalize_term <- function(x) {
  gsub("\\s+", "", x)
}

term_type <- function(term) {
  if (startsWith(term, "ti(")) {
    return("interaction")
  } else if (startsWith(term, "s(")) {
    return("smooth")
  } else if (startsWith(term, "factor(")) {
    return("factor")
  } else {
    return("other")
  }
}

# 构建公式
build_formula_full <- function(selected_smooth, selected_factor, selected_ti,
                               k_smooth, spatial_term) {
  terms <- c()
  
  for (var in selected_smooth) {
    if (var == "Slope_m") {
      terms <- c(terms, sprintf("s(%s, bs = 'ad', k = %d)", var, k_smooth))
    } else {
      terms <- c(terms, sprintf("s(%s, k = %d)", var, k_smooth))
    }
  }
  
  for (ti_item in selected_ti) {
    terms <- c(
      terms,
      sprintf(
        "ti(%s, %s, k = %d)",
        ti_item$vars[1],
        ti_item$vars[2],
        ti_item$k
      )
    )
  }
  
  for (var in selected_factor) {
    terms <- c(terms, sprintf("as.factor(%s)", var))
  }
  
  terms <- c(terms, spatial_term)
  
  formula_str <- paste("log_Mobility ~", paste(terms, collapse = " + "))
  as.formula(formula_str)
}

# 安全拟合
safe_fit <- function(formula, data) {
  tryCatch({
    model <- mgcv::bam(
      formula,
      family   = gaussian(),
      data     = data,
      method   = "fREML",
      select   = TRUE,
      nthreads = 4
    )
    return(model)
  }, error = function(e) {
    return(NULL)
  })
}

# 提取 smooth table
extract_smooth_table <- function(model) {
  sm <- summary(model)
  
  if (is.null(sm$s.table) || nrow(sm$s.table) == 0) {
    return(data.frame())
  }
  
  st <- as.data.frame(sm$s.table)
  st$Term_raw  <- rownames(st)
  st$Term_norm <- normalize_term(st$Term_raw)
  
  p_col <- NULL
  if ("p-value" %in% colnames(st)) {
    p_col <- "p-value"
  } else if ("p.value" %in% colnames(st)) {
    p_col <- "p.value"
  } else {
    candidate_cols <- grep("^p", colnames(st), value = TRUE)
    if (length(candidate_cols) > 0) p_col <- candidate_cols[1]
  }
  
  if (!is.null(p_col)) {
    st$p_value <- as.numeric(st[[p_col]])
  } else {
    st$p_value <- NA_real_
  }
  
  st$edf <- as.numeric(st$edf)
  
  st
}

# 根据 EDF + 显著性进行二次精简
post_selection_pruning <- function(model,
                                   selected_smooth,
                                   selected_factor,
                                   selected_ti,
                                   data,
                                   p_threshold = 0.50,
                                   edf_quantile = 0.25,
                                   refit_after_pruning = TRUE) {
  
  selected_smooth_terms <- smooth_to_string(selected_smooth)
  selected_ti_terms     <- sapply(selected_ti, ti_to_string)
  
  selected_terms <- c(selected_smooth_terms, selected_ti_terms)
  
  st <- extract_smooth_table(model)
  
  if (nrow(st) == 0 || length(selected_terms) == 0) {
    final_terms <- c(
      smooth_to_string(selected_smooth),
      factor_to_string(selected_factor),
      sapply(selected_ti, ti_to_string)
    )
    
    return(list(
      selected_smooth = selected_smooth,
      selected_factor = selected_factor,
      selected_ti     = selected_ti,
      model           = model,
      diagnostics     = data.frame(),
      selected_terms  = final_terms
    ))
  }
  
  selected_norm <- normalize_term(selected_terms)
  
  diag <- st[st$Term_norm %in% selected_norm, , drop = FALSE]
  
  if (nrow(diag) == 0) {
    final_terms <- c(
      smooth_to_string(selected_smooth),
      factor_to_string(selected_factor),
      sapply(selected_ti, ti_to_string)
    )
    
    return(list(
      selected_smooth = selected_smooth,
      selected_factor = selected_factor,
      selected_ti     = selected_ti,
      model           = model,
      diagnostics     = data.frame(),
      selected_terms  = final_terms
    ))
  }
  
  term_lookup <- data.frame(
    Term      = selected_terms,
    Term_norm = normalize_term(selected_terms),
    stringsAsFactors = FALSE
  )
  
  diag$Term <- term_lookup$Term[match(diag$Term_norm, term_lookup$Term_norm)]
  
  diag <- diag[order(diag$edf, decreasing = TRUE), ]
  diag$EDF_rank_desc <- seq_len(nrow(diag))
  
  if (edf_quantile >= 1) {
    edf_cutoff <- Inf
  } else {
    edf_cutoff <- as.numeric(
      quantile(diag$edf, probs = edf_quantile, na.rm = TRUE)
    )
  }
  
  diag$Low_EDF <- diag$edf <= edf_cutoff
  diag$Weak_support <- !is.na(diag$p_value) & diag$p_value > p_threshold
  diag$Remove <- diag$Low_EDF & diag$Weak_support
  
  remove_terms_norm <- diag$Term_norm[diag$Remove]
  
  if (length(remove_terms_norm) > 0) {
    selected_smooth <- selected_smooth[
      !(normalize_term(smooth_to_string(selected_smooth)) %in% remove_terms_norm)
    ]
    
    if (length(selected_ti) > 0) {
      selected_ti <- selected_ti[
        !(normalize_term(sapply(selected_ti, ti_to_string)) %in% remove_terms_norm)
      ]
    }
  }
  
  pruned_model <- model
  
  if (refit_after_pruning && length(remove_terms_norm) > 0) {
    pruned_formula <- build_formula_full(
      selected_smooth,
      selected_factor,
      selected_ti,
      num.baseFun,
      spatial_term
    )
    
    tmp_model <- safe_fit(pruned_formula, data)
    if (!is.null(tmp_model)) {
      pruned_model <- tmp_model
    }
  }
  
  final_terms <- c(
    smooth_to_string(selected_smooth),
    factor_to_string(selected_factor),
    sapply(selected_ti, ti_to_string)
  )
  
  final_terms <- final_terms[!is.na(final_terms)]
  
  list(
    selected_smooth = selected_smooth,
    selected_factor = selected_factor,
    selected_ti     = selected_ti,
    model           = pruned_model,
    diagnostics     = diag,
    selected_terms  = final_terms
  )
}

# 检查当前训练集中可用变量
get_available_candidates <- function(data) {
  available_smooth <- candidate_smooth_vars[candidate_smooth_vars %in% names(data)]
  
  available_smooth <- available_smooth[
    sapply(available_smooth, function(v) {
      length(unique(data[[v]])) > 3 && sd(data[[v]], na.rm = TRUE) > 0
    })
  ]
  
  available_factor <- candidate_factor_vars[candidate_factor_vars %in% names(data)]
  
  available_factor <- available_factor[
    sapply(available_factor, function(v) {
      length(unique(data[[v]])) > 1
    })
  ]
  
  available_ti_terms <- list()
  
  for (ti_item in candidate_ti_terms) {
    v1 <- ti_item$vars[1]
    v2 <- ti_item$vars[2]
    
    if (
      v1 %in% names(data) &&
      v2 %in% names(data) &&
      length(unique(data[[v1]])) > 3 &&
      length(unique(data[[v2]])) > 3 &&
      sd(data[[v1]], na.rm = TRUE) > 0 &&
      sd(data[[v2]], na.rm = TRUE) > 0
    ) {
      available_ti_terms <- c(available_ti_terms, list(ti_item))
    }
  }
  
  list(
    smooth = available_smooth,
    factor = available_factor,
    ti     = available_ti_terms
  )
}

# ========================================
# 5. 单次 Stepwise AIC 筛选函数
# ========================================

run_stepwise_selection <- function(data,
                                   fold_label = "Full",
                                   verbose = TRUE) {
  
  cand <- get_available_candidates(data)
  
  available_smooth <- cand$smooth
  available_factor <- cand$factor
  available_ti     <- cand$ti
  
  selected_smooth <- c()
  selected_factor <- c()
  selected_ti     <- list()
  
  remaining_smooth <- available_smooth
  remaining_factor <- available_factor
  
  selection_history <- data.frame(
    Fold      = character(),
    Step      = integer(),
    Stage     = character(),
    Variable  = character(),
    Type      = character(),
    AIC       = numeric(),
    Delta_AIC = numeric(),
    R2        = numeric(),
    Dev_Expl  = numeric(),
    stringsAsFactors = FALSE
  )
  
  if (verbose) {
    cat("\n", paste(rep("=", 80), collapse = ""), "\n", sep = "")
    cat("Fold:", fold_label, "\n")
    cat("训练样本数:", nrow(data), "\n")
    cat("可用连续变量:", length(available_smooth), "\n")
    cat("可用分类变量:", length(available_factor), "\n")
    cat("可用交互项:", length(available_ti), "\n")
    cat(paste(rep("=", 80), collapse = ""), "\n", sep = "")
  }
  
  # baseline model
  base_formula <- as.formula(paste("log_Mobility ~", spatial_term))
  base_model <- safe_fit(base_formula, data)
  
  if (is.null(base_model)) {
    warning(sprintf("Fold %s: baseline model failed.", fold_label))
    return(NULL)
  }
  
  current_model <- base_model
  current_AIC   <- AIC(base_model)
  
  selection_history <- rbind(selection_history, data.frame(
    Fold      = fold_label,
    Step      = 0,
    Stage     = "Baseline",
    Variable  = "Spatial_Only",
    Type      = "baseline",
    AIC       = current_AIC,
    Delta_AIC = NA,
    R2        = summary(base_model)$r.sq,
    Dev_Expl  = summary(base_model)$dev.expl * 100,
    stringsAsFactors = FALSE
  ))
  
  # ------------------------------
  # 第一阶段：主效应选择
  # ------------------------------
  
  step <- 0
  continue_main <- TRUE
  
  while (
    continue_main &&
    (length(remaining_smooth) > 0 || length(remaining_factor) > 0)
  ) {
    
    step <- step + 1
    
    best_AIC   <- Inf
    best_var   <- NULL
    best_type  <- NULL
    best_model <- NULL
    
    # 连续变量
    for (var in remaining_smooth) {
      test_smooth <- c(selected_smooth, var)
      
      test_formula <- build_formula_full(
        test_smooth,
        selected_factor,
        selected_ti,
        num.baseFun,
        spatial_term
      )
      
      test_model <- safe_fit(test_formula, data)
      
      if (!is.null(test_model)) {
        test_AIC <- AIC(test_model)
        
        if (test_AIC < best_AIC) {
          best_AIC   <- test_AIC
          best_var   <- var
          best_type  <- "smooth"
          best_model <- test_model
        }
      }
    }
    
    # 分类变量
    for (var in remaining_factor) {
      test_factor <- c(selected_factor, var)
      
      test_formula <- build_formula_full(
        selected_smooth,
        test_factor,
        selected_ti,
        num.baseFun,
        spatial_term
      )
      
      test_model <- safe_fit(test_formula, data)
      
      if (!is.null(test_model)) {
        test_AIC <- AIC(test_model)
        
        if (test_AIC < best_AIC) {
          best_AIC   <- test_AIC
          best_var   <- var
          best_type  <- "factor"
          best_model <- test_model
        }
      }
    }
    
    if (!is.null(best_var)) {
      delta_AIC <- current_AIC - best_AIC
      
      if (verbose) {
        cat(sprintf(
          "[%s] Main Step %d: best = %s, AIC %.2f -> %.2f, Delta = %.2f\n",
          fold_label, step, best_var, current_AIC, best_AIC, delta_AIC
        ))
      }
      
      if (delta_AIC > AIC_threshold_main) {
        if (best_type == "smooth") {
          selected_smooth <- c(selected_smooth, best_var)
          remaining_smooth <- setdiff(remaining_smooth, best_var)
        } else if (best_type == "factor") {
          selected_factor <- c(selected_factor, best_var)
          remaining_factor <- setdiff(remaining_factor, best_var)
        }
        
        current_model <- best_model
        current_AIC   <- best_AIC
        
        selection_history <- rbind(selection_history, data.frame(
          Fold      = fold_label,
          Step      = step,
          Stage     = "Main",
          Variable  = best_var,
          Type      = best_type,
          AIC       = best_AIC,
          Delta_AIC = delta_AIC,
          R2        = summary(best_model)$r.sq,
          Dev_Expl  = summary(best_model)$dev.expl * 100,
          stringsAsFactors = FALSE
        ))
        
      } else {
        continue_main <- FALSE
      }
      
    } else {
      continue_main <- FALSE
    }
  }
  
  # ------------------------------
  # 第二阶段：交互项选择
  # ------------------------------
  
  if (ti_candidate_rule == "all_predefined") {
    remaining_ti <- available_ti
  } else if (ti_candidate_rule == "at_least_one_selected_main") {
    remaining_ti <- list()
    for (ti_item in available_ti) {
      if (any(ti_item$vars %in% selected_smooth)) {
        remaining_ti <- c(remaining_ti, list(ti_item))
      }
    }
  } else if (ti_candidate_rule == "both_selected_main") {
    remaining_ti <- list()
    for (ti_item in available_ti) {
      if (all(ti_item$vars %in% selected_smooth)) {
        remaining_ti <- c(remaining_ti, list(ti_item))
      }
    }
  } else {
    stop("ti_candidate_rule 设置错误。")
  }
  
  continue_ti <- TRUE
  
  while (continue_ti && length(remaining_ti) > 0) {
    
    step <- step + 1
    
    best_AIC    <- Inf
    best_ti     <- NULL
    best_ti_idx <- NULL
    best_model  <- NULL
    
    for (i in seq_along(remaining_ti)) {
      ti_item <- remaining_ti[[i]]
      test_ti <- c(selected_ti, list(ti_item))
      
      test_formula <- build_formula_full(
        selected_smooth,
        selected_factor,
        test_ti,
        num.baseFun,
        spatial_term
      )
      
      test_model <- safe_fit(test_formula, data)
      
      if (!is.null(test_model)) {
        test_AIC <- AIC(test_model)
        
        if (test_AIC < best_AIC) {
          best_AIC    <- test_AIC
          best_ti     <- ti_item
          best_ti_idx <- i
          best_model  <- test_model
        }
      }
    }
    
    if (!is.null(best_ti)) {
      delta_AIC <- current_AIC - best_AIC
      
      if (verbose) {
        cat(sprintf(
          "[%s] Interaction Step %d: best = %s, AIC %.2f -> %.2f, Delta = %.2f\n",
          fold_label, step, ti_to_string(best_ti), current_AIC, best_AIC, delta_AIC
        ))
      }
      
      if (delta_AIC > AIC_threshold_ti) {
        selected_ti <- c(selected_ti, list(best_ti))
        remaining_ti <- remaining_ti[-best_ti_idx]
        
        current_model <- best_model
        current_AIC   <- best_AIC
        
        selection_history <- rbind(selection_history, data.frame(
          Fold      = fold_label,
          Step      = step,
          Stage     = "Interaction",
          Variable  = ti_to_string(best_ti),
          Type      = "ti",
          AIC       = best_AIC,
          Delta_AIC = delta_AIC,
          R2        = summary(best_model)$r.sq,
          Dev_Expl  = summary(best_model)$dev.expl * 100,
          stringsAsFactors = FALSE
        ))
        
      } else {
        continue_ti <- FALSE
      }
      
    } else {
      continue_ti <- FALSE
    }
  }
  
  # ------------------------------
  # 第三阶段：EDF + 显著性二次精简
  # ------------------------------
  
  prune_result <- post_selection_pruning(
    model           = current_model,
    selected_smooth = selected_smooth,
    selected_factor = selected_factor,
    selected_ti     = selected_ti,
    data            = data,
    p_threshold     = p_remove_threshold,
    edf_quantile    = edf_low_quantile,
    refit_after_pruning = TRUE
  )
  
  list(
    fold_label       = fold_label,
    selected_smooth  = prune_result$selected_smooth,
    selected_factor  = prune_result$selected_factor,
    selected_ti      = prune_result$selected_ti,
    selected_terms   = prune_result$selected_terms,
    model            = prune_result$model,
    selection_history = selection_history,
    pruning_diagnostics = prune_result$diagnostics
  )
}

# ========================================
# 6. 构建或读取空间交叉验证折叠
# ========================================

make_or_read_spatial_folds <- function(data, K = 10, fold_column = NULL) {
  
  if (!is.null(fold_column)) {
    if (!fold_column %in% names(data)) {
      stop(sprintf("指定的 fold_column = '%s' 不存在于数据中。", fold_column))
    }
    
    fold_id <- data[[fold_column]]
    return(as.factor(fold_id))
  }
  
  candidate_fold_cols <- c(
    "SCV_fold", "scv_fold", "Fold", "fold",
    "CV_fold", "cv_fold", "spatial_fold", "SpatialFold"
  )
  
  detected <- candidate_fold_cols[candidate_fold_cols %in% names(data)]
  
  if (length(detected) > 0) {
    message(sprintf("检测到已有折叠列: %s，将直接使用。", detected[1]))
    return(as.factor(data[[detected[1]]]))
  }
  
  message("未检测到已有 SCV 折叠列，将基于 CoordX/CoordY 使用 k-means 生成空间折叠。")
  
  coords <- scale(data[, c("CoordX", "CoordY")])
  
  km <- kmeans(
    coords,
    centers = K,
    nstart  = 50,
    iter.max = 100
  )
  
  as.factor(km$cluster)
}

data4fit$SCV_Stability_Fold <- make_or_read_spatial_folds(
  data4fit,
  K = K_fold,
  fold_column = fold_column
)

fold_ids <- levels(data4fit$SCV_Stability_Fold)
K_actual <- length(fold_ids)

cat("\n空间折叠数量:", K_actual, "\n")
print(table(data4fit$SCV_Stability_Fold))

# ========================================
# 7. 全数据筛选：作为参考最终模型变量组合
# ========================================

cat("\n开始全数据变量筛选，作为参考变量组合...\n")

full_selection <- run_stepwise_selection(
  data       = data4fit,
  fold_label = "Full_Data",
  verbose    = TRUE
)

if (is.null(full_selection)) {
  stop("全数据变量筛选失败。")
}

full_terms <- full_selection$selected_terms

cat("\n全数据最终保留项:\n")
print(full_terms)

# ========================================
# 8. 每个 SCV 训练折内重复变量筛选
# ========================================

fold_results <- list()

for (fold in fold_ids) {
  
  cat("\n", paste(rep("#", 80), collapse = ""), "\n", sep = "")
  cat(sprintf("稳定性分析：当前 validation fold = %s\n", fold))
  cat("仅使用其余 folds 作为训练集重新执行变量筛选。\n")
  cat(paste(rep("#", 80), collapse = ""), "\n", sep = "")
  
  train_data <- data4fit[data4fit$SCV_Stability_Fold != fold, ]
  
  res <- run_stepwise_selection(
    data       = train_data,
    fold_label = paste0("Train_without_Fold_", fold),
    verbose    = TRUE
  )
  
  fold_results[[as.character(fold)]] <- res
}

# 移除失败的 fold
valid_fold_names <- names(fold_results)[
  sapply(fold_results, function(x) !is.null(x))
]

fold_results <- fold_results[valid_fold_names]

if (length(fold_results) == 0) {
  stop("所有 folds 的变量筛选均失败，无法进行稳定性分析。")
}

# ========================================
# 9. 整理稳定性矩阵
# ========================================

fold_term_list <- lapply(fold_results, function(x) x$selected_terms)

all_terms <- sort(unique(c(full_terms, unlist(fold_term_list))))

if (length(all_terms) == 0) {
  stop("没有任何变量被选中，无法进行稳定性分析。")
}

stability_matrix <- matrix(
  0,
  nrow = length(fold_term_list),
  ncol = length(all_terms)
)

rownames(stability_matrix) <- paste0("Train_without_Fold_", names(fold_term_list))
colnames(stability_matrix) <- all_terms

for (i in seq_along(fold_term_list)) {
  selected <- fold_term_list[[i]]
  stability_matrix[i, colnames(stability_matrix) %in% selected] <- 1
}

# 选择频率
selection_count <- colSums(stability_matrix)
selection_freq  <- selection_count / nrow(stability_matrix)

frequency_df <- data.frame(
  Term = all_terms,
  Type = sapply(all_terms, term_type),
  Full_Data_Selected = all_terms %in% full_terms,
  Selected_Folds = as.integer(selection_count),
  Total_Folds = nrow(stability_matrix),
  Selection_Frequency = selection_freq,
  Selection_Frequency_Percent = selection_freq * 100,
  Stable_80 = selection_freq >= stability_threshold,
  stringsAsFactors = FALSE
)

frequency_df <- frequency_df[
  order(
    -frequency_df$Selection_Frequency,
    frequency_df$Type,
    frequency_df$Term
  ),
]

# 长格式 fold-selected terms
fold_selected_long <- do.call(
  rbind,
  lapply(names(fold_term_list), function(fold) {
    terms <- fold_term_list[[fold]]
    
    if (length(terms) == 0) {
      return(data.frame(
        Fold = paste0("Train_without_Fold_", fold),
        Term = NA_character_,
        Type = NA_character_,
        stringsAsFactors = FALSE
      ))
    }
    
    data.frame(
      Fold = paste0("Train_without_Fold_", fold),
      Term = terms,
      Type = sapply(terms, term_type),
      stringsAsFactors = FALSE
    )
  })
)

# ========================================
# 10. Jaccard 相似度矩阵
# ========================================

calc_jaccard <- function(a, b) {
  union_ab <- union(a, b)
  
  if (length(union_ab) == 0) {
    return(NA_real_)
  }
  
  length(intersect(a, b)) / length(union_ab)
}

jaccard_matrix <- matrix(
  NA_real_,
  nrow = length(fold_term_list),
  ncol = length(fold_term_list)
)

rownames(jaccard_matrix) <- paste0("Train_without_Fold_", names(fold_term_list))
colnames(jaccard_matrix) <- paste0("Train_without_Fold_", names(fold_term_list))

for (i in seq_along(fold_term_list)) {
  for (j in seq_along(fold_term_list)) {
    jaccard_matrix[i, j] <- calc_jaccard(
      fold_term_list[[i]],
      fold_term_list[[j]]
    )
  }
}

# ========================================
# 11. 保存稳定性分析表格
# ========================================

write.csv(
  frequency_df,
  file = file.path(OUT_DIR, "stability_selection_frequency.csv"),
  row.names = FALSE
)

write.csv(
  stability_matrix,
  file = file.path(OUT_DIR, "stability_fold_term_matrix.csv"),
  row.names = TRUE
)

write.csv(
  fold_selected_long,
  file = file.path(OUT_DIR, "stability_fold_selected_terms_long.csv"),
  row.names = FALSE
)

write.csv(
  jaccard_matrix,
  file = file.path(OUT_DIR, "stability_jaccard_matrix.csv"),
  row.names = TRUE
)

# 保存每个 fold 的 pruning diagnostics
pruning_diagnostics_all <- do.call(
  rbind,
  lapply(names(fold_results), function(fold) {
    diag <- fold_results[[fold]]$pruning_diagnostics
    
    if (is.null(diag) || nrow(diag) == 0) {
      return(NULL)
    }
    
    diag$Fold <- paste0("Train_without_Fold_", fold)
    diag
  })
)

if (!is.null(pruning_diagnostics_all) && nrow(pruning_diagnostics_all) > 0) {
  write.csv(
    pruning_diagnostics_all,
    file = file.path(OUT_DIR, "stability_pruning_diagnostics.csv"),
    row.names = FALSE
  )
}

# ========================================
# 12. 稳定性图 1：选择频率柱状图
# ========================================

frequency_plot_df <- frequency_df
frequency_plot_df$Term <- factor(
  frequency_plot_df$Term,
  levels = rev(frequency_plot_df$Term)
)

p_freq <- ggplot(
  frequency_plot_df,
  aes(x = Term, y = Selection_Frequency_Percent)
) +
  geom_col() +
  coord_flip() +
  geom_hline(
    yintercept = stability_threshold * 100,
    linetype = "dashed",
    linewidth = 0.6
  ) +
  labs(
    title = "Selection Stability of GAM Terms Across Spatial CV Training Folds",
    subtitle = sprintf(
      "Dashed line indicates %.0f%% selection-frequency threshold",
      stability_threshold * 100
    ),
    x = "Selected terms",
    y = "Selection frequency (%)"
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold"),
    axis.text.y = element_text(size = 8)
  )

ggsave(
  filename = file.path(OUT_DIR, "stability_selection_frequency.pdf"),
  plot = p_freq,
  width = 12,
  height = max(6, 0.25 * nrow(frequency_plot_df)),
  units = "in"
)

ggsave(
  filename = file.path(OUT_DIR, "stability_selection_frequency.png"),
  plot = p_freq,
  width = 12,
  height = max(6, 0.25 * nrow(frequency_plot_df)),
  units = "in",
  dpi = 300
)

# ========================================
# 13. 稳定性图 2：fold × term 热图
# ========================================

stability_long <- as.data.frame(as.table(stability_matrix))
colnames(stability_long) <- c("Fold", "Term", "Selected")

term_order <- frequency_df$Term

stability_long$Term <- factor(
  stability_long$Term,
  levels = rev(term_order)
)

p_heat <- ggplot(
  stability_long,
  aes(x = Fold, y = Term, fill = factor(Selected))
) +
  geom_tile(color = "grey80", linewidth = 0.2) +
  scale_fill_manual(
    values = c("0" = "white", "1" = "black"),
    name = "Selected",
    labels = c("No", "Yes")
  ) +
  labs(
    title = "Fold-wise Stability Matrix of Selected GAM Terms",
    x = "Spatial CV training fold",
    y = "Terms"
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.y = element_text(size = 8)
  )

ggsave(
  filename = file.path(OUT_DIR, "stability_fold_term_heatmap.pdf"),
  plot = p_heat,
  width = 10,
  height = max(6, 0.25 * length(all_terms)),
  units = "in"
)

ggsave(
  filename = file.path(OUT_DIR, "stability_fold_term_heatmap.png"),
  plot = p_heat,
  width = 10,
  height = max(6, 0.25 * length(all_terms)),
  units = "in",
  dpi = 300
)

# ========================================
# 14. 稳定性图 3：Jaccard 相似度热图
# ========================================

jaccard_long <- as.data.frame(as.table(jaccard_matrix))
colnames(jaccard_long) <- c("Fold_1", "Fold_2", "Jaccard")

p_jaccard <- ggplot(
  jaccard_long,
  aes(x = Fold_1, y = Fold_2, fill = Jaccard)
) +
  geom_tile(color = "grey80", linewidth = 0.2) +
  geom_text(aes(label = sprintf("%.2f", Jaccard)), size = 3) +
  scale_fill_gradient(
    low = "white",
    high = "black",
    limits = c(0, 1),
    name = "Jaccard"
  ) +
  labs(
    title = "Jaccard Similarity of Selected Term Sets Across Spatial CV Folds",
    x = "Training fold",
    y = "Training fold"
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

ggsave(
  filename = file.path(OUT_DIR, "stability_jaccard_heatmap.pdf"),
  plot = p_jaccard,
  width = 8,
  height = 7,
  units = "in"
)

ggsave(
  filename = file.path(OUT_DIR, "stability_jaccard_heatmap.png"),
  plot = p_jaccard,
  width = 8,
  height = 7,
  units = "in",
  dpi = 300
)

# ========================================
# 15. 控制台仅输出稳定性分析摘要
# ========================================

cat("\n", paste(rep("=", 80), collapse = ""), "\n", sep = "")
cat("变量选择稳定性分析完成\n")
cat(paste(rep("=", 80), collapse = ""), "\n", sep = "")

cat("\n稳定性频率表：\n")
print(frequency_df)

cat("\nFold 间 Jaccard 相似度矩阵：\n")
print(round(jaccard_matrix, 3))

cat("\n稳定项，选择频率 >= ", stability_threshold * 100, "%：\n", sep = "")
print(frequency_df[frequency_df$Stable_80, ])

cat("\n输出文件位置：\n")
cat(normalizePath(file.path(getwd(), OUT_DIR)), "\n")

cat("\n已输出文件：\n")
print(list.files(OUT_DIR, full.names = FALSE))

cat(paste(rep("=", 80), collapse = ""), "\n", sep = "")

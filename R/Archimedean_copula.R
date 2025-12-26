#' Create a simdata design using an Archimedean copula as joint generator
#'
#' @param copula An object inheriting from class 'copula' (e.g., claytonCopula, gumbelCopula, frankCopula, joeCopula).
#' @param dist List of marginal quantile functions. Each function maps u in (0,1) to x.
#' @param names_final Optional character vector of final variable names.
#' @param process_final Optional simdata post-processing list passed to simdata::simdesign.
#' @param name Optional design name.
#' @param eps Numeric in (0,1). Clipping applied to U to avoid 0/1 in quantiles.
#' @param ... Further arguments stored in the simdata design object.
#'
#' @return A 'simdesign' object usable with simdata::simulate_data().
#' @export
simdesign_archimedean_copula <- function( # 用 Archimedean copula 作为“联合分布/依赖结构”的 simdata 设计入口（本质复用通用 copula 设计器）
  copula, # copula 对象（例如 claytonCopula / gumbelCopula / frankCopula / joeCopula），决定依赖结构/尾部依赖形态
  dist, # 边缘分布分位数函数列表：dist[[j]](u) 把 u∈(0,1) 映射到第 j 个变量的实际尺度
  names_final = NULL, # 最终输出列名（可选）
  process_final = list(), # simdata 的后处理步骤（可选），在 transform_initial 之后执行
  name = "Archimedean copula design", # 该设计对象的名称
  eps = 1e-6, # 对 U 做裁剪的数值稳定参数，避免 0/1 进入分位数函数导致 ±Inf
  ... # 其余参数透传到 simdata::simdesign（通过 simdesign_elliptical_copula 统一封装）
) { # 函数体开始
  simdesign_elliptical_copula( # 直接复用通用 copula→U→边缘分位数映射的实现（此处不区分椭圆/阿基米德，只要是 copula 类即可）
    copula = copula, # 依赖结构
    dist = dist, # 边缘分布
    names_final = names_final, # 最终列名
    process_final = process_final, # 最终后处理
    name = name, # 设计名称
    eps = eps, # U 裁剪参数
    ... # 继续透传其他元数据/参数
  ) # 返回 simdesign
} # simdesign_archimedean_copula 结束

#' Convenience constructor for Clayton copula simdata designs
#'
#' @param dist List of marginal quantile functions.
#' @param theta Copula parameter.
#' @param dim Dimension (>= 2).
#' @param names_final Optional names for final variables.
#' @param ... Passed to simdesign_archimedean_copula().
#'
#' @return A 'simdesign' object.
#' @export
simdesign_clayton_copula <- function( # 便捷构造器：用 Clayton copula（常见下尾依赖）创建 simdata 设计
  dist, # 边缘分位数函数列表
  theta, # Clayton copula 参数（控制依赖强度；具体范围由 copula 包定义）
  dim, # 维度 d（变量个数，>=2）
  names_final = NULL, # 最终列名（可选）
  ... # 透传给 simdesign_archimedean_copula
) { # 函数体开始
  if (is.null(dim) || !is.numeric(dim) || length(dim) != 1 || dim < 2) { # 检查 dim 合法性
    stop("dim must be a single numeric >= 2.", call. = FALSE) # 不合法则报错
  }
  cop <- copula::claytonCopula(param = theta, dim = as.integer(dim)) # 构造 Clayton copula 对象（依赖结构）

  simdesign_archimedean_copula( # 用通用 Archimedean 入口封装为 simdesign
    copula = cop, # Clayton copula
    dist = dist, # 边缘分布
    names_final = names_final, # 最终列名
    ... # 其余参数透传
  ) # 返回 simdesign
} # simdesign_clayton_copula 结束

#' Convenience constructor for Gumbel copula simdata designs
#'
#' @param dist List of marginal quantile functions.
#' @param theta Copula parameter.
#' @param dim Dimension (>= 2).
#' @param names_final Optional names for final variables.
#' @param ... Passed to simdesign_archimedean_copula().
#'
#' @return A 'simdesign' object.
#' @export
simdesign_gumbel_copula <- function( # 便捷构造器：用 Gumbel copula（常见上尾依赖）创建 simdata 设计
  dist, # 边缘分位数函数列表
  theta, # Gumbel copula 参数（控制依赖强度；通常 theta>=1）
  dim, # 维度 d（变量个数，>=2）
  names_final = NULL, # 最终列名（可选）
  ... # 透传给 simdesign_archimedean_copula
) { # 函数体开始
  if (is.null(dim) || !is.numeric(dim) || length(dim) != 1 || dim < 2) { # 检查 dim 合法性
    stop("dim must be a single numeric >= 2.", call. = FALSE) # 不合法则报错
  }
  cop <- copula::gumbelCopula(param = theta, dim = as.integer(dim)) # 构造 Gumbel copula 对象

  simdesign_archimedean_copula( # 封装为 simdesign
    copula = cop, # Gumbel copula
    dist = dist, # 边缘分布
    names_final = names_final, # 最终列名
    ... # 其余参数透传
  ) # 返回 simdesign
} # simdesign_gumbel_copula 结束

#' Convenience constructor for Frank copula simdata designs
#'
#' @param dist List of marginal quantile functions.
#' @param theta Copula parameter (cannot be 0).
#' @param dim Dimension (>= 2).
#' @param names_final Optional names for final variables.
#' @param ... Passed to simdesign_archimedean_copula().
#'
#' @return A 'simdesign' object.
#' @export
simdesign_frank_copula <- function( # 便捷构造器：用 Frank copula（对称依赖、无尾依赖）创建 simdata 设计
  dist, # 边缘分位数函数列表
  theta, # Frank copula 参数（theta=0 退化为独立；具体范围由 copula 包定义）
  dim, # 维度 d（变量个数，>=2）
  names_final = NULL, # 最终列名（可选）
  ... # 透传给 simdesign_archimedean_copula
) { # 函数体开始
  if (is.null(dim) || !is.numeric(dim) || length(dim) != 1 || dim < 2) { # 检查 dim 合法性
    stop("dim must be a single numeric >= 2.", call. = FALSE) # 不合法则报错
  }
  cop <- copula::frankCopula(param = theta, dim = as.integer(dim)) # 构造 Frank copula 对象

  simdesign_archimedean_copula( # 封装为 simdesign
    copula = cop, # Frank copula
    dist = dist, # 边缘分布
    names_final = names_final, # 最终列名
    ... # 其余参数透传
  ) # 返回 simdesign
} # simdesign_frank_copula 结束

#' Convenience constructor for Joe copula simdata designs
#'
#' @param dist List of marginal quantile functions.
#' @param theta Copula parameter.
#' @param dim Dimension (>= 2).
#' @param names_final Optional names for final variables.
#' @param ... Passed to simdesign_archimedean_copula().
#'
#' @return A 'simdesign' object.
#' @export
simdesign_joe_copula <- function( # 便捷构造器：用 Joe copula（常见上尾依赖且更“尖锐”）创建 simdata 设计
  dist, # 边缘分位数函数列表
  theta, # Joe copula 参数（控制依赖强度；具体范围由 copula 包定义）
  dim, # 维度 d（变量个数，>=2）
  names_final = NULL, # 最终列名（可选）
  ... # 透传给 simdesign_archimedean_copula
) { # 函数体开始
  if (is.null(dim) || !is.numeric(dim) || length(dim) != 1 || dim < 2) { # 检查 dim 合法性
    stop("dim must be a single numeric >= 2.", call. = FALSE) # 不合法则报错
  }
  cop <- copula::joeCopula(param = theta, dim = as.integer(dim)) # 构造 Joe copula 对象

  simdesign_archimedean_copula( # 封装为 simdesign
    copula = cop, # Joe copula
    dist = dist, # 边缘分布
    names_final = names_final, # 最终列名
    ... # 其余参数透传
  ) # 返回 simdesign
} # simdesign_joe_copula 结束

#' Fit an Archimedean copula design from data (empirical margins)
#'
#' @param data data.frame containing variables.
#' @param vars character vector of variable names to model.
#' @param family One of "clayton", "gumbel", "frank", "joe".
#' @param fit_method "itau" or "itau_mpl". "itau_mpl" tries MPL after ITAU init.
#' @param qtype Quantile type for empirical margins.
#' @param eps Clipping for pseudo-observations and simulation U.
#' @param name Optional simdesign name.
#' @param ... Stored in the returned simdesign (e.g., for metadata).
#'
#' @return A 'simdesign' object usable with simdata::simulate_data().
#' @export
simdesign_archimedean_copula_from_data <- function( # 从数据拟合 Archimedean copula 的依赖结构，并配合经验边缘生成可模拟的 simdesign
  data, # 原始数据（data.frame）
  vars, # 需要建模的变量名向量（长度 >= 2）
  family = c("clayton", "gumbel", "frank", "joe"), # 指定 Archimedean 家族（不同家族对应不同尾依赖形态）
  fit_method = c("itau", "itau_mpl"), # 拟合方法：itau（稳健）或 itau 初始化后尝试 mpl（可能失败则回退）
  qtype = 8, # 经验分位数算法（stats::quantile 的 type）
  eps = 1e-6, # 数值稳定裁剪：用于 pobs 与后续模拟的 U
  name = "Archimedean copula design (fit from data)", # 设计对象名称
  ... # 其余参数作为元数据保存在 simdesign 中
) { # 函数体开始
  family <- match.arg(family) # 规范化 family
  fit_method <- match.arg(fit_method) # 规范化 fit_method

  if (!is.data.frame(data)) { # 检查 data 类型
    stop("data must be a data.frame.", call. = FALSE) # 不是 data.frame 则报错
  }
  if (!is.character(vars) || length(vars) < 2) { # 检查 vars 合法性
    stop("vars must be a character vector of length >= 2.", call. = FALSE) # 必须是长度>=2的字符向量
  }
  if (!all(vars %in% names(data))) { # 检查 vars 是否都存在于 data
    miss <- setdiff(vars, names(data)) # 找出缺失列
    stop("vars missing from data: ", paste(miss, collapse = ", "), call. = FALSE) # 报错提示缺失列名
  }

  df_in <- data[, vars, drop = FALSE] # 提取待建模变量子集（保持 data.frame）
  for (nm in vars) { # MVP：仅支持数值型变量
    if (!is.numeric(df_in[[nm]])) { # 若某列非数值
      stop("All vars must be numeric for MVP. Non-numeric: ", nm, call. = FALSE) # 报错并指出变量名
    }
  }
  xmat <- as.matrix(df_in) # 转为矩阵便于后续计算
  if (!all(is.finite(xmat))) { # 检查是否含 NA/NaN/Inf
    stop("Data contains non-finite values in selected vars.", call. = FALSE) # 若存在非有限值则报错
  }

  dist <- lapply(df_in, .make_empirical_q, qtype = qtype, eps = eps) # 为每列构造经验边缘分位数函数（模拟时保持原边缘形状）

  u_hat <- copula::pobs(xmat) # 将样本转换为伪观测 U（每列近似 U(0,1)，用于拟合 copula）
  u_hat <- .clip_unit(u_hat, eps = eps) # 裁剪到 (eps,1-eps) 避免 0/1 边界影响拟合

  dim <- ncol(u_hat) # 变量维度 d

  cop0 <- switch( # 按家族选择一个“初始 copula 对象”（用于 fitCopula 初始化）
    family, # 分支选择器
    clayton = copula::claytonCopula(param = 0.5, dim = dim), # Clayton 初值（示例：中等依赖强度）
    gumbel  = copula::gumbelCopula(param = 1.1, dim = dim), # Gumbel 初值（略强于独立：theta=1 为独立）
    frank   = copula::frankCopula(param = 1, dim = dim), # Frank 初值（theta=0 为独立，取 1 作为温和初值）
    joe     = copula::joeCopula(param = 1.1, dim = dim), # Joe 初值（theta=1 为独立，取 1.1 作为温和初值）
    stop("Unsupported family.", call. = FALSE) # 不支持的家族直接报错
  )

  fit_itau <- copula::fitCopula(cop0, data = u_hat, method = "itau") # 使用 ITAU（基于 Kendall's tau 的反演）进行稳健拟合/初始化
  fit_final <- fit_itau # 默认将 itau 结果作为最终拟合

  if (fit_method == "itau_mpl") { # 若选择 itau 后再尝试 MPL（最大伪似然）
    fit_final <- tryCatch( # MPL 可能在某些数据/初值下不稳定；失败则回退
      copula::fitCopula(fit_itau@copula, data = u_hat, method = "mpl"), # 以 itau 拟合的 copula 作为 MPL 起点
      error = function(e) fit_itau # 若 MPL 抛错则回退到 itau
    )
  }

  simdesign_archimedean_copula( # 用“拟合后的依赖结构 + 经验边缘”封装为可模拟的 simdesign
    copula = fit_final@copula, # 拟合后的 copula 对象（含估计参数）
    dist = dist, # 经验边缘分位数函数列表
    names_final = vars, # 输出列名与建模变量名一致
    name = name, # 设计名称
    eps = eps, # U 裁剪参数
    copula_fit = fit_final, # 附带保存拟合对象，便于用户检查参数/拟合信息
    copula_family = family, # 元数据：copula 家族
    margins = "empirical", # 元数据：边缘采用经验分布
    ... # 其余元数据透传
  ) # 返回 simdesign
} # simdesign_archimedean_copula_from_data 结束
.make_empirical_q <- function(x, qtype = 8, eps = 1e-6) { # 构造“经验分位数函数”：把 u∈(0,1) 映射到样本 x 的分位数
  force(x) # 强制捕获 x，避免闭包里 x 发生惰性求值/被外部同名变量影响
  function(u) { # 返回一个函数：输入 u（概率），输出对应分位数
    u <- pmin(pmax(u, eps), 1 - eps) # 将 u 截断到 (eps, 1-eps)，避免 quantile 在 0/1 处产生边界问题
    as.numeric(stats::quantile(x, probs = u, type = qtype, names = FALSE)) # 计算经验分位数并转为纯数值向量
  } # 内部函数结束
} # .make_empirical_q 结束

.assert_is_function_list <- function(x, expected_len, name) { # 断言：x 必须是指定长度的“函数列表”（用于边缘分布 dist）
  if (!is.list(x) || length(x) != expected_len) { # 检查是否为 list 且长度符合预期
    stop(name, " must be a list of length ", expected_len, ".", call. = FALSE) # 不满足则抛错，并关闭调用栈提示
  }
  ok <- vapply(x, is.function, logical(1)) # 逐个元素检查是否为函数，返回逻辑向量
  if (!all(ok)) { # 任意一个不是函数就失败
    stop(name, " must be a list of functions.", call. = FALSE) # 抛错提示必须是函数列表
  }
  invisible(TRUE) # 通过校验：返回不可见 TRUE，便于在管道/函数内部使用
} # .assert_is_function_list 结束

.clip_unit <- function(u, eps = 1e-6) { # 将数值截断到 (eps, 1-eps)，用于 pseudo-observations 或 U(0,1) 模拟值
  pmin(pmax(u, eps), 1 - eps) # 先下界截断再上界截断（向量化）
} # .clip_unit 结束

.get_copula_dim <- function(copula) { # 获取 copula 对象的维度（替代 copula::dim 非导出导致的问题）
  d <- tryCatch(base::dim(copula), error = function(e) NULL) # 优先尝试 base::dim（若类实现了 dim 方法会返回维度）
  if (!is.null(d) && length(d) == 1 && is.finite(d)) { # 若成功拿到单个有限数值
    return(as.integer(d)) # 转为整数并返回
  }

  d <- tryCatch(methods::slot(copula, "dimension"), error = function(e) NULL) # 再尝试从 S4 slot "dimension" 中读取
  if (!is.null(d) && length(d) == 1 && is.finite(d)) { # 同样校验其有效性
    return(as.integer(d)) # 转为整数并返回
  }

  stop("Unable to determine copula dimension.", call. = FALSE) # 两种方式都失败：终止并提示无法确定维度
} # .get_copula_dim 结束

#' Create a simdata design using an elliptical copula as joint generator
#'
#' This is the core building block: a copula generates dependent U(0,1)
#' samples, then marginal quantile functions map each column to the final scale.
#'
#' @param copula An object inheriting from class 'copula' (e.g., normalCopula, tCopula).
#' @param dist List of marginal quantile functions. Each function maps u in (0,1) to x.
#' @param names_final Optional character vector of final variable names.
#' @param process_final Optional simdata post-processing list passed to simdata::simdesign.
#' @param name Optional design name.
#' @param eps Numeric in (0,1). Clipping applied to U to avoid 0/1 in quantiles.
#' @param ... Further arguments stored in the simdata design object.
#'
#' @return A 'simdesign' object usable with simdata::simulate_data().
#' @export
simdesign_elliptical_copula <- function( # 使用“椭圆族 copula”构造一个 simdata 设计对象：copula 生成相关的 U(0,1)，再用边缘分位数函数映射到目标尺度
  copula, # copula 对象（例如 normalCopula / tCopula），决定变量间依赖结构
  dist, # 边缘分布的分位数函数列表：dist[[j]](u) 把 u∈(0,1) 映射到第 j 个变量
  names_final = NULL, # 最终输出数据列名（可选）
  process_final = list(), # simdata 的后处理步骤（可选），在 transform_initial 后执行
  name = "Elliptical copula design", # 该设计对象的名称
  eps = 1e-6, # 对 U 做裁剪的数值稳定参数，避免 0/1 进入分位数函数
  ... # 其余参数透传给 simdata::simdesign（例如元数据等）
) { # 函数体开始
  if (!inherits(copula, "copula")) { # 校验 copula 必须继承自类 "copula"
    stop("copula must inherit from class 'copula'.", call. = FALSE) # 不满足则直接报错
  }

  dim <- .get_copula_dim(copula) # 从 copula 对象推断维度 d（变量个数）
  .assert_is_function_list(dist, dim, "dist") # 校验 dist 是长度为 d 的函数列表

  generator <- function(n_obs, ...) { # simdata 的 generator：只负责产生“初始态”的随机数（这里是 U 矩阵）
    u <- copula::rCopula(n_obs, copula) # 从 copula 抽样得到 n_obs×d 的相关 U(0,1)
    u <- .clip_unit(u, eps = eps) # 将 U 裁剪到 (eps, 1-eps)，避免极端值导致分位数为 ±Inf
    u # 返回 U（作为 transform_initial 的输入）
  } # generator 结束

  transform_initial <- function(u) { # simdata 的 transform_initial：把 U 映射到最终变量尺度
    if (!is.matrix(u) && !is.data.frame(u)) { # 防御式检查：generator 应该返回二维对象
      stop("Internal error: generator did not return a 2D object.", call. = FALSE) # 若不是二维则报内部错误
    }
    u <- as.matrix(u) # 统一转成矩阵，便于按列处理
    if (ncol(u) != dim) { # 检查列数是否与 copula 维度一致
      stop("Internal error: U has unexpected number of columns.", call. = FALSE) # 不一致则报内部错误
    }

    x <- lapply(seq_len(dim), function(j) dist[[j]](u[, j])) # 对每一列 U 应用对应边缘分位数函数，得到各变量的样本
    x <- as.data.frame(x, optional = TRUE, stringsAsFactors = FALSE) # 组装为 data.frame（不强制设定因子）

    x # 返回映射后的数据框（后续再由 simdata 负责命名与 post-process）
  } # transform_initial 结束

  simdata::simdesign( # 构建并返回 simdata 设计对象
    generator = generator, # 初始随机生成器
    transform_initial = transform_initial, # 初始到最终变量的变换
    names_final = names_final, # 最终列名
    process_final = process_final, # 最终后处理流程
    name = name, # 设计对象名称
    ... # 其余参数继续透传
  ) # 返回 simdesign
} # simdesign_elliptical_copula 结束

#' Convenience constructor for Gaussian copula simdata designs
#'
#' @param dist List of marginal quantile functions.
#' @param rho Scalar correlation used when structure = "ex".
#' @param Sigma Correlation matrix used when structure = "un".
#' @param dim Dimension (required for structure = "ex"; inferred from Sigma for "un").
#' @param structure Correlation structure: "ex" (exchangeable) or "un" (unstructured).
#' @param names_final Optional names for final variables.
#' @param ... Passed to simdesign_elliptical_copula().
#'
#' @return A 'simdesign' object.
#' @export
simdesign_gaussian_copula <- function( # 便捷构造器：用“高斯 copula（normalCopula）”快速生成 simdata 设计
  dist, # 边缘分位数函数列表（长度等于维度）
  rho = 0, # 当 structure = "ex"（交换相关）时使用的单一相关系数
  Sigma = NULL, # 当 structure = "un"（非结构化相关）时使用的相关矩阵
  dim = NULL, # 当 structure = "ex" 时需要显式给出维度；structure = "un" 时从 Sigma 推断
  structure = c("ex", "un"), # 相关结构："ex"=exchangeable（同相关）；"un"=unstructured（任意相关）
  names_final = NULL, # 最终列名（可选）
  ... # 透传给 simdesign_elliptical_copula
) { # 函数体开始
  structure <- match.arg(structure) # 将 structure 规范化为允许值之一

  if (structure == "ex") { # 交换相关结构：只用一个 rho 描述所有非对角相关
    if (is.null(dim) || !is.numeric(dim) || length(dim) != 1 || dim < 2) { # 检查 dim 合法性
      stop("For structure='ex', dim must be a single numeric >= 2.", call. = FALSE) # dim 不合法则报错
    }
    cop <- copula::normalCopula(param = rho, dim = as.integer(dim), dispstr = "ex") # 构造 exchangeable 的高斯 copula
  } else { # 非结构化相关：由完整相关矩阵 Sigma 决定
    if (is.null(Sigma) || !is.matrix(Sigma) || nrow(Sigma) != ncol(Sigma)) { # 检查 Sigma 是否为方阵
      stop("For structure='un', Sigma must be a square matrix.", call. = FALSE) # 不是方阵则报错
    }
    dim <- ncol(Sigma) # 维度由 Sigma 的列数决定
    if (dim < 2) { # 维度至少要 2
      stop("Sigma must have dimension >= 2.", call. = FALSE) # 否则报错
    }
    if (!all(is.finite(Sigma))) { # 检查相关矩阵是否包含非有限值
      stop("Sigma must be finite.", call. = FALSE) # 否则报错
    }
    if (max(abs(diag(Sigma) - 1)) > 1e-10) { # 检查对角线是否为 1（相关矩阵要求）
      stop("Sigma must be a correlation matrix with diag = 1.", call. = FALSE) # 不满足则报错
    }
    cop <- copula::normalCopula(param = copula::P2p(Sigma), dim = dim, dispstr = "un") # 将 Sigma 转为参数向量并构造 unstructured 高斯 copula
  }

  simdesign_elliptical_copula( # 复用通用椭圆 copula 设计构造函数
    copula = cop, # 刚构造的高斯 copula
    dist = dist, # 边缘分布
    names_final = names_final, # 最终列名
    ... # 其余参数透传
  ) # 返回 simdesign
} # simdesign_gaussian_copula 结束

#' Fit an elliptical copula design from data (empirical margins)
#'
#' MVP supports Gaussian copula + exchangeable/unstructured structure with robust defaults.
#'
#' @param data data.frame containing variables.
#' @param vars character vector of variable names to model.
#' @param family "gaussian" or "t" (t currently requires df provided).
#' @param structure "ex" or "un".
#' @param fit_method "itau" or "itau_mpl". "itau_mpl" tries MPL after ITAU init.
#' @param df Degrees of freedom if family="t".
#' @param qtype Quantile type for empirical margins.
#' @param eps Clipping for pseudo-observations and simulation U.
#' @param name Optional simdesign name.
#' @param ... Stored in the returned simdesign (e.g., for metadata).
#'
#' @return A 'simdesign' object usable with simdata::simulate_data().
#' @export
simdesign_elliptical_copula_from_data <- function( # 从已有数据拟合椭圆 copula（依赖结构）并用经验边缘（分位数函数）构造 simdata 设计
  data, # 原始数据（data.frame）
  vars, # 需要建模的变量名向量（长度 >= 2）
  family = c("gaussian", "t"), # copula 家族：高斯或 t；t 族需要提供 df
  structure = c("ex", "un"), # 相关结构：交换相关 or 非结构化
  fit_method = c("itau", "itau_mpl"), # 拟合方法：仅 itau（稳健）或 itau 初始化后尝试 mpl
  df = NULL, # t copula 的自由度（family="t" 时必填且 > 2）
  qtype = 8, # 经验分位数的 type（stats::quantile 的算法选择）
  eps = 1e-6, # 数值稳定裁剪：用于 pobs 与后续模拟的 U
  name = "Elliptical copula design (fit from data)", # 设计对象名称
  ... # 其余参数作为元数据保存在 simdesign 中
) { # 函数体开始
  family <- match.arg(family) # 规范化 family
  structure <- match.arg(structure) # 规范化 structure
  fit_method <- match.arg(fit_method) # 规范化 fit_method

  if (!is.data.frame(data)) { # 检查输入 data 类型
    stop("data must be a data.frame.", call. = FALSE) # 不是 data.frame 则报错
  }
  if (!is.character(vars) || length(vars) < 2) { # 检查 vars 合法性
    stop("vars must be a character vector of length >= 2.", call. = FALSE) # vars 必须是长度>=2的字符向量
  }
  if (!all(vars %in% names(data))) { # 检查 vars 是否都在 data 中
    miss <- setdiff(vars, names(data)) # 找出缺失的列名
    stop("vars missing from data: ", paste(miss, collapse = ", "), call. = FALSE) # 报错提示缺失列
  }

  df_in <- data[, vars, drop = FALSE] # 只取出待建模变量子集（保持 data.frame）
  for (nm in vars) { # 逐个变量做最小可行（MVP）检查
    if (!is.numeric(df_in[[nm]])) { # 当前实现仅支持数值型变量（便于 pobs 与相关计算）
      stop("All vars must be numeric for MVP. Non-numeric: ", nm, call. = FALSE) # 非数值则报错并指出变量名
    }
  }
  xmat <- as.matrix(df_in) # 转成数值矩阵便于后续计算
  if (!all(is.finite(xmat))) { # 检查是否包含 NA/NaN/Inf
    stop("Data contains non-finite values in selected vars.", call. = FALSE) # 存在非有限值则报错
  }

  dist <- lapply(df_in, .make_empirical_q, qtype = qtype, eps = eps) # 为每个变量构造经验边缘分位数函数（用于模拟时还原到原尺度）

  u_hat <- copula::pobs(xmat) # 将样本转换为伪观测（pseudo-observations），近似得到每列的 U(0,1)
  u_hat <- .clip_unit(u_hat, eps = eps) # 裁剪到 (eps,1-eps)，避免边界 0/1 影响拟合与后续分位数

  dim <- ncol(u_hat) # copula 维度 d（变量个数）

  if (family == "gaussian") { # 拟合高斯 copula
    if (structure == "ex") { # 交换相关：只需要一个参数
      cop0 <- copula::normalCopula(param = 0, dim = dim, dispstr = "ex") # 初始值设为 0 相关
    } else { # 非结构化相关：需要 d(d-1)/2 个参数
      cop0 <- copula::normalCopula(param = rep(0, dim * (dim - 1) / 2), dim = dim, dispstr = "un") # 初始相关全为 0
    }
  } else { # 拟合 t copula（尾部相关更强）
    if (is.null(df) || !is.numeric(df) || length(df) != 1 || df <= 2) { # t copula 的 df 约束
      stop("For family='t', df must be a single numeric > 2.", call. = FALSE) # df 不合法则报错
    }
    if (structure == "ex") { # 交换相关的 t copula
      cop0 <- copula::tCopula(param = 0, dim = dim, dispstr = "ex", df = df, df.fixed = TRUE) # df.fixed=TRUE 表示 df 不参与估计
    } else { # 非结构化相关的 t copula
      cop0 <- copula::tCopula(param = rep(0, dim * (dim - 1) / 2), dim = dim, dispstr = "un", df = df, df.fixed = TRUE) # 相关参数初值为 0
    }
  }

  fit_itau <- copula::fitCopula(cop0, data = u_hat, method = "itau") # 先用 Kendall's tau 的反演（ITAU）做稳健拟合/初始化
  fit_final <- fit_itau # 默认把 itau 结果作为最终结果

  if (fit_method == "itau_mpl") { # 若选择 itau 后再尝试 MPL（最大伪似然）
    fit_final <- tryCatch( # MPL 在某些数据/初值下可能失败；失败则回退到 itau
      copula::fitCopula(fit_itau@copula, data = u_hat, method = "mpl"), # 用 itau 拟合得到的 copula 作为 MPL 的起点
      error = function(e) fit_itau # 若 MPL 抛错，则返回 itau 拟合结果
    )
  }

  simdesign_elliptical_copula( # 用拟合出的依赖结构 + 经验边缘，构造可直接 simulate_data 的 simdesign
    copula = fit_final@copula, # 拟合后的 copula 对象（包含相关参数/结构等）
    dist = dist, # 经验边缘分位数函数列表（保持原数据的边缘分布形状）
    names_final = vars, # 输出列名与建模变量名一致
    name = name, # 设计名称
    eps = eps, # U 裁剪参数（模拟时也会使用）
    copula_fit = fit_final, # 附带保存拟合对象，便于用户检查参数/收敛信息
    copula_family = family, # 元数据：copula 家族
    copula_structure = structure, # 元数据：相关结构
    margins = "empirical", # 元数据：边缘采用经验分布
    ... # 继续透传用户的其他元数据字段
  ) # 返回 simdesign
} # simdesign_elliptical_copula_from_data 结束
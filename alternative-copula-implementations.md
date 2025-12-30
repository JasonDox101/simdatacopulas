# 不依赖 `copula` 包实现 Copula 的方法与可选替代包

本文面向 `simdatacopulas` 的使用与扩展场景：当你不想依赖 R 包 `copula`，仍希望实现“copula 依赖结构 + 任意边缘分布（quantile 映射）”时，可采用哪些数学构造与替代包。

`simdatacopulas` 的核心工作流可以概括为：

1. 生成相关的 \(U\in(0,1)^d\)（依赖结构）
2. 用边缘分位数函数 `dist[[j]](u)` 把每一列 \(U_j\) 映射到目标尺度（边缘分布）

也就是说，只要你能生成一个合理的 `U` 矩阵（并且可控相关/尾部行为），就能完成 copula 级别的联合生成。

---

## 1) 直接用数学构造（Base R 即可）

### 1.1 Gaussian copula（高斯 copula）

定义：

\[
Z \sim N(0,\Sigma),\quad U = \Phi(Z)
\]

其中 \(\Sigma\) 是相关矩阵（对角为 1），\(\Phi\) 是标准正态 CDF。生成步骤：

1. 对 \(\Sigma\) 做 Cholesky 分解：\(\Sigma = LL^\top\)
2. 采样 \(E\sim N(0,I)\)，令 \(Z = E L^\top\)
3. \(U = \Phi(Z)\)

Base R 示例（只依赖 `stats`）：

```r
r_gaussian_copula_base <- function(n, Sigma, eps = 1e-6) {
  if (!is.matrix(Sigma) || nrow(Sigma) != ncol(Sigma)) stop("Sigma must be square.")
  d <- ncol(Sigma)
  if (d < 2) stop("Need d >= 2.")
  if (max(abs(diag(Sigma) - 1)) > 1e-10) stop("Sigma diag must be 1.")

  L <- chol(Sigma)
  E <- matrix(stats::rnorm(n * d), nrow = n, ncol = d)
  Z <- E %*% t(L)
  U <- stats::pnorm(Z)
  pmin(pmax(U, eps), 1 - eps)
}
```

适用性与限制：

- 优点：实现简单、稳定、速度快，常作为多变量模拟的“基线依赖结构”。
- 限制：尾部依赖为“弱尾/对称”，难以表达极端事件同步增强的现象（这通常需要 t copula 或 Archimedean/vine）。

### 1.2 t copula

定义：

\[
T \sim t_\nu(0,\Sigma),\quad U = F_{t,\nu}(T)
\]

一个常用的构造方式是“正态 + 卡方缩放”：

\[
Y\sim N(0,\Sigma),\quad W\sim \chi^2_\nu,\quad T = \frac{Y}{\sqrt{W/\nu}}
\]

Base R 示例：

```r
r_t_copula_base <- function(n, Sigma, df, eps = 1e-6) {
  if (!is.numeric(df) || length(df) != 1 || !is.finite(df) || df <= 2) stop("df must be > 2.")
  if (!is.matrix(Sigma) || nrow(Sigma) != ncol(Sigma)) stop("Sigma must be square.")
  d <- ncol(Sigma)
  if (d < 2) stop("Need d >= 2.")
  if (max(abs(diag(Sigma) - 1)) > 1e-10) stop("Sigma diag must be 1.")

  L <- chol(Sigma)
  E <- matrix(stats::rnorm(n * d), nrow = n, ncol = d)
  Y <- E %*% t(L)
  W <- stats::rchisq(n, df = df)
  T <- Y / sqrt(W / df)
  U <- stats::pt(T, df = df)
  pmin(pmax(U, eps), 1 - eps)
}
```

适用性与限制：

- 优点：比 Gaussian copula 更能表达“尾部相关/极端同步”（df 越小尾越厚）。
- 限制：df 的估计/选择需要额外策略（网格 profile、AIC/BIC、或固定 df）。

### 1.3 直接实现 Archimedean copula（可行但不建议作为首选）

Archimedean copula 的数学形式是：

\[
C(u_1,\dots,u_d)=\varphi^{-1}\left(\varphi(u_1)+\cdots+\varphi(u_d)\right)
\]

其中 \(\varphi\) 为生成元（generator）。从“数学上”确实可以不靠任何 copula 包实现 Clayton/Gumbel/Frank/Joe 等家族，但工程实践上难点在于：

- 每个家族的采样算法不同（常用的是 mixture representation 或特定的采样构造）。
- 参数边界与数值稳定性处理复杂（高维、极端参数、尾部区域更明显）。
- 拟合（ITAU/MPL/MLE）需要更多优化与约束处理。

建议：如果你的目标是稳定模拟与可解释建模，优先使用成熟包（`copula`、`nacopula`、`VineCopula` 等）或采用 vine / 因子模型替代。

---

## 2) 使用其他 R 包实现（不使用 `copula` 包本身）

下面列的是“实现同类能力”的常见方向，是否在你的环境中可用取决于你安装的包与许可约束。

### 2.1 仅用多元分布包实现 Gaussian/t copula

如果你只需要 Gaussian/t copula，通常不需要专用 copula 包：

- `mvtnorm`：提供多元正态/多元 t 的采样与密度等工具，可用来生成 \(Z\) 或 \(T\)，然后做 CDF 映射得到 \(U\)。
- `MASS`：`mvrnorm()` 可采样多元正态（t 仍需你自行做卡方缩放）。

特点：

- 优点：依赖更“底层”，通常比完整 copula 包更轻。
- 限制：只覆盖椭圆族（Gaussian/t），不覆盖 Archimedean、vine 的丰富结构。

### 2.2 Vine copula：`VineCopula` / `rvinecopulib`

vine（pair-copula construction）在高维相关建模中非常常用，特点是用大量二维 copula 拼出高维结构。

- `VineCopula`：经典实现，支持结构选择、拟合、模拟。
- `rvinecopulib`：常见的高性能实现（C++），也提供结构选择与模拟。

特点：

- 优点：高维能力强、家族选择灵活、拟合工具成熟。
- 限制：接口与对象结构与 `copula` 不同，需要单独适配（`simdatacopulas` 已提供 `VineCopula` 入口）。

### 2.3 Archimedean / nested Archimedean：`nacopula`

`nacopula` 主要面向 Archimedean 与嵌套 Archimedean 的建模与采样（具体能力依版本而定）。

特点：

- 优点：在 Archimedean 方向更专门。
- 限制：生态与 API 与 `copula` 不同，拟合策略可能需要重新梳理。

---

## 3) “不用 copula 包但实现同一目标”的替代建模路线

如果你的目标是“生成多终点联合分布”，copula 只是实现路径之一。以下方法不一定以“copula”命名，但能实现同样的两件事：依赖结构 + 边缘控制。

### 3.1 经验 copula / rank resampling（重采样联合秩结构）

做法：

1. 从历史数据得到伪观测 \(U\)（秩变换）
2. 对 \(U\) 行重采样（可选 jitter）
3. 用边缘分位数函数把 \(U\) 映射回目标尺度

特点：

- 优点：不需要选 copula 家族，实施简单，能“复刻”已观测到的依赖。
- 限制：外推能力弱，尤其是尾部结构（极端事件）难以凭空增强。

### 3.2 因子高斯模型（低秩相关结构）

做法：

\[
Z = F\Lambda^\top + E,\quad U=\Phi(\tilde Z)
\]

特点：

- 优点：高维更稳定，参数更少（低秩 + 对角），易解释（公共因子）。
- 限制：依赖形态仍偏“高斯化”，尾部结构不如 vine/Archimedean 灵活。

### 3.3 潜变量阈值模型（离散/有序结局的联合生成）

当变量是二值/多分类/有序等级时，直接对原始变量做 `pobs()` 往往会遇到 ties 多、信息损失与拟合不稳定。潜变量阈值模型提供了更自然的联合生成路径：

1. 在潜在连续空间生成相关的 \(Z\)
2. 用阈值把 \(Z\) 切分成离散类别

特点：

- 优点：对离散/有序终点更贴近生成机制，ties 不是问题。
- 限制：模型假设是“潜变量高斯结构”，不等价于任意离散 copula。

---

## 4) 如何对接到 `simdatacopulas` 的设计模式（只要你能生成 `U`）

`simdatacopulas` 的核心设计是 `generator` 产生 `U`，再用 `transform_initial` 应用边缘分位数函数 `dist`。

即使完全不使用 `copula` 包，你也可以用 `simdata::simdesign()` 包一层：

```r
simdesign_from_u_generator <- function(rU, dist, names_final = NULL, eps = 1e-6, ...) {
  d <- length(dist)
  generator <- function(n_obs, ...) {
    u <- rU(n_obs)
    u <- pmin(pmax(u, eps), 1 - eps)
    u
  }
  transform_initial <- function(u) {
    u <- as.matrix(u)
    x <- lapply(seq_len(d), function(j) dist[[j]](u[, j]))
    as.data.frame(x, optional = TRUE, stringsAsFactors = FALSE)
  }
  simdata::simdesign(
    generator = generator,
    transform_initial = transform_initial,
    names_final = names_final,
    ...
  )
}
```

其中 `rU` 是你自己实现的函数，例如：

- `function(n) r_gaussian_copula_base(n, Sigma)`
- `function(n) r_t_copula_base(n, Sigma, df)`
- `function(n) { u[idx, ] }`（经验 copula 重采样）

这样你就能保持 `simdatacopulas` 的“统一接口体验”，同时把底层依赖结构的实现替换成你需要的数学构造或其他包。

---

## 5) 选择建议（按使用目标）

- 主要是连续变量 + 相关结构为主：优先 Gaussian copula（Base R / `mvtnorm` 即可）。
- 明显关心极端同步（尾部相关）：优先 t copula（Base R 构造可行；df 用网格 profile 选择更稳）。
- 高维（d 大）且希望灵活依赖：优先 vine（`VineCopula` / `rvinecopulib`）。
- 不想选家族、只想复刻真实依赖：经验 copula（rank resampling）。
- 变量包含二值/有序/多分类：优先潜变量阈值模型或其他离散生成机制，而不是强行对原变量 `pobs()`。


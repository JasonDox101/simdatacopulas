# simdatacopulas

`simdatacopulas` 是对 `simdata` 的扩展：把不同类型的 copula（依赖结构）封装为 `simdata::simdesign()` 模板对象，从而可以用统一的方式生成“边缘分布任意 + 相关结构可控”的多变量数据。

本包的核心思想是两步：

1. 用某个 copula 生成相关的 \(U\in(0,1)^d\)（联合依赖结构）
2. 对每一列 \(U_j\) 施加边缘分位数函数 `dist[[j]](u)`，把 \(U\) 映射到你真正想要的变量尺度（边缘分布）

因此 copula 负责“相关/尾部依赖/高维结构”，边缘分位数函数负责“单变量分布形状（连续/离散/截断等）”。

---

## 安装与依赖

- 运行依赖：`simdata`、`copula`
- 可选依赖：`VineCopula`（仅 vine copula 系列需要）

---

## 快速上手（通用工作流）

```r
library(simdata)
library(simdatacopulas)

dist <- list(
  function(u) qnorm(u, mean = 0, sd = 1),
  function(u) qexp(u, rate = 0.2),
  function(u) qbeta(u, 2, 8)
)

design <- simdesign_gaussian_copula(
  dist = dist,
  dim = 3,
  rho = 0.4,
  names_final = c("x1", "t", "b")
)

x <- simulate_data(design, n_obs = 1000, seed = 1)
head(x)
```

---

## 模型总览（每类 copula 的“具体作用”）

本包目前提供 6 组联合生成模块：

- 椭圆族 copula：`simdesign_elliptical_copula()` / `simdesign_gaussian_copula()` / `simdesign_elliptical_copula_from_data()`
- 阿基米德 copula：`simdesign_archimedean_copula()` 及 Clayton/Gumbel/Frank/Joe 的便捷构造器与 `*_from_data()`
- Vine copula（pair-copula construction）：`simdesign_vine_copula()` / `simdesign_vine_copula_from_data()`
- 经验 copula（重采样联合秩结构）：`simdesign_empirical_copula()` / `simdesign_empirical_copula_from_data()`
- 因子 copula（低秩依赖，适合高维）：`simdesign_factor_copula()` / `simdesign_factor_copula_from_data()`
- 潜变量阈值高斯模型（离散/等级终点联合生成）：`simdesign_latent_threshold_gaussian()` / `simdesign_latent_threshold_gaussian_from_data()`

下面逐一说明每个模型的作用、适用场景与用法示例。

---

## 1) 椭圆族 copula（Gaussian / t）

### 作用

- 用“相关矩阵（或交换相关系数）”定义依赖结构，适合模拟**线性相关/对称依赖**的多变量终点。
- `Gaussian copula` 适用于“相关结构主要由相关系数刻画”的场景。
- `t copula` 在高斯基础上增加自由度（`df`），能表达更强的**尾部相关**（极端事件同步出现的概率更大）。

### 何时用

- 你有目标相关矩阵（或想要交换相关结构），希望把多终点（如 PFS、OS、ORR 连续评分等）关联起来。
- 维度中等、结构不复杂，或希望先用一个稳健的“基线依赖结构”。

### 关键接口

- `simdesign_elliptical_copula(copula, dist, ...)`：传入 `copula::normalCopula()` 或 `copula::tCopula()` 等对象
- `simdesign_gaussian_copula(dist, structure = "ex"/"un", ...)`：高斯 copula 便捷入口
- `simdesign_elliptical_copula_from_data(data, vars, family = "gaussian"/"t", ...)`：从数据拟合依赖结构，边缘默认用经验分位数

### 示例：从数据拟合 Gaussian copula

```r
library(simdata)
library(simdatacopulas)

set.seed(1)
pilot <- data.frame(
  pfs = rexp(500, 0.15),
  os  = rexp(500, 0.10),
  biomarker = rnorm(500, 0, 1)
)

design_fit <- simdesign_elliptical_copula_from_data(
  data = pilot,
  vars = c("pfs", "os", "biomarker"),
  family = "gaussian",
  structure = "un"
)

sim <- simulate_data(design_fit, n_obs = 1000, seed = 2)
head(sim)
```

---

## 2) 阿基米德 copula（Clayton / Gumbel / Frank / Joe）

### 作用

阿基米德 copula 通过一个参数（通常记为 `theta`）刻画依赖强度，常用于模拟**非对称尾部依赖**：

- Clayton：**下尾依赖**更强（低值/差结局更容易一起出现）
- Gumbel：**上尾依赖**更强（高值/好结局或极端高值更容易一起出现）
- Frank：对称依赖、通常**无尾依赖**（中间区域依赖更明显）
- Joe：上尾依赖（常比 Gumbel 更“尖锐”）

### 何时用

- 你关心“极端事件/尾部事件的联动”，例如严重不良事件与疗效极端反应的同步性。
- 你希望用少量参数表达依赖（模型更轻、更好控），但仍保留尾部行为差异。

### 关键接口

- `simdesign_archimedean_copula(copula, dist, ...)`
- 便捷构造器：`simdesign_clayton_copula()`、`simdesign_gumbel_copula()`、`simdesign_frank_copula()`、`simdesign_joe_copula()`
- `simdesign_archimedean_copula_from_data(data, vars, family = ..., ...)`：从数据拟合 `theta`（并用经验边缘）

### 示例：用 Gumbel copula 模拟“上尾依赖”

```r
library(simdata)
library(simdatacopulas)

dist <- list(
  function(u) qlnorm(u, meanlog = 0, sdlog = 0.6),
  function(u) qlnorm(u, meanlog = 0.2, sdlog = 0.6)
)

design <- simdesign_gumbel_copula(
  dist = dist,
  theta = 2,
  dim = 2,
  names_final = c("endpoint1", "endpoint2")
)

x <- simulate_data(design, n_obs = 2000, seed = 1)
head(x)
```

---

## 3) Vine copula（C-vine / D-vine / R-vine）

### 作用

Vine copula 用“成对 copula 的树结构（pair-copula construction）”来拼接高维依赖结构：

- 在高维下比“一个整体 copula”更灵活：不同变量对可以选择不同家族、不同尾部行为
- 能表达复杂的条件依赖关系：例如 \(X_1\) 与 \(X_3\) 的相关依赖可能主要通过 \(X_2\) 传递

### 何时用

- 维度较高（例如 6+ 个终点/协变量）且你不希望依赖结构被一个单一相关矩阵/单参数限制。
- 你希望数据驱动地选择 pair-copula 家族（AIC/BIC 等准则）。

### 依赖说明

Vine 模块依赖 `VineCopula` 包；未安装时会报错提示安装。

### 关键接口

- `simdesign_vine_copula(vine, dist, ...)`：传入 `VineCopula` 的 vine 对象（如 `RVineMatrix`）
- `simdesign_vine_copula_from_data(data, vars, vine_type = "rvine"/"cvine"/"dvine", ...)`：从数据拟合 vine 结构与 pair-copula

### 示例：从数据拟合 R-vine（自动结构选择）

```r
library(simdata)
library(simdatacopulas)

set.seed(1)
pilot <- data.frame(
  x1 = rnorm(600),
  x2 = rnorm(600),
  x3 = rnorm(600),
  x4 = rnorm(600)
)

design_fit <- simdesign_vine_copula_from_data(
  data = pilot,
  vars = c("x1", "x2", "x3", "x4"),
  vine_type = "rvine",
  selectioncrit = "AIC",
  method = "itau"
)

sim <- simulate_data(design_fit, n_obs = 1000, seed = 2)
head(sim)
```

---

## 4) 经验 copula（Empirical copula / rank resampling）

### 作用

经验 copula 不显式假设某个参数化 copula 家族，而是：

- 从已有样本的 `U`（伪观测）中**重采样**联合秩结构
- 结合你指定的边缘分位数函数，生成新样本

它的核心价值是：在样本量足够、你相信样本联合秩结构具有代表性时，可以非常直接地“复刻”依赖结构。

### 何时用

- 你有可靠的历史/真实数据，想尽量保留其中的联合依赖结构（尤其是复杂、难以参数化的依赖）。
- 你想避免 copula 家族选择带来的建模偏差。

### 注意点

- 经验 copula 的外推能力有限：它擅长“重现已观察到的依赖”，不擅长“推断未观察到的尾部结构”。
- 可选 `jitter` 用于缓解大量 ties 或离散化带来的重复点问题。

### 关键接口

- `simdesign_empirical_copula(u_data, dist, replace = TRUE, jitter = 0, ...)`
- `simdesign_empirical_copula_from_data(data, vars, ...)`：自动对数据做 `pobs()` 得到 `u_data`

### 示例：直接从数据生成经验 copula 设计

```r
library(simdata)
library(simdatacopulas)

set.seed(1)
pilot <- data.frame(
  a = rexp(800, 0.3),
  b = rlnorm(800, 0, 0.5),
  c = rnorm(800)
)

design_fit <- simdesign_empirical_copula_from_data(
  data = pilot,
  vars = c("a", "b", "c"),
  jitter = 0
)

sim <- simulate_data(design_fit, n_obs = 1000, seed = 2)
head(sim)
```

---

## 5) 因子 copula（Factor Gaussian copula）

### 作用

因子 copula 用“少数公共因子 + 特异噪声”的结构生成高维相关：

\[
Z = F \Lambda^\top + E,\quad U = \Phi(\tilde Z)
\]

- \(\Lambda\)（`loadings`）控制每个变量对公共因子的敏感度
- `uniq_var` 控制每个变量的特异方差

它相当于把协方差/相关结构限制为“低秩 + 对角”，从而在高维下更稳定、参数更少。

### 何时用

- 维度较高（例如多终点 + 多协变量），你希望用少量参数刻画“整体相关驱动因素”（例如疾病严重程度、暴露量、治疗敏感性等潜因子）。
- 你需要一个比完整相关矩阵更稳定、更可解释的相关结构。

### 关键接口

- `simdesign_factor_copula(loadings, uniq_var, dist, ...)`：直接给定因子结构与边缘
- `simdesign_factor_copula_from_data(data, vars, n_factors = 1, ...)`：先把数据转换为伪观测并映射到高斯空间，再做因子拟合

### 示例：从数据拟合 1 因子 copula

```r
library(simdata)
library(simdatacopulas)

set.seed(1)
pilot <- data.frame(
  x1 = rexp(2000, 0.2),
  x2 = rlnorm(2000, 0, 0.6),
  x3 = rnorm(2000)
)

design_fit <- simdesign_factor_copula_from_data(
  data = pilot,
  vars = c("x1", "x2", "x3"),
  n_factors = 1
)

sim <- simulate_data(design_fit, n_obs = 1000, seed = 2)
head(sim)
```

---

## 6) 潜变量阈值高斯模型（Latent threshold Gaussian）

### 作用

该模块用于联合生成**离散/等级终点**（如 AE 分级、ORR 二分类、毒性等级等）：

- 先用因子结构生成连续潜变量 \(Z\)
- 每个离散变量通过阈值向量 `thresholds[[j]]` 把 \(Z_j\) 切分为 1/2/… 多个等级
- 可选 `levels` 为每列指定输出的有序因子标签

这类模型的价值是：你能用高斯因子结构（或多元高斯）表达“离散终点之间的相关”，同时保持各离散变量的边缘分布（由阈值决定）。

### 何时用

- 终点包含等级/二分类变量，且你希望它们与其他终点共享潜在相关来源（公共因子）。
- 你希望显式控制各等级比例（阈值）与等级间相关结构。

### 关键接口

- `simdesign_latent_threshold_gaussian(loadings, uniq_var, thresholds, levels = NULL, ...)`
- `simdesign_latent_threshold_gaussian_from_data(data, vars, n_factors = 1, ...)`：
  - 自动从每列数据估计阈值（按类别累计概率映射到 `qnorm` 切点）
  - 将数据转为伪观测并映射到高斯空间后拟合因子结构

### 示例：手工指定阈值生成二分类 + 三分类联合分布

```r
library(simdata)
library(simdatacopulas)

loadings <- matrix(c(0.9, 0.6, 0.7), nrow = 3, ncol = 1)
uniq <- c(0.4, 0.7, 0.6)

thresholds <- list(
  qnorm(0.30),
  qnorm(0.55),
  qnorm(c(0.20, 0.70))
)

design <- simdesign_latent_threshold_gaussian(
  loadings = loadings,
  uniq_var = uniq,
  thresholds = thresholds,
  names_final = c("y1", "y2", "y3")
)

sim <- simulate_data(design, n_obs = 2000, seed = 1)
head(sim)
```

---

## 与 TrialSimulator 的连接方式（建议）

`TrialSimulator` 对数据生成器的要求是“一个 `generator(n)` 返回终点列的数据框”。本包返回的是 `simdata` 的 `simdesign` 对象，因此你通常会这样桥接：

1. 用 `simdatacopulas::*` 构造 `design`
2. 用 `simdata::simulate_data(design, n_obs = n, seed = ...)` 生成一批“联合相关”的基础变量
3. 在 `TrialSimulator` 的 endpoint generator 内，把这些基础变量映射成 TrialSimulator 需要的终点格式（例如 `*_event`、删失规则等）

这样做的好处是：联合依赖结构与边缘分布建模与试验流程（入组/删失/分析）解耦，便于复用与迭代。


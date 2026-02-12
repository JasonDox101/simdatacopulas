# simdesign_factor_latent_threshold

## 相关函数

- simdesign_factor_copula
- simdesign_factor_copula_from_data
- simdesign_latent_threshold_gaussian
- simdesign_latent_threshold_gaussian_from_data

## 用法

```
simdesign_factor_copula(loadings, uniq_var = NULL, dist, names_final = NULL, process_final = list(), name = "Factor Gaussian copula design", eps = 1e-6, ...)
simdesign_factor_copula_from_data(data, vars, n_factors = 1, rotation = "none", qtype = 8, eps = 1e-6, name = "Factor Gaussian copula design (fit from data)", ...)
simdesign_latent_threshold_gaussian(loadings, uniq_var = NULL, thresholds, levels = NULL, names_final = NULL, process_final = list(), name = "Latent Gaussian threshold design", eps = 1e-6, ...)
simdesign_latent_threshold_gaussian_from_data(data, vars, n_factors = 1, rotation = "none", tie_jitter = 1e-8, eps = 1e-6, name = "Latent Gaussian threshold design (fit from data)", ...)
```

## 参数

- loadings：因子载荷矩阵
- uniq_var：唯一性方差向量
- dist：边缘分布的分位数函数列表
- names_final：输出变量名
- process_final：simdata 的后处理列表
- name：设计名称
- eps：数值截断参数
- data：输入数据框
- vars：建模变量名
- n_factors：因子数
- rotation：因子旋转方法
- qtype：经验分位类型
- thresholds：阈值列表
- levels：有序水平列表
- tie_jitter：抖动幅度
- ...：附加元数据

## 返回

- simdesign 对象，可用于 simdata::simulate_data()

## 说明

用于构建基于因子结构或阈值模型的 simdesign，支持从数据拟合参数。

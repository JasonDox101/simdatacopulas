# simdesign_elliptical_copula

## 相关函数

- simdesign_elliptical_copula
- simdesign_gaussian_copula
- simdesign_elliptical_copula_from_data

## 用法

```
simdesign_elliptical_copula(copula, dist, names_final = NULL, process_final = list(), name = "Elliptical copula design", eps = 1e-6, ...)
simdesign_gaussian_copula(dist, rho = 0, Sigma = NULL, dim = NULL, structure = c("ex", "un"), names_final = NULL, ...)
simdesign_elliptical_copula_from_data(data, vars, family = c("gaussian", "t"), structure = c("ex", "un"), fit_method = c("itau", "itau_mpl"), df = NULL, qtype = 8, eps = 1e-6, name = "Elliptical copula design (fit from data)", ...)
```

## 参数

- copula：copula 对象
- dist：边缘分布的分位数函数列表
- names_final：输出变量名
- process_final：simdata 的后处理列表
- name：设计名称
- eps：数值截断参数
- rho：交换型相关系数
- Sigma：相关矩阵
- dim：维度
- structure：相关结构 ex 或 un
- data：输入数据框
- vars：建模变量名
- family：gaussian 或 t
- fit_method：itau 或 itau_mpl
- df：t 分布自由度
- qtype：经验分位类型
- ...：附加元数据

## 返回

- simdesign 对象，可用于 simdata::simulate_data()

## 说明

用于构建椭圆 copula 相关结构的 simdesign，支持直接指定 copula 或从数据拟合。

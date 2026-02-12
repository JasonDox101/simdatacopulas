# simdesign_archimedean_copula

## 相关函数

- simdesign_archimedean_copula
- simdesign_clayton_copula
- simdesign_gumbel_copula
- simdesign_frank_copula
- simdesign_joe_copula
- simdesign_archimedean_copula_from_data

## 用法

```
simdesign_archimedean_copula(copula, dist, names_final = NULL, process_final = list(), name = "Archimedean copula design", eps = 1e-6, ...)
simdesign_clayton_copula(dist, theta, dim, names_final = NULL, ...)
simdesign_gumbel_copula(dist, theta, dim, names_final = NULL, ...)
simdesign_frank_copula(dist, theta, dim, names_final = NULL, ...)
simdesign_joe_copula(dist, theta, dim, names_final = NULL, ...)
simdesign_archimedean_copula_from_data(data, vars, family = c("clayton", "gumbel", "frank", "joe"), fit_method = c("itau", "itau_mpl"), qtype = 8, eps = 1e-6, name = "Archimedean copula design (fit from data)", ...)
```

## 参数

- copula：copula 对象
- dist：边缘分布的分位数函数列表
- names_final：输出变量名
- process_final：simdata 的后处理列表
- name：设计名称
- eps：数值截断参数
- theta：copula 参数
- dim：维度
- data：输入数据框
- vars：建模变量名
- family：clayton、gumbel、frank 或 joe
- fit_method：itau 或 itau_mpl
- qtype：经验分位类型
- ...：附加元数据

## 返回

- simdesign 对象，可用于 simdata::simulate_data()

## 说明

用于构建阿基米德 copula 的 simdesign，支持参数化构造与数据拟合。

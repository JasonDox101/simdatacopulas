# simdesign_empirical_copula

## 相关函数

- simdesign_empirical_copula
- simdesign_empirical_copula_from_data

## 用法

```
simdesign_empirical_copula(u_data, dist, names_final = NULL, process_final = list(), replace = TRUE, jitter = 0, name = "Empirical copula design", eps = 1e-6, ...)
simdesign_empirical_copula_from_data(data, vars, qtype = 8, replace = TRUE, jitter = 0, eps = 1e-6, name = "Empirical copula design (fit from data)", ...)
```

## 参数

- u_data：伪观测矩阵或数据框
- dist：边缘分布的分位数函数列表
- names_final：输出变量名
- process_final：simdata 的后处理列表
- replace：是否有放回抽样
- jitter：抖动幅度
- name：设计名称
- eps：数值截断参数
- data：输入数据框
- vars：建模变量名
- qtype：经验分位类型
- ...：附加元数据

## 返回

- simdesign 对象，可用于 simdata::simulate_data()

## 说明

用于构建经验 copula 的 simdesign，可直接提供伪观测或从数据拟合。

# BYC 复现差异与风险说明

## 已完成

-  from src.byc_model import BYCInput, run_model​inp = BYCInput(    z_ref_u=10.0,    z_ref_tq=10.0,    u_ref=7.0,    t_air_c=27.0,    rh=0.82,    sst_c=28.5,    p_hpa=1010.0,)​out = run_model(inp, z_max=50.0, dz=0.2)print(out["evaporation_duct_height_m"])python
- 完成符号映射（`spec/symbols.yaml`）。
- 完成最小测试集（`tests/test_byc_model.py`）。

## 差异点

1. 稳定函数
- 当前使用 Businger-Dyer 形式。
- 不同 BYC 文献版本可能对稳定条件分段与系数有修订。

2. 粗糙度参数化
- 当前使用 Charnock + 粘性项。
- 若原文存在特定海况修正项（浪龄/海浪状态），需替换。

3. 低风速与阵风项
- 部分工程论文在 BYC 外加低风速修正和阵风项。
- 当前版本未引入该扩展项，避免与“原始 BYC”混淆。

## 风险

- 若目标是“逐字逐式同版复现”，必须拿到目标论文原文（至少公式页 + 参数表）。

## 下一步

- 接入你提供的 BYC 原文 PDF。
- 按页码逐式核对：稳定函数、常数、边界条件、迭代准则。
- 追加基于同一观测样本的指标对齐实验（RMSE/相关系数/EDH 分布对比）。

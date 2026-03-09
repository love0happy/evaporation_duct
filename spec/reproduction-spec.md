# BYC 蒸发波导模型复现规范（Python）

## 1. 目标与范围

- 目标：复现 Babin-Young-Carton (BYC) 体系下的海面蒸发波导高度计算流程。
- 范围：输入海气常规观测（风速、气温、湿度、海温、气压），输出修正折射率廓线 `M(z)` 与蒸发波导高度 `EDH`。
- 高度范围：默认 `0.1 m ~ 50 m`。

## 2. 核心方程

> 说明：为了兼容 Typora，这里同时给出 **纯文本公式** + **LaTeX 公式**。

1) 折射率与修正折射率

- 纯文本：`N(z) = 77.6 * P / T(z) + 3.73e5 * e(z) / T(z)^2`
- 纯文本：`M(z) = N(z) + 0.157 * z`

$$
N(z)=77.6\frac{P}{T(z)}+3.73\times 10^5\frac{e(z)}{T(z)^2}
$$
$$
M(z)=N(z)+0.157z
$$

2) MO 相似理论剖面

- 纯文本：`U(z) = (u*/kappa) * [ln(z/z0m) - psi_m(z/L) + psi_m(z0m/L)]`
- 纯文本：`theta(z)-theta_s = (theta*/kappa) * [ln(z/z0h) - psi_h(z/L) + psi_h(z0h/L)]`
- 纯文本：`q(z)-q_s = (q*/kappa) * [ln(z/z0q) - psi_h(z/L) + psi_h(z0q/L)]`

$$
U(z)=\frac{u_*}{\kappa}\left[\ln\frac{z}{z_{0m}}-\psi_m\left(\frac{z}{L}\right)+\psi_m\left(\frac{z_{0m}}{L}\right)\right]
$$
$$
\theta(z)-\theta_s=\frac{\theta_*}{\kappa}\left[\ln\frac{z}{z_{0h}}-\psi_h\left(\frac{z}{L}\right)+\psi_h\left(\frac{z_{0h}}{L}\right)\right]
$$
$$
q(z)-q_s=\frac{q_*}{\kappa}\left[\ln\frac{z}{z_{0q}}-\psi_h\left(\frac{z}{L}\right)+\psi_h\left(\frac{z_{0q}}{L}\right)\right]
$$

3) 粗糙度参数化

- 纯文本：`z0m = alpha * u*^2 / g + beta * nu / u*`

$$
z_{0m}=\alpha\frac{u_*^2}{g}+\beta\frac{\nu}{u_*}
$$

4) 波导高度判据

- 在离散高度网格中寻找 `M(z)` 的首个局部极小值；对应高度记为 `EDH`。

## 3. 输入输出定义

### 输入

- `z_ref_u`：风速参考高度（m）
- `z_ref_tq`：温湿参考高度（m）
- `u_ref`：参考高度风速（m/s）
- `t_air_c`：气温（°C）
- `rh`：相对湿度（0~1）
- `sst_c`：海表温度（°C）
- `p_hpa`：气压（hPa）

### 输出

- `state`: `u_*`, `t_*`, `q_*`, `z0m`, `z0h`, `z0q`, `L`
- `profile`: `z`, `N(z)`, `M(z)`
- `evaporation_duct_height_m`: 蒸发波导高度（m, 可空）

## 4. 数值策略

- 当前实现已接入 `pycoare`，由 `pycoare.coare_35` 计算通量状态量：`u_*`, `t_*`, `q_*`, `z0m`, `z0h`, `z0q`, `L`。
- 输入映射：本项目 `rh` 使用 0~1，传入 pycoare 时转换为百分比（`rh * 100`）。
- 近地层剖面重建时，稳定函数调用 `pycoare.util.psit_26`。
- 这样做的目的：将通量算法交给可追溯的 COARE 实现，减少手写迭代误差。

## 5. 量纲与一致性检查

- `N`、`M` 量纲：N-unit / M-unit。
- `u_*`：m/s；`z0*`：m；`L`：m。
- 单元测试检查：
  - `N` 在典型海面环境下落在合理范围；
  - `M-N=0.157z` 恒成立；
  - 粗糙度和摩擦速度保持正值。

## 6. 参考文献线索（检索结果）

- Babin et al., 1997, *A New Model of the Oceanic Evaporation Duct*, DOI: 10.1175/1520-0450(1997)036<0193:ANMOTO>2.0.CO;2
- Levy, 2002, *LKB-Based Evaporation Duct Model Comparison with Buoy Data*, DOI: 10.1175/1520-0450(2002)041<0434:LBEDMC>2.0.CO;2
- Zhao et al., 2021, IEEE TAP, DOI: 10.1109/TAP.2021.3076478（将 BYC 作为 baseline）

## 7. 当前复现缺口

- 由于原文受访问限制，BYC 原始文内某些经验常数与稳定函数细节可能存在版本差异。
- 若提供原文 PDF（尤其公式页），可做逐项常数校准，完成“文内同版”严格复现。

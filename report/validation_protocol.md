# BYC 模型正确性验证协议（可证伪）

> 目标：验证当前实现是否符合 Babin 模型真实计算步骤、输入输出定义和数值行为。

---

## 0. 验证结论分级

- **A 级（严格同版）**：与目标 Babin 论文版本逐式一致，算例误差在阈值内。
- **B 级（工程同族）**：流程与物理机制一致，但通量内核版本存在差异（例如使用 pycoare/coare3.5）。
- **C 级（不可接受）**：流程、参数或数值存在关键错误。

当前目标：先完成 B 级可审计验证，再冲刺 A 级逐式验证。

---

## 1. 输入输出参数验证（参数正确性）

## 1.1 输入参数（必须项）

| 参数名 | 含义 | 单位 | 当前实现字段 | 检查规则 |
|---|---|---:|---|---|
| z_ref_u | 风速观测高度 | m | `BYCInput.z_ref_u` | > 0 |
| z_ref_tq | 温湿观测高度 | m | `BYCInput.z_ref_tq` | > 0 |
| u_ref | 参考高度风速 | m/s | `BYCInput.u_ref` | >= 0 |
| t_air_c | 气温 | °C | `BYCInput.t_air_c` | 合理区间 |
| rh | 相对湿度 | 0~1 | `BYCInput.rh` | 0 <= rh <= 1 |
| sst_c | 海表温度 | °C | `BYCInput.sst_c` | 合理区间 |
| p_hpa | 气压 | hPa | `BYCInput.p_hpa` | > 800 且 < 1100（常规） |
| lat_deg | 纬度 | deg | `BYCInput.lat_deg` | -90~90 |

关键映射确认：`rh` 传入 pycoare 前必须转换为百分比 `rh*100`。

## 1.2 输出参数（必须项）

| 输出项 | 含义 | 单位 | 代码位置 |
|---|---|---:|---|
| state.u_star | 摩擦速度 | m/s | `solve_state` |
| state.t_star | 温度尺度 | K | `solve_state` |
| state.q_star | 湿度尺度 | kg/kg | `solve_state` |
| state.z0m/z0h/z0q | 粗糙度长度 | m | `solve_state` |
| state.L | Obukhov 长度 | m | `solve_state` |
| profile.z | 高度网格 | m | `profile` |
| profile.N | 折射率 | N-unit | `profile` |
| profile.M | 修正折射率 | M-unit | `profile` |
| evaporation_duct_height_m | 蒸发波导高度 | m | `run_model` |

---

## 2. 计算流程验证（步骤正确性）

按步骤检查（任一步失败即流程失败）：

1) 读取观测输入并完成单位一致化。  
2) 调用 COARE 通量算法求状态量（u*、t*、q*、z0、L）。  
3) 基于相似理论重建温湿垂直廓线。  
4) 计算水汽压 e(z)。  
5) 计算 N(z)：`N = 77.6*P/T + 3.73e5*e/T^2`。  
6) 计算 M(z)：`M = N + 0.157*z`。  
7) 从 M(z) 提取首个局部极小点作为 EDH。  

代码对应关系：
- 步骤2：`solve_state()`（当前来源 `pycoare.coare_35`）
- 步骤3~6：`profile()`
- 步骤7：`duct_height()`
- 全流程：`run_model()`

---

## 3. 数值与物理行为验证（结果正确性）

## 3.1 基础数值检查

- `M - N = 0.157*z` 必须逐点成立。  
- `u_star > 0`, `z0m/z0h/z0q > 0`。  
- `profile.z/N/M` 三数组长度一致。

## 3.2 物理趋势检查（定性）

- 风速上升时，状态量变化应平滑，无非物理突变。  
- 极端低风速场景不应出现明显数值爆炸。  
- 若无局部极小点，EDH 允许为 `None`（需记录条件）。

## 3.3 对标误差检查（定量）

当获得论文算例后，使用如下指标：

- 绝对误差：`|EDH_model - EDH_ref|`
- 相对误差：`|EDH_model - EDH_ref| / max(EDH_ref, eps)`
- 廓线误差：`RMSE(M_model(z), M_ref(z))`

建议阈值（可协商）：
- EDH 相对误差 <= 10%
- M(z) RMSE <= 2 M-unit

---

## 4. A/B 双轨验证策略

- **A轨（严格 Babin 轨）**：以目标论文公式、常数、分段函数为唯一标准。  
- **B轨（工程实现轨）**：以 pycoare 作为通量核心，验证工程稳定性与可复现性。  

说明：当前代码属于 B轨。若要达到 A轨，需要论文原文公式页与参数细节。

---

## 5. 验证执行清单（打勾）

- [ ] 输入参数范围与单位校验通过
- [ ] 输出字段完整性校验通过
- [ ] 流程 1~7 步骤逐项核对通过
- [ ] 基础数值检查通过
- [ ] 物理趋势检查通过
- [ ] 与基准算例定量误差检查通过
- [ ] 生成结论等级（A/B/C）

---

## 6. 下一步你我协作输入

你提供任意一项即可进入定量验证：

1) Babin 目标论文 PDF（优先：公式页、参数页、算例页）  
2) 你已有的观测样本 + 可信结果（EDH 或 M(z)）  
3) 你认可的对照实现输出（MATLAB/Fortran/历史脚本）

我将据此输出：
- `validation_report.md`（通过/失败项）
- `error_budget.csv`（误差分解）
- `step_mapping_table.md`（论文步骤 -> 代码行）

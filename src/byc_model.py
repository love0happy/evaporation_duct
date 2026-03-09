"""BYC 蒸发波导模型（基于 pycoare 通量核心）。

说明：
1) 本模块使用 pycoare.coare_35 计算海气通量相关状态量（u*、t*、q*、粗糙度、Obukhov 长度）。
2) 在此基础上，按近地层相似理论重建温湿廓线，再计算 N(z)、M(z) 与蒸发波导高度。
3) 公式采用人类可读表达：
   - N = 77.6*P/T + 3.73e5*e/T^2
   - M = N + 0.157*z
"""

from __future__ import annotations

from dataclasses import dataclass
import math
from typing import Dict, List, Optional

# 冯·卡门常数
KAPPA = 0.4


@dataclass
class BYCInput:
    """BYC 模型输入参数。

    字段说明：
    - z_ref_u: 风速观测高度（m）
    - z_ref_tq: 温湿观测高度（m）
    - u_ref: 参考高度风速（m/s）
    - t_air_c: 气温（摄氏度）
    - rh: 相对湿度（0~1）
    - sst_c: 海表温度（摄氏度）
    - p_hpa: 气压（hPa）
    - lat_deg: 纬度（度），供 pycoare 使用
    """

    z_ref_u: float
    z_ref_tq: float
    u_ref: float
    t_air_c: float
    rh: float
    sst_c: float
    p_hpa: float
    lat_deg: float = 20.0


@dataclass
class BYCState:
    """BYC 通量状态量。

    字段说明：
    - u_star: 摩擦速度 u*（m/s）
    - t_star: 温度尺度 t*（K）
    - q_star: 比湿尺度 q*（kg/kg）
    - z0m: 动量粗糙度（m）
    - z0h: 温度粗糙度（m）
    - z0q: 湿度粗糙度（m）
    - L: Monin-Obukhov 长度（m）
    - source: 状态量来源标识
    """

    u_star: float
    t_star: float
    q_star: float
    z0m: float
    z0h: float
    z0q: float
    L: float
    source: str


def saturation_vapor_pressure_hpa(t_c: float) -> float:
    """计算饱和水汽压（hPa）。

    作用：根据摄氏温度计算饱和水汽压，供比湿与折射率计算使用。
    方法：Tetens 经验公式。
    """
    return 6.112 * math.exp((17.67 * t_c) / (t_c + 243.5))


def specific_humidity_from_e(e_hpa: float, p_hpa: float) -> float:
    """根据水汽压和总压计算比湿（kg/kg）。

    作用：将水汽压 e 与气压 p 转为比湿 q。
    """
    return 0.622 * e_hpa / (p_hpa - 0.378 * e_hpa)


def refractivity_n(p_hpa: float, t_k: float, e_hpa: float) -> float:
    """计算大气折射率 N。

    作用：根据气压、温度、水汽压计算折射率 N（N-unit）。
    公式：N = 77.6*P/T + 3.73e5*e/T^2
    """
    return 77.6 * p_hpa / t_k + 3.73e5 * e_hpa / (t_k * t_k)


def modified_refractivity_m(n: float, z_m: float) -> float:
    """计算修正折射率 M。

    作用：由折射率 N 和高度 z 得到修正折射率 M（M-unit）。
    公式：M = N + 0.157*z
    """
    return n + 0.157 * z_m


def _as_scalar(x) -> float:
    """将标量/单元素数组统一转换为 float。

    作用：兼容 pycoare 返回的 ndarray，便于下游统一处理。
    """
    try:
        return float(x[0])
    except Exception:
        return float(x)


def solve_state(inp: BYCInput) -> BYCState:
    """调用 pycoare 计算通量状态量。

    作用：
    - 将输入观测量映射到 pycoare.coare_35。
    - 返回 u*、t*、q*、z0m/z0h/z0q、L 等关键状态量。

    关键映射：
    - 本项目 rh 为 0~1；pycoare 需要百分比，因此传入 rh*100。
    - p_hpa 与 pycoare 的 p（mb/hPa）一致。
    """
    try:
        from pycoare import coare_35
    except Exception as exc:  # pragma: no cover
        raise RuntimeError("缺少 pycoare 依赖，请先执行: pip install pycoare") from exc

    c = coare_35(
        [inp.u_ref],
        t=[inp.t_air_c],
        rh=[inp.rh * 100.0],
        zu=[inp.z_ref_u],
        zt=[inp.z_ref_tq],
        zq=[inp.z_ref_tq],
        zrf=[inp.z_ref_u],
        ts=[inp.sst_c],
        p=[inp.p_hpa],
        lat=[inp.lat_deg],
        nits=10,
        jcool=1,
    )

    return BYCState(
        u_star=_as_scalar(c.velocities.usr),
        t_star=_as_scalar(c.stability_parameters.tsr),
        q_star=_as_scalar(c.stability_parameters.qsr),
        z0m=_as_scalar(c.stability_parameters.zo),
        z0h=_as_scalar(c.stability_parameters.zot),
        z0q=_as_scalar(c.stability_parameters.zoq),
        L=_as_scalar(c.stability_parameters.obukL),
        source="pycoare.coare_35",
    )


def _psi_h(z_over_L: float) -> float:
    """计算热量/湿度稳定函数 psi_h。

    作用：调用 pycoare 的稳定函数实现，避免手写版本偏差。
    """
    from pycoare.util import psit_26

    return _as_scalar(psit_26([z_over_L]))


def profile(inp: BYCInput, state: BYCState, z_max: float = 50.0, dz: float = 0.1) -> Dict[str, List[float]]:
    """重建垂直廓线并计算 N(z)、M(z)。

    作用：
    - 使用状态量（t*、q*、L）将参考高度温湿外推到不同高度 z。
    - 逐层计算 e(z)、N(z)、M(z)。

    参数：
    - z_max: 最大高度（m）
    - dz: 高度步长（m）

    返回：
    - 字典包含 z、N、M 三个等长数组。
    """
    if dz <= 0:
        raise ValueError("dz 必须大于 0")

    t_air_k = inp.t_air_c + 273.15
    e_air = inp.rh * saturation_vapor_pressure_hpa(inp.t_air_c)
    q_air_ref = specific_humidity_from_e(e_air, inp.p_hpa)

    z_values: List[float] = []
    n_values: List[float] = []
    m_values: List[float] = []

    # 防止 L 过小导致数值不稳定
    L_eff = state.L if abs(state.L) > 1e-9 else 1e9

    # 从近海面最小有效高度开始积分
    z = max(0.1, min(inp.z_ref_tq, inp.z_ref_u) * 0.1)
    while z <= z_max + 1e-12:
        psi_ref = _psi_h(inp.z_ref_tq / L_eff)
        psi_z = _psi_h(z / L_eff)

        # 相似理论重建温湿剖面（人类可读）
        # q(z) = q_ref - (q*/kappa) * [ln(z_ref/z) - psi_h(z_ref/L) + psi_h(z/L)]
        # T(z) = T_ref - (t*/kappa) * [ln(z_ref/z) - psi_h(z_ref/L) + psi_h(z/L)]
        qz = q_air_ref - (state.q_star / KAPPA) * (math.log(inp.z_ref_tq / z) - psi_ref + psi_z)
        tz = t_air_k - (state.t_star / KAPPA) * (math.log(inp.z_ref_tq / z) - psi_ref + psi_z)

        # 由比湿反推水汽压，再计算折射率
        ez = (qz * inp.p_hpa) / (0.622 + 0.378 * qz)
        nz = refractivity_n(inp.p_hpa, tz, ez)
        mz = modified_refractivity_m(nz, z)

        z_values.append(z)
        n_values.append(nz)
        m_values.append(mz)
        z += dz

    return {"z": z_values, "N": n_values, "M": m_values}


def duct_height(z: List[float], m: List[float]) -> Optional[float]:
    """提取蒸发波导高度。

    作用：在离散的 M(z) 中寻找首个局部极小值位置，作为蒸发波导高度。
    返回：
    - 找到则返回高度（m）
    - 未找到返回 None
    """
    if len(z) < 3 or len(z) != len(m):
        return None
    for i in range(1, len(z) - 1):
        if m[i] <= m[i - 1] and m[i] <= m[i + 1]:
            return z[i]
    return None


def run_model(inp: BYCInput, z_max: float = 50.0, dz: float = 0.1) -> Dict[str, object]:
    """执行 BYC 计算主流程。

    作用：
    1) 先计算通量状态量；
    2) 再重建廓线；
    3) 最后提取蒸发波导高度。

    返回：
    - state: 状态量
    - profile: 廓线数据（z、N、M）
    - evaporation_duct_height_m: 波导高度
    """
    state = solve_state(inp)
    prof = profile(inp, state, z_max=z_max, dz=dz)
    edh = duct_height(prof["z"], prof["M"])
    return {
        "state": state,
        "profile": prof,
        "evaporation_duct_height_m": edh,
    }

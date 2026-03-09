from __future__ import annotations

"""全球海洋蒸发波导高度（EDH）格点计算模块。

设计目标：
1. 读取 ECMWF C1D 单文件 GRIB 数据；
2. 在不改动现有单点 `src/byc_model.py` 的前提下，新增格点批量计算能力；
3. 输出全球海洋 EDH 格点场（NetCDF）。

方法说明（与单点代码保持一致的核心思路）：
- 通量状态量：使用 `pycoare.coare_35` 计算 `t*`、`q*`、`L`；
- 廓线重建：按相似理论重建 `T(z)`、`q(z)`；
- 折射率计算：`N = 77.6*P/T + 3.73e5*e/T^2`；
- 修正折射率：`M = N + 0.157*z`；
- EDH 定义：取 `M(z)` 最小值对应高度。
"""

from dataclasses import dataclass
from pathlib import Path
from typing import Dict

import numpy as np
import xarray as xr


@dataclass
class GridConfig:
    """格点计算配置参数。"""

    z_max: float = 40.0
    dz: float = 0.5
    chunk_size: int = 200_000
    use_ocean_mask: bool = True
    sst_min_c: float = -2.5
    sst_max_c: float = 45.0


def _open_grib_var(grib_path: str, short_name: str) -> xr.Dataset:
    """按 shortName 读取 GRIB 变量。

    参数：
    - grib_path: C1D 文件路径
    - short_name: 变量 shortName（例如 10u/10v/2t/2d/sst/msl）
    """

    return xr.open_dataset(
        grib_path,
        engine="cfgrib",
        backend_kwargs={"filter_by_keys": {"shortName": short_name}, "indexpath": ""},
    )


def load_c1d_required_fields(grib_path: str) -> Dict[str, np.ndarray]:
    """读取计算 EDH 所需的 C1D 近地面变量。"""

    ds_10u = _open_grib_var(grib_path, "10u")
    ds_10v = _open_grib_var(grib_path, "10v")
    ds_2t = _open_grib_var(grib_path, "2t")
    ds_2d = _open_grib_var(grib_path, "2d")
    ds_sst = _open_grib_var(grib_path, "sst")
    ds_msl = _open_grib_var(grib_path, "msl")

    return {
        "lat": ds_10u["latitude"].values,
        "lon": ds_10u["longitude"].values,
        "u10": ds_10u["u10"].values.astype(np.float32),
        "v10": ds_10v["v10"].values.astype(np.float32),
        "t2m_k": ds_2t["t2m"].values.astype(np.float32),
        "d2m_k": ds_2d["d2m"].values.astype(np.float32),
        "sst_k": ds_sst["sst"].values.astype(np.float32),
        "msl_pa": ds_msl["msl"].values.astype(np.float32),
    }


def rh_from_t_td_percent(t_c: np.ndarray, td_c: np.ndarray) -> np.ndarray:
    """由气温与露点温度计算相对湿度（百分数）。"""

    es = np.exp((17.625 * t_c) / (243.04 + t_c))
    ed = np.exp((17.625 * td_c) / (243.04 + td_c))
    rh = 100.0 * ed / np.maximum(es, 1e-12)
    return np.clip(rh, 1.0, 100.0).astype(np.float32)


def build_ocean_mask(lat_1d: np.ndarray, lon_1d: np.ndarray, enable: bool) -> np.ndarray:
    """构建海洋掩膜（True=海洋，False=陆地）。"""

    lat2d = np.repeat(lat_1d[:, None], len(lon_1d), axis=1)
    lon2d = np.repeat(lon_1d[None, :], len(lat_1d), axis=0)

    if not enable:
        return np.ones_like(lat2d, dtype=bool)

    from global_land_mask import globe

    # global_land_mask 使用 -180~180 经度
    lon_wrapped = ((lon2d + 180.0) % 360.0) - 180.0
    return globe.is_ocean(lat2d, lon_wrapped)


def _compute_chunk_edh(
    wind10: np.ndarray,
    t2m_c: np.ndarray,
    rh_pct: np.ndarray,
    sst_c: np.ndarray,
    p_hpa: np.ndarray,
    lat_deg: np.ndarray,
    z_max: float,
    dz: float,
) -> np.ndarray:
    """对单个扁平分块执行 EDH 计算。"""

    from pycoare import coare_35
    from pycoare.util import psit_26, qair

    kappa = 0.4

    # 1) 通量核心：coare_35
    c = coare_35(
        wind10,
        t=t2m_c,
        rh=rh_pct,
        zu=np.full_like(wind10, 10.0),
        zt=np.full_like(wind10, 2.0),
        zq=np.full_like(wind10, 2.0),
        zrf=np.full_like(wind10, 10.0),
        ts=sst_c,
        p=p_hpa,
        lat=lat_deg,
        nits=10,
        jcool=1,
    )

    t_star = np.asarray(c.stability_parameters.tsr, dtype=np.float32)
    q_star = np.asarray(c.stability_parameters.qsr, dtype=np.float32)
    L = np.asarray(c.stability_parameters.obukL, dtype=np.float32)

    # qair 输出 g/kg，这里统一成 kg/kg
    q_ref = np.asarray(qair(t2m_c, p_hpa, rh_pct), dtype=np.float32) / 1000.0
    t_ref_k = (t2m_c + 273.15).astype(np.float32)

    z_ref = 2.0
    L_eff = np.where(np.abs(L) < 1e-8, np.where(L >= 0, 1e8, -1e8), L).astype(np.float32)

    z_levels = np.arange(0.1, z_max + 1e-12, dz, dtype=np.float32)
    m_min = np.full_like(wind10, np.inf, dtype=np.float32)
    z_at_min = np.full_like(wind10, np.nan, dtype=np.float32)

    psi_ref = np.asarray(psit_26(z_ref / L_eff), dtype=np.float32)

    # 2) 沿垂直方向扫描，逐层更新最小 M
    for z in z_levels:
        psi_z = np.asarray(psit_26(z / L_eff), dtype=np.float32)
        delta = np.log(z_ref / z) - psi_ref + psi_z

        qz = q_ref - (q_star / kappa) * delta
        tz = t_ref_k - (t_star / kappa) * delta

        ez = (qz * p_hpa) / np.maximum(0.622 + 0.378 * qz, 1e-8)
        n = 77.6 * p_hpa / np.maximum(tz, 150.0) + 3.73e5 * ez / np.maximum(tz * tz, 1.0)
        m = n + 0.157 * z

        better = m < m_min
        m_min[better] = m[better]
        z_at_min[better] = z

    return z_at_min


def compute_global_ocean_edh(grib_path: str, cfg: GridConfig) -> xr.Dataset:
    """基于 C1D 单个 GRIB 文件计算全球海洋 EDH。"""

    f = load_c1d_required_fields(grib_path)

    lat = f["lat"]
    lon = f["lon"]

    wind10 = np.sqrt(f["u10"] * f["u10"] + f["v10"] * f["v10"]).astype(np.float32)
    t2m_c = (f["t2m_k"] - 273.15).astype(np.float32)
    td2m_c = (f["d2m_k"] - 273.15).astype(np.float32)
    sst_c = (f["sst_k"] - 273.15).astype(np.float32)
    p_hpa = (f["msl_pa"] / 100.0).astype(np.float32)

    rh_pct = rh_from_t_td_percent(t2m_c, td2m_c)

    lat2d = np.repeat(lat[:, None], len(lon), axis=1).astype(np.float32)
    ocean = build_ocean_mask(lat, lon, cfg.use_ocean_mask)

    qc = (
        np.isfinite(wind10)
        & np.isfinite(t2m_c)
        & np.isfinite(sst_c)
        & np.isfinite(p_hpa)
        & (sst_c >= cfg.sst_min_c)
        & (sst_c <= cfg.sst_max_c)
        & (p_hpa > 800.0)
        & (p_hpa < 1100.0)
    )
    valid = ocean & qc

    edh_flat = np.full(wind10.size, np.nan, dtype=np.float32)
    valid_idx = np.flatnonzero(valid.ravel())

    wind_flat = wind10.ravel()
    t_flat = t2m_c.ravel()
    rh_flat = rh_pct.ravel()
    sst_flat = sst_c.ravel()
    p_flat = p_hpa.ravel()
    lat_flat = lat2d.ravel()

    for start in range(0, len(valid_idx), cfg.chunk_size):
        idx = valid_idx[start : start + cfg.chunk_size]
        edh_part = _compute_chunk_edh(
            wind10=wind_flat[idx],
            t2m_c=t_flat[idx],
            rh_pct=rh_flat[idx],
            sst_c=sst_flat[idx],
            p_hpa=p_flat[idx],
            lat_deg=lat_flat[idx],
            z_max=cfg.z_max,
            dz=cfg.dz,
        )
        edh_flat[idx] = edh_part

    edh = edh_flat.reshape(wind10.shape)

    ds = xr.Dataset(
        data_vars={"edh": (("latitude", "longitude"), edh)},
        coords={"latitude": lat, "longitude": lon},
        attrs={
            "title": "Global ocean evaporation duct height from ECMWF C1D",
            "method": "pycoare.coare_35 + MO profile + argmin(M)",
            "source_grib": str(grib_path),
            "z_max_m": cfg.z_max,
            "dz_m": cfg.dz,
        },
    )
    ds["edh"].attrs["units"] = "m"
    return ds


def save_to_netcdf(ds: xr.Dataset, out_path: str) -> None:
    """保存结果到 NetCDF 文件。"""

    out = Path(out_path)
    out.parent.mkdir(parents=True, exist_ok=True)
    ds.to_netcdf(out)

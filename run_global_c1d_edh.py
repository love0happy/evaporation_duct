from __future__ import annotations

"""C1D 全球海洋 EDH 格点计算入口脚本。"""

import argparse
import time

from src.grid_global_edh import GridConfig, compute_global_ocean_edh, save_to_netcdf


def build_parser() -> argparse.ArgumentParser:
    """构建命令行参数。"""

    p = argparse.ArgumentParser(description="使用 ECMWF C1D GRIB 计算全球海洋蒸发波导高度")
    p.add_argument("--grib", required=True, help="C1D GRIB 文件路径")
    p.add_argument("--out", required=True, help="输出 NetCDF 文件路径")
    p.add_argument("--z-max", type=float, default=40.0, help="垂直最大高度（m）")
    p.add_argument("--dz", type=float, default=0.5, help="垂直步长（m）")
    p.add_argument("--chunk-size", type=int, default=200000, help="分块大小")
    p.add_argument("--no-ocean-mask", action="store_true", help="关闭海洋掩膜")
    return p


def main() -> None:
    """主流程：读取参数、计算 EDH、保存文件。"""

    args = build_parser().parse_args()
    cfg = GridConfig(
        z_max=args.z_max,
        dz=args.dz,
        chunk_size=args.chunk_size,
        use_ocean_mask=not args.no_ocean_mask,
    )

    t0 = time.time()
    ds = compute_global_ocean_edh(args.grib, cfg)
    save_to_netcdf(ds, args.out)
    t1 = time.time()

    print(f"输出文件: {args.out}")
    print(f"网格大小: {ds['edh'].shape}")
    print(f"耗时: {t1 - t0:.1f} 秒")


if __name__ == "__main__":
    main()

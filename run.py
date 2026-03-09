from src.byc_model import BYCInput, run_model

# 1) 构造输入参数
inp = BYCInput(
z_ref_u=10.0, # 风速观测高度(m)
z_ref_tq=10.0, # 温湿观测高度(m)
u_ref=7.0, # 10m风速(m/s)
t_air_c=27.0, # 气温(°C)
rh=0.82, # 相对湿度(0~1)，注意不是82
sst_c=28.5, # 海温(°C)
p_hpa=1010.0, # 气压(hPa)
lat_deg=20.0 # 纬度(度)
)

# 2) 运行模型
out = run_model(inp, z_max=50.0, dz=0.2)

# 3) 读取结果
print("蒸发波导高度 EDH(m):", out["evaporation_duct_height_m"])
print("状态量:", out["state"]) # u*、t*、q*、z0m/z0h/z0q、L
print("廓线点数:", len(out["profile"]["z"]))
print("前5个 M(z):", out["profile"]["M"][:5])
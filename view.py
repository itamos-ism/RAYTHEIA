import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.cm as cm
import matplotlib.colors as mcolors
import numpy as np

# -----------------------------------------------------------
# 配置参数
# -----------------------------------------------------------
FILENAME = 'octree_check_000.txt'  # 您的输出文件名
LOG_SCALE = True                   # 使用对数色标
COLORMAP = 'viridis'               # 颜色风格
SLICE_Z_INDEX = 0                  # 切片位置 (Z=0)
DENSITY_THRESHOLD = 1e-30          # 密度阈值
# -----------------------------------------------------------

# 1. 读取数据
try:
    df = pd.read_csv(FILENAME, skipinitialspace=True)
except FileNotFoundError:
    print(f"错误: 找不到文件 {FILENAME}")
    exit()

# [新增步骤] 处理 Fortran 的 T/F 布尔值
# Fortran 输出的 L1 格式是 T 或 F，Pandas 读进来可能是字符串，需要转换
# 如果读进来已经是 bool 则会自动忽略这一步，如果是字符串则转为 bool
if df['Is_Padding'].dtype == 'object':
    df['Is_Padding'] = df['Is_Padding'].astype(str).str.strip() == 'T'

# 2. 准备绘图
fig, ax = plt.subplots(figsize=(12, 10))

# 3. 设置颜色映射
valid_densities = df[df['Density'] > DENSITY_THRESHOLD]['Density']

if len(valid_densities) == 0:
    print("警告: 数据密度过低")
    vmax = 1.0; vmin = 1e-5
else:
    vmax = valid_densities.max()
    vmin = valid_densities.min()

print(f"Plotting Density Range: [{vmin:.3e}, {vmax:.3e}]")

if LOG_SCALE:
    norm = mcolors.LogNorm(vmin=vmin, vmax=vmax)
else:
    norm = mcolors.Normalize(vmin=vmin, vmax=vmax)

cmap = cm.get_cmap(COLORMAP)

# 4. 绘制矩形
count = 0
for index, row in df.iterrows():
    
    # 切片过滤
    if row['Global_Z'] <= SLICE_Z_INDEX < (row['Global_Z'] + row['Size_Phys']):
        count += 1
        
        dens = row['Density']
        is_padding = row['Is_Padding'] # 读取 Padding 标志
        
        # 逻辑：只要不是 Padding 且密度大于阈值，就上色
        # 如果是 Padding，或者密度极低，就画斜线
        is_vacuum = (dens <= DENSITY_THRESHOLD)
        
        # 注意：这里我们依然画出 Padding 区域（作为背景参考），
        # 但稍后我们会限制坐标轴只显示物理区域。
        
        if is_vacuum or is_padding:
            face_color = 'none' 
            edge_color = 'lightgray'
            line_style = '--'
            alpha = 0.3
            hatch = '///' 
        else:
            face_color = cmap(norm(dens))
            edge_color = 'black'
            line_style = '-'
            alpha = 1.0
            hatch = None

        rect = patches.Rectangle(
            (row['Global_X'], row['Global_Y']), 
            row['Size_Phys'], row['Size_Phys'], 
            linewidth=0.5, 
            edgecolor=edge_color, 
            facecolor=face_color, 
            linestyle=line_style,
            alpha=alpha,
            hatch=hatch
        )
        ax.add_patch(rect)

print(f"Drawn {count} boxes on Z={SLICE_Z_INDEX} plane.")

# 5. 添加 Colorbar
sm = cm.ScalarMappable(cmap=cmap, norm=norm)
sm.set_array([]) 
cbar = fig.colorbar(sm, ax=ax, fraction=0.046, pad=0.04)
cbar.set_label(f"Density ({'Log' if LOG_SCALE else 'Linear'})", rotation=270, labelpad=20)

# ===========================================================
# [核心修改] 6. 自动计算并设置物理网格范围
# ===========================================================

# 筛选出所有 "非 Padding" 的节点 (即物理计算域内的节点)
phys_df = df[~df['Is_Padding']]

if not phys_df.empty:
    # 物理域的最小坐标
    min_x = phys_df['Global_X'].min()
    min_y = phys_df['Global_Y'].min()
    
    # 物理域的最大坐标 = 节点起点 + 节点大小
    max_x = (phys_df['Global_X'] + phys_df['Size_Phys']).max()
    max_y = (phys_df['Global_Y'] + phys_df['Size_Phys']).max()
    
    print(f"Auto-detected Physical Domain: X[{min_x:.2f}, {max_x:.2f}], Y[{min_y:.2f}, {max_y:.2f}]")
    
    # 设置坐标轴范围
    ax.set_xlim(min_x, max_x)
    ax.set_ylim(min_y, max_y)
else:
    print("警告: 没有检测到物理节点 (全都是 Padding?)，使用默认全范围")
    # 回退方案
    max_x_all = df['Global_X'].max() + df['Size_Phys'].max()
    max_y_all = df['Global_Y'].max() + df['Size_Phys'].max()
    ax.set_xlim(df['Global_X'].min(), max_x_all)
    ax.set_ylim(df['Global_Y'].min(), max_y_all)

ax.set_aspect('equal')

plt.title(f"Linear Octree Density Map (Rank {FILENAME.split('_')[-1].split('.')[0]}, Z=0 slice)")
plt.xlabel("Global X")
plt.ylabel("Global Y")
plt.tight_layout()
plt.show()
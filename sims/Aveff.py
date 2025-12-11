import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.legend_handler import HandlerTuple
from matplotlib import gridspec
from scipy.ndimage import gaussian_filter1d
import sys

import numpy as np

def plot_bins(ax, rho, values, color, ls, alpha, label=None, num_bins=100, sigma=1.0):
    """
    对单个 level 的数据，按密度分箱，绘制物种丰度的均值线与范围阴影。

    Parameters:
        ax: matplotlib axis
        rho: array-like, 氢密度 (nH)
        values: array-like, 物种相对丰度
        color: str, 线条/阴影颜色
        ls: str, 线型
        alpha: float, 阴影透明度
        label: str, 图例标签
        num_bins: int, 分箱数量
        sigma: float, 高斯平滑参数
    """
    # 过滤有效数据
    valid = (rho > 0) & (values > 0)
    if not np.any(valid):
        return
    
    log_rho = np.log10(rho[valid])
    values_valid = values[valid]

    # 创建分箱
    bins = np.linspace(log_rho.min(), log_rho.max(), num_bins)
    bin_indices = np.digitize(log_rho, bins)

    mean_vals, max_vals, min_vals, centers = [], [], [], []

    for i in range(1, len(bins)):
        mask = (bin_indices == i)
        if np.sum(mask) == 0:
            continue
        subset = values_valid[mask]
        mean_vals.append(np.mean(subset))
        max_vals.append(np.max(subset))
        min_vals.append(np.min(subset))
        # 计算对数空间 bin 中心，转回线性
        center_log = 0.5 * (bins[i-1] + bins[i])
        centers.append(10**center_log)

    if len(centers) == 0:
        return

    # 平滑
    min_smooth = gaussian_filter1d(min_vals, sigma)
    max_smooth = gaussian_filter1d(max_vals, sigma)
    mean_smooth = gaussian_filter1d(mean_vals, sigma)

    # 绘图
    ax.fill_between(centers, min_smooth, max_smooth, color=color, alpha=alpha, linewidth=0)

    plot_kwargs = {
        'color': color,
        'linestyle': ls,
        'linewidth': 2
    }
    if label is not None:
        plot_kwargs['label'] = label

    ax.plot(centers, mean_smooth, **plot_kwargs)

# 输入参数
num_bins = 100
sigma = 1.0

# 设置全局样式
plt.style.use('default')
plt.rcParams.update({
    'font.size': 24,
    'axes.labelsize': 16,
    'xtick.labelsize': 16,
    'ytick.labelsize': 16,
    'legend.fontsize': 16,
    'axes.linewidth': 1.5
})

# 读取数据
data0 = np.loadtxt('Aveff.dat')
data1 = np.loadtxt('base.dat')
data6 = np.loadtxt('Aveff_level6.dat')


rho = data0[:,0]
Aveff0 = data0[:,1]
Aveff1 = data1[:,1]
Aveff6 = data6[:,1]

error0 = np.abs(Aveff0-Aveff6)/Aveff6
error1 = np.abs(Aveff1-Aveff6)/Aveff6

n_H = np.logspace(0, 4, 400)
Av_fun = 0.05 * np.exp(1.6 * n_H**0.12)

# 设置颜色和样式配置
level_styles = {
    'level0': {'color': '#2E86AB', 'ls': '--', 'alpha': 0.2, 'label': r'$AMR$'},
    'level1': {'color': '#A23B72', 'ls': '--', 'alpha': 0.2, 'label': r'$Base$'},
    'level2': {'color': '#F18F01', 'ls': '--', 'alpha': 0.2, 'label': r'${\cal N}_{\mathrm {rays}}=192$'},
    'level3': {'color': '#C73E1D', 'ls': '--', 'alpha': 0.2, 'label': r'${\cal N}_{\mathrm {rays}}=768$'}
}

fig, ax1 = plt.subplots(figsize=(8, 6))

# 处理所有level数据
for level, data in zip(['level0', 'level1'], [error0, error1]):
    style = level_styles[level]
    plot_bins(
        ax=ax1,
        rho=rho,
        values=data,
        color=style['color'],
        ls=style['ls'],
        alpha=style['alpha'],
        label=style['label'],
        num_bins=num_bins,
        sigma=sigma
    )

ax1.set_ylabel(r'$e_{A_\mathrm {V,eff}}$', fontsize=22)
ax1.set_xlabel(r'$n_\mathrm{H}$ [cm$^{-3}$]')
ax1.set_xscale('log')
ax1.set_yscale('log')
ax1.set_xlim(1,10000)
ax1.set_ylim(2E-8,3)
ax1.legend(loc='lower right', frameon=False, ncol=1)

plt.tight_layout()
plt.savefig('Aveff.png', dpi=300, bbox_inches='tight')

# 多站协作 NLoS 感知：非稀疏参数化贝叶斯模型与实验路线

> 本文档根据现阶段实验思路重新整理技术路线。接收信号模型、多接收站几何关系以及发射站位置和定时信息不完美的场景保持不变；主要变化是取消位置网格稀疏建模，不再引入支撑变量或伯努利—高斯稀疏先验，而是直接对每条路径的散射系数、角度、测量时延、目标位置、UE 位置和时钟偏差建立参数化分层贝叶斯模型。

## 1. 研究问题与技术路线

### 1.1 研究场景

考虑发射站位置与定时信息不完美条件下的多站协作上行 ISAC 环境感知：

- 单天线 UE 发送基于 OFDM 的上行 ISAC 导频信号；
- $N_{\mathrm{Rx}}\geq2$ 个配置 ULA 的接收站共享观测，并由融合中心联合处理；
- UE 与接收站之间存在未知公共定时偏差 $\Delta t$；为使几何状态各分量的尺度一致，定义其对应的距离偏差 $\beta:=c\Delta t$，后续 LM 与贝叶斯几何推断均以 $\beta$ 为待估变量；
- 融合中心仅掌握 UE 的粗位置；
- 场景为 NLoS 探测，不存在可直接用于同步和位置校准的 LoS 路径；
- 每条反射路径同时受到目标位置、UE 位置与公共时钟偏差的影响。

设环境中共有 $K$ 个目标。对于第 $n$ 个接收站，散射系数、归一化空间相位和归一化时延相位分别写为

$$
\boldsymbol h_n=[h_{n,1},\ldots,h_{n,K}]^T,
$$

$$
\boldsymbol\theta_n=[\theta_{n,1},\ldots,\theta_{n,K}]^T,
\qquad
\boldsymbol\tau_n=[\tau_{n,1},\ldots,\tau_{n,K}]^T.
$$

记物理入射角和秒单位测量时延分别为 $\phi_{n,k}$ 和 $t_{n,k}$。定义

$$
\theta_{n,k}=\kappa_a\cos\phi_{n,k},
\qquad
\tau_{n,k}=\kappa_\tau t_{n,k},
$$

其中 $\kappa_a=2\pi d/\lambda$，$d$ 为阵元间距，$\lambda$ 为载波波长；当导频子载波按间隔 $D$ 均匀抽取时，$\kappa_\tau=2\pi D\Delta f$。因此 $\theta_{n,k}$ 和 $\tau_{n,k}$ 均为无量纲相位参数，而 $t_{n,k}$ 是包含公共时钟偏差的物理测量时延。距离偏差 $\beta$ 的单位为米，最终再通过 $\Delta t=\beta/c$ 恢复秒单位定时偏差。

### 1.2 技术路线

1. 保留原天线域—子载波域接收信号模型，但直接以 $K$ 条物理路径为参数，不再将感兴趣区域离散成位置网格，也不构造位置域稀疏反射向量。
2. 分别对 $N_{\mathrm{Rx}}$ 个接收站的观测执行过采样 2D-FFT，获得每条路径归一化空间相位 $\theta_{n,k}$ 和归一化时延相位 $\tau_{n,k}$ 的粗估计，并反归一化为物理入射角 $\phi_{n,k}$ 和测量时延 $t_{n,k}$。
3. 利用多接收站物理角度射线的联合交汇残差和接收站间测量距离差完成跨站路径关联，再以全部已关联射线的加权最小二乘交点作为目标位置粗估计；双接收站方法是该过程的特例。
4. 将物理入射角和测量时延关于目标位置、UE 位置及公共距离偏差 $\beta=c\Delta t$ 的非线性几何方程组成残差向量，使用 LM 算法获得 $\{\boldsymbol p_k^{(0)}\}_{k=1}^{K}$、$\boldsymbol p_{\mathrm{UE}}^{(0)}$ 和 $\beta^{(0)}$。
5. 以上述粗估计为几何变量高斯先验的均值，建立关于原始观测 $\boldsymbol y_n$ 的分层贝叶斯模型：
   - $\boldsymbol h_n$ 采用零均值复高斯先验，协方差为对角矩阵；
   - 噪声精度 $\gamma_n$ 采用 Gamma 先验；
   - 归一化空间相位和时延相位通过 Dirac delta 因子与目标位置、UE 位置和距离偏差保持确定性几何关系；
   - 目标位置、UE 位置和距离偏差采用以粗估计为均值的高斯先验。
6. 第一部分执行连续线谱贝叶斯估计：首轮以 2D-FFT 粗估计为展开点，后续轮次以几何模块反馈的 $\boldsymbol g_n(\boldsymbol\mu_z^{(l-1)})$ 为展开点，对参数化路径矩阵作一阶 Taylor 近似，并依次更新角度—时延相位、散射系数和噪声精度的替代后验分布。
7. 第二部分将 Dirac 因子内部的几何映射在上一轮 $\boldsymbol z^{(l-1)}$ 处作一阶 Taylor 展开，把线谱模块输出的角度—时延联合消息等价转换成关于几何参数 $\boldsymbol z$ 的高斯信息消息，并与几何先验融合得到后验分布；随后在更新后的几何后验均值处重新线性化 Dirac 因子，分别按 $\boldsymbol j\mathbf C_p\boldsymbol j^T$ 计算各时延和角度分量的标量方差，将几何后验转换成传回线谱模块的对角 covariance Gaussian 消息，并令下一轮线谱 Taylor 展开点等于 $\boldsymbol g_n(\boldsymbol\mu_z^{(l)})$。

整体信息流为

```text
多 Rx 导频观测 Y_1, ..., Y_NRx
    ↓ 去除已知导频
过采样 2D-FFT
    ↓
归一化相位粗估计 theta 和 tau
    ↓ 反归一化
物理入射角 phi 和测量时延 t
    ↓ 多站射线联合交汇 + 距离差关联
目标位置粗估计
    ↓ 几何非线性方程组
LM：目标位置 + UE 位置 + 距离偏差 beta 粗估计
    ↓ 作为几何高斯先验均值
非稀疏参数化分层贝叶斯模型
    ├─ 第一部分：Taylor 线性化连续线谱贝叶斯估计
    │    └─ 输出 psi_n=[theta_n^T,tau_n^T]^T 的联合概率消息
    └─ 第二部分：线性化联合 Dirac 几何因子并计算 z 的高斯后验
         └─ 将几何后验转换为传回线谱模块的高斯消息，并更新 Taylor 展开点
    ↺ 两部分之间迭代传递消息
    ↓ 收敛
目标位置、UE 位置、距离偏差及其后验不确定性（最终换算为时钟偏差）
```

## 2. 接收信号模型

### 2.1 归一化阵列响应与时延响应

第 $n$ 个接收站配置 $N_a$ 阵元 ULA。将物理入射角 $\phi$ 映射为归一化空间相位

$$
\theta=\kappa_a\cos\phi,
\qquad
\kappa_a=\frac{2\pi d}{\lambda}.
$$

采用正指数阵列约定时，单位范数阵列响应写为

$$
\boldsymbol a(\theta)
=\frac{1}{\sqrt{N_a}}
\left[
1,e^{j\theta},\ldots,
e^{j(N_a-1)\theta}
\right]^T
=\frac{1}{\sqrt{N_a}}
\left[e^{jm\theta}\right]_{m=0}^{N_a-1}.
$$

当 $d=\lambda/2$ 时，有 $\theta=\pi\cos\phi$。类似地，设第 $q$ 个导频对应的子载波索引为 $qD$，忽略可吸收到复散射系数中的公共起始频率相位，并定义归一化时延相位

$$
\tau=\kappa_\tau t,
\qquad
\kappa_\tau=2\pi D\Delta f,
$$

则单位范数时延响应为

$$
\boldsymbol b(\tau)
=\frac{1}{\sqrt{N_p}}
\left[
1,e^{-j\tau},\ldots,
e^{-j(N_p-1)\tau}
\right]^T
=\frac{1}{\sqrt{N_p}}
\left[e^{-jq\tau}\right]_{q=0}^{N_p-1}.
$$

因此，$\boldsymbol a(\theta)$ 和 $\boldsymbol b(\tau)$ 的指数中不再显式出现阵元间距、波长、子载波间隔或物理时延，所有物理尺度均包含在无量纲参数 $\theta$ 和 $\tau$ 中。单位范数因子用于固定 $\boldsymbol h_n$ 的幅度尺度；若代码中省略该因子，则必须相应调整散射系数及其先验协方差。

当前 MATLAB 数据生成代码使用负指数阵列约定。若概率模型采用上述 $e^{jm\theta}$ 形式，则实现时应使用 $\theta_{\mathrm{code}}=-\theta$，或将角度维 IFFT 改为 FFT；两种处理物理等价，但整条链路必须固定一种符号约定。

### 2.2 目标、UE 与接收站几何

记第 $k$ 个目标的二维位置为 $\boldsymbol p_k\in\mathbb R^2$，UE 位置为 $\boldsymbol p_{\mathrm{UE}}\in\mathbb R^2$，第 $n$ 个接收站的位置和 ULA 轴向单位向量分别为 $\boldsymbol p_{\mathrm{Rx},n}$ 和 $\boldsymbol e_n$。若融合中心掌握的 UE 粗位置为 $\widetilde{\boldsymbol p}_{\mathrm{UE}}$，则发射站位置偏差可写成 $\Delta\boldsymbol p_{\mathrm{UE}}=\boldsymbol p_{\mathrm{UE}}-\widetilde{\boldsymbol p}_{\mathrm{UE}}$；估计 $\boldsymbol p_{\mathrm{UE}}$ 与估计该位置偏差是等价参数化。

归一化空间频率与归一化时延频率共同构成二维线谱频率。它们与目标位置、UE 位置、接收站位置及公共距离偏差 $\beta=c\Delta t$ 的关系直接写为

$$
\theta_{n,k}
=
g_{\theta,n,k}(\boldsymbol p_k)
:=
\kappa_a
\frac{
\boldsymbol e_n^T
(\boldsymbol p_k-\boldsymbol p_{\mathrm{Rx},n})
}{
\|\boldsymbol p_k-\boldsymbol p_{\mathrm{Rx},n}\|
},
$$

$$
\begin{aligned}
\tau_{n,k}
&=
g_{\tau,n,k}(
\boldsymbol p_k,\boldsymbol p_{\mathrm{UE}},\beta)
\\
&:=
\frac{\kappa_\tau}{c}
\left[
\|\boldsymbol p_k-\boldsymbol p_{\mathrm{UE}}\|
+
\|\boldsymbol p_k-\boldsymbol p_{\mathrm{Rx},n}\|
+
\beta
\right].
\end{aligned}
$$

其中，\(g_{\theta,n,k}\) 和 \(g_{\tau,n,k}\) 仅作为上述直接几何映射的简写，以便后续书写 Dirac 因子及其 Jacobian；不再引入物理入射角或物理测量时延作为中间变量。

为直接构造 LM 和第 5.2 节中的几何 Jacobian，定义目标到第 $n$ 个接收站以及目标到 UE 的位移、距离和单位方向向量

$$
\begin{aligned}
\boldsymbol d_{n,k}^{\mathrm R}
&:=\boldsymbol p_k-\boldsymbol p_{\mathrm{Rx},n},
&\rho_{n,k}^{\mathrm R}
&:=\|\boldsymbol d_{n,k}^{\mathrm R}\|,
&\boldsymbol u_{n,k}^{\mathrm R}
&:=\frac{\boldsymbol d_{n,k}^{\mathrm R}}{\rho_{n,k}^{\mathrm R}},
\\
\boldsymbol d_k^{\mathrm U}
&:=\boldsymbol p_k-\boldsymbol p_{\mathrm{UE}},
&\rho_k^{\mathrm U}
&:=\|\boldsymbol d_k^{\mathrm U}\|,
&\boldsymbol u_k^{\mathrm U}
&:=\frac{\boldsymbol d_k^{\mathrm U}}{\rho_k^{\mathrm U}}.
\end{aligned}
$$

则角度映射的非零导数为

$$
\frac{\partial g_{\theta,n,k}}
{\partial\boldsymbol p_k^T}
=
\frac{\kappa_a}{\rho_{n,k}^{\mathrm R}}
\boldsymbol e_n^T
\left(
\mathbf I_2-
\boldsymbol u_{n,k}^{\mathrm R}
\left(\boldsymbol u_{n,k}^{\mathrm R}\right)^T
\right),
$$

并且 $\partial g_{\theta,n,k}/\partial\boldsymbol p_j^T=\boldsymbol 0^T$（$j\ne k$）、$\partial g_{\theta,n,k}/\partial\boldsymbol p_{\mathrm{UE}}^T=\boldsymbol 0^T$、$\partial g_{\theta,n,k}/\partial\beta=0$。时延映射的非零导数为

$$
\frac{\partial g_{\tau,n,k}}
{\partial\boldsymbol p_k^T}
=
\frac{\kappa_\tau}{c}
\left(
\boldsymbol u_k^{\mathrm U}
+\boldsymbol u_{n,k}^{\mathrm R}
\right)^T,
$$

$$
\frac{\partial g_{\tau,n,k}}
{\partial\boldsymbol p_{\mathrm{UE}}^T}
=
-\frac{\kappa_\tau}{c}
\left(\boldsymbol u_k^{\mathrm U}\right)^T,
\qquad
\frac{\partial g_{\tau,n,k}}
{\partial\beta}
=\frac{\kappa_\tau}{c},
$$

且 $\partial g_{\tau,n,k}/\partial\boldsymbol p_j^T=\boldsymbol 0^T$（$j\ne k$）。后续所有几何 Jacobian 均由这些行向量按照 $\boldsymbol z=[\boldsymbol p_1^T,\ldots,\boldsymbol p_K^T,\boldsymbol p_{\mathrm{UE}}^T,\beta]^T$ 的列顺序组装，不再使用数值差分。这里对距离偏差 $\beta$ 求导，而不对秒单位的 $\Delta t$ 求导，因此 $\boldsymbol z$ 中的位置坐标与定时偏差分量均采用米作为单位。

### 2.3 非稀疏参数化观测模型

定义第 $n$ 个接收站的参数化路径矩阵

$$
\boldsymbol\Phi_n(\boldsymbol\theta_n,\boldsymbol\tau_n)
=\left[
\boldsymbol a(\theta_{n,1})\otimes\boldsymbol b(\tau_{n,1}),
\ldots,
\boldsymbol a(\theta_{n,K})\otimes\boldsymbol b(\tau_{n,K})
\right].
$$

按“导频子载波索引 $q$ 为行、阵元索引 $m$ 为列”定义去除已知导频后的观测矩阵 $\mathbf Y_n\in\mathbb C^{N_p\times N_a}$。若原始接收样本为 $R_n[q,m]$、已知非零导频为 $X[q]$，则 $Y_n[q,m]=R_n[q,m]/X[q]$，并且

$$
\mathbf Y_n
=
\sum_{k=1}^{K}
h_{n,k}
\boldsymbol b(\tau_{n,k})
\boldsymbol a^T(\theta_{n,k})
+\mathbf W_n.
$$

采用 MATLAB 列优先顺序进行向量化，即子载波索引 $q$ 变化最快：

$$
\boldsymbol y_n
:=\operatorname{vec}(\mathbf Y_n)
=\operatorname{reshape}(\mathbf Y_n,M,1),
\qquad
M=N_pN_a.
$$

利用 $\operatorname{vec}(\boldsymbol b\boldsymbol a^T)=\boldsymbol a\otimes\boldsymbol b$，第 $n$ 个接收站的向量观测为

$$
\boldsymbol y_n
=\boldsymbol\Phi_n(\boldsymbol\theta_n,\boldsymbol\tau_n)
\boldsymbol h_n+\boldsymbol w_n,
$$

其中

$$
\boldsymbol w_n\sim
\mathcal{CN}(\boldsymbol 0,\gamma_n^{-1}\mathbf I_M),
\qquad
\boldsymbol w_n=\operatorname{vec}(\mathbf W_n).
$$

实现中各变量的固定维度为

$$
\mathbf Y_n\in\mathbb C^{N_p\times N_a},
\quad
\boldsymbol y_n\in\mathbb C^M,
\quad
\boldsymbol\Phi_n\in\mathbb C^{M\times K},
\quad
\boldsymbol h_n\in\mathbb C^K,
\quad
\boldsymbol\psi_n\in\mathbb R^{2K},
\quad
\boldsymbol z\in\mathbb R^{2K+3}.
$$

该模型的未知向量 $\boldsymbol h_n$ 只有 $K$ 个元素，每一列直接对应一条物理路径；这里不再引入大小为位置网格数 $N_g$ 的稀疏向量。

## 3. 分层贝叶斯概率模型

### 3.1 观测似然

给定 $\boldsymbol h_n$、$\boldsymbol\theta_n$、$\boldsymbol\tau_n$ 和噪声精度 $\gamma_n$，观测概率分布为

$$
p(\boldsymbol y_n\mid
\boldsymbol h_n,\boldsymbol\theta_n,
\boldsymbol\tau_n,\gamma_n)
=\mathcal{CN}\!\left(
\boldsymbol y_n;
\boldsymbol\Phi_n(\boldsymbol\theta_n,\boldsymbol\tau_n)
\boldsymbol h_n,
\gamma_n^{-1}\mathbf I_M
\right).
$$

本文采用的圆对称复高斯密度约定为

$$
\mathcal{CN}\!\left(
\boldsymbol y;\boldsymbol\mu,
\gamma^{-1}\mathbf I_M
\right)
=
\frac{\gamma^M}{\pi^M}
\exp\!\left(
-\gamma\|\boldsymbol y-\boldsymbol\mu\|^2
\right),
$$

该约定决定第 5.1.2 节中 Gamma 后验的 shape 增量为 $M$，而不是 $M/2$。

不同接收站噪声条件独立时

$$
p(\mathcal Y\mid\mathcal H,\boldsymbol\Theta,
\boldsymbol T,\boldsymbol\gamma)
=\prod_{n=1}^{N_{\mathrm{Rx}}}
p(\boldsymbol y_n\mid
\boldsymbol h_n,\boldsymbol\theta_n,
\boldsymbol\tau_n,\gamma_n),
$$

其中 $\mathcal Y=\{\boldsymbol y_n\}_n$、$\mathcal H=\{\boldsymbol h_n\}_n$、$\boldsymbol\Theta=\{\boldsymbol\theta_n\}_n$、$\boldsymbol T=\{\boldsymbol\tau_n\}_n$。

### 3.2 散射系数的复高斯先验

对每个接收站的路径散射系数采用零均值复高斯先验

$$
p(\boldsymbol h_n)
=\mathcal{CN}(\boldsymbol h_n;
\boldsymbol 0,\boldsymbol\Sigma_{h,n}),
$$

其中

$$
\boldsymbol\Sigma_{h,n}
=\operatorname{diag}
\left(\sigma_{h,n,1}^2,\ldots,
\sigma_{h,n,K}^2\right).
$$

对角协方差表示同一接收站内各路径散射系数先验独立，但允许不同路径具有不同先验功率。不同接收站的散射系数可条件独立建模，因为同一目标在不同接收站上的复反射系数不要求相干一致：

$$
p(\mathcal H)=\prod_n p(\boldsymbol h_n).
$$

当前模型不采用伯努利—高斯先验，也不引入路径支撑变量。若后续需要学习 $\boldsymbol\Sigma_{h,n}$，应另行指定其超先验或确定的经验贝叶斯更新规则。

### 3.3 噪声精度的 Gamma 先验

采用 shape-rate 参数化：

$$
p(\gamma_n)
=\operatorname{Ga}(\gamma_n;a_{\gamma,0},b_{\gamma,0})
=\frac{b_{\gamma,0}^{a_{\gamma,0}}}
{\Gamma(a_{\gamma,0})}
\gamma_n^{a_{\gamma,0}-1}
e^{-b_{\gamma,0}\gamma_n}.
$$

于是

$$
p(\boldsymbol\gamma)=\prod_n p(\gamma_n).
$$

这里的 $\gamma_n$ 是复高斯噪声的精度，即单个复观测维度的方差为 $\gamma_n^{-1}$。实验实现中必须固定 Gamma 分布采用 shape-rate 而非 shape-scale 参数化。

### 3.4 归一化空间相位与位置之间的确定性概率因子

归一化空间相位由目标位置和接收站几何唯一确定，因此写为 Dirac delta 条件分布

$$
p(\boldsymbol\Theta\mid\mathbf P)
=\prod_{n=1}^{N_{\mathrm{Rx}}}
\prod_{k=1}^{K}
\delta\!\left(
\theta_{n,k}-g_{\theta,n,k}(\boldsymbol p_k)
\right),
$$

其中

$$
\mathbf P=[\boldsymbol p_1,\ldots,\boldsymbol p_K].
$$

该 delta 因子表示严格的几何映射，而不是“角度测量误差为零”的独立观测模型。$\theta_{n,k}$ 的不确定性最终来自目标位置后验的不确定性，其对应的物理入射角由 $\phi_{n,k}=\arccos(\theta_{n,k}/\kappa_a)$ 恢复。

### 3.5 归一化时延相位与几何参数之间的确定性概率因子

归一化时延相位由目标位置、UE 位置和公共距离偏差共同决定：

$$
p(\boldsymbol T\mid
\mathbf P,\boldsymbol p_{\mathrm{UE}},\beta)
=\prod_{n=1}^{N_{\mathrm{Rx}}}
\prod_{k=1}^{K}
\delta\!\left(
\tau_{n,k}
-g_{\tau,n,k}(\boldsymbol p_k,
\boldsymbol p_{\mathrm{UE}},\beta)
\right).
$$

因为 $\tau_{n,k}=\kappa_\tau t_{n,k}$ 已包含公共距离偏差对应的归一化相位，所以在 $\boldsymbol b(\tau_{n,k})$ 中不能再次叠加 $\kappa_\tau\beta/c$。秒单位测量时延由 $t_{n,k}=\tau_{n,k}/\kappa_\tau$ 恢复。

### 3.6 几何变量的高斯先验

令 2D-FFT、多站关联和 LM 给出的粗估计为

$$
\boldsymbol p_k^{(0)},\qquad
\boldsymbol p_{\mathrm{UE}}^{(0)},\qquad
\beta^{(0)}=c\Delta t^{(0)}.
$$

目标位置先验为

$$
p(\mathbf P)
=\prod_{k=1}^{K}
\mathcal N\!\left(
\boldsymbol p_k;
\boldsymbol p_k^{(0)},\mathbf C_{p,k}^{(0)}
\right).
$$

UE 位置与距离偏差先验分别为

$$
p(\boldsymbol p_{\mathrm{UE}})
=\mathcal N\!\left(
\boldsymbol p_{\mathrm{UE}};
\boldsymbol p_{\mathrm{UE}}^{(0)},
\mathbf C_{\mathrm{UE}}^{(0)}
\right),
$$

$$
p(\beta)
=\mathcal N\!\left(
\beta;\beta^{(0)},
\sigma_{\beta,0}^{2}
\right).
$$

若原始配置给出的是秒单位时钟偏差先验 $\Delta t\sim\mathcal N(\Delta t^{(0)},\sigma_{\Delta t,0}^2)$，则变量变换 $\beta=c\Delta t$ 给出 $\beta^{(0)}=c\Delta t^{(0)}$ 和 $\sigma_{\beta,0}^2=c^2\sigma_{\Delta t,0}^2$。因此该重参数化不改变物理先验，只改变几何状态的单位。

将全部几何变量堆叠为

$$
\boldsymbol z=
\left[
\boldsymbol p_1^T,\ldots,\boldsymbol p_K^T,
\boldsymbol p_{\mathrm{UE}}^T,\beta
\right]^T,
$$

则可以把全部几何先验合并为

$$
p(\boldsymbol z)
=\mathcal N(\boldsymbol z;
\boldsymbol z^{(0)},\mathbf C_{z,0}),
$$

其中

$$
\boldsymbol z^{(0)}
=\left[
(\boldsymbol p_1^{(0)})^T,\ldots,
(\boldsymbol p_K^{(0)})^T,
(\boldsymbol p_{\mathrm{UE}}^{(0)})^T,
\beta^{(0)}
\right]^T,
$$

$$
\mathbf C_{z,0}
=\operatorname{blkdiag}\!\left(
\mathbf C_{p,1}^{(0)},\ldots,
\mathbf C_{p,K}^{(0)},
\mathbf C_{\mathrm{UE}}^{(0)},
\sigma_{\beta,0}^{2}
\right).
$$

这些协方差决定贝叶斯精化对粗估计的信任程度。协方差过小会使后验难以修正粗估计偏差；协方差过大则会削弱 2D-FFT 和 LM 初始化的稳定作用。

### 3.7 完整联合概率分布

完整联合概率分布可写为

$$
\begin{aligned}
&p(\mathcal Y,\mathcal H,\boldsymbol\Theta,
\boldsymbol T,\boldsymbol\gamma,
\mathbf P,\boldsymbol p_{\mathrm{UE}},\beta)
\\
&=\left\{\prod_{n=1}^{N_{\mathrm{Rx}}}
p(\boldsymbol y_n\mid
\boldsymbol h_n,\boldsymbol\theta_n,
\boldsymbol\tau_n,\gamma_n)
p(\boldsymbol h_n)p(\gamma_n)
\right\}
\\
&\quad\times
p(\boldsymbol\Theta\mid\mathbf P)
p(\boldsymbol T\mid
\mathbf P,\boldsymbol p_{\mathrm{UE}},\beta)
p(\mathbf P)
p(\boldsymbol p_{\mathrm{UE}})
p(\beta).
\end{aligned}
$$

目标是计算或近似后验分布

$$
p(\mathcal H,\boldsymbol\gamma,
\mathbf P,\boldsymbol p_{\mathrm{UE}},\beta
\mid\mathcal Y).
$$

## 4. 2D-FFT 与 LM 粗初始化

### 4.1 过采样 2D-FFT

按第 2.3 节约定将 $\boldsymbol y_n$ 恢复为 $\mathbf Y_n=\operatorname{reshape}(\boldsymbol y_n,N_p,N_a)$。设时延维和角度维的过采样倍数分别为 $\rho_\tau\geq1$ 和 $\rho_\theta\geq1$，FFT 点数固定为

$$
N_\tau^{\mathrm{FFT}}
=\max\!\left(N_p,\left\lceil\rho_\tau N_p\right\rceil\right),
\qquad
N_\theta^{\mathrm{FFT}}
=\max\!\left(N_a,\left\lceil\rho_\theta N_a\right\rceil\right).
$$

在 $\boldsymbol a(\theta)\propto[e^{jm\theta}]_m$ 和 $\boldsymbol b(\tau)\propto[e^{-jq\tau}]_q$ 的符号约定下，过采样二维变换定义为

$$
Z_n[p,\ell]
:=
\sum_{q=0}^{N_p-1}
\sum_{m=0}^{N_a-1}
Y_n[q,m]
\exp\!\left(
j\frac{2\pi qp}{N_\tau^{\mathrm{FFT}}}
\right)
\exp\!\left(
-j\frac{2\pi m\ell}{N_\theta^{\mathrm{FFT}}}
\right),
$$

其中 $p=0,\ldots,N_\tau^{\mathrm{FFT}}-1$，$\ell=-\lfloor N_\theta^{\mathrm{FFT}}/2\rfloor,\ldots,N_\theta^{\mathrm{FFT}}-1-\lfloor N_\theta^{\mathrm{FFT}}/2\rfloor$。因此代码中对子载波维执行 IFFT、对阵元维执行 FFT，并仅对阵元维结果执行 `fftshift`；常数归一化不影响峰值位置。搜索功率谱为

$$
P_n[p,\ell]=|Z_n[p,\ell]|^2.
$$

若允许的方向余弦和秒单位时延区间分别为 $[u_{\min},u_{\max}]$ 和 $[t_{\min},t_{\max}]$，则搜索集合明确为

$$
\Omega_{\mathrm{search}}
:=
\left\{
(p,\ell):
u_{\min}
\leq
\frac{2\pi\ell}{\kappa_aN_\theta^{\mathrm{FFT}}}
\leq u_{\max},
\quad
t_{\min}
\leq
\frac{2\pi p}{\kappa_\tau N_\tau^{\mathrm{FFT}}}
\leq t_{\max}
\right\}.
$$

本版本将路径数 $K$ 作为已知输入，并使用二维非极大值抑制递归提取 $K$ 个峰值。默认抑制半径为

$$
R_\tau
=
\max\!\left(
1,\left\lceil
\frac{N_\tau^{\mathrm{FFT}}}{N_p}
\right\rceil
\right),
\qquad
R_\theta
=
\max\!\left(
1,\left\lceil
\frac{N_\theta^{\mathrm{FFT}}}{N_a}
\right\rceil
\right).
$$

定义角度循环距离

$$
d_{\mathrm{circ}}(\ell_1,\ell_2)
:=
\min_{v\in\mathbb Z}
\left|
\ell_1-\ell_2+vN_\theta^{\mathrm{FFT}}
\right|,
$$

则第 $r$ 个已选峰周围的抑制邻域为

$$
\mathcal N_r
:=
\left\{
(p,\ell):
|p-p_{n,r}|\leq R_\tau,
\ d_{\mathrm{circ}}(\ell,\ell_{n,r})\leq R_\theta
\right\}.
$$

第 $k$ 个整数峰递归定义为

$$
(p_{n,k},\ell_{n,k})
=
\arg\max_{(p,\ell)\in
\Omega_{\mathrm{search}}\setminus
\bigcup_{r=1}^{k-1}\mathcal N_r}
P_n[p,\ell],
$$

其中 $\Omega_{\mathrm{search}}$ 由给定的方向余弦范围和时延范围确定。为进行亚栅格修正，令 $L_n[p,\ell]=\log\max(P_n[p,\ell],\varepsilon)$。沿某一维取相邻三个值 $L_{-1},L_0,L_{+1}$ 时，抛物线峰值偏移为

$$
\delta
=
\operatorname{clip}_{[-1/2,1/2]}
\left[
\frac{1}{2}
\frac{L_{-1}-L_{+1}}
{L_{-1}-2L_0+L_{+1}}
\right].
$$

当抛物线公式的分母为零或非有限值时固定令 $\delta=0$。分别沿时延维和角度维计算 $\delta_{p,n,k}$ 和 $\delta_{\ell,n,k}$；角度维邻点按周期索引，时延维边界峰不插值并令 $\delta_{p,n,k}=0$。最终的归一化相位粗估计为

$$
\widehat\theta_{n,k}^{(0)}
=
\frac{2\pi(\ell_{n,k}+\delta_{\ell,n,k})}
{N_\theta^{\mathrm{FFT}}},
\qquad
\widehat\tau_{n,k}^{(0)}
=
\frac{2\pi(p_{n,k}+\delta_{p,n,k})}
{N_\tau^{\mathrm{FFT}}}.
$$

对应的方向余弦和秒单位测量时延为

$$
\widehat u_{n,k}^{(0)}
=\frac{\widehat\theta_{n,k}^{(0)}}{\kappa_a},
\qquad
\widehat t_{n,k}^{(0)}
=\frac{\widehat\tau_{n,k}^{(0)}}{\kappa_\tau}.
$$

实现前将 $\widehat u_{n,k}^{(0)}$ 裁剪到 $[-1,1]$。当 $d\leq\lambda/2$ 时，主值区间内不存在空间混叠；均匀抽取导频时，无模糊时延范围为 $0\leq t<1/(D\Delta f)$。本版本假设所有路径均位于上述无模糊区间。

当前实验固定 $\rho_\tau=\rho_\theta=8$。零填充降低 FFT 参数栅格的量化误差，但不会突破阵列孔径和有效带宽决定的物理分辨率。

现有 `CoarseEstimation2DFFT.m` 针对代码中的负指数阵列约定沿两个维度均使用 IFFT，并直接输出方向余弦和秒单位时延。生成新代码时必须统一到本文的正指数阵列约定：阵元维改用 FFT，函数内部首先输出 $\widehat\theta_{n,k}^{(0)}$ 和 $\widehat\tau_{n,k}^{(0)}$，再按上式提供方向余弦和秒单位时延。

### 4.2 多接收站路径关联与目标位置粗估计

设 $N_{\mathrm{Rx}}\geq2$，每个接收站均从 2D-FFT 结果中提取 $K$ 个峰。仅由 ULA 方向余弦存在前后向二义性，因此必须给定每个接收站的可见半平面。令 $\boldsymbol e_n^\perp$ 为与阵列轴向 $\boldsymbol e_n$ 正交的单位向量，$s_{n,i}\in\{-1,+1\}$ 表示由已知可见半平面确定的分支，则第 $n$ 个接收站第 $i$ 个峰对应的全局单位射线方向为

$$
\widehat{\boldsymbol u}_{n,i}
=
\widehat u_{n,i}^{(0)}\boldsymbol e_n
+s_{n,i}
\sqrt{1-(\widehat u_{n,i}^{(0)})^2}
\boldsymbol e_n^\perp.
$$

定义与该射线正交的投影矩阵

$$
\boldsymbol\Pi_{n,i}^{\perp}
:=
\mathbf I_2-
\widehat{\boldsymbol u}_{n,i}
\widehat{\boldsymbol u}_{n,i}^T.
$$

对包含每个接收站一个峰的候选多站路径组合 $\boldsymbol\iota=(i_1,\ldots,i_{N_{\mathrm{Rx}}})$，其目标位置候选值定义为全部射线的加权最小二乘交点

$$
\boldsymbol p_{\boldsymbol\iota}^{\mathrm{ray}}
=
\arg\min_{\boldsymbol p\in\mathbb R^2}
\sum_{n=1}^{N_{\mathrm{Rx}}}
w_{\theta,n,i_n}
\left\|
\boldsymbol\Pi_{n,i_n}^{\perp}
(\boldsymbol p-\boldsymbol p_{\mathrm{Rx},n})
\right\|^2,
$$

其中 $w_{\theta,n,i_n}>0$ 为由 FFT 峰值强度或峰曲率给出的方向置信权重；尚未标定权重时固定取 $w_{\theta,n,i_n}=1$。定义

$$
\mathbf Q_{\boldsymbol\iota}
=
\sum_{n=1}^{N_{\mathrm{Rx}}}
w_{\theta,n,i_n}\boldsymbol\Pi_{n,i_n}^{\perp},
\qquad
\boldsymbol q_{\boldsymbol\iota}
=
\sum_{n=1}^{N_{\mathrm{Rx}}}
w_{\theta,n,i_n}\boldsymbol\Pi_{n,i_n}^{\perp}
\boldsymbol p_{\mathrm{Rx},n},
$$

则 $\boldsymbol p_{\boldsymbol\iota}^{\mathrm{ray}}$ 通过求解 $\mathbf Q_{\boldsymbol\iota}\boldsymbol p_{\boldsymbol\iota}^{\mathrm{ray}}=\boldsymbol q_{\boldsymbol\iota}$ 得到，不显式计算矩阵逆。若 $\mathbf Q_{\boldsymbol\iota}$ 不满秩或条件数超过阈值 $T_{\mathrm{cond}}$，则该组合因射线几何退化而无效。有效组合还必须满足全部前向射线约束

$$
\alpha_{n,\boldsymbol\iota}
:=
\widehat{\boldsymbol u}_{n,i_n}^T
(\boldsymbol p_{\boldsymbol\iota}^{\mathrm{ray}}-\boldsymbol p_{\mathrm{Rx},n})
>0,\qquad n=1,\ldots,N_{\mathrm{Rx}}.
$$

定义第 $n$ 个接收站第 $i$ 个峰对应的测量总路程 $\widehat d_{n,i}:=c\widehat t_{n,i}^{(0)}$。为了使关联代价不依赖参考接收站的选择，对所有接收站对 $1\leq n<m\leq N_{\mathrm{Rx}}$ 定义距离差残差

$$
\begin{aligned}
r_{\tau,nm}(\boldsymbol\iota)
&:=
\left(\widehat d_{n,i_n}-\widehat d_{m,i_m}\right)
\\
&\quad-
\left(
\|\boldsymbol p_{\boldsymbol\iota}^{\mathrm{ray}}-\boldsymbol p_{\mathrm{Rx},n}\|
-
\|\boldsymbol p_{\boldsymbol\iota}^{\mathrm{ray}}-\boldsymbol p_{\mathrm{Rx},m}\|
\right).
\end{aligned}
$$

由于同一目标在各接收站处的测量总路程满足 $d_{n,k}=\|\boldsymbol p_k-\boldsymbol p_{\mathrm{UE}}\|+\|\boldsymbol p_k-\boldsymbol p_{\mathrm{Rx},n}\|+\beta$，上述距离差同时消除了 UE—目标距离和公共距离偏差 $\beta$，因此可以在尚未知晓 UE 精确位置和定时偏差时用于路径关联。按照固定的字典序将所有 $r_{\tau,nm}(\boldsymbol\iota)$ 堆叠成 $\boldsymbol r_{\tau,\boldsymbol\iota}\in\mathbb R^{L_\tau}$，其中 $L_\tau=N_{\mathrm{Rx}}(N_{\mathrm{Rx}}-1)/2$，并定义多站组合代价

$$
\begin{aligned}
C_{\boldsymbol\iota}
&:=
\frac{1}{T_{\perp}^{2}}
\sum_{n=1}^{N_{\mathrm{Rx}}}
w_{\theta,n,i_n}
\left\|
\boldsymbol\Pi_{n,i_n}^{\perp}
(\boldsymbol p_{\boldsymbol\iota}^{\mathrm{ray}}-\boldsymbol p_{\mathrm{Rx},n})
\right\|^2
\\
&\quad+
\boldsymbol r_{\tau,\boldsymbol\iota}^T
\mathbf W_{\tau,\boldsymbol\iota}
\boldsymbol r_{\tau,\boldsymbol\iota},
\end{aligned}
$$

其中 $T_{\perp}>0$ 为射线垂距归一化尺度，$\mathbf W_{\tau,\boldsymbol\iota}\succeq0$ 为距离差残差权重。若没有可靠的 FFT 时延误差协方差，则固定 $\mathbf W_{\tau,\boldsymbol\iota}=T_d^{-2}\mathbf I_{L_\tau}$，其中 $T_d$ 为距离差归一化尺度。全部成对距离差存在代数冗余且统计相关；当可获得各峰的时延方差时，应按固定堆叠顺序构造距离差协方差，并以其广义逆或正则化逆作为 $\mathbf W_{\tau,\boldsymbol\iota}$。

设 $\mathfrak S$ 为选出的 $K$ 个多站路径组合的集合。严格的一对一多站关联定义为

$$
\begin{aligned}
\mathfrak S^\star
&=
\arg\min_{\mathfrak S}
\sum_{\boldsymbol\iota\in\mathfrak S}
C_{\boldsymbol\iota}
\\
\text{s.t.}\quad
&|\mathfrak S|=K,
\qquad
\{i_n:\boldsymbol\iota\in\mathfrak S\}
=\{1,\ldots,K\},
\quad \forall n.
\end{aligned}
$$

该约束保证每个接收站的每个峰恰好分配给一个目标。$N_{\mathrm{Rx}}=2$ 时，上式退化为标准二分图指派，可直接使用 Hungarian 算法；$N_{\mathrm{Rx}}>2$ 时属于多维指派问题。对于较小的 $K$ 和 $N_{\mathrm{Rx}}$，可以使用整数规划或分支定界获得全局解并作为算法验证基准。

为避免正式实验中枚举 $K^{N_{\mathrm{Rx}}}$ 个组合，本文固定采用锚接收站对驱动的可扩展初始化。若 $r$ 为锚参考站，第 $k$ 个当前目标假设在该站对应峰索引为 $i_{r,k}$，则将待加入接收站 $n$ 的第 $j$ 个峰分配给目标 $k$ 的代价明确写为

$$
\begin{aligned}
C_{k,j}^{(n)}
&:=
\frac{w_{\theta,n,j}}{T_{\perp}^{2}}
\left\|
\boldsymbol\Pi_{n,j}^{\perp}
(\boldsymbol p_k^{\mathrm{ray}}-\boldsymbol p_{\mathrm{Rx},n})
\right\|^2
\\
&\quad+
\frac{1}{T_d^2}
\left[
(\widehat d_{n,j}-\widehat d_{r,i_{r,k}})
-
(\|\boldsymbol p_k^{\mathrm{ray}}-\boldsymbol p_{\mathrm{Rx},n}\|
-\|\boldsymbol p_k^{\mathrm{ray}}-\boldsymbol p_{\mathrm{Rx},r}\|)
\right]^2.
\end{aligned}
$$

该接收站的一对一指派为

$$
\pi_n^\star
=
\arg\min_{\pi\in\mathcal P_K}
\sum_{k=1}^{K}C_{k,\pi(k)}^{(n)}.
$$

具体流程为：

1. 对每一对接收站构造二站组合代价，排除射线退化、反向射线或残差超过门限的组合，并使用 Hungarian 算法求得完整二站指派；选择归一化总代价最小的可行接收站对作为锚接收站对，总代价相同时依次按站间基线更长和接收站固定标识的字典序确定，以避免输入排列改变锚站选择。
2. 以锚接收站对的 $K$ 个关联结果建立目标标签和初始位置。其余接收站按照与锚参考站的基线长度降序加入，基线相同时按固定标识排序。对于每个待加入接收站，构造其 $K$ 个峰与当前 $K$ 个目标假设之间的 $K\times K$ 代价矩阵，代价由射线垂距和相对于锚参考站的距离差残差共同组成，再使用 Hungarian 算法完成一对一关联。
3. 每加入一个接收站后，使用当前已关联的全部射线重新求解上述加权最小二乘交点，并按射线垂距、距离差残差、前向约束和矩阵条件数重新检查可行性。
4. 全部接收站加入后，保持锚接收站对的标签不变，依次固定其他接收站的关联并重新指派每个非锚接收站，随后更新全部目标位置；重复该逐站重指派过程，直至所有路径标签不再变化或达到最大关联迭代次数 $I_{\mathrm{assoc}}$。

若任一阶段不存在完整有限指派，任一目标少于两个非平行有效射线，或最终残差超过 $T_{\perp}$ 和 $T_d$ 对应门限，则当前固定 $K$ 模型的初始化失败并返回诊断，不静默补造路径。将该可扩展流程得到的关联记为 $\widehat{\mathfrak S}$；当 $N_{\mathrm{Rx}}>2$ 时，它是严格多维指派解 $\mathfrak S^\star$ 的近似，不宣称具有全局最优性。最终令 $\boldsymbol p_k^{\mathrm{ray}}$ 等于第 $k$ 个已关联多站射线组的加权最小二乘交点，作为第 4.3 节联合 LM 的目标位置初值。

### 4.3 LM 几何粗估计

将多站关联后的归一化相位观测堆叠为

$$
\widehat{\boldsymbol\psi}
:=
\begin{bmatrix}
\widehat{\boldsymbol\psi}_1\\
\vdots\\
\widehat{\boldsymbol\psi}_{N_{\mathrm{Rx}}}
\end{bmatrix},
\qquad
\widehat{\boldsymbol\psi}_n
=
\begin{bmatrix}
\widehat{\boldsymbol\theta}_n^{(0)}\\
\widehat{\boldsymbol\tau}_n^{(0)}
\end{bmatrix},
$$

并定义完整几何映射与残差

$$
\boldsymbol g(\boldsymbol z)
:=
\begin{bmatrix}
\boldsymbol g_1(\boldsymbol z)\\
\vdots\\
\boldsymbol g_{N_{\mathrm{Rx}}}(\boldsymbol z)
\end{bmatrix},
\qquad
\boldsymbol r(\boldsymbol z)
:=
\widehat{\boldsymbol\psi}
-\boldsymbol g(\boldsymbol z).
$$

以 $\boldsymbol z_{\mathrm{init}}=[(\boldsymbol p_1^{\mathrm{ray}})^T,\ldots,(\boldsymbol p_K^{\mathrm{ray}})^T,\widetilde{\boldsymbol p}_{\mathrm{UE}}^T,0]^T$ 为初值，其中最后一个零表示米单位距离偏差的初值 $\beta_{\mathrm{init}}=0$。令 $\mathbf W_0\succeq0$ 为粗相位观测权重矩阵；若尚未根据 FFT 峰曲率标定测量协方差，则固定 $\mathbf W_0=\mathbf I_{2N_{\mathrm{Rx}}K}$。LM 粗估计定义为

$$
\boldsymbol z^{(0)}
=
\arg\min_{\boldsymbol z}
\frac{1}{2}
\boldsymbol r^T(\boldsymbol z)
\mathbf W_0
\boldsymbol r(\boldsymbol z).
$$

在第 $s$ 次 LM 迭代中，使用第 2.2 节的解析导数组装

$$
\mathbf H^{(s)}
:=
\left.
\frac{\partial\boldsymbol g(\boldsymbol z)}
{\partial\boldsymbol z^T}
\right|_{\boldsymbol z=\boldsymbol z^{(s)}},
\qquad
\mathbf B^{(s)}
:=\mathbf H^{(s)T}\mathbf W_0\mathbf H^{(s)}.
$$

由于几何状态中的位置坐标和距离偏差 $\beta$ 均以米为单位，LM 不再面对位置与秒单位时钟偏差之间的基本量纲失配。为改善不同几何方向上的数值条件，仍采用对角阻尼缩放矩阵

$$
\mathbf D^{(s)}
=
\operatorname{diag}\!\left(
\max\!(
\operatorname{diag}(\mathbf B^{(s)}),
\epsilon_{\mathrm{LM}}\mathbf 1_{2K+3}
)
\right),
$$

其中最大值逐元素计算。LM 试探步由线性方程

$$
\left(
\mathbf B^{(s)}
+\lambda_s\mathbf D^{(s)}
\right)
\Delta\boldsymbol z^{(s)}
=
\mathbf H^{(s)T}
\mathbf W_0
\boldsymbol r(\boldsymbol z^{(s)})
$$

求得，并令 $\boldsymbol z_{\mathrm{trial}}=\boldsymbol z^{(s)}+\Delta\boldsymbol z^{(s)}$。定义 $J(\boldsymbol z)=\boldsymbol r^T(\boldsymbol z)\mathbf W_0\boldsymbol r(\boldsymbol z)/2$；若 $J(\boldsymbol z_{\mathrm{trial}})<J(\boldsymbol z^{(s)})$，则接受试探点并令 $\lambda_{s+1}=\max(\lambda_s/\nu_\downarrow,\lambda_{\min})$，否则拒绝试探点并令 $\lambda_s\leftarrow\min(\nu_\uparrow\lambda_s,\lambda_{\max})$ 后重新求解，其中 $\nu_\downarrow>1$ 且 $\nu_\uparrow>1$。

LM 在满足下列任一条件时停止：

$$
\frac{\|\Delta\boldsymbol z^{(s)}\|}
{\|\boldsymbol z^{(s)}\|+\epsilon_{\mathrm{LM}}}
<\varepsilon_{z,\mathrm{LM}},
\qquad
\frac{|J(\boldsymbol z^{(s+1)})-J(\boldsymbol z^{(s)})|}
{J(\boldsymbol z^{(s)})+\epsilon_{\mathrm{LM}}}
<\varepsilon_{J,\mathrm{LM}},
$$

或达到最大迭代次数 $S_{\max}$。只有当最终 $\mathbf H^{(s)}$ 对待估方向满列秩时，才将 LM 输出解释为局部唯一粗估计。对于 $N_{\mathrm{Rx}}$ 个接收站，观测数必须满足必要条件 $2N_{\mathrm{Rx}}K\geq2K+3$，但该计数条件不能替代 Jacobian 秩检查；接收站或目标几何退化时，即使观测数充足仍可能不可辨识。

现有 `estimate_bias_by_coherence.m` 仍采用“固定射线交点，只更新 UE 位置和距离偏差”的简化实现。在该简化模型中，每个目标只提供一个独立的 UE—偏差约束，因此至少需要 3 个几何分散的目标。生成与本文概率模型一致的新代码时，应实现上述同时包含 $\boldsymbol\theta$ 和 $\boldsymbol\tau$ 残差、联合更新全部 $\boldsymbol z$ 的 LM 公式，并令其最后一个状态分量为 $\beta$。

## 5. 两部分贝叶斯推断与消息传递

第 $l$ 轮推断由连续线谱更新和几何消息传递两部分组成。第 5.1 节在当前局部观测模型固定的条件下，按照 $\boldsymbol\psi_n\rightarrow\boldsymbol h_n\rightarrow\gamma_n$ 的顺序进行 KL 坐标最小化，并向几何模块输出角度—时延联合外信息；第 5.2 节利用 Dirac 几何因子更新几何后验，再把几何不确定性转换为下一轮线谱输入消息，同时通过完整几何后验均值的非线性映射确定下一轮局部模型的参考点。因此，一轮中的替代后验坐标更新与两轮之间的局部模型重构属于两个不同层次。

### 5.1 第一部分：连续线谱贝叶斯估计

连续线谱模块直接推断第 \(n\) 个接收站对应的散射系数 \(\boldsymbol h_n\)、联合角度—时延向量 \(\boldsymbol\psi_n\) 以及噪声精度 \(\gamma_n\)，其中

$$
\boldsymbol\psi_n
:=
\begin{bmatrix}
\boldsymbol\theta_n\\
\boldsymbol\tau_n
\end{bmatrix}
\in\mathbb R^{2K}.
$$

由于 $\boldsymbol\Phi_n(\boldsymbol\theta_n,\boldsymbol\tau_n)$ 的各列由阵列指数响应和时延指数响应构成，观测均值 $\boldsymbol\Phi_n(\boldsymbol\theta_n,\boldsymbol\tau_n)\boldsymbol h_n$ 关于 $\boldsymbol\psi_n$ 是非线性的，因而不能直接得到 $\boldsymbol\psi_n$ 的高斯共轭更新。为此，第 $l$ 轮线谱推断在参考点

$$
\boldsymbol\psi_n^{(l-1)}
:=
\begin{bmatrix}
\boldsymbol\theta_n^{(l-1)}\\
\boldsymbol\tau_n^{(l-1)}
\end{bmatrix}
$$

附近对路径矩阵作一阶 Taylor 近似。该近似只处理路径矩阵关于角度—时延参数的非线性：推断 $\boldsymbol\psi_n$ 时，将近似观测模型整理成关于 $\boldsymbol\psi_n$ 的仿射形式；推断 $\boldsymbol h_n$ 时，则将同一近似整理成关于 $\boldsymbol h_n$ 的线性形式。两种等价整理及其期望项分别在第 5.1.1 节和第 5.1.2 节中给出。

首轮参考点由过采样 2D-FFT 粗估计给出，即 $\boldsymbol\psi_n^{(0)}=[\widehat{\boldsymbol\theta}_{n,\mathrm{FFT}}^T,\widehat{\boldsymbol\tau}_{n,\mathrm{FFT}}^T]^T$；从第 2 轮开始，$\boldsymbol\psi_n^{(l-1)}$ 由上一轮几何后验均值的非线性映射 $\boldsymbol g_n(\boldsymbol\mu_z^{(l-1)})$ 给出，具体反馈过程见第 5.2 节。这里的参考点与 $q^{(l-1)}(\boldsymbol\psi_n)$ 的后验均值不是同一概念。

第 $l$ 轮迭代结束时的替代后验分解为

$$
q_n^{(l)}(\boldsymbol h_n,\boldsymbol\psi_n,\gamma_n)
=
q^{(l)}(\boldsymbol h_n)
q^{(l)}(\boldsymbol\psi_n)
q^{(l)}(\gamma_n),
$$

其中 \(q^{(l)}(\boldsymbol\psi_n)\) 保留角度与时延之间的后验相关性。

为使第 1 轮 $\boldsymbol\psi_n$ 更新具有确定的散射系数和噪声精度初值，先令

$$
q^{(0)}(\gamma_n)
=p(\gamma_n),
\qquad
\bar\gamma_n^{(0)}
=\frac{a_{\gamma,0}}{b_{\gamma,0}},
$$

并在 2D-FFT 展开点构造

$$
\boldsymbol\Phi_n^{(0)}
=
\boldsymbol\Phi_n\!\left(
\widehat{\boldsymbol\theta}_n^{(0)},
\widehat{\boldsymbol\tau}_n^{(0)}
\right).
$$

散射系数初始替代后验取为条件线性高斯更新

$$
q^{(0)}(\boldsymbol h_n)
=
\mathcal{CN}\!\left(
\boldsymbol h_n;
\boldsymbol m_{h,n}^{(0)},
\mathbf C_{h,n}^{(0)}
\right),
$$

其中

$$
\left(\mathbf C_{h,n}^{(0)}\right)^{-1}
=
\bar\gamma_n^{(0)}
\boldsymbol\Phi_n^{(0)H}
\boldsymbol\Phi_n^{(0)}
+\boldsymbol\Sigma_{h,n}^{-1},
$$

$$
\boldsymbol m_{h,n}^{(0)}
=
\mathbf C_{h,n}^{(0)}
\bar\gamma_n^{(0)}
\boldsymbol\Phi_n^{(0)H}
\boldsymbol y_n.
$$



第 $l$ 轮中来自几何模块的输入消息写成信息形式

$$
m_{\mathrm{geo}\rightarrow\psi_n}^{(l)}(\boldsymbol\psi_n)
\propto
\exp\!\left(
-\frac{1}{2}\boldsymbol\psi_n^T
\boldsymbol\Lambda_{\psi,n}^{\mathrm{in},(l)}
\boldsymbol\psi_n
+\boldsymbol\eta_{\psi,n}^{\mathrm{in},(l)\,T}
\boldsymbol\psi_n
\right).
$$

首次执行线谱估计时，尚无上一轮几何后验，但已经具有第 3.6 节的几何先验。令

$$
\boldsymbol g_n^{(0)}
=\boldsymbol g_n(\boldsymbol z^{(0)}),
$$

并按第 5.2 节分别计算每条路径在 $\boldsymbol z^{(0)}$ 处的 $\boldsymbol j_{\tau,n,k}^{(0)}$ 和 $\boldsymbol j_{\theta,n,k}^{(0)}$。定义初始局部 ext 协方差

$$
\mathbf C_{p,n,k}^{\mathrm{ext},(0)}
:=
\mathbf S_k\mathbf C_{z,0}\mathbf S_k^T,
$$

则初始时延和角度消息方差分别为标量

$$
\sigma_{\tau,n,k}^{2,\mathrm{in},(1)}
=
\boldsymbol j_{\tau,n,k}^{(0)}
\mathbf C_{p,n,k}^{\mathrm{ext},(0)}
\boldsymbol j_{\tau,n,k}^{(0)T},
$$

$$
\sigma_{\theta,n,k}^{2,\mathrm{in},(1)}
=
\boldsymbol j_{\theta,n,k}^{(0)}
\mathbf C_{p,n,k}^{\mathrm{ext},(0)}
\boldsymbol j_{\theta,n,k}^{(0)T}.
$$

在角度、时延及不同路径相互独立的假设下，将 $p(\boldsymbol z)=\mathcal N(\boldsymbol z;\boldsymbol z^{(0)},\mathbf C_{z,0})$ 通过初始线性化 Dirac 因子推前，得到第 1 轮几何输入消息

$$
m_{\mathrm{geo}\rightarrow\psi_n}^{(1)}(\boldsymbol\psi_n)
\approx
\mathcal N\!\left(
\boldsymbol\psi_n;
\boldsymbol g_n^{(0)},
\mathbf V_{\psi,n}^{\mathrm{in},(1)}
\right).
$$

因此

$$
\mathbf V_{\psi,n}^{\mathrm{in},(1)}
=
\operatorname{diag}\!\left(
\sigma_{\theta,n,1}^{2,\mathrm{in},(1)},\ldots,
\sigma_{\theta,n,K}^{2,\mathrm{in},(1)},
\sigma_{\tau,n,1}^{2,\mathrm{in},(1)},\ldots,
\sigma_{\tau,n,K}^{2,\mathrm{in},(1)}
\right)
+\epsilon_{\mathrm{geo}}\mathbf I_{2K},
\qquad
\boldsymbol\mu_{\psi,n}^{\mathrm{in},(1)}
=\boldsymbol g_n^{(0)},
$$

$$
\boldsymbol\Lambda_{\psi,n}^{\mathrm{in},(1)}
=
\left(\mathbf V_{\psi,n}^{\mathrm{in},(1)}\right)^{-1},
\qquad
\boldsymbol\eta_{\psi,n}^{\mathrm{in},(1)}
=
\boldsymbol\Lambda_{\psi,n}^{\mathrm{in},(1)}
\boldsymbol\mu_{\psi,n}^{\mathrm{in},(1)}.
$$

其中 $\epsilon_{\mathrm{geo}}>0$ 为实验配置中固定记录的小正则量；不再默认把首轮几何消息设为零。

将下面两个小节构造的一阶近似似然记为 $\widetilde p^{(l)}(\boldsymbol y_n\mid\boldsymbol h_n,\boldsymbol\psi_n,\gamma_n)$。在第 $l$ 轮固定展开点 $\boldsymbol\psi_n^{(l-1)}$ 后，局部目标后验为

$$
\begin{aligned}
&\widetilde p_n^{(l)}(
\boldsymbol h_n,\boldsymbol\psi_n,\gamma_n
\mid\boldsymbol y_n)
\\
&\quad\propto
\widetilde p^{(l)}(
\boldsymbol y_n\mid
\boldsymbol h_n,\boldsymbol\psi_n,\gamma_n)
p(\boldsymbol h_n)p(\gamma_n)
m_{\mathrm{geo}\rightarrow\psi_n}^{(l)}(\boldsymbol\psi_n),
\end{aligned}
$$

该局部目标后验仅在当前 Taylor 邻域内近似原始非线性后验。替代后验通过求解

$$
q_n^{(l),\star}
=
\arg\min_{q_n^{(l)}}
\operatorname{KL}\!\left[
q_n^{(l)}(\boldsymbol h_n,\boldsymbol\psi_n,\gamma_n)
\middle\|
\widetilde p_n^{(l)}(
\boldsymbol h_n,\boldsymbol\psi_n,\gamma_n
\mid\boldsymbol y_n)
\right]
$$

获得。由于

$$
\ln\widetilde p_n^{(l)}(\boldsymbol y_n)
=
\mathcal L_n^{(l)}(q_n^{(l)})
+
\operatorname{KL}\!\left[
q_n^{(l)}
\middle\|
\widetilde p_n^{(l)}(\cdot\mid\boldsymbol y_n)
\right],
$$

最小化 KL 散度等价于最大化局部证据下界

$$
\begin{aligned}
\mathcal L_n^{(l)}(q_n^{(l)})
={}&
\mathbb E_{q_n^{(l)}}\!\left[
\ln\widetilde p_n^{(l)}(
\boldsymbol y_n,\boldsymbol h_n,
\boldsymbol\psi_n,\gamma_n)
\right]
\\
&-
\mathbb E_{q_n^{(l)}}\!\left[
\ln q_n^{(l)}(
\boldsymbol h_n,\boldsymbol\psi_n,\gamma_n)
\right].
\end{aligned}
$$

根据一轮内“先更新 $\boldsymbol\psi_n$，再更新 $\boldsymbol h_n$，最后更新 $\gamma_n$”的顺序，三个坐标最优条件分别为

$$
\ln q^{(l)}(\boldsymbol\psi_n)
=
\mathbb E_{q^{(l-1)}(\boldsymbol h_n)
q^{(l-1)}(\gamma_n)}\!\left[
\ln\widetilde p_n^{(l)}(
\boldsymbol y_n,\boldsymbol h_n,
\boldsymbol\psi_n,\gamma_n)
\right]
+\mathrm{const}.
$$

$$
\ln q^{(l)}(\boldsymbol h_n)
=
\mathbb E_{q^{(l)}(\boldsymbol\psi_n)
q^{(l-1)}(\gamma_n)}\!\left[
\ln\widetilde p_n^{(l)}(
\boldsymbol y_n,\boldsymbol h_n,
\boldsymbol\psi_n,\gamma_n)
\right]
+\mathrm{const}.
$$

$$
\ln q^{(l)}(\gamma_n)
=
\mathbb E_{q^{(l)}(\boldsymbol\psi_n)
q^{(l)}(\boldsymbol h_n)}\!\left[
\ln\widetilde p_n^{(l)}(
\boldsymbol y_n,\boldsymbol h_n,
\boldsymbol\psi_n,\gamma_n)
\right]
+\mathrm{const}.
$$

记上一轮噪声精度的后验均值为

$$
\bar\gamma_n^{(l-1)}
:=
\mathbb E_{q^{(l-1)}(\gamma_n)}[\gamma_n].
$$

#### 5.1.1 关于 \(\boldsymbol\psi_n\) 的一阶近似及其推断

为了获得关于 \(\boldsymbol\psi_n\) 的条件线性模型，沿用第 2.3 节定义的路径矩阵 \(\boldsymbol\Phi_n(\boldsymbol\theta_n,\boldsymbol\tau_n)\)。

第 $l$ 轮采用本节开头定义的 $\boldsymbol\psi_n^{(l-1)}$ 作为展开点。由于 $\boldsymbol\Phi_n$ 的第 $k$ 列只依赖 $(\theta_{n,k},\tau_{n,k})$，定义两个 $M\times K$ 的列导数矩阵

$$
\boldsymbol\Phi_{\theta,n}^{(l-1)}
:=
\left[
\left.
\frac{\partial\boldsymbol\phi_{n,1}}
{\partial\theta_{n,1}}\right|_{l-1},
\ldots,
\left.
\frac{\partial\boldsymbol\phi_{n,K}}
{\partial\theta_{n,K}}\right|_{l-1}
\right],
$$

$$
\boldsymbol\Phi_{\tau,n}^{(l-1)}
:=
\left[
\left.
\frac{\partial\boldsymbol\phi_{n,1}}
{\partial\tau_{n,1}}\right|_{l-1},
\ldots,
\left.
\frac{\partial\boldsymbol\phi_{n,K}}
{\partial\tau_{n,K}}\right|_{l-1}
\right],
$$

其中

$$
\boldsymbol\phi_{n,k}
=
\boldsymbol a(\theta_{n,k})
\otimes
\boldsymbol b(\tau_{n,k}).
$$

于是路径矩阵的一阶 Taylor 近似为

$$
\begin{aligned}
\boldsymbol\Phi_n(\boldsymbol\theta_n,\boldsymbol\tau_n)
\approx{}&
\boldsymbol\Phi_n^{(l-1)}
+
\boldsymbol\Phi_{\theta,n}^{(l-1)}
\operatorname{diag}\!\left(
\boldsymbol\theta_n-\boldsymbol\theta_n^{(l-1)}
\right)
\\
&+
\boldsymbol\Phi_{\tau,n}^{(l-1)}
\operatorname{diag}\!\left(
\boldsymbol\tau_n-\boldsymbol\tau_n^{(l-1)}
\right),
\end{aligned}
$$

其中

$$
\boldsymbol\Phi_n^{(l-1)}
=
\boldsymbol\Phi_n\!\left(
\boldsymbol\theta_n^{(l-1)},
\boldsymbol\tau_n^{(l-1)}
\right).
$$

若按照一般矩阵对向量求导，\(\partial\boldsymbol\Phi_n/\partial\boldsymbol\theta_n\)会形成三阶张量。上述列导数定义利用了各路径列之间的参数独立性，使一阶近似中的两个右乘对角矩阵具有明确维度。

令

$$
\mathbf D_a
=
\operatorname{diag}(0,1,\ldots,N_a-1),
\qquad
\mathbf D_b
=
\operatorname{diag}(0,1,\ldots,N_p-1),
$$

则归一化导向向量的解析导数为

$$
\frac{\partial\boldsymbol a(\theta)}{\partial\theta}
=
j\mathbf D_a\boldsymbol a(\theta),
\qquad
\frac{\partial\boldsymbol b(\tau)}{\partial\tau}
=
-j\mathbf D_b\boldsymbol b(\tau),
$$

从而

$$
\frac{\partial\boldsymbol\phi_{n,k}}
{\partial\theta_{n,k}}
=
\left(
j\mathbf D_a\boldsymbol a(\theta_{n,k})
\right)
\otimes
\boldsymbol b(\tau_{n,k}),
$$

$$
\frac{\partial\boldsymbol\phi_{n,k}}
{\partial\tau_{n,k}}
=
\boldsymbol a(\theta_{n,k})
\otimes
\left(
-j\mathbf D_b\boldsymbol b(\tau_{n,k})
\right).
$$

因此代码中两个列导数矩阵应直接按列构造为

$$
\left[
\boldsymbol\Phi_{\theta,n}^{(l-1)}
\right]_{:,k}
=
\left(
j\mathbf D_a
\boldsymbol a(\theta_{n,k}^{(l-1)})
\right)
\otimes
\boldsymbol b(\tau_{n,k}^{(l-1)}),
$$

$$
\left[
\boldsymbol\Phi_{\tau,n}^{(l-1)}
\right]_{:,k}
=
\boldsymbol a(\theta_{n,k}^{(l-1)})
\otimes
\left(
-j\mathbf D_b
\boldsymbol b(\tau_{n,k}^{(l-1)})
\right),
\qquad k=1,\ldots,K.
$$

给定散射系数，定义联合线性化矩阵

$$
\mathbf G_n^{(l-1)}(\boldsymbol h_n)
:=
\left[
\boldsymbol\Phi_{\theta,n}^{(l-1)}
\operatorname{diag}(\boldsymbol h_n),
\;
\boldsymbol\Phi_{\tau,n}^{(l-1)}
\operatorname{diag}(\boldsymbol h_n)
\right].
$$

因此 $\mathbf G_n^{(l-1)}(\boldsymbol h_n)\in\mathbb C^{M\times2K}$。

观测模型关于 \(\boldsymbol\psi_n\) 的一阶近似为

$$
\boldsymbol y_n
\approx
\boldsymbol\Phi_n^{(l-1)}\boldsymbol h_n
+
\mathbf G_n^{(l-1)}(\boldsymbol h_n)
\left(
\boldsymbol\psi_n-\boldsymbol\psi_n^{(l-1)}
\right)
+
\boldsymbol w_n.
$$

因此，条件于 \(\boldsymbol h_n\) 时，局部似然可写为

$$
\widetilde p^{(l)}\!\left(
\boldsymbol y_n\mid
\boldsymbol h_n,\boldsymbol\psi_n,\gamma_n
\right)
=
\mathcal{CN}\!\left(
\boldsymbol y_n;
\boldsymbol\Phi_n^{(l-1)}\boldsymbol h_n
+
\mathbf G_n^{(l-1)}(\boldsymbol h_n)
\left(
\boldsymbol\psi_n-\boldsymbol\psi_n^{(l-1)}
\right),
\gamma_n^{-1}\mathbf I_M
\right).
$$

首先固定上一轮的 \(q^{(l-1)}(\boldsymbol h_n)\) 和 \(q^{(l-1)}(\gamma_n)\)。为推导 \(q^{(l)}(\boldsymbol\psi_n)\)，定义

$$
\boldsymbol c_n^{(l-1)}(\boldsymbol h_n)
:=
\boldsymbol\Phi_n^{(l-1)}\boldsymbol h_n
-
\mathbf G_n^{(l-1)}(\boldsymbol h_n)
\boldsymbol\psi_n^{(l-1)}.
$$

其中 $\boldsymbol c_n^{(l-1)}(\boldsymbol h_n)\in\mathbb C^M$。

于是关于 \(\boldsymbol\psi_n\) 的局部观测模型可以写成

$$
\boldsymbol y_n
\approx
\mathbf G_n^{(l-1)}(\boldsymbol h_n)
\boldsymbol\psi_n
+
\boldsymbol c_n^{(l-1)}(\boldsymbol h_n)
+
\boldsymbol w_n.
$$

根据 KL 坐标最优条件，保留所有与 \(\boldsymbol\psi_n\) 有关的项，有

$$
\begin{aligned}
\ln q^{(l)}(\boldsymbol\psi_n)
={}&
\mathbb E_{
q^{(l-1)}(\boldsymbol h_n)
q^{(l-1)}(\gamma_n)}
\left[
\ln\widetilde p^{(l)}(
\boldsymbol y_n\mid
\boldsymbol h_n,\boldsymbol\psi_n,\gamma_n)
\right]
\\
&+
\ln m_{\mathrm{geo}\rightarrow\psi_n}^{(l)}
(\boldsymbol\psi_n)
+
\mathrm{const}
\\
={}&
-\bar\gamma_n^{(l-1)}
\mathbb E_{q^{(l-1)}(\boldsymbol h_n)}
\left[
\left\|
\boldsymbol y_n
-
\boldsymbol c_n^{(l-1)}(\boldsymbol h_n)
-
\mathbf G_n^{(l-1)}(\boldsymbol h_n)
\boldsymbol\psi_n
\right\|^2
\right]
\\
&-
\frac{1}{2}\boldsymbol\psi_n^T
\boldsymbol\Lambda_{\psi,n}^{\mathrm{in},(l)}
\boldsymbol\psi_n
+
\boldsymbol\eta_{\psi,n}^{\mathrm{in},(l)\,T}
\boldsymbol\psi_n
+
\mathrm{const}.
\end{aligned}
$$

定义关于 \(q^{(l-1)}(\boldsymbol h_n)\) 的两个期望量

$$
\mathbf R_{G,n}^{(l-1)}
:=
\mathbb E_{q^{(l-1)}(\boldsymbol h_n)}
\left[
\mathbf G_n^{(l-1)H}(\boldsymbol h_n)
\mathbf G_n^{(l-1)}(\boldsymbol h_n)
\right],
$$

$$
\begin{aligned}
\boldsymbol d_{G,n}^{(l-1)}
:={}&
\mathbb E_{q^{(l-1)}(\boldsymbol h_n)}
\Big[
\mathbf G_n^{(l-1)H}(\boldsymbol h_n)
\\
&\qquad\cdot
\big(
\boldsymbol y_n
-
\boldsymbol\Phi_n^{(l-1)}\boldsymbol h_n
+
\mathbf G_n^{(l-1)}(\boldsymbol h_n)
\boldsymbol\psi_n^{(l-1)}
\big)
\Big].
\end{aligned}
$$

其中，$\mathbf R_{G,n}^{(l-1)}\in\mathbb C^{2K\times2K}$， $\boldsymbol d_{G,n}^{(l-1)}\in\mathbb C^{2K}$。

由于 \(\boldsymbol\psi_n\) 为实向量，上述期望残差关于 \(\boldsymbol\psi_n\) 的二次项可展开为

$$
\begin{aligned}
&\mathbb E_{q^{(l-1)}(\boldsymbol h_n)}
\left[
\left\|
\boldsymbol y_n
-
\boldsymbol c_n^{(l-1)}(\boldsymbol h_n)
-
\mathbf G_n^{(l-1)}(\boldsymbol h_n)
\boldsymbol\psi_n
\right\|^2
\right]
\\
&\quad=
\boldsymbol\psi_n^T
\operatorname{Re}\!\left\{
\mathbf R_{G,n}^{(l-1)}
\right\}
\boldsymbol\psi_n
-
2\boldsymbol\psi_n^T
\operatorname{Re}\!\left\{
\boldsymbol d_{G,n}^{(l-1)}
\right\}
+
\mathrm{const}.
\end{aligned}
$$

下面进一步给出 \(\mathbf R_{G,n}^{(l-1)}\) 和 \(\boldsymbol d_{G,n}^{(l-1)}\) 的显式表达式。由上一轮散射系数替代后验

$$
q^{(l-1)}(\boldsymbol h_n)
=
\mathcal{CN}\!\left(
\boldsymbol h_n;
\boldsymbol m_{h,n}^{(l-1)},
\mathbf C_{h,n}^{(l-1)}
\right),
$$

定义

$$
\begin{aligned}
\mathbf Q_{h,n}^{(l-1)}
&:=
\mathbb E_{q^{(l-1)}(\boldsymbol h_n)}
\left[
\boldsymbol h_n^*\boldsymbol h_n^T
\right]
\\
&=
\left(\mathbf C_{h,n}^{(l-1)}\right)^T
+
\left(\boldsymbol m_{h,n}^{(l-1)}\right)^*
\left(\boldsymbol m_{h,n}^{(l-1)}\right)^T.
\end{aligned}
$$

则第一个期望项为

$$
\mathbf R_{G,n}^{(l-1)}
=
\begin{bmatrix}
\left[
\left(\boldsymbol\Phi_{\theta,n}^{(l-1)}\right)^H
\boldsymbol\Phi_{\theta,n}^{(l-1)}
\right]\odot\mathbf Q_{h,n}^{(l-1)}
&
\left[
\left(\boldsymbol\Phi_{\theta,n}^{(l-1)}\right)^H
\boldsymbol\Phi_{\tau,n}^{(l-1)}
\right]\odot\mathbf Q_{h,n}^{(l-1)}
\\
\left[
\left(\boldsymbol\Phi_{\tau,n}^{(l-1)}\right)^H
\boldsymbol\Phi_{\theta,n}^{(l-1)}
\right]\odot\mathbf Q_{h,n}^{(l-1)}
&
\left[
\left(\boldsymbol\Phi_{\tau,n}^{(l-1)}\right)^H
\boldsymbol\Phi_{\tau,n}^{(l-1)}
\right]\odot\mathbf Q_{h,n}^{(l-1)}
\end{bmatrix},
$$

其中 \(\odot\) 表示 Hadamard 积。第二个期望项为

$$
\begin{aligned}
\boldsymbol d_{G,n}^{(l-1)}
={}&
\begin{bmatrix}
\operatorname{diag}\!\left[
\left(\boldsymbol m_{h,n}^{(l-1)}\right)^*
\right]
\left(\boldsymbol\Phi_{\theta,n}^{(l-1)}\right)^H
\boldsymbol y_n
\\
\operatorname{diag}\!\left[
\left(\boldsymbol m_{h,n}^{(l-1)}\right)^*
\right]
\left(\boldsymbol\Phi_{\tau,n}^{(l-1)}\right)^H
\boldsymbol y_n
\end{bmatrix}
\\
&-
\begin{bmatrix}
\left(
\left[
\left(\boldsymbol\Phi_{\theta,n}^{(l-1)}\right)^H
\boldsymbol\Phi_n^{(l-1)}
\right]
\odot
\mathbf Q_{h,n}^{(l-1)}
\right)\mathbf 1_K
\\
\left(
\left[
\left(\boldsymbol\Phi_{\tau,n}^{(l-1)}\right)^H
\boldsymbol\Phi_n^{(l-1)}
\right]
\odot
\mathbf Q_{h,n}^{(l-1)}
\right)\mathbf 1_K
\end{bmatrix}
+
\mathbf R_{G,n}^{(l-1)}
\boldsymbol\psi_n^{(l-1)},
\end{aligned}
$$

其中 \(\mathbf 1_K\) 为 \(K\) 维全一向量。该表达式同时包含 \(\boldsymbol h_n\) 的后验均值和协方差；若忽略 \(\mathbf C_{h,n}^{(l-1)}\)，则上述矩不再等于关于 \(q^{(l-1)}(\boldsymbol h_n)\) 的完整期望。

将二次型代回 \(\ln q^{(l)}(\boldsymbol\psi_n)\)，可得

$$
\begin{aligned}
\ln q^{(l)}(\boldsymbol\psi_n)
={}&
-\frac{1}{2}
\boldsymbol\psi_n^T
\Big[
\boldsymbol\Lambda_{\psi,n}^{\mathrm{in},(l)}
+
2\bar\gamma_n^{(l-1)}
\operatorname{Re}\!\left\{
\mathbf R_{G,n}^{(l-1)}
\right\}
\Big]
\boldsymbol\psi_n
\\
&+
\Big[
\boldsymbol\eta_{\psi,n}^{\mathrm{in},(l)}
+
2\bar\gamma_n^{(l-1)}
\operatorname{Re}\!\left\{
\boldsymbol d_{G,n}^{(l-1)}
\right\}
\Big]^T
\boldsymbol\psi_n
+
\mathrm{const}.
\end{aligned}
$$

因此，联合角度—时延向量的替代后验为

$$
q^{(l)}(\boldsymbol\psi_n)
=
\mathcal N\!\left(
\boldsymbol\psi_n;
\boldsymbol m_{\psi,n}^{(l)},
\mathbf C_{\psi,n}^{(l)}
\right),
$$

其中

$$
\left(\mathbf C_{\psi,n}^{(l)}\right)^{-1}
=
\boldsymbol\Lambda_{\psi,n}^{\mathrm{in},(l)}
+
2\bar\gamma_n^{(l-1)}
\operatorname{Re}\!\left\{
\mathbf R_{G,n}^{(l-1)}
\right\},
$$

$$
\boldsymbol m_{\psi,n}^{(l)}
=
\mathbf C_{\psi,n}^{(l)}
\left[
\boldsymbol\eta_{\psi,n}^{\mathrm{in},(l)}
+
2\bar\gamma_n^{(l-1)}
\operatorname{Re}\!\left\{
\boldsymbol d_{G,n}^{(l-1)}
\right\}
\right].
$$

上式以 \(2K\) 维向量形式联合更新 \(\boldsymbol\psi_n\)，从而保留角度与时延之间的后验交叉协方差。因子 2 来自圆对称复高斯似然与实值参数 \(\boldsymbol\psi_n\) 的组合。角度与时延的后验均值分别为

$$
\boldsymbol\mu_{\theta,n}^{\mathrm{post},(l)}
=
[\boldsymbol m_{\psi,n}^{(l)}]_{1:K},
\qquad
\boldsymbol\mu_{\tau,n}^{\mathrm{post},(l)}
=
[\boldsymbol m_{\psi,n}^{(l)}]_{K+1:2K}.
$$

为避免重复使用几何模块输入的信息，发送给几何模块的联合外信息为

$$
\boldsymbol\Lambda_{\psi,n}^{\mathrm{ext},(l)}
=
\left(\mathbf C_{\psi,n}^{(l)}\right)^{-1}
-
\boldsymbol\Lambda_{\psi,n}^{\mathrm{in},(l)},
\qquad
\boldsymbol\eta_{\psi,n}^{\mathrm{ext},(l)}
=
\left(\mathbf C_{\psi,n}^{(l)}\right)^{-1}
\boldsymbol m_{\psi,n}^{(l)}
-
\boldsymbol\eta_{\psi,n}^{\mathrm{in},(l)}.
$$

当外信息精度矩阵非奇异时，线谱模块输出的联合高斯消息参数为

$$
\mathbf V_{\psi,n}^{\mathrm{LS},(l)}
=
\left(
\boldsymbol\Lambda_{\psi,n}^{\mathrm{ext},(l)}
\right)^{-1},
\qquad
\boldsymbol\mu_{\psi,n}^{\mathrm{LS},(l)}
=
\mathbf V_{\psi,n}^{\mathrm{LS},(l)}
\boldsymbol\eta_{\psi,n}^{\mathrm{ext},(l)}.
$$

若外信息精度矩阵奇异，则固定保留 $\boldsymbol\Lambda_{\psi,n}^{\mathrm{ext},(l)}$ 和 $\boldsymbol\eta_{\psi,n}^{\mathrm{ext},(l)}$ 的信息形式，不构造 $\mathbf V_{\psi,n}^{\mathrm{LS},(l)}$ 与 $\boldsymbol\mu_{\psi,n}^{\mathrm{LS},(l)}$；第 5.2 节的前向几何消息可直接使用该信息形式。

#### 5.1.2 关于 \(\boldsymbol h_n\) 的一阶近似及其推断

为了获得关于 \(\boldsymbol h_n\) 的条件线性模型，将上一小节的路径矩阵一阶近似记为

$$
\begin{aligned}
\widehat{\boldsymbol\Phi}_n^{(l)}
(\boldsymbol\psi_n)
:={}&
\boldsymbol\Phi_n^{(l-1)}
+
\boldsymbol\Phi_{\theta,n}^{(l-1)}
\operatorname{diag}\!\left(
\boldsymbol\theta_n-\boldsymbol\theta_n^{(l-1)}
\right)
\\
&+
\boldsymbol\Phi_{\tau,n}^{(l-1)}
\operatorname{diag}\!\left(
\boldsymbol\tau_n-\boldsymbol\tau_n^{(l-1)}
\right).
\end{aligned}
$$

于是，给定 \(\boldsymbol\psi_n\) 时，观测模型的一阶近似为

$$
\boldsymbol y_n
\approx
\widehat{\boldsymbol\Phi}_n^{(l)}
(\boldsymbol\psi_n)
\boldsymbol h_n
+
\boldsymbol w_n,
$$

同一局部似然等价写为

$$
\widetilde p^{(l)}\!\left(
\boldsymbol y_n\mid
\boldsymbol h_n,\boldsymbol\psi_n,\gamma_n
\right)
=
\mathcal{CN}\!\left(
\boldsymbol y_n;
\widehat{\boldsymbol\Phi}_n^{(l)}
(\boldsymbol\psi_n)\boldsymbol h_n,
\gamma_n^{-1}\mathbf I_M
\right).
$$

然后固定本轮新的 \(q^{(l)}(\boldsymbol\psi_n)\) 和上一轮的 \(q^{(l-1)}(\gamma_n)\)。根据 KL 坐标最优条件，保留所有与 \(\boldsymbol h_n\) 有关的项，有

$$
\begin{aligned}
\ln q^{(l)}(\boldsymbol h_n)
={}&
\mathbb E_{
q^{(l)}(\boldsymbol\psi_n)
q^{(l-1)}(\gamma_n)}
\left[
\ln\widetilde p^{(l)}(
\boldsymbol y_n\mid
\boldsymbol h_n,\boldsymbol\psi_n,\gamma_n)
\right]
+
\ln p(\boldsymbol h_n)
+
\mathrm{const}
\\
={}&
-\boldsymbol h_n^H
\left[
\bar\gamma_n^{(l-1)}
\mathbf R_{\Phi,n}^{(l)}
+
\boldsymbol\Sigma_{h,n}^{-1}
\right]
\boldsymbol h_n
\\
&+
2\operatorname{Re}\!\left\{
\boldsymbol h_n^H
\bar\gamma_n^{(l-1)}
\left(
\overline{\boldsymbol\Phi}_n^{(l)}
\right)^H
\boldsymbol y_n
\right\}
+
\mathrm{const},
\end{aligned}
$$

其中定义了关于 \(q^{(l)}(\boldsymbol\psi_n)\) 的两个矩阵期望

$$
\overline{\boldsymbol\Phi}_n^{(l)}
:=
\mathbb E_{q^{(l)}(\boldsymbol\psi_n)}
\left[
\widehat{\boldsymbol\Phi}_n^{(l)}
(\boldsymbol\psi_n)
\right],
$$

$$
\mathbf R_{\Phi,n}^{(l)}
:=
\mathbb E_{q^{(l)}(\boldsymbol\psi_n)}
\left[
\left[
\widehat{\boldsymbol\Phi}_n^{(l)}
(\boldsymbol\psi_n)
\right]^H
\widehat{\boldsymbol\Phi}_n^{(l)}
(\boldsymbol\psi_n)
\right].
$$

因此 $\overline{\boldsymbol\Phi}_n^{(l)}\in\mathbb C^{M\times K}$，$\mathbf R_{\Phi,n}^{(l)}\in\mathbb C^{K\times K}$。

下面给出这两个期望的显式表达式。将 \(\boldsymbol\psi_n\) 的后验均值和协方差按照角度与时延分块为

$$
\boldsymbol m_{\psi,n}^{(l)}
=
\begin{bmatrix}
\boldsymbol m_{\theta,n}^{(l)}
\\
\boldsymbol m_{\tau,n}^{(l)}
\end{bmatrix},
\qquad
\mathbf C_{\psi,n}^{(l)}
=
\begin{bmatrix}
\mathbf C_{\theta\theta,n}^{(l)}
&
\mathbf C_{\theta\tau,n}^{(l)}
\\
\mathbf C_{\tau\theta,n}^{(l)}
&
\mathbf C_{\tau\tau,n}^{(l)}
\end{bmatrix}.
$$

定义相对于当前 Taylor 展开点的后验均值增量及其对角矩阵

$$
\overline{\boldsymbol\delta}_{\theta,n}^{(l)}
:=
\boldsymbol m_{\theta,n}^{(l)}
-
\boldsymbol\theta_n^{(l-1)},
\qquad
\overline{\boldsymbol\delta}_{\tau,n}^{(l)}
:=
\boldsymbol m_{\tau,n}^{(l)}
-
\boldsymbol\tau_n^{(l-1)},
$$

$$
\overline{\mathbf D}_{\theta,n}^{(l)}
:=
\operatorname{diag}\!\left(
\overline{\boldsymbol\delta}_{\theta,n}^{(l)}
\right),
\qquad
\overline{\mathbf D}_{\tau,n}^{(l)}
:=
\operatorname{diag}\!\left(
\overline{\boldsymbol\delta}_{\tau,n}^{(l)}
\right).
$$

进一步定义四个增量二阶矩

$$
\begin{aligned}
\mathbf S_{\theta\theta,n}^{(l)}
&:=
\mathbf C_{\theta\theta,n}^{(l)}
+
\overline{\boldsymbol\delta}_{\theta,n}^{(l)}
\overline{\boldsymbol\delta}_{\theta,n}^{(l)T},
\\
\mathbf S_{\theta\tau,n}^{(l)}
&:=
\mathbf C_{\theta\tau,n}^{(l)}
+
\overline{\boldsymbol\delta}_{\theta,n}^{(l)}
\overline{\boldsymbol\delta}_{\tau,n}^{(l)T},
\\
\mathbf S_{\tau\theta,n}^{(l)}
&:=
\mathbf C_{\tau\theta,n}^{(l)}
+
\overline{\boldsymbol\delta}_{\tau,n}^{(l)}
\overline{\boldsymbol\delta}_{\theta,n}^{(l)T},
\\
\mathbf S_{\tau\tau,n}^{(l)}
&:=
\mathbf C_{\tau\tau,n}^{(l)}
+
\overline{\boldsymbol\delta}_{\tau,n}^{(l)}
\overline{\boldsymbol\delta}_{\tau,n}^{(l)T}.
\end{aligned}
$$

于是第一个矩阵期望为

$$
\begin{aligned}
\overline{\boldsymbol\Phi}_n^{(l)}
={}&
\boldsymbol\Phi_n^{(l-1)}
+
\boldsymbol\Phi_{\theta,n}^{(l-1)}
\overline{\mathbf D}_{\theta,n}^{(l)}
+
\boldsymbol\Phi_{\tau,n}^{(l-1)}
\overline{\mathbf D}_{\tau,n}^{(l)}.
\end{aligned}
$$

第二个矩阵期望为

$$
\begin{aligned}
\mathbf R_{\Phi,n}^{(l)}
={}&
\boldsymbol\Phi_n^{(l-1)H}
\boldsymbol\Phi_n^{(l-1)}
\\
&+
\boldsymbol\Phi_n^{(l-1)H}
\boldsymbol\Phi_{\theta,n}^{(l-1)}
\overline{\mathbf D}_{\theta,n}^{(l)}
+
\boldsymbol\Phi_n^{(l-1)H}
\boldsymbol\Phi_{\tau,n}^{(l-1)}
\overline{\mathbf D}_{\tau,n}^{(l)}
\\
&+
\overline{\mathbf D}_{\theta,n}^{(l)}
\boldsymbol\Phi_{\theta,n}^{(l-1)H}
\boldsymbol\Phi_n^{(l-1)}
+
\overline{\mathbf D}_{\tau,n}^{(l)}
\boldsymbol\Phi_{\tau,n}^{(l-1)H}
\boldsymbol\Phi_n^{(l-1)}
\\
&+
\left[
\boldsymbol\Phi_{\theta,n}^{(l-1)H}
\boldsymbol\Phi_{\theta,n}^{(l-1)}
\right]
\odot
\mathbf S_{\theta\theta,n}^{(l)}
\\
&+
\left[
\boldsymbol\Phi_{\theta,n}^{(l-1)H}
\boldsymbol\Phi_{\tau,n}^{(l-1)}
\right]
\odot
\mathbf S_{\theta\tau,n}^{(l)}
\\
&+
\left[
\boldsymbol\Phi_{\tau,n}^{(l-1)H}
\boldsymbol\Phi_{\theta,n}^{(l-1)}
\right]
\odot
\mathbf S_{\tau\theta,n}^{(l)}
\\
&+
\left[
\boldsymbol\Phi_{\tau,n}^{(l-1)H}
\boldsymbol\Phi_{\tau,n}^{(l-1)}
\right]
\odot
\mathbf S_{\tau\tau,n}^{(l)}.
\end{aligned}
$$

上述表达式使用了恒等式 \(\mathbb E[\operatorname{diag}(\boldsymbol x)\mathbf B\operatorname{diag}(\boldsymbol y)]=\mathbf B\odot\mathbb E[\boldsymbol x\boldsymbol y^T]\)，因此完整保留了 \(q^{(l)}(\boldsymbol\psi_n)\) 的角度方差、时延方差以及角度—时延交叉协方差。

将关于 \(\boldsymbol h_n\) 的二次型配方，可得散射系数的替代后验分布

$$
q^{(l)}(\boldsymbol h_n)
=
\mathcal{CN}\!\left(
\boldsymbol h_n;
\boldsymbol m_{h,n}^{(l)},
\mathbf C_{h,n}^{(l)}
\right),
$$

其中

$$
\left(\mathbf C_{h,n}^{(l)}\right)^{-1}
=
\bar\gamma_n^{(l-1)}
\mathbf R_{\Phi,n}^{(l)}
+
\boldsymbol\Sigma_{h,n}^{-1},
$$

$$
\boldsymbol m_{h,n}^{(l)}
=
\mathbf C_{h,n}^{(l)}
\bar\gamma_n^{(l-1)}
\left(
\overline{\boldsymbol\Phi}_n^{(l)}
\right)^H
\boldsymbol y_n.
$$

完成本轮 \(q^{(l)}(\boldsymbol\psi_n)\) 和 \(q^{(l)}(\boldsymbol h_n)\) 更新后，再更新噪声精度。根据 KL 坐标最优条件和 Gamma 先验，有

$$
\begin{aligned}
\ln q^{(l)}(\gamma_n)
={}&
\mathbb E_{
q^{(l)}(\boldsymbol\psi_n)
q^{(l)}(\boldsymbol h_n)}
\left[
\ln\widetilde p^{(l)}(
\boldsymbol y_n\mid
\boldsymbol h_n,\boldsymbol\psi_n,\gamma_n)
\right]
+
\ln p(\gamma_n)
+
\mathrm{const}
\\
={}&
\left(a_{\gamma,0}+M-1\right)\ln\gamma_n
-
\left(
b_{\gamma,0}
+
\varepsilon_n^{(l)}
\right)\gamma_n
+
\mathrm{const},
\end{aligned}
$$

其中

$$
\begin{aligned}
\varepsilon_n^{(l)}
:={}&
\mathbb E_{
q^{(l)}(\boldsymbol h_n)
q^{(l)}(\boldsymbol\psi_n)}
\left[
\left\|
\boldsymbol y_n
-
\widehat{\boldsymbol\Phi}_n^{(l)}
(\boldsymbol\psi_n)\boldsymbol h_n
\right\|^2
\right]
\\
={}&
\boldsymbol y_n^H\boldsymbol y_n
-
2\operatorname{Re}\!\left\{
\boldsymbol y_n^H
\overline{\boldsymbol\Phi}_n^{(l)}
\boldsymbol m_{h,n}^{(l)}
\right\}
\\
&+
\operatorname{tr}\!\left\{
\mathbf R_{\Phi,n}^{(l)}
\left[
\mathbf C_{h,n}^{(l)}
+
\boldsymbol m_{h,n}^{(l)}
\boldsymbol m_{h,n}^{(l)H}
\right]
\right\}.
\end{aligned}
$$

上式使用了散射系数的二阶矩

$$
\mathbb E_{q^{(l)}(\boldsymbol h_n)}
\left[
\boldsymbol h_n\boldsymbol h_n^H
\right]
=
\mathbf C_{h,n}^{(l)}
+
\boldsymbol m_{h,n}^{(l)}
\boldsymbol m_{h,n}^{(l)H},
$$

以及前述关于 \(q^{(l)}(\boldsymbol\psi_n)\) 的两个显式矩阵期望 \(\overline{\boldsymbol\Phi}_n^{(l)}\) 和 \(\mathbf R_{\Phi,n}^{(l)}\)。因此，期望残差同时包含散射系数不确定性、角度—时延不确定性以及二者均值对应的拟合误差。

由 Gamma 分布的共轭形式可得

$$
q^{(l)}(\gamma_n)
=
\operatorname{Ga}\!\left(
\gamma_n;
a_{\gamma,n}^{(l)},
b_{\gamma,n}^{(l)}
\right),
$$

其中

$$
a_{\gamma,n}^{(l)}
=
a_{\gamma,0}+M,
\qquad
b_{\gamma,n}^{(l)}
=
b_{\gamma,0}
+
\varepsilon_n^{(l)},
$$

且

$$
\bar\gamma_n^{(l)}
=
\mathbb E_{q^{(l)}(\gamma_n)}[\gamma_n]
=
\frac{a_{\gamma,n}^{(l)}}{b_{\gamma,n}^{(l)}}.
$$

在第 \(l\) 轮局部推断中，固定以 \(\boldsymbol\psi_n^{(l-1)}\) 为中心构造的 Taylor 模型，并按照

$$
q^{(l)}(\boldsymbol\psi_n)
\longrightarrow
q^{(l)}(\boldsymbol h_n)
\longrightarrow
q^{(l)}(\gamma_n)
$$

的顺序完成一轮向量坐标更新。每一步均对当前局部目标后验执行 KL 坐标最小化，等价于不降低固定 Taylor 模型对应的局部 ELBO。

第 $l$ 轮线谱更新得到的 $\boldsymbol m_{\psi,n}^{(l)}$ 是 $q^{(l)}(\boldsymbol\psi_n)$ 的后验均值，用于构造传向几何模块的外信息，但不直接作为下一轮线谱 Taylor 展开点。几何模块更新 $b^{(l)}(\boldsymbol z)$ 后，由第 5.2 节给出的非线性几何映射生成几何一致的 $\boldsymbol\psi_n^{(l)}$，再在该点重新计算 $\boldsymbol\Phi_n^{(l)}$、$\boldsymbol\Phi_{\theta,n}^{(l)}$ 和 $\boldsymbol\Phi_{\tau,n}^{(l)}$，作为第 $l+1$ 轮的展开量。Taylor 展开点改变后，局部目标后验也随之改变，因此固定展开点下的 ELBO 单调性不能直接推广到外层重线性化；外层迭代统一采用第 5.3 节的几何均值—消息双重停止判据，并记录原始非线性残差作为近似有效性的诊断量。

### 5.2 第二部分：Dirac 几何因子的线性化、后验更新与反馈消息

沿用第 5.1 节定义的联合绝对相位向量 $\boldsymbol\psi_n=[\boldsymbol\theta_n^T,\boldsymbol\tau_n^T]^T$，并定义

$$
\boldsymbol g_n(\boldsymbol z)
:=
\begin{bmatrix}
\boldsymbol g_{\theta,n}(\mathbf P)\\
\boldsymbol g_{\tau,n}(
\mathbf P,\boldsymbol p_{\mathrm{UE}},\beta)
\end{bmatrix},
$$

其中两个向量映射按路径标签堆叠：

$$
\boldsymbol g_{\theta,n}(\mathbf P)
=
\left[
g_{\theta,n,1}(\boldsymbol p_1),
\ldots,
g_{\theta,n,K}(\boldsymbol p_K)
\right]^T,
$$

$$
\boldsymbol g_{\tau,n}(
\mathbf P,\boldsymbol p_{\mathrm{UE}},\beta)
=
\left[
g_{\tau,n,1}(\boldsymbol p_1,\boldsymbol p_{\mathrm{UE}},\beta),
\ldots,
g_{\tau,n,K}(\boldsymbol p_K,\boldsymbol p_{\mathrm{UE}},\beta)
\right]^T.
$$

以及联合 Dirac 几何因子

$$
f_n(\boldsymbol\psi_n,\boldsymbol z)
=\delta\!\left(
\boldsymbol\psi_n-\boldsymbol g_n(\boldsymbol z)
\right).
$$

线谱模块输出采用第 5.1.1 节给出的联合高斯外信息参数

$$
m_{\mathrm{LS},n}^{(l)}(\boldsymbol\psi_n)
\propto
\exp\!\left(
-\frac{1}{2}\boldsymbol\psi_n^T
\boldsymbol\Lambda_{\psi,n}^{\mathrm{ext},(l)}
\boldsymbol\psi_n
+\boldsymbol\eta_{\psi,n}^{\mathrm{ext},(l)T}
\boldsymbol\psi_n
\right).
$$

当 $\boldsymbol\Lambda_{\psi,n}^{\mathrm{ext},(l)}$ 正定时，其矩参数为

$$
\mathbf V_{\psi,n}^{\mathrm{LS},(l)}
=
\left(
\boldsymbol\Lambda_{\psi,n}^{\mathrm{ext},(l)}
\right)^{-1},
\qquad
\boldsymbol\mu_{\psi,n}^{\mathrm{LS},(l)}
=
\mathbf V_{\psi,n}^{\mathrm{LS},(l)}
\boldsymbol\eta_{\psi,n}^{\mathrm{ext},(l)}.
$$

若外信息精度半正定或奇异，前向几何更新仍直接使用其信息形式，不计算伪逆。该消息定义在当前无模糊局部分支上；构造消息前应以 $\boldsymbol\psi_n^{(l-1)}$ 为参考分别对角度和时延分量解缠，避免在 $\pm\pi$ 或 $2\pi$ 边界产生虚假的大增量。

在第 $l$ 次几何更新中，以上一次几何估计 $\boldsymbol z^{(l-1)}$ 为展开点，并记

$$
\boldsymbol g_n^{(l-1)}
=\boldsymbol g_n\!\left(\boldsymbol z^{(l-1)}\right),
\qquad
\mathbf J_n^{(l-1)}
=\left.
\frac{\partial\boldsymbol g_n(\boldsymbol z)}
{\partial\boldsymbol z^T}
\right|_{\boldsymbol z=\boldsymbol z^{(l-1)}}.
$$

按照 $\boldsymbol g_n=[\boldsymbol g_{\theta,n}^T,\boldsymbol g_{\tau,n}^T]^T$ 的行顺序和 $\boldsymbol z=[\boldsymbol p_1^T,\ldots,\boldsymbol p_K^T,\boldsymbol p_{\mathrm{UE}}^T,\beta]^T$ 的列顺序，$\mathbf J_n^{(l-1)}\in\mathbb R^{2K\times(2K+3)}$；其每个非零分块直接使用第 2.2 节给出的解析导数，最后一列固定为对 $\beta$ 的导数，而不是对 $\Delta t$ 的导数。

令 $\Delta\boldsymbol z^{(l)}=\boldsymbol z-\boldsymbol z^{(l-1)}$，则联合 Dirac 因子内部的几何映射近似为

$$
\boldsymbol g_n(\boldsymbol z)
\approx
\boldsymbol g_n^{(l-1)}
+\mathbf J_n^{(l-1)}
\Delta\boldsymbol z^{(l)}.
$$

定义该仿射近似的常数项

$$
\boldsymbol a_{z,n}^{(l-1)}
:=
\boldsymbol g_n^{(l-1)}
-\mathbf J_n^{(l-1)}
\boldsymbol z^{(l-1)},
$$

则 $\boldsymbol g_n(\boldsymbol z)\approx\boldsymbol a_{z,n}^{(l-1)}+\mathbf J_n^{(l-1)}\boldsymbol z$。

这里的一阶近似作用于 delta 函数内部的非线性映射，而不是对 delta 广义函数本身直接求 Taylor 展开。

Dirac 因子传递到几何变量的消息为

$$
\begin{aligned}
m_{f_n\rightarrow z}^{(l)}(\boldsymbol z)
&=
\int
f_n(\boldsymbol\psi_n,\boldsymbol z)
m_{\mathrm{LS},n}^{(l)}(\boldsymbol\psi_n)
\mathrm d\boldsymbol\psi_n
\\
&\propto
m_{\mathrm{LS},n}^{(l)}\!\left(
\boldsymbol a_{z,n}^{(l-1)}
+\mathbf J_n^{(l-1)}\boldsymbol z
\right).
\end{aligned}
$$

联合消息对绝对几何变量的信息矩阵为

$$
\boldsymbol\Lambda_{z,n}^{(l)}
=\mathbf J_n^{(l-1)T}
\boldsymbol\Lambda_{\psi,n}^{\mathrm{ext},(l)}
\mathbf J_n^{(l-1)},
$$

对应的信息向量为

$$
\boldsymbol\eta_{z,n}^{(l)}
=\mathbf J_n^{(l-1)T}
\left[
\boldsymbol\eta_{\psi,n}^{\mathrm{ext},(l)}
-\boldsymbol\Lambda_{\psi,n}^{\mathrm{ext},(l)}
\boldsymbol a_{z,n}^{(l-1)}
\right].
$$

其中 $\boldsymbol\Lambda_{z,n}^{(l)}\in\mathbb R^{(2K+3)\times(2K+3)}$，$\boldsymbol\eta_{z,n}^{(l)}\in\mathbb R^{2K+3}$。

因此

$$
m_{f_n\rightarrow z}^{(l)}(\boldsymbol z)
\propto
\exp\!\left(
-\frac{1}{2}\boldsymbol z^T
\boldsymbol\Lambda_{z,n}^{(l)}
\boldsymbol z
+\boldsymbol\eta_{z,n}^{(l)T}
\boldsymbol z
\right).
$$

单个接收站的 $\boldsymbol\Lambda_{z,n}^{(l)}$ 可能秩亏，因此相应消息不一定能单独归一化为完整维度的高斯分布，但其信息形式仍然有效。融合全部接收站消息与第 3.6 节的固定几何先验 $p(\boldsymbol z)=\mathcal N(\boldsymbol z;\boldsymbol z^{(0)},\mathbf C_{z,0})$，得到第 $l$ 次更新的几何后验

$$
b^{(l)}(\boldsymbol z)
\propto
p(\boldsymbol z)
\prod_{n=1}^{N_{\mathrm{Rx}}}
m_{f_n\rightarrow z}^{(l)}(\boldsymbol z)
=
\mathcal N\!\left(
\boldsymbol z;
\boldsymbol\mu_z^{(l)},
\mathbf C_z^{(l)}
\right),
$$

其中

$$
\left(\mathbf C_z^{(l)}\right)^{-1}
=\mathbf C_{z,0}^{-1}
+\sum_{n=1}^{N_{\mathrm{Rx}}}
\boldsymbol\Lambda_{z,n}^{(l)},
$$

$$
\boldsymbol\mu_z^{(l)}
=\mathbf C_z^{(l)}
\left[
\mathbf C_{z,0}^{-1}\boldsymbol z^{(0)}
+\sum_{n=1}^{N_{\mathrm{Rx}}}
\boldsymbol\eta_{z,n}^{(l)}
\right].
$$

得到 $b^{(l)}(\boldsymbol z)$ 后，令 $\boldsymbol z^{(l)}=\boldsymbol\mu_z^{(l)}$，并在更新后的几何后验均值处重新计算

$$
\boldsymbol g_n^{(l)}
=\boldsymbol g_n\!\left(\boldsymbol\mu_z^{(l)}\right),
\qquad
\mathbf J_n^{(l)}
=\left.
\frac{\partial\boldsymbol g_n(\boldsymbol z)}
{\partial\boldsymbol z^T}
\right|_{\boldsymbol z=\boldsymbol\mu_z^{(l)}}.
$$

据此，供反向消息传递使用的线性化 Dirac 因子为

$$
\widetilde f_n^{(l+1)}(\boldsymbol\psi_n,\boldsymbol z)
=
\delta\!\left(
\boldsymbol\psi_n
-\boldsymbol g_n^{(l)}
-\mathbf J_n^{(l)}
\left(\boldsymbol z-\boldsymbol\mu_z^{(l)}\right)
\right).
$$

若直接对完整几何后验进行一阶投影并采用本版本的角度—时延独立假设，则先定义仅保留边缘方差的对角 covariance

$$
\mathbf V_{\psi,n}^{\mathrm{bel},(l+1)}
:=
\operatorname{Diag}\!\left(
\operatorname{diag}\!\left(
\mathbf J_n^{(l)}
\mathbf C_z^{(l)}
\mathbf J_n^{(l)T}
\right)
\right),
$$

再得到

$$
\begin{aligned}
\widetilde m_{\mathrm{geo}\rightarrow\psi_n}^{(l+1)}(\boldsymbol\psi_n)
&=
\int
\widetilde f_n^{(l+1)}(\boldsymbol\psi_n,\boldsymbol z)
b^{(l)}(\boldsymbol z)
\mathrm d\boldsymbol z
\\
&\approx
\mathcal N\!\left(
\boldsymbol\psi_n;
\boldsymbol g_n^{(l)},
\mathbf V_{\psi,n}^{\mathrm{bel},(l+1)}
\right).
\end{aligned}
$$

其中 $\operatorname{diag}(\cdot)$ 提取矩阵对角元素，$\operatorname{Diag}(\cdot)$ 将向量组成对角矩阵。该式只保留各相位分量的边缘方差，不保留互协方差。由于 $b^{(l)}(\boldsymbol z)$ 已包含第 $n$ 个 Dirac 因子传入的线谱信息，若将上述完整后验投影直接作为下一轮输入消息，会在迭代中重复使用本站信息。严格的因子图消息应首先构造排除当前因子 $f_n$ 的几何 cavity 分布

$$
\begin{aligned}
b_{\setminus n}^{(l)}(\boldsymbol z)
&\propto
p(\boldsymbol z)
\prod_{m\ne n}
m_{f_m\rightarrow z}^{(l)}(\boldsymbol z)
\\
&=
\mathcal N\!\left(
\boldsymbol z;
\boldsymbol\mu_{z\setminus n}^{(l)},
\mathbf C_{z\setminus n}^{(l)}
\right).
\end{aligned}
$$

利用已经计算的完整后验信息参数，可写为

$$
\left(\mathbf C_{z\setminus n}^{(l)}\right)^{-1}
=
\left(\mathbf C_z^{(l)}\right)^{-1}
-\boldsymbol\Lambda_{z,n}^{(l)},
$$

$$
\boldsymbol\mu_{z\setminus n}^{(l)}
=
\mathbf C_{z\setminus n}^{(l)}
\left[
\left(\mathbf C_z^{(l)}\right)^{-1}
\boldsymbol\mu_z^{(l)}
-\boldsymbol\eta_{z,n}^{(l)}
\right].
$$

下面将 geo 到 $\boldsymbol\psi_n$ 的高斯消息按路径展开，以明确 Jacobian 和协方差的具体构造。对第 $k$ 条路径，定义只包含该路径相关几何变量的局部状态及其选择矩阵

$$
\boldsymbol z_k
:=
\begin{bmatrix}
\boldsymbol p_k\\
\boldsymbol p_{\mathrm{UE}}\\
\beta
\end{bmatrix}
=\mathbf S_k\boldsymbol z,
\qquad
\mathbf S_k\in\mathbb R^{5\times(2K+3)}.
$$

在本小节中，$g_{\tau,n,k}(\boldsymbol z_k)$ 简记第 2.2 节的 $g_{\tau,n,k}(\boldsymbol p_k,\boldsymbol p_{\mathrm{UE}},\beta)$，而 $g_{\theta,n,k}(\boldsymbol z_k)$ 简记 $g_{\theta,n,k}(\boldsymbol p_k)$。

记第 $l$ 轮的局部线性化点为 $\overline{\boldsymbol z}_k^{(l)}:=\mathbf S_k\boldsymbol\mu_z^{(l)}$。分别计算归一化时延相位和归一化空间相位关于 $[\boldsymbol p_k^T,\boldsymbol p_{\mathrm{UE}}^T,\beta]^T$ 的 Jacobian 行向量：

$$
\begin{aligned}
\boldsymbol j_{\tau,n,k}^{(l)}
&:=
\left.
\frac{\partial g_{\tau,n,k}(\boldsymbol z_k)}
{\partial\boldsymbol z_k^T}
\right|_{\boldsymbol z_k=\overline{\boldsymbol z}_k^{(l)}}
\\
&=
\begin{bmatrix}
\dfrac{\kappa_\tau}{c}
\left(\boldsymbol u_k^{\mathrm U}+\boldsymbol u_{n,k}^{\mathrm R}\right)^T
&
-\dfrac{\kappa_\tau}{c}
\left(\boldsymbol u_k^{\mathrm U}\right)^T
&
\dfrac{\kappa_\tau}{c}
\end{bmatrix}
\in\mathbb R^{1\times5},
\end{aligned}
$$

$$
\begin{aligned}
\boldsymbol j_{\theta,n,k}^{(l)}
&:=
\left.
\frac{\partial g_{\theta,n,k}(\boldsymbol z_k)}
{\partial\boldsymbol z_k^T}
\right|_{\boldsymbol z_k=\overline{\boldsymbol z}_k^{(l)}}
\\
&=
\begin{bmatrix}
\dfrac{\kappa_a}{\rho_{n,k}^{\mathrm R}}
\boldsymbol e_n^T
\left(
\mathbf I_2-
\boldsymbol u_{n,k}^{\mathrm R}
\left(\boldsymbol u_{n,k}^{\mathrm R}\right)^T
\right)
&
\boldsymbol 0_{1\times2}
&
0
\end{bmatrix}
\in\mathbb R^{1\times5}.
\end{aligned}
$$

其中所有距离和单位方向向量均在 $\overline{\boldsymbol z}_k^{(l)}$ 处计算。若需要在完整几何状态上实现矩阵乘法，则对应的宽 Jacobian 行向量为

$$
\widetilde{\boldsymbol j}_{\tau,n,k}^{(l)}
:=
\boldsymbol j_{\tau,n,k}^{(l)}\mathbf S_k,
\qquad
\widetilde{\boldsymbol j}_{\theta,n,k}^{(l)}
:=
\boldsymbol j_{\theta,n,k}^{(l)}\mathbf S_k,
$$

其中 $\widetilde{\boldsymbol j}_{\tau,n,k}^{(l)},\widetilde{\boldsymbol j}_{\theta,n,k}^{(l)}\in\mathbb R^{1\times(2K+3)}$。

将传递到第 $n$ 个接收站 Dirac 因子的几何 ext 消息在局部变量 $\boldsymbol z_k$ 上的均值和协方差记为

$$
\boldsymbol\mu_{p,n,k}^{\mathrm{ext},(l)}
:=
\mathbf S_k\boldsymbol\mu_{z\setminus n}^{(l)},
\qquad
\mathbf C_{p,n,k}^{\mathrm{ext},(l)}
:=
\mathbf S_k\mathbf C_{z\setminus n}^{(l)}\mathbf S_k^T
\in\mathbb R^{5\times5}.
$$

下文将 $\mathbf C_{p,n,k}^{\mathrm{ext},(l)}$ 简记为 $\mathbf C_p$。这里的 $\mathbf C_p$ 是传入线性化 Dirac 因子的几何 ext 消息协方差，包含 $\boldsymbol p_k$、$\boldsymbol p_{\mathrm{UE}}$ 和 $\beta$ 之间的协方差，不是第 3.6 节中的目标位置先验协方差 $\mathbf C_{p,k}^{(0)}$。

由一阶近似和标量 Dirac 映射，时延相位消息为

$$
\begin{aligned}
m_{\mathrm{geo}\rightarrow\tau_{n,k}}^{(l+1)}(\tau_{n,k})
&\approx
\int
\delta\!\left(
\tau_{n,k}-g_{\tau,n,k}(\overline{\boldsymbol z}_k^{(l)})
-\boldsymbol j_{\tau,n,k}^{(l)}
(\boldsymbol z_k-\overline{\boldsymbol z}_k^{(l)})
\right)
\\
&\qquad\times
\mathcal N\!\left(
\boldsymbol z_k;
\boldsymbol\mu_{p,n,k}^{\mathrm{ext},(l)},
\mathbf C_{p,n,k}^{\mathrm{ext},(l)}
\right)
\mathrm d\boldsymbol z_k
\\
&=
\mathcal N\!\left(
\tau_{n,k};
\mu_{\tau,n,k}^{\mathrm{in},(l+1)},
\sigma_{\tau,n,k}^{2,\mathrm{in},(l+1)}
\right),
\end{aligned}
$$

其中

$$
\mu_{\tau,n,k}^{\mathrm{in},(l+1)}
=
g_{\tau,n,k}(\overline{\boldsymbol z}_k^{(l)})
+\boldsymbol j_{\tau,n,k}^{(l)}
\left(
\boldsymbol\mu_{p,n,k}^{\mathrm{ext},(l)}
-\overline{\boldsymbol z}_k^{(l)}
\right),
$$

$$
\begin{aligned}
\sigma_{\tau,n,k}^{2,\mathrm{in},(l+1)}
&=
\boldsymbol j_{\tau,n,k}^{(l)}
\mathbf C_{p,n,k}^{\mathrm{ext},(l)}
\boldsymbol j_{\tau,n,k}^{(l)T}
\\
&=
\widetilde{\boldsymbol j}_{\tau,n,k}^{(l)}
\mathbf C_{z\setminus n}^{(l)}
\widetilde{\boldsymbol j}_{\tau,n,k}^{(l)T}.
\end{aligned}
$$

同理，角度相位消息为

$$
\begin{aligned}
m_{\mathrm{geo}\rightarrow\theta_{n,k}}^{(l+1)}(\theta_{n,k})
&\approx
\mathcal N\!\left(
\theta_{n,k};
\mu_{\theta,n,k}^{\mathrm{in},(l+1)},
\sigma_{\theta,n,k}^{2,\mathrm{in},(l+1)}
\right),
\end{aligned}
$$

其中

$$
\mu_{\theta,n,k}^{\mathrm{in},(l+1)}
=
g_{\theta,n,k}(\overline{\boldsymbol z}_k^{(l)})
+\boldsymbol j_{\theta,n,k}^{(l)}
\left(
\boldsymbol\mu_{p,n,k}^{\mathrm{ext},(l)}
-\overline{\boldsymbol z}_k^{(l)}
\right),
$$

$$
\begin{aligned}
\sigma_{\theta,n,k}^{2,\mathrm{in},(l+1)}
&=
\boldsymbol j_{\theta,n,k}^{(l)}
\mathbf C_{p,n,k}^{\mathrm{ext},(l)}
\boldsymbol j_{\theta,n,k}^{(l)T}
\\
&=
\widetilde{\boldsymbol j}_{\theta,n,k}^{(l)}
\mathbf C_{z\setminus n}^{(l)}
\widetilde{\boldsymbol j}_{\theta,n,k}^{(l)T}.
\end{aligned}
$$

因此，$\sigma_{\tau,n,k}^{2,\mathrm{in},(l+1)}$ 和 $\sigma_{\theta,n,k}^{2,\mathrm{in},(l+1)}$ 都是由各自 Jacobian 行向量按照 $\boldsymbol j\mathbf C_p\boldsymbol j^T$ 得到的标量。本版本假设所有 $\{\tau_{n,k},\theta_{n,k}\}_{k=1}^{K}$ 在 geo→线谱消息中相互独立，不计算 $\boldsymbol j_{\tau,n,k}^{(l)}\mathbf C_p\boldsymbol j_{\theta,n,k}^{(l)T}$，也不保留不同路径之间的交叉协方差。

按照全文固定的 $\boldsymbol\psi_n=[\boldsymbol\theta_n^T,\boldsymbol\tau_n^T]^T$ 排列，定义

$$
\boldsymbol\mu_{\psi,n}^{\mathrm{in},(l+1)}
=
\begin{bmatrix}
\mu_{\theta,n,1}^{\mathrm{in},(l+1)}\\
\vdots\\
\mu_{\theta,n,K}^{\mathrm{in},(l+1)}\\
\mu_{\tau,n,1}^{\mathrm{in},(l+1)}\\
\vdots\\
\mu_{\tau,n,K}^{\mathrm{in},(l+1)}
\end{bmatrix},
$$

$$
\begin{aligned}
\mathbf V_{\psi,n}^{\mathrm{in},(l+1)}
=
\operatorname{diag}\!\big(&
\sigma_{\theta,n,1}^{2,\mathrm{in},(l+1)},\ldots,
\sigma_{\theta,n,K}^{2,\mathrm{in},(l+1)},
\\
&\sigma_{\tau,n,1}^{2,\mathrm{in},(l+1)},\ldots,
\sigma_{\tau,n,K}^{2,\mathrm{in},(l+1)}
\big).
\end{aligned}
$$

于是 geo→线谱消息分解为

$$
\begin{aligned}
m_{\mathrm{geo}\rightarrow\psi_n}^{(l+1)}(\boldsymbol\psi_n)
&=
\prod_{k=1}^{K}
m_{\mathrm{geo}\rightarrow\theta_{n,k}}^{(l+1)}(\theta_{n,k})
m_{\mathrm{geo}\rightarrow\tau_{n,k}}^{(l+1)}(\tau_{n,k})
\\
&=
\mathcal N\!\left(
\boldsymbol\psi_n;
\boldsymbol\mu_{\psi,n}^{\mathrm{in},(l+1)},
\mathbf V_{\psi,n}^{\mathrm{in},(l+1)}
\right).
\end{aligned}
$$

这里的独立性是对线性化 Dirac 推前分布作出的对角 Gaussian 近似：它保留每个分量由 $\boldsymbol j\mathbf C_p\boldsymbol j^T$ 得到的边缘方差，但主动舍弃底层几何变量共享所诱导的非对角协方差。该假设只约束 geo→线谱输入消息；线谱模块更新得到的替代后验 $q^{(l)}(\boldsymbol\psi_n)$ 仍可因观测似然产生非对角协方差。若使用完整几何后验而非 cavity ext 消息作单次分布投影，则只需将 $\boldsymbol\mu_{z\setminus n}^{(l)}$ 和 $\mathbf C_{z\setminus n}^{(l)}$ 替换为 $\boldsymbol\mu_z^{(l)}$ 和 $\mathbf C_z^{(l)}$，但不能将其当作无重复信息的外信息消息循环反馈。

为得到唯一且数值稳定的信息参数，代码中固定使用正则化传播协方差

$$
\widetilde{\mathbf V}_{\psi,n}^{\mathrm{in},(l+1)}
:=
\mathbf V_{\psi,n}^{\mathrm{in},(l+1)}
+\epsilon_{\mathrm{geo}}\mathbf I_{2K}.
$$

因此，第 $l+1$ 轮线谱推断所需的输入信息参数更新为

$$
\boldsymbol\Lambda_{\psi,n}^{\mathrm{in},(l+1)}
=
\left(\widetilde{\mathbf V}_{\psi,n}^{\mathrm{in},(l+1)}\right)^{-1},
\qquad
\boldsymbol\eta_{\psi,n}^{\mathrm{in},(l+1)}
=
\boldsymbol\Lambda_{\psi,n}^{\mathrm{in},(l+1)}
\boldsymbol\mu_{\psi,n}^{\mathrm{in},(l+1)},
$$

若任一标量传播方差 $\boldsymbol j\mathbf C_p\boldsymbol j^T$ 为零或数值上过小，则对应分量的未正则化推前分布退化；$\epsilon_{\mathrm{geo}}\mathbf I_{2K}$ 为每个独立分量提供统一的方差下限，同时必须记录相应 Jacobian 范数和几何协方差秩以提示不可辨识方向。本文算法固定执行外信息消除，不采用把 $\boldsymbol\mu_{z\setminus n}^{(l)}$ 和 $\mathbf C_{z\setminus n}^{(l)}$ 替换为完整后验参数的反馈方式。

线谱模块下一轮的 Taylor 展开点不取上述高斯消息的均值，而取完整几何后验均值的非线性映射

$$
\boldsymbol\psi_n^{(l)}
:=
\boldsymbol g_n\!\left(\boldsymbol\mu_z^{(l)}\right)
=
\begin{bmatrix}
\boldsymbol g_{\theta,n}\!\left(\boldsymbol\mu_P^{(l)}\right)
\\
\boldsymbol g_{\tau,n}\!\left(
\boldsymbol\mu_P^{(l)},
\boldsymbol\mu_{\mathrm{UE}}^{(l)},
\mu_{\beta}^{(l)}
\right)
\end{bmatrix},
$$

其中 $\boldsymbol\mu_P^{(l)}$、$\boldsymbol\mu_{\mathrm{UE}}^{(l)}$ 和 $\mu_{\beta}^{(l)}$ 分别为 $\boldsymbol\mu_z^{(l)}$ 中目标位置、UE 位置和公共距离偏差的分块。随后在 $\boldsymbol\psi_n^{(l)}$ 处重新计算 $\boldsymbol\Phi_n^{(l)}$、$\boldsymbol\Phi_{\theta,n}^{(l)}$ 和 $\boldsymbol\Phi_{\tau,n}^{(l)}$，并将 $m_{\mathrm{geo}\rightarrow\psi_n}^{(l+1)}(\boldsymbol\psi_n)$ 作为第 $l+1$ 轮连续线谱估计的输入消息。这里使用非线性映射 $\boldsymbol g_n(\boldsymbol\mu_z^{(l)})$ 而不是一阶消息均值 $\boldsymbol\mu_{\psi,n}^{\mathrm{in},(l+1)}$，可保证展开点严格满足当前后验均值对应的几何关系；后者仅承担不确定性与外信息传递的作用。

每次重新线性化都应重新使用固定先验 $p(\boldsymbol z)$ 与当前线谱外信息构造后验，不能把上一次已经融合过相同消息的 $b^{(l-1)}(\boldsymbol z)$ 再当作新先验，否则会重复计算信息。两部分持续迭代，直至几何后验均值及返回线谱模块的消息同时满足第 5.3 节的收敛条件。

最终从 $\boldsymbol\mu_z^{(l)}$ 的对应分块读取各目标位置、UE 位置和公共距离偏差的后验均值，从 $\mathbf C_z^{(l)}$ 的对应对角分块读取其后验协方差。若 $\mu_\beta^{(l)}$ 和 $[\mathbf C_z^{(l)}]_{\beta\beta}$ 分别表示距离偏差的后验均值和方差，则秒单位时钟偏差的输出为 $\widehat{\Delta t}=\mu_\beta^{(l)}/c$，$\operatorname{var}(\Delta t\mid\mathcal Y)=[\mathbf C_z^{(l)}]_{\beta\beta}/c^2$；其与其他状态的交叉协方差也应按 $\operatorname{cov}(\Delta t,\boldsymbol x\mid\mathcal Y)=\operatorname{cov}(\beta,\boldsymbol x\mid\mathcal Y)/c$ 换算。

### 5.3 统一推断流程

```text
输入：N_Rx 个 Rx 的原始观测、Rx 几何、OFDM/阵列参数、UE 粗位置

1. 去除已知导频。
2. 对每个 Rx 执行过采样 2D-FFT，提取 K 个归一化空间相位—时延相位峰值。
3. 反归一化得到物理角度和秒单位测量时延，再使用多站射线交汇残差和接收站间距离差完成路径关联，并以全部已关联射线的加权最小二乘交点作为目标位置初值。
4. 用 LM 求得目标、UE 和距离偏差 beta 的粗估计 z^(0)，最后一个状态分量不使用秒单位 Delta t。
5. 以 z^(0) 为几何高斯先验均值，初始化 q^(0)(gamma_n) 和 q^(0)(h_n)，并将几何先验经初始线性化 Dirac 因子转换为 m_geo->psi_n^(1)。

第一部分：连续线谱贝叶斯估计
6. 以 2D-FFT 的 theta_n^(0)、tau_n^(0) 为首次 Taylor 展开中心。
7. 首先对 psi_n 进行推断：计算 Phi_theta,n^(l-1)、Phi_tau,n^(l-1) 和 G_n^(l-1)(h_n)，使用 q^(l-1)(h_n) 与 q^(l-1)(gamma_n) 联合更新 q^(l)(psi_n)。
8. 然后对 h_n 进行推断：构造给定 psi_n 时的等价一阶观测矩阵 Phihat_n^(l)(psi_n)，使用 q^(l)(psi_n) 与 q^(l-1)(gamma_n) 更新 q^(l)(h_n)。
9. 使用 q^(l)(psi_n) 和 q^(l)(h_n) 更新 q^(l)(gamma_n)，从而按 psi_n→h_n→gamma_n 的顺序完成第 l 轮 KL 坐标最小化（等价于 ELBO 坐标上升）。
10. 输出保留角度—时延交叉协方差的联合外信息；本轮线谱后验均值不直接作为下一轮 Taylor 展开中心。

第二部分：Dirac 几何因子的线性化消息传递
11. 以 z^(l-1) 为展开点，对各 Dirac 因子内部的几何映射 g_n(z) 作一阶 Taylor 展开。
12. 将线谱联合高斯外信息通过线性化后的 Dirac 因子，转换成传向 z 的高斯信息消息。
13. 融合各 Rx 的信息消息与固定几何先验 p(z)，解析计算 b^(l)(z)=N(z;mu_z^(l),C_z^(l))。
14. 令 z^(l)=mu_z^(l)，在更新后的后验均值处，分别计算每条路径 tau_(n,k) 和 theta_(n,k) 关于 [p_k^T,p_UE^T,beta]^T 的 1×5 Jacobian 行向量，并通过选择矩阵形成两个 1×(2K+3) 宽 Jacobian 行向量。
15. 构造排除当前因子信息的几何 cavity ext 消息及其局部协方差 C_p；分别使用 tau 和 theta 的 Jacobian 按 J*C_p*J^T 计算两个标量方差，并在独立性假设下组装对角 covariance，形成传回第 l+1 轮线谱估计的 Gaussian 消息。
16. 令下一轮线谱 Taylor 展开点 psi_n^(l)=g_n(mu_z^(l))，在该点重新计算 Phi_n^(l)、Phi_theta,n^(l) 和 Phi_tau,n^(l)。
17. 重复两部分的消息传递直至收敛，再从最终几何后验的相应分块输出目标位置、UE 位置和距离偏差的均值与协方差，并通过 Delta t=beta/c 换算时钟偏差结果。
```

为使循环出口可直接实现，定义几何均值和返回线谱模块的信息参数的相对变化量

$$
\delta_z^{(l)}
:=
\frac{
\|\boldsymbol\mu_z^{(l)}-\boldsymbol\mu_z^{(l-1)}\|
}{
\|\boldsymbol\mu_z^{(l-1)}\|+\epsilon_{\mathrm{stop}}
},
$$

$$
\delta_{\mathrm{msg}}^{(l)}
:=
\max_n
\max\!\left{
\frac{
\|\boldsymbol\Lambda_{\psi,n}^{\mathrm{in},(l+1)}-
\boldsymbol\Lambda_{\psi,n}^{\mathrm{in},(l)}\|_{\mathrm F}
}{
\|\boldsymbol\Lambda_{\psi,n}^{\mathrm{in},(l)}\|_{\mathrm F}
+\epsilon_{\mathrm{stop}}
},
\frac{
\|\boldsymbol\eta_{\psi,n}^{\mathrm{in},(l+1)}-
\boldsymbol\eta_{\psi,n}^{\mathrm{in},(l)}\|
}{
\|\boldsymbol\eta_{\psi,n}^{\mathrm{in},(l)}\|
+\epsilon_{\mathrm{stop}}
}
\right}.
$$

当 $\delta_z^{(l)}<\varepsilon_z$ 且 $\delta_{\mathrm{msg}}^{(l)}<\varepsilon_{\mathrm{msg}}$ 时停止外层迭代；否则继续至最大迭代次数 $L_{\max}$。每轮还应计算原始非线性负对数似然

$$
\mathcal J_{\mathrm{NL}}^{(l)}
=
\sum_{n=1}^{N_{\mathrm{Rx}}}
\bar\gamma_n^{(l)}
\left\|
\boldsymbol y_n
-\boldsymbol\Phi_n\!\left(
\boldsymbol g_{\theta,n}(\boldsymbol\mu_P^{(l)}),
\boldsymbol g_{\tau,n}(\boldsymbol\mu_P^{(l)},
\boldsymbol\mu_{\mathrm{UE}}^{(l)},
\mu_{\beta}^{(l)})
\right)
\boldsymbol m_{h,n}^{(l)}
\right\|^2,
$$

用于诊断一阶近似是否离开有效邻域，但不以该诊断量替代上述双重收敛判据。

## 6. 可辨识性与数值稳定性

### 6.1 路径标签与数据关联

参数化模型要求所有接收站的第 $k$ 条路径对应同一目标。路径排列本身具有置换不变性，因此必须通过粗估计阶段的多站数据关联固定标签，或在贝叶斯推断中显式处理路径置换。当前方案采用前者。

### 6.2 几何可辨识性

- 目标数量只是可辨识性的必要条件之一；
- 目标、UE 与接收站的退化几何会导致雅可比秩亏或条件数过大；
- 两个目标角度或时延过近时，2D-FFT 无法稳定分离对应路径；
- UE 位置与公共距离偏差可能存在强相关后验，不能只报告边缘点估计而忽略协方差。

建议记录几何可辨识性矩阵的最小奇异值、条件数以及 $b(\boldsymbol z)$ 中 UE 位置与 $\beta$ 的后验相关系数；若需要报告秒单位时钟偏差的交叉协方差，再按 $1/c$ 进行换算。

### 6.3 数值检查

- 构造 $\arccos$ 前将方向余弦裁剪到 $[-1,1]$；
- 检查所有目标到 UE 和接收站的距离分母不为零；
- 保证 $\gamma_n>0$，并统一 Gamma 分布的 rate/scale 约定；
- 保证 $\boldsymbol\Sigma_{h,n}$ 和几何先验协方差为 Hermitian/实对称正定矩阵；
- 用线性方程求解代替显式矩阵求逆；
- 对 $\mathbf C_{h,n}$、$\mathbf C_{\psi,n}$ 和几何后验协方差 $\mathbf C_z^{(l)}$ 使用 Cholesky 分解或对称化处理；
- 以当前 Taylor 中心为参考对周期相位增量解缠；
- 每次线谱相位中心更新后重新计算 $\boldsymbol\Phi_n$、$\boldsymbol\Phi_{\theta,n}$ 和 $\boldsymbol\Phi_{\tau,n}$；
- 记录 $\|\boldsymbol\theta_n-\boldsymbol\theta_n^{(l-1)}\|$、$\|\boldsymbol\tau_n-\boldsymbol\tau_n^{(l-1)}\|$ 和原始非线性残差，用于诊断 Taylor 近似是否有效；
- 检查几何 Jacobian 的秩和条件数，并保留单站秩亏消息的信息形式；
- 每次重新线性化都由固定几何先验和当前线谱外信息重新构造后验，避免重复累计同一信息；
- 对解析矩阵矩与数值采样结果做小规模交叉验证，确认协方差项实现正确。

## 7. 实验整理

### 7.1 固定配置

- [ ] UE、$N_{\mathrm{Rx}}\geq2$ 个接收站和 $K$ 个目标的二维坐标；
- [ ] 两个 ULA 的朝向、阵元数、阵元间距和载频；
- [ ] 导频子载波索引、子载波间隔、均匀抽取间隔和导频符号；
- [ ] 散射系数协方差 $\boldsymbol\Sigma_{h,n}$；
- [ ] Gamma 先验参数 $a_{\gamma,0},b_{\gamma,0}$；
- [ ] 已知路径数 $K$、2D-FFT 过采样倍数和峰值抑制半径；
- [ ] 多站路径关联的射线垂距尺度 $T_{\perp}$、距离差尺度 $T_d$、条件数阈值 $T_{\mathrm{cond}}$、最大关联迭代次数 $I_{\mathrm{assoc}}$、权重构造方式与冲突消解规则；
- [ ] LM 初始阻尼、停止门限和最大迭代次数；
- [ ] 三类几何先验协方差 $\mathbf C_{p,k}^{(0)}$、$\mathbf C_{\mathrm{UE}}^{(0)}$、$\sigma_{\beta,0}^2$；若输入为秒单位方差，则先换算为 $\sigma_{\beta,0}^2=c^2\sigma_{\Delta t,0}^2$；
- [ ] 外层均值—消息停止门限、最大迭代次数和非线性残差诊断阈值；
- [ ] 几何 Jacobian 的秩判定阈值和病态矩阵处理规则；Jacobian 本身固定使用第 2.2 节的解析式；

### 7.2 建议对比方法

1. 仅 2D-FFT + 几何 LM，不进行贝叶斯精化；
2. Taylor 线谱估计后直接将角度和时延后验均值代入几何方程；
3. Taylor 线谱估计 + 单次线性化 Dirac 几何消息传递；
4. Taylor 线谱估计 + 迭代重线性化的几何后验更新；
5. 固定噪声精度、固定 UE 位置或固定时钟偏差的消融方法；
6. 原位置网格稀疏方法可作为历史基线，但不属于当前概率模型。

### 7.3 评价指标

- 目标位置 RMSE 或多目标 Chamfer distance；
- UE 位置 RMSE；
- 时钟偏差 RMSE；
- 归一化空间相位和时延相位的粗估计误差，以及反归一化后的物理角度和测量时延误差；
- 散射系数 NMSE；
- 噪声精度估计误差；
- 负对数似然、ELBO 和收敛迭代次数；
- 后验可信区间覆盖率；
- 运行时间与峰值内存。

设独立 Monte Carlo 实验次数为 $R$，第 $r$ 次实验的估计量加上上标 $(r)$。目标位置采用标签置换不变的 RMSE

$$
\operatorname{RMSE}_{P}
=
\sqrt{
\frac{1}{RK}
\sum_{r=1}^{R}
\min_{\pi\in\mathcal P_K}
\sum_{k=1}^{K}
\left\|
\widehat{\boldsymbol p}_{\pi(k)}^{(r)}
-\boldsymbol p_k^{(r)}
\right\|^2
}.
$$

UE 位置和时钟偏差 RMSE 分别为，其中每次实验先由几何估计量换算 $\widehat{\Delta t}^{(r)}=\widehat{\beta}^{(r)}/c$，

$$
\operatorname{RMSE}_{\mathrm{UE}}
=
\sqrt{
\frac{1}{R}
\sum_{r=1}^{R}
\left\|
\widehat{\boldsymbol p}_{\mathrm{UE}}^{(r)}
-\boldsymbol p_{\mathrm{UE}}^{(r)}
\right\|^2
},
$$

$$
\operatorname{RMSE}_{\Delta t}
=
\sqrt{
\frac{1}{R}
\sum_{r=1}^{R}
\left(
\widehat{\Delta t}^{(r)}-\Delta t^{(r)}
\right)^2
}.
$$

散射系数 NMSE 和噪声精度相对 RMSE 定义为

$$
\operatorname{NMSE}_{h}
=
\frac{
\sum_{r=1}^{R}
\sum_{n=1}^{N_{\mathrm{Rx}}}
\left\|
\widehat{\boldsymbol h}_n^{(r)}
-\boldsymbol h_n^{(r)}
\right\|^2
}{
\sum_{r=1}^{R}
\sum_{n=1}^{N_{\mathrm{Rx}}}
\left\|\boldsymbol h_n^{(r)}\right\|^2
},
$$

$$
\operatorname{RRMSE}_{\gamma}
=
\sqrt{
\frac{1}{RN_{\mathrm{Rx}}}
\sum_{r=1}^{R}
\sum_{n=1}^{N_{\mathrm{Rx}}}
\left(
\frac{\widehat\gamma_n^{(r)}-\gamma_n^{(r)}}
{\gamma_n^{(r)}}
\right)^2
}.
$$

归一化相位误差使用主值差 $\operatorname{wrap}_{(-\pi,\pi]}(\widehat\theta-\theta)$ 和 $\operatorname{wrap}_{(-\pi,\pi]}(\widehat\tau-\tau)$ 后再计算 RMSE，避免将跨越周期边界的等价相位计为大误差。

### 7.4 正确性测试

- [ ] 无噪声单径情况下，2D-FFT 峰值映射到正确的归一化空间相位和时延相位，并能反归一化为正确的方向余弦和秒单位测量时延；
- [ ] 使用真实角度和时延时，LM 能恢复目标、UE 与距离偏差 $\beta$，并能通过 $\Delta t=\beta/c$ 恢复时钟偏差；
- [ ] 数值差分验证 $\boldsymbol\Phi_{\theta,n}$ 和 $\boldsymbol\Phi_{\tau,n}$；
- [ ] 验证 Taylor 近似误差在小增量范围内随增量范数呈二阶衰减；
- [ ] 固定几何时，$q(\boldsymbol h_n)$ 与标准线性复高斯后验一致；
- [ ] Gamma 更新中的期望残差包含均值和协方差两部分；
- [ ] 线谱模块每轮替代后验更新后局部 ELBO 不降低，或满足规定的接受准则；
- [ ] 在线性化误差可忽略时，Dirac 因子的第 $l$ 轮消息满足 $m_{f_n\rightarrow z}^{(l)}(\boldsymbol z)\approx m_{\mathrm{LS},n}^{(l)}(\boldsymbol g_n^{(l-1)}+\mathbf J_n^{(l-1)}\Delta\boldsymbol z^{(l)}))$；
- [ ] 由信息参数解析计算的几何后验与直接高斯乘积结果一致；
- [ ] 数值差分分别验证每条路径的 $1\times5$ 时延 Jacobian 和角度 Jacobian，并检查选择矩阵嵌入后的两个宽 Jacobian 行向量位置正确；
- [ ] 对给定几何 ext 协方差 $\mathbf C_p$，检查每个时延和角度消息方差均等于各自的标量 $\boldsymbol j\mathbf C_p\boldsymbol j^T$，且组装后的 geo→线谱 covariance 为对角矩阵；
- [ ] 对 $N_{\mathrm{Rx}}$ 个接收站任意重排输入顺序后，多站关联结果和目标粗位置在统一目标标签下保持不变；
- [ ] 路径排列改变后，经数据关联得到的几何结果保持一致；
- [ ] 几何秩亏时能够给出诊断，而不是输出虚假的高置信度结果。

## 8. 当前代码接口与后续实现边界

当前仓库已经具备：

- `CoarseEstimation2DFFT.m`：过采样 2D-FFT、二维峰值抑制和亚栅格插值；
- `SpecEstimation2D.m`：兼容旧实验调用的 FFT 包装接口；
- `estimate_bias_by_coherence.m`：双站路径关联以及固定目标位置条件下的 UE—距离偏差简化 LM 粗估计。

这三个函数仍采用旧接口或简化模型：`CoarseEstimation2DFFT.m` 的阵元维符号和输出参数与本文归一化模型不同，`estimate_bias_by_coherence.m` 没有联合更新目标位置。因此生成新代码时只能复用其峰值搜索、指派和阻尼控制思路，不能直接把其输入输出当作本文算法接口。

与本文公式一一对应的目标模块接口如下；函数名是建议名称，输入输出的数学含义为固定约定。

| 模块 | 关键输入 | 必须输出 | 对应公式位置 |
|---|---|---|---|
| `prepare_observation` | 原始样本、非零导频、$N_p,N_a$ | $\mathbf Y_n,\boldsymbol y_n$ | 第 2.3 节 |
| `coarse_estimation_2d_fft` | $\mathbf Y_n,K,\rho_\tau,\rho_\theta$、搜索区间 | $\widehat{\boldsymbol\theta}_n^{(0)},\widehat{\boldsymbol\tau}_n^{(0)}$、峰值状态 | 第 4.1 节 |
| `associate_multi_rx_paths` | 全部接收站的粗相位、位置与朝向、$T_{\perp},T_d,T_{\mathrm{cond}},I_{\mathrm{assoc}}$ | $\widehat{\mathfrak S},\{\boldsymbol p_k^{\mathrm{ray}}\}$、关联诊断 | 第 4.2 节 |
| `lm_geometry_initialization` | 关联后的粗相位、最后一维为 $\beta$ 的 $\boldsymbol z_{\mathrm{init}}$、$\mathbf W_0$ | $\boldsymbol z^{(0)}$、最终 Jacobian 秩、代价历史 | 第 4.3 节 |
| `build_parametric_path_matrix` | $\boldsymbol\theta_n,\boldsymbol\tau_n,N_a,N_p$ | $\boldsymbol\Phi_n$ | 第 2.1、2.3 节 |
| `build_path_derivative_matrices` | 当前展开点 | $\boldsymbol\Phi_{\theta,n},\boldsymbol\Phi_{\tau,n}$ | 第 5.1.1 节 |
| `initialize_bayesian_model` | FFT 展开点、$\boldsymbol z^{(0)},\mathbf C_{z,0}$、先验参数 | $q^{(0)}(\boldsymbol h_n),q^{(0)}(\gamma_n),m_{\mathrm{geo}\rightarrow\psi_n}^{(1)}$ | 第 5.1 节开头 |
| `update_angle_delay_posterior` | 当前 Taylor 模型、$q^{(l-1)}(\boldsymbol h_n),q^{(l-1)}(\gamma_n)$、几何输入消息 | $\boldsymbol m_{\psi,n}^{(l)},\mathbf C_{\psi,n}^{(l)}$ | 第 5.1.1 节 |
| `update_scattering_posterior` | $q^{(l)}(\boldsymbol\psi_n),q^{(l-1)}(\gamma_n)$ | $\boldsymbol m_{h,n}^{(l)},\mathbf C_{h,n}^{(l)}$ | 第 5.1.2 节 |
| `update_noise_precision_posterior` | $q^{(l)}(\boldsymbol\psi_n),q^{(l)}(\boldsymbol h_n)$ | $a_{\gamma,n}^{(l)},b_{\gamma,n}^{(l)},\bar\gamma_n^{(l)}$ | 第 5.1.2 节 |
| `form_linespectral_extrinsic_message` | $q^{(l)}(\boldsymbol\psi_n)$、本轮几何输入信息 | $\boldsymbol\Lambda_{\psi,n}^{\mathrm{ext},(l)},\boldsymbol\eta_{\psi,n}^{\mathrm{ext},(l)}$ | 第 5.1.1 节末尾 |
| `build_geometry_jacobian` | $\boldsymbol z$、接收站几何、路径索引 $k$ | $\boldsymbol j_{\tau,n,k},\boldsymbol j_{\theta,n,k}\in\mathbb R^{1\times5}$、对应宽 Jacobian 行向量及按 $\boldsymbol\psi_n$ 顺序组装的 $\mathbf J_n$ | 第 2.2、5.2 节 |
| `propagate_linespectral_to_geometry` | 线谱外信息、$\boldsymbol z^{(l-1)},\mathbf J_n^{(l-1)}$ | $\boldsymbol\Lambda_{z,n}^{(l)},\boldsymbol\eta_{z,n}^{(l)}$ | 第 5.2 节 |
| `update_geometry_posterior` | 固定先验、全部接收站几何消息 | $\boldsymbol\mu_z^{(l)},\mathbf C_z^{(l)}$ | 第 5.2 节 |
| `propagate_geometry_to_linespectral` | 几何后验、本站 cavity ext 参数、逐路径时延和角度 Jacobian | 各 $\sigma_{\tau,n,k}^2$、$\sigma_{\theta,n,k}^2$，以及由对角 covariance 得到的 $\boldsymbol\Lambda_{\psi,n}^{\mathrm{in},(l+1)}$ 和 $\boldsymbol\eta_{\psi,n}^{\mathrm{in},(l+1)}$ | 第 5.2 节 |
| `check_message_convergence` | 相邻两轮几何均值和输入消息 | $\delta_z^{(l)},\delta_{\mathrm{msg}}^{(l)}$、停止标志 | 第 5.3 节 |

初始化调用顺序固定为

```text
prepare_observation
coarse_estimation_2d_fft
associate_multi_rx_paths
lm_geometry_initialization
initialize_bayesian_model
```

第 l 轮外层调用顺序固定为

```text
build_parametric_path_matrix(at psi_n^(l-1))
build_path_derivative_matrices(at psi_n^(l-1))
update_angle_delay_posterior
update_scattering_posterior
update_noise_precision_posterior
form_linespectral_extrinsic_message
build_geometry_jacobian(at z^(l-1))
propagate_linespectral_to_geometry
update_geometry_posterior
build_geometry_jacobian(at mu_z^(l))
propagate_geometry_to_linespectral
set psi_n^(l) = g_n(mu_z^(l))
check_message_convergence
```

现有 `JSLES*.m`、`OGVAMP.m` 及位置网格感知矩阵函数实现的是原稀疏网格模型，不能直接视为本概率模型的变分实现。若保留，应在实验中标记为历史基线；当前技术路线的新推断模块应围绕 $K$ 列参数化路径矩阵 $\boldsymbol\Phi_n$、其角度/时延列导数矩阵和局部 Taylor 观测模型重新实现。

## 9. 本版本固定的实现约定

为避免依据本文档生成代码时出现多种不兼容实现，本版本固定采用以下约定：

1. $K$ 为已知场景参数，每个接收站均提取并关联恰好 $K$ 条路径；路径数估计、漏检和虚警处理不属于当前模型。
2. 阵列响应固定使用 $e^{jm\theta}$，时延响应固定使用 $e^{-jq\tau}$，观测向量固定采用 MATLAB 列优先的 $\operatorname{vec}(\mathbf Y_n)$ 顺序。
3. $\boldsymbol\Sigma_{h,n}$ 为用户输入的固定正定对角协方差，本版本不学习该超参数。
4. $\mathbf C_{p,k}^{(0)}$、$\mathbf C_{\mathrm{UE}}^{(0)}$ 和 $\sigma_{\beta,0}^2$ 为实验配置输入，并按第 3.6 节组成固定的 $\mathbf C_{z,0}$；若配置文件给出 $\sigma_{\Delta t,0}^2$，必须先乘以 $c^2$。几何推断保留完整 $\mathbf C_z^{(l)}$，不强制投影成目标间块对角矩阵。
5. 几何内部变量固定使用米单位的距离偏差 $\beta=c\Delta t$，LM 与贝叶斯消息传递中的状态、先验和 Jacobian 均关于 $\beta$ 定义，不对 $\Delta t$ 求导；仅在最终报告结果时使用 $\widehat{\Delta t}=\widehat{\beta}/c$ 及相应协方差变换。
6. 路径矩阵导数和几何 Jacobian 均使用本文给出的解析式；数值差分只用于单元测试，不用于正式迭代。
7. 正则量 $\epsilon_{\mathrm{geo}}$、多站关联参数 $T_{\perp},T_d,T_{\mathrm{cond}},I_{\mathrm{assoc}}$、LM 阻尼参数、$\varepsilon_z$、$\varepsilon_{\mathrm{msg}}$ 和最大迭代次数均为数值配置，不改变概率模型；这些数值必须在实验脚本中显式记录。
8. 角度和时延均假设处于第 4.1 节给出的无模糊区间，ULA 前后向分支由已知可见半平面确定；本版本不执行整数周相位解模糊。
9. 线谱模块输出到几何模块、几何模块返回线谱模块时均使用外信息或 cavity 信息，不能把包含接收方自身旧消息的完整后验直接循环反馈。
10. $N_{\mathrm{Rx}}=2$ 时路径关联使用一次 Hungarian 指派；$N_{\mathrm{Rx}}>2$ 时正式实验默认采用第 4.2 节的锚接收站对初始化、逐站 Hungarian 指派和多射线位置重估流程，并将严格多维指派仅作为小规模验证基准。
11. geo→线谱消息固定假设全部 $\{\tau_{n,k},\theta_{n,k}\}_{k=1}^{K}$ 相互独立；每个分量的标量方差分别由其 Jacobian 行向量按照 $\boldsymbol j\mathbf C_p\boldsymbol j^T$ 计算，最终 $\mathbf V_{\psi,n}^{\mathrm{in}}$ 为对角矩阵，不传播时延—角度或跨路径互协方差。

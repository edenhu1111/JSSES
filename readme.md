# 多站协作 NLoS 感知：非稀疏参数化贝叶斯模型与实验路线

> 本文档根据现阶段实验思路重新整理技术路线。接收信号模型、双接收站几何关系以及发射站位置和定时信息不完美的场景保持不变；主要变化是取消位置网格稀疏建模，不再引入支撑变量或伯努利—高斯稀疏先验，而是直接对每条路径的散射系数、角度、测量时延、目标位置、UE 位置和时钟偏差建立参数化分层贝叶斯模型。

## 1. 研究问题与技术路线

### 1.1 研究场景

考虑发射站位置与定时信息不完美条件下的多站协作上行 ISAC 环境感知：

- 单天线 UE 发送基于 OFDM 的上行 ISAC 导频信号；
- 两个配置 ULA 的接收站共享观测，并由融合中心联合处理；
- UE 与接收站之间存在未知公共定时偏差 $\Delta t$；
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

其中 $\kappa_a=2\pi d/\lambda$，$d$ 为阵元间距，$\lambda$ 为载波波长；当导频子载波按间隔 $D$ 均匀抽取时，$\kappa_\tau=2\pi D\Delta f$。因此 $\theta_{n,k}$ 和 $\tau_{n,k}$ 均为无量纲相位参数，而 $t_{n,k}$ 是包含公共时钟偏差的物理测量时延。

### 1.2 修改后的技术路线

1. 保留原天线域—子载波域接收信号模型，但直接以 $K$ 条物理路径为参数，不再将感兴趣区域离散成位置网格，也不构造位置域稀疏反射向量。
2. 分别对两个接收站的观测执行过采样 2D-FFT，获得每条路径归一化空间相位 $\theta_{n,k}$ 和归一化时延相位 $\tau_{n,k}$ 的粗估计，并反归一化为物理入射角 $\phi_{n,k}$ 和测量时延 $t_{n,k}$。
3. 根据两站物理角度射线的交点和测量时延差完成跨接收站路径关联，得到目标位置粗估计。
4. 将物理入射角和测量时延关于目标位置、UE 位置及公共时钟偏差的非线性几何方程组成残差向量，使用 LM 算法获得
   $\{\boldsymbol p_k^{(0)}\}_{k=1}^{K}$、$\boldsymbol p_{\mathrm{UE}}^{(0)}$ 和 $\Delta t^{(0)}$。
5. 以上述粗估计为几何变量高斯先验的均值，建立关于原始观测 $\boldsymbol y_n$ 的分层贝叶斯模型：
   - $\boldsymbol h_n$ 采用零均值复高斯先验，协方差为对角矩阵；
   - 噪声精度 $\gamma_n$ 采用 Gamma 先验；
   - 归一化空间相位和时延相位通过 Dirac delta 因子与目标位置、UE 位置和时钟偏差保持确定性几何关系；
   - 目标位置、UE 位置和时钟偏差采用以粗估计为均值的高斯先验。
6. 第一部分执行连续线谱贝叶斯估计：以上一轮归一化角度和时延相位为展开点，对参数化路径矩阵作一阶泰勒展开，并更新散射系数、噪声精度以及角度—时延相位的后验消息。
7. 第二部分执行非线性几何消息传递：将线谱模块输出的角度—时延联合消息输入 Dirac 几何因子，推断目标位置、UE 位置和时钟偏差；必要时将几何外信息反馈给线谱模块继续迭代。

整体信息流为

```text
双 Rx 导频观测 Y_1, Y_2
    ↓ 去除已知导频
过采样 2D-FFT
    ↓
归一化相位粗估计 theta 和 tau
    ↓ 反归一化
物理入射角 phi 和测量时延 t
    ↓ 双站射线交汇 + 时延差关联
目标位置粗估计
    ↓ 几何非线性方程组
LM：目标位置 + UE 位置 + 时钟偏差粗估计
    ↓ 作为几何高斯先验均值
非稀疏参数化分层贝叶斯模型
    ├─ 第一部分：Taylor 线性化连续线谱贝叶斯估计
    │    └─ 输出 (theta_n,tau_n) 联合概率消息
    └─ 第二部分：Dirac 几何因子上的非线性消息传递
    ↓
目标位置、UE 位置、时钟偏差及其后验不确定性
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

定义物理入射角和物理测量时延的几何函数

$$
g_{\phi,n,k}(\boldsymbol p_k)
=\arccos\!\left(
\frac{\boldsymbol e_n^T
(\boldsymbol p_k-\boldsymbol p_{\mathrm{Rx},n})}
{\|\boldsymbol p_k-\boldsymbol p_{\mathrm{Rx},n}\|}
\right).
$$

$$
g_{t,n,k}(\boldsymbol p_k,
\boldsymbol p_{\mathrm{UE}},\Delta t)
=\frac{
\|\boldsymbol p_k-\boldsymbol p_{\mathrm{UE}}\|
+\|\boldsymbol p_k-\boldsymbol p_{\mathrm{Rx},n}\|}{c}
+\Delta t.
$$

将物理尺度吸收到归一化相位后，定义

$$
g_{\theta,n,k}(\boldsymbol p_k)
:=\kappa_a\cos\!\left(g_{\phi,n,k}(\boldsymbol p_k)\right)
=\kappa_a
\frac{\boldsymbol e_n^T
(\boldsymbol p_k-\boldsymbol p_{\mathrm{Rx},n})}
{\|\boldsymbol p_k-\boldsymbol p_{\mathrm{Rx},n}\|},
$$

$$
g_{\tau,n,k}(\boldsymbol p_k,
\boldsymbol p_{\mathrm{UE}},\Delta t)
:=\kappa_\tau
g_{t,n,k}(\boldsymbol p_k,
\boldsymbol p_{\mathrm{UE}},\Delta t).
$$

因此，进入导向向量的两个无量纲参数满足

$$
\theta_{n,k}=g_{\theta,n,k}(\boldsymbol p_k),
\qquad
\tau_{n,k}=g_{\tau,n,k}(
\boldsymbol p_k,\boldsymbol p_{\mathrm{UE}},\Delta t).
$$

### 2.3 非稀疏参数化观测模型

定义第 $n$ 个接收站的参数化路径矩阵

$$
\mathbf A_n(\boldsymbol\theta_n,\boldsymbol\tau_n)
=\left[
\boldsymbol a(\theta_{n,1})\otimes\boldsymbol b(\tau_{n,1}),
\ldots,
\boldsymbol a(\theta_{n,K})\otimes\boldsymbol b(\tau_{n,K})
\right].
$$

去除已知恒模导频后，第 $n$ 个接收站的观测为

$$
\boldsymbol y_n
=\mathbf A_n(\boldsymbol\theta_n,\boldsymbol\tau_n)
\boldsymbol h_n+\boldsymbol w_n,
$$

其中

$$
\boldsymbol w_n\sim
\mathcal{CN}(\boldsymbol 0,\gamma_n^{-1}\mathbf I_M),
\qquad M=N_aN_p.
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
\mathbf A_n(\boldsymbol\theta_n,\boldsymbol\tau_n)
\boldsymbol h_n,
\gamma_n^{-1}\mathbf I_M
\right).
$$

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

对角协方差表示同一接收站内各路径散射系数先验独立，但允许不同路径具有不同先验功率。两个接收站的散射系数可条件独立建模，因为同一目标在不同接收站上的复反射系数不要求相干一致：

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

归一化时延相位由目标位置、UE 位置和公共时钟偏差共同决定：

$$
p(\boldsymbol T\mid
\mathbf P,\boldsymbol p_{\mathrm{UE}},\Delta t)
=\prod_{n=1}^{N_{\mathrm{Rx}}}
\prod_{k=1}^{K}
\delta\!\left(
\tau_{n,k}
-g_{\tau,n,k}(\boldsymbol p_k,
\boldsymbol p_{\mathrm{UE}},\Delta t)
\right).
$$

因为 $\tau_{n,k}=\kappa_\tau t_{n,k}$ 已包含公共时钟偏差对应的归一化相位，所以在 $\boldsymbol b(\tau_{n,k})$ 中不能再次叠加 $\kappa_\tau\Delta t$。秒单位测量时延由 $t_{n,k}=\tau_{n,k}/\kappa_\tau$ 恢复。

### 3.6 几何变量的高斯先验

令 2D-FFT、双站关联和 LM 给出的粗估计为

$$
\boldsymbol p_k^{(0)},\qquad
\boldsymbol p_{\mathrm{UE}}^{(0)},\qquad
\Delta t^{(0)}.
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

UE 位置与时钟偏差先验分别为

$$
p(\boldsymbol p_{\mathrm{UE}})
=\mathcal N\!\left(
\boldsymbol p_{\mathrm{UE}};
\boldsymbol p_{\mathrm{UE}}^{(0)},
\mathbf C_{\mathrm{UE}}^{(0)}
\right),
$$

$$
p(\Delta t)
=\mathcal N\!\left(
\Delta t;\Delta t^{(0)},
\sigma_{\Delta t,0}^{2}
\right).
$$

将全部几何变量堆叠为

$$
\boldsymbol z=
\left[
\boldsymbol p_1^T,\ldots,\boldsymbol p_K^T,
\boldsymbol p_{\mathrm{UE}}^T,\Delta t
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
\Delta t^{(0)}
\right]^T,
$$

$$
\mathbf C_{z,0}
=\operatorname{blkdiag}\!\left(
\mathbf C_{p,1}^{(0)},\ldots,
\mathbf C_{p,K}^{(0)},
\mathbf C_{\mathrm{UE}}^{(0)},
\sigma_{\Delta t,0}^{2}
\right).
$$

这些协方差决定贝叶斯精化对粗估计的信任程度。协方差过小会使后验难以修正粗估计偏差；协方差过大则会削弱 2D-FFT 和 LM 初始化的稳定作用。

### 3.7 完整联合概率分布

在几何先验相互独立的假设下，记

$$
p(\boldsymbol z)
=p(\mathbf P)
p(\boldsymbol p_{\mathrm{UE}})
p(\Delta t).
$$

完整联合概率分布可写为

$$
\begin{aligned}
&p(\mathcal Y,\mathcal H,\boldsymbol\Theta,
\boldsymbol T,\boldsymbol\gamma,
\mathbf P,\boldsymbol p_{\mathrm{UE}},\Delta t)
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
\mathbf P,\boldsymbol p_{\mathrm{UE}},\Delta t)
p(\mathbf P)
p(\boldsymbol p_{\mathrm{UE}})
p(\Delta t).
\end{aligned}
$$

目标是计算或近似后验分布

$$
p(\mathcal H,\boldsymbol\gamma,
\mathbf P,\boldsymbol p_{\mathrm{UE}},\Delta t
\mid\mathcal Y).
$$

## 4. Dirac 因子图及其等价后验

Dirac delta 因子将 $\boldsymbol\theta_n$ 和 $\boldsymbol\tau_n$ 严格限制为 $\boldsymbol z$ 的确定性函数。从概率分布角度，可将其直接代入路径矩阵：

$$
\widetilde{\mathbf A}_n(\boldsymbol z)
=\mathbf A_n\!\left(
\boldsymbol g_{\theta,n}(\mathbf P),
\boldsymbol g_{\tau,n}(
\mathbf P,\boldsymbol p_{\mathrm{UE}},\Delta t)
\right).
$$

于是等价观测模型为

$$
\boldsymbol y_n
=\widetilde{\mathbf A}_n(\boldsymbol z)
\boldsymbol h_n+\boldsymbol w_n,
$$

后验可简化为

$$
\begin{aligned}
p(\mathcal H,\boldsymbol\gamma,\boldsymbol z\mid\mathcal Y)
\propto{}&
\left\{\prod_{n=1}^{N_{\mathrm{Rx}}}
\mathcal{CN}\!\left(
\boldsymbol y_n;
\widetilde{\mathbf A}_n(\boldsymbol z)\boldsymbol h_n,
\gamma_n^{-1}\mathbf I_M
\right)
p(\boldsymbol h_n)p(\gamma_n)
\right\}
p(\boldsymbol z).
\end{aligned}
$$

这是完整模型对应的等价后验，可用于核对两模块消息传递是否与原概率模型一致。实际算法保留
$\boldsymbol\theta_n$、$\boldsymbol\tau_n$ 和 Dirac 因子作为线谱模块与几何模块之间的接口。线谱模块可以维护
$q(\boldsymbol\theta_n,\boldsymbol\tau_n)$ 的局部近似并将其作为消息输出，但该分布不是与
$\boldsymbol z$ 无关的全局均值场因子；它必须通过 Dirac 几何因子与
$p(\boldsymbol z)$ 交换消息。

2D-FFT 输出的 $\widehat{\boldsymbol\theta}_n$ 和 $\widehat{\boldsymbol\tau}_n$ 是由同一原始观测 $\boldsymbol y_n$ 得到的初始化量，在当前模型中不再作为独立观测似然与 $p(\boldsymbol y_n\mid\cdot)$ 重复相乘，否则会重复使用同一份数据。若后续确实希望在贝叶斯阶段继续使用 FFT 测量值，则需要另外引入带有限协方差的 $p(\widehat{\boldsymbol\theta}_n\mid\boldsymbol\theta_n)$ 和 $p(\widehat{\boldsymbol\tau}_n\mid\boldsymbol\tau_n)$，此时它们不应写成几何 Dirac 因子。

## 5. 2D-FFT 与 LM 粗初始化

### 5.1 过采样 2D-FFT

按“子载波索引为行、阵元索引为列”将第 $n$ 个接收站的观测恢复成矩阵。在归一化参数下，采用 $\boldsymbol a(\theta)\propto[e^{jm\theta}]_m$ 和 $\boldsymbol b(\tau)\propto[e^{-jq\tau}]_q$ 时，单条路径对应

$$
Y_n[q,m]\propto
e^{jm\theta_{n,k}}
e^{-jq\tau_{n,k}}.
$$

因此，阵元维使用 FFT，子载波维使用 IFFT，并对有正负号的空间相位维执行 fftshift。峰值索引首先映射到无量纲相位参数：

$$
\widehat\theta_l
=\frac{2\pi l}{N_{\theta}^{\mathrm{FFT}}},
\qquad
\widehat\tau_p
=\frac{2\pi p}{N_{\tau}^{\mathrm{FFT}}}.
$$

其中 $l$ 为 fftshift 后的有符号空间频率索引，$p$ 为时延维索引。对应的物理角度和秒单位测量时延由反归一化得到：

$$
\widehat\phi_{n,k}
=\arccos\!\left(
\frac{\widehat\theta_{n,k}}{\kappa_a}
\right),
\qquad
\widehat t_{n,k}
=\frac{\widehat\tau_{n,k}}{\kappa_\tau}.
$$

当 $d\leq\lambda/2$ 时，主值区间 $\theta\in[-\pi,\pi]$ 内不存在空间混叠；均匀抽取导频时，无模糊时延范围为 $0\leq t<1/(D\Delta f)$。超出上述范围时，2D-FFT 只能得到相位折叠后的参数，需要另行进行解模糊处理。

当前实现默认采用 8 倍过采样、二维非极大值抑制和局部对数谱抛物线插值。零填充降低 FFT 参数栅格的量化误差，但不会突破阵列孔径和有效带宽决定的物理分辨率。

现有 `CoarseEstimation2DFFT.m` 针对代码中的负指数阵列约定沿阵元维使用 IFFT，并直接输出方向余弦和秒单位时延。若后续将代码接口也统一为本节的正指数归一化参数，应同步修改变换方向和输出映射。

### 5.2 双站关联与目标位置粗估计

1. 将归一化空间相位除以 $\kappa_a$，恢复方向余弦和全局物理角度；
2. 由每个角度从对应接收站发出射线；
3. 两站射线两两求交，排除平行射线和反向交点；
4. 对候选交点比较测量距离差与接收站—目标几何距离差；
5. 使用一对一匹配得到跨接收站路径对应关系和目标位置粗估计 $\boldsymbol p_k^{(0)}$。

距离差中 UE—目标距离和公共时钟偏差均被抵消，因此适合作为两站路径关联的物理一致性量。

### 5.3 LM 几何粗估计

以关联后的目标位置、UE 粗位置和零时钟偏差为初值，构造

$$
r_{n,k}(\boldsymbol z)
=\widehat t_{n,k}
-\frac{
\|\boldsymbol p_k-\boldsymbol p_{\mathrm{UE}}\|
+\|\boldsymbol p_k-\boldsymbol p_{\mathrm{Rx},n}\|}{c}
-\Delta t.
$$

LM 求解

$$
\boldsymbol z^{(0)}
=\arg\min_{\boldsymbol z}
\sum_{n=1}^{N_{\mathrm{Rx}}}
\sum_{k=1}^{K}r_{n,k}^2(\boldsymbol z).
$$

代码中可改用距离偏差 $b=c\Delta t$ 作为优化变量，使其与位置变量具有相近尺度。每轮仅在残差下降时接受试探步，并相应减小阻尼；否则增大阻尼后重新计算。

在当前分阶段实现中，目标位置先由 AoA 射线固定。此时同一目标的两个接收站观测对 UE 位置和公共时钟偏差只提供一个独立和式约束，因此至少需要 3 个几何分散的目标，并应检查 LM 雅可比矩阵秩。若秩亏，不能将 LM 输出解释为唯一的 UE—时钟解。

## 6. 两部分贝叶斯推断与消息传递

### 6.1 路径矩阵关于角度和时延的一阶展开

为与泰勒展开的记号一致，以下记

$$
\boldsymbol\Phi_n(\boldsymbol\theta_n,\boldsymbol\tau_n)
:=\mathbf A_n(\boldsymbol\theta_n,\boldsymbol\tau_n).
$$

第 $l$ 轮推断以上一轮后验均值对应的
$\boldsymbol\theta_n^{(l-1)}$ 和
$\boldsymbol\tau_n^{(l-1)}$ 为展开点，并定义

$$
\Delta\boldsymbol\theta_n^{(l)}
=\boldsymbol\theta_n-\boldsymbol\theta_n^{(l-1)},
\qquad
\Delta\boldsymbol\tau_n^{(l)}
=\boldsymbol\tau_n-\boldsymbol\tau_n^{(l-1)}.
$$

由于 $\boldsymbol\Phi_n$ 的第 $k$ 列只依赖
$(\theta_{n,k},\tau_{n,k})$，定义两个 $M\times K$ 的列导数矩阵

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
=\boldsymbol a(\theta_{n,k})
\otimes\boldsymbol b(\tau_{n,k}).
$$

于是相应的一阶泰勒近似可严格写为

$$
\begin{aligned}
\boldsymbol\Phi_n(\boldsymbol\theta_n,\boldsymbol\tau_n)
\approx{}&
\boldsymbol\Phi_n^{(l-1)}
+\boldsymbol\Phi_{\theta,n}^{(l-1)}
\operatorname{diag}\!\left(
\Delta\boldsymbol\theta_n^{(l)}
\right)
\\
&+\boldsymbol\Phi_{\tau,n}^{(l-1)}
\operatorname{diag}\!\left(
\Delta\boldsymbol\tau_n^{(l)}
\right),
\end{aligned}
$$

其中

$$
\boldsymbol\Phi_n^{(l-1)}
=\boldsymbol\Phi_n\!\left(
\boldsymbol\theta_n^{(l-1)},
\boldsymbol\tau_n^{(l-1)}
\right).
$$

注意，$\partial\boldsymbol\Phi_n/\partial\boldsymbol\theta_n$
若按一般矩阵对向量求导会形成三阶张量；上述列导数定义利用了各路径列之间的参数独立性，使泰勒式中的两个右乘对角矩阵具有明确维度。

令

$$
\mathbf D_a=\operatorname{diag}(0,1,\ldots,N_a-1),
\qquad
\mathbf D_b=\operatorname{diag}(0,1,\ldots,N_p-1),
$$

则归一化导向向量的解析导数为

$$
\frac{\partial\boldsymbol a(\theta)}{\partial\theta}
=j\mathbf D_a\boldsymbol a(\theta),
\qquad
\frac{\partial\boldsymbol b(\tau)}{\partial\tau}
=-j\mathbf D_b\boldsymbol b(\tau),
$$

从而

$$
\frac{\partial\boldsymbol\phi_{n,k}}{\partial\theta_{n,k}}
=
\left(j\mathbf D_a\boldsymbol a(\theta_{n,k})\right)
\otimes\boldsymbol b(\tau_{n,k}),
$$

$$
\frac{\partial\boldsymbol\phi_{n,k}}{\partial\tau_{n,k}}
=
\boldsymbol a(\theta_{n,k})\otimes
\left(-j\mathbf D_b\boldsymbol b(\tau_{n,k})\right).
$$

### 6.2 第一部分：连续线谱贝叶斯估计

第一部分只处理观测模型中的连续线谱参数
$\boldsymbol h_n$、$\boldsymbol\theta_n$、
$\boldsymbol\tau_n$ 和 $\gamma_n$，不在该模块中对位置或时钟变量作 Taylor 展开。定义联合相位增量

$$
\boldsymbol\xi_n^{(l)}
:=
\begin{bmatrix}
\Delta\boldsymbol\theta_n^{(l)}\\
\Delta\boldsymbol\tau_n^{(l)}
\end{bmatrix}
\in\mathbb R^{2K},
$$

以及给定散射系数时的线性化矩阵

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

由第 6.1 节的一阶展开，观测可近似写为

$$
\boldsymbol y_n
\approx
\boldsymbol\Phi_n^{(l-1)}\boldsymbol h_n
+\mathbf G_n^{(l-1)}(\boldsymbol h_n)
\boldsymbol\xi_n^{(l)}
+\boldsymbol w_n.
$$

等价地，定义

$$
\widehat{\boldsymbol\Phi}_n^{(l)}
(\boldsymbol\xi_n)
:=
\boldsymbol\Phi_n^{(l-1)}
+\boldsymbol\Phi_{\theta,n}^{(l-1)}
\operatorname{diag}(\Delta\boldsymbol\theta_n)
+\boldsymbol\Phi_{\tau,n}^{(l-1)}
\operatorname{diag}(\Delta\boldsymbol\tau_n),
$$

则

$$
\boldsymbol y_n
\approx
\widehat{\boldsymbol\Phi}_n^{(l)}
(\boldsymbol\xi_n)\boldsymbol h_n
+\boldsymbol w_n.
$$

线谱模块采用

$$
q(\boldsymbol h_n,\gamma_n,\boldsymbol\xi_n)
=q(\boldsymbol h_n)q(\gamma_n)q(\boldsymbol\xi_n)
$$

的局部均值场分解。来自几何模块的输入消息用信息形式表示为

$$
m_{\mathrm{geo}\rightarrow\xi_n}(\boldsymbol\xi_n)
\propto
\exp\!\left(
-\frac{1}{2}\boldsymbol\xi_n^T
\boldsymbol\Lambda_{\xi,n}^{\mathrm{in}}\boldsymbol\xi_n
+\boldsymbol\eta_{\xi,n}^{\mathrm{in}\,T}
\boldsymbol\xi_n
\right).
$$

首次执行线谱估计时，可令
$\boldsymbol\Lambda_{\xi,n}^{\mathrm{in}}=\mathbf 0$，
仅将 2D-FFT 结果作为 Taylor 展开中心；也可将几何高斯先验通过 Dirac 因子传播后形成初始输入消息。不能把由同一
$\boldsymbol y_n$ 得到的 2D-FFT 峰值再次作为独立概率观测相乘。

### 6.3 线谱模块的变分更新与输出消息

散射系数和相位增量的变分因子分别取为

$$
q(\boldsymbol h_n)
=\mathcal{CN}(\boldsymbol h_n;
\boldsymbol m_{h,n},\mathbf C_{h,n}),
\qquad
q(\boldsymbol\xi_n)
=\mathcal N(\boldsymbol\xi_n;
\boldsymbol m_{\xi,n},\mathbf C_{\xi,n}).
$$

在第 $l$ 轮局部线性模型下，散射系数更新为

$$
\widehat{\boldsymbol\Phi}_n^{(l)}
=\boldsymbol\Phi_n^{(l-1)}
+\sum_{r=1}^{2K}
\mathbf F_{n,r}^{(l-1)}\xi_{n,r},
$$

其中，令 $\boldsymbol e_r$ 表示第 $r$ 个 $K$ 维标准基向量，则

$$
\mathbf F_{n,r}^{(l-1)}
=
\begin{cases}
\boldsymbol\Phi_{\theta,n}^{(l-1)}
\operatorname{diag}(\boldsymbol e_r),
&1\leq r\leq K,\\
\boldsymbol\Phi_{\tau,n}^{(l-1)}
\operatorname{diag}(\boldsymbol e_{r-K}),
&K<r\leq 2K.
\end{cases}
$$

若当前
$q(\boldsymbol\xi_n)
=\mathcal N(\boldsymbol m_{\xi,n},\mathbf C_{\xi,n})$，
则局部路径矩阵的一阶矩为

$$
\mathbb E[
\widehat{\boldsymbol\Phi}_n^{(l)}
]
=
\boldsymbol\Phi_n^{(l-1)}
+\sum_r
\mathbf F_{n,r}^{(l-1)}m_{\xi,n,r},
$$

二阶矩为

$$
\begin{aligned}
\mathbb E[
\widehat{\boldsymbol\Phi}_n^{(l)H}
\widehat{\boldsymbol\Phi}_n^{(l)}
]
={}&
\boldsymbol\Phi_n^{(l-1)H}
\boldsymbol\Phi_n^{(l-1)}
\\
&+\sum_r m_{\xi,n,r}
\left(
\mathbf F_{n,r}^{(l-1)H}\boldsymbol\Phi_n^{(l-1)}
+\boldsymbol\Phi_n^{(l-1)H}\mathbf F_{n,r}^{(l-1)}
\right)
\\
&+\sum_{r,s}
\left(
[\mathbf C_{\xi,n}]_{r,s}
+m_{\xi,n,r}m_{\xi,n,s}
\right)
\mathbf F_{n,r}^{(l-1)H}
\mathbf F_{n,s}^{(l-1)}.
\end{aligned}
$$

因此

$$
\mathbf C_{h,n}^{-1}
=\mathbb E[\gamma_n]\,
\mathbb E_{q(\boldsymbol\xi_n)}
\left[
\widehat{\boldsymbol\Phi}_n^{(l)H}
\widehat{\boldsymbol\Phi}_n^{(l)}
\right]
+\boldsymbol\Sigma_{h,n}^{-1},
$$

$$
\boldsymbol m_{h,n}
=\mathbf C_{h,n}
\mathbb E[\gamma_n]\,
\mathbb E_{q(\boldsymbol\xi_n)}
\left[
\widehat{\boldsymbol\Phi}_n^{(l)H}
\right]\boldsymbol y_n.
$$

由于 $\widehat{\boldsymbol\Phi}_n^{(l)}$ 关于
$\boldsymbol\xi_n$ 是仿射函数，以上一阶和二阶矩阵期望只依赖
$q(\boldsymbol\xi_n)$ 的均值与协方差，可以解析计算。

噪声精度仍为 Gamma 分布：

$$
q(\gamma_n)
=\operatorname{Ga}(\gamma_n;a_{\gamma,n},b_{\gamma,n}),
\qquad
a_{\gamma,n}=a_{\gamma,0}+M,
$$

$$
b_{\gamma,n}
=b_{\gamma,0}
+\mathbb E_q\!\left[
\left\|
\boldsymbol y_n-
\widehat{\boldsymbol\Phi}_n^{(l)}
(\boldsymbol\xi_n)
\boldsymbol h_n
\right\|^2
\right],
\qquad
\mathbb E[\gamma_n]
=\frac{a_{\gamma,n}}{b_{\gamma,n}}.
$$

期望残差应同时保留散射系数协方差和相位增量协方差产生的二阶矩，不能只将
$\boldsymbol h_n$ 与 $\boldsymbol\xi_n$ 替换为后验均值。具体地，

$$
\begin{aligned}
&\mathbb E_q\!\left[
\left\|
\boldsymbol y_n-
\widehat{\boldsymbol\Phi}_n^{(l)}
\boldsymbol h_n
\right\|^2
\right]
\\
={}&
\boldsymbol y_n^H\boldsymbol y_n
-2\operatorname{Re}\!\left\{
\boldsymbol y_n^H
\mathbb E[\widehat{\boldsymbol\Phi}_n^{(l)}]
\boldsymbol m_{h,n}
\right\}
\\
&+
\operatorname{tr}\!\left(
\mathbb E[
\widehat{\boldsymbol\Phi}_n^{(l)H}
\widehat{\boldsymbol\Phi}_n^{(l)}
]
\left[
\mathbf C_{h,n}
+\boldsymbol m_{h,n}\boldsymbol m_{h,n}^H
\right]
\right).
\end{aligned}
$$

给定当前 $q(\boldsymbol h_n)$ 和 $q(\gamma_n)$，相位增量的局部高斯更新为

$$
\begin{aligned}
\mathbf C_{\xi,n}^{-1}
={}&
\boldsymbol\Lambda_{\xi,n}^{\mathrm{in}}
\\
&+2\mathbb E[\gamma_n]\,
\operatorname{Re}\!\left\{
\mathbb E_{q(\boldsymbol h_n)}
\left[
\mathbf G_n^{(l-1)H}\mathbf G_n^{(l-1)}
\right]
\right\},
\end{aligned}
$$

$$
\begin{aligned}
\boldsymbol m_{\xi,n}
=\mathbf C_{\xi,n}
\Bigg[
&\boldsymbol\eta_{\xi,n}^{\mathrm{in}}
\\
&+2\mathbb E[\gamma_n]\,
\operatorname{Re}\!\left\{
\mathbb E_{q(\boldsymbol h_n)}
\left[
\mathbf G_n^{(l-1)H}
\left(
\boldsymbol y_n-
\boldsymbol\Phi_n^{(l-1)}\boldsymbol h_n
\right)
\right]
\right\}
\Bigg].
\end{aligned}
$$

这里的因子 2 来自圆对称复高斯似然和实值
$\boldsymbol\xi_n$ 的组合。线谱模块得到的绝对相位均值为

$$
\boldsymbol\mu_{\theta,n}^{\mathrm{post}}
=\boldsymbol\theta_n^{(l-1)}
+[\boldsymbol m_{\xi,n}]_{1:K},
\qquad
\boldsymbol\mu_{\tau,n}^{\mathrm{post}}
=\boldsymbol\tau_n^{(l-1)}
+[\boldsymbol m_{\xi,n}]_{K+1:2K},
$$

并保留 $\mathbf C_{\xi,n}$ 中角度、时延及不同路径之间的相关性。为避免在两个模块之间重复使用先验信息，发送给几何模块的外信息应从后验信息参数中扣除输入消息：

$$
\boldsymbol\Lambda_{\xi,n}^{\mathrm{ext}}
=\mathbf C_{\xi,n}^{-1}
-\boldsymbol\Lambda_{\xi,n}^{\mathrm{in}},
\qquad
\boldsymbol\eta_{\xi,n}^{\mathrm{ext}}
=\mathbf C_{\xi,n}^{-1}\boldsymbol m_{\xi,n}
-\boldsymbol\eta_{\xi,n}^{\mathrm{in}}.
$$

记当前绝对相位展开中心为

$$
\boldsymbol x_n^{(l-1)}
=
\begin{bmatrix}
\boldsymbol\theta_n^{(l-1)}\\
\boldsymbol\tau_n^{(l-1)}
\end{bmatrix}.
$$

当 $\boldsymbol\Lambda_{\xi,n}^{\mathrm{ext}}$ 非奇异时，发送给几何模块的高斯外信息参数为

$$
\mathbf V_{x,n}^{\mathrm{LS}}
=\left(
\boldsymbol\Lambda_{\xi,n}^{\mathrm{ext}}
\right)^{-1},
\qquad
\boldsymbol\mu_{x,n}^{\mathrm{LS}}
=\boldsymbol x_n^{(l-1)}
+\mathbf V_{x,n}^{\mathrm{LS}}
\boldsymbol\eta_{\xi,n}^{\mathrm{ext}}.
$$

若外信息精度矩阵奇异，则应保留其信息形式或使用定义在可辨识子空间上的伪逆，不能人为加入高精度正则项制造虚假的确定性。

在接受当前线谱增量后，更新
$\boldsymbol\theta_n^{(l)}=
\boldsymbol\theta_n^{(l-1)}
+[\boldsymbol m_{\xi,n}]_{1:K}$ 和
$\boldsymbol\tau_n^{(l)}=
\boldsymbol\tau_n^{(l-1)}
+[\boldsymbol m_{\xi,n}]_{K+1:2K}$，
并在新的相位中心重新计算
$\boldsymbol\Phi_n$、$\boldsymbol\Phi_{\theta,n}$ 和
$\boldsymbol\Phi_{\tau,n}$。若增量过大，应采用阻尼或信赖域保证一阶近似有效。

### 6.4 第二部分：Dirac 几何因子上的非线性消息传递

第二部分不再对几何映射作一阶 Taylor 展开，而是保留第 3.4 和 3.5 节中的非线性 Dirac 因子。定义

$$
\boldsymbol x_n
:=
\begin{bmatrix}
\boldsymbol\theta_n\\
\boldsymbol\tau_n
\end{bmatrix},
\qquad
\boldsymbol g_n(\boldsymbol z)
:=
\begin{bmatrix}
\boldsymbol g_{\theta,n}(\mathbf P)\\
\boldsymbol g_{\tau,n}(
\mathbf P,\boldsymbol p_{\mathrm{UE}},\Delta t)
\end{bmatrix},
$$

以及联合几何因子

$$
f_n(\boldsymbol x_n,\boldsymbol z)
=\delta\!\left(
\boldsymbol x_n-\boldsymbol g_n(\boldsymbol z)
\right).
$$

将线谱模块输出的外信息转换成关于绝对相位的联合消息

$$
m_{\mathrm{LS},n}(\boldsymbol x_n)
\propto
\mathcal N(
\boldsymbol x_n;
\boldsymbol\mu_{x,n}^{\mathrm{LS}},
\mathbf V_{x,n}^{\mathrm{LS}}).
$$

该高斯消息定义在当前无模糊局部分支上；构造消息前应以 Taylor 中心为参考对
$\theta$ 和 $\tau$ 解缠，避免在 $\pm\pi$ 或 $2\pi$ 边界产生虚假的大增量。

这里应尽量保留 $\boldsymbol\theta_n$ 与
$\boldsymbol\tau_n$ 的交叉协方差，而不是将其拆成互相独立的两个标量似然。Dirac 因子向几何变量发送的精确消息为

$$
\begin{aligned}
m_{f_n\rightarrow z}(\boldsymbol z)
&=
\int
f_n(\boldsymbol x_n,\boldsymbol z)
m_{\mathrm{LS},n}(\boldsymbol x_n)
\mathrm d\boldsymbol x_n
\\
&=
m_{\mathrm{LS},n}\!\left(
\boldsymbol g_n(\boldsymbol z)
\right).
\end{aligned}
$$

因此几何变量的后验信念满足

$$
b(\boldsymbol z)
\propto
p(\boldsymbol z)
\prod_{n=1}^{N_{\mathrm{Rx}}}
m_{\mathrm{LS},n}\!\left(
\boldsymbol g_n(\boldsymbol z)
\right).
$$

若线谱外信息采用高斯表示，则

$$
\begin{aligned}
\log b(\boldsymbol z)
= {}&
\log p(\boldsymbol z)
\\
&-\frac{1}{2}
\sum_n
\left[
\boldsymbol g_n(\boldsymbol z)
-\boldsymbol\mu_{x,n}^{\mathrm{LS}}
\right]^T
\left(\mathbf V_{x,n}^{\mathrm{LS}}\right)^{-1}
\left[
\boldsymbol g_n(\boldsymbol z)
-\boldsymbol\mu_{x,n}^{\mathrm{LS}}
\right]
+\mathrm{const}.
\end{aligned}
$$

由于 $\boldsymbol g_n(\boldsymbol z)$ 包含方向余弦、欧氏距离和公共时钟偏差，上述消息关于
$\boldsymbol z$ 通常不是高斯分布。可采用粒子消息传递、sigma-point 消息传递或局部 Laplace 高斯投影计算其均值与协方差。此处非线性几何关系仍由 Dirac 因子精确表达；近似只发生在消息的数值表示和投影步骤。

若只执行一次两阶段估计，则直接由 $b(\boldsymbol z)$ 输出目标位置、UE 位置和时钟偏差的后验。若需要在线谱域与几何域之间迭代，则几何因子返回给相位变量的消息为

$$
m_{f_n\rightarrow x_n}(\boldsymbol x_n)
\propto
\int
\delta\!\left(
\boldsymbol x_n-\boldsymbol g_n(\boldsymbol z)
\right)
m_{z\rightarrow f_n}(\boldsymbol z)
\mathrm d\boldsymbol z,
$$

即几何空腔信念经过非线性映射
$\boldsymbol g_n$ 后的推前分布。这里的
$m_{z\rightarrow f_n}$ 必须排除接收站 $n$ 自身传入的
$m_{f_n\rightarrow z}$，以免信息自反馈。

在以 $\boldsymbol z$ 为联合几何变量节点的表示下，该空腔消息为

$$
m_{z\rightarrow f_n}(\boldsymbol z)
\propto
p(\boldsymbol z)
\prod_{\substack{m=1\\m\neq n}}^{N_{\mathrm{Rx}}}
m_{f_m\rightarrow z}(\boldsymbol z).
$$

将该推前分布近似为
$\mathcal N(\boldsymbol\mu_{x,n}^{\mathrm{geo}},
\mathbf V_{x,n}^{\mathrm{geo}})$ 后，相对于当前 Taylor 中心

$$
\boldsymbol x_n^{(l-1)}
=
\begin{bmatrix}
\boldsymbol\theta_n^{(l-1)}\\
\boldsymbol\tau_n^{(l-1)}
\end{bmatrix}
$$

的线谱输入消息参数为

$$
\boldsymbol\Lambda_{\xi,n}^{\mathrm{in}}
=\left(\mathbf V_{x,n}^{\mathrm{geo}}\right)^{-1},
\qquad
\boldsymbol\eta_{\xi,n}^{\mathrm{in}}
=\boldsymbol\Lambda_{\xi,n}^{\mathrm{in}}
\left(
\boldsymbol\mu_{x,n}^{\mathrm{geo}}
-\boldsymbol x_n^{(l-1)}
\right).
$$

这样，Taylor 展开只服务于第一部分的连续线谱估计，而目标位置、UE 位置和时钟偏差始终通过非线性 Dirac 几何因子进行消息传递。

### 6.5 统一推断流程

```text
输入：双 Rx 原始观测、Rx 几何、OFDM/阵列参数、UE 粗位置

1. 去除已知导频。
2. 对每个 Rx 执行过采样 2D-FFT，提取 K 个归一化空间相位—时延相位峰值。
3. 反归一化得到物理角度和秒单位测量时延，再使用射线交点和距离差完成双 Rx 路径关联。
4. 用 LM 求得目标、UE 和时钟偏差的粗估计 z^(0)。
5. 以 z^(0) 为几何高斯先验均值，并初始化 q(h_n) 和 q(gamma_n)。

第一部分：连续线谱贝叶斯估计
6. 以 2D-FFT 的 theta_n^(0)、tau_n^(0) 为首次 Taylor 展开中心。
7. 计算 Phi_n、Phi_theta,n 和 Phi_tau,n，构造关于 xi_n=[Delta theta_n;Delta tau_n] 的局部线性模型。
8. 交替更新 q(h_n)、q(gamma_n) 和 q(xi_n)，直至线谱模块收敛。
9. 更新 Taylor 中心并按需重新线性化，输出保留角度—时延交叉协方差的联合外信息。

第二部分：非线性几何消息传递
10. 将线谱联合外信息送入 p(theta|P) 和 p(tau|P,p_UE,Delta t) 对应的 Dirac 因子。
11. 在非线性几何因子图上融合各 Rx 消息与几何高斯先验，得到 b(z)。
12. 使用粒子、sigma-point 或局部 Laplace 投影表示非高斯几何消息。
13. 默认由 b(z) 输出目标位置、UE 位置和时钟偏差的后验。
14. 若采用迭代式两模块推断，则由几何空腔信念生成外信息并反馈给线谱模块，直至消息或后验均值收敛。
```

## 7. 可辨识性与数值稳定性

### 7.1 路径标签与数据关联

参数化模型要求两个接收站的第 $k$ 条路径对应同一目标。路径排列本身具有置换不变性，因此必须通过粗估计阶段的数据关联固定标签，或在贝叶斯推断中显式处理路径置换。当前方案采用前者。

### 7.2 几何可辨识性

- 目标数量只是可辨识性的必要条件之一；
- 目标、UE 与接收站的退化几何会导致雅可比秩亏或条件数过大；
- 两个目标角度或时延过近时，2D-FFT 无法稳定分离对应路径；
- UE 位置与公共时钟偏差可能存在强相关后验，不能只报告边缘点估计而忽略协方差。

建议记录几何可辨识性矩阵的最小奇异值、条件数以及 $b(\boldsymbol z)$ 中 UE 位置与 $\Delta t$ 的后验相关系数。

### 7.3 数值检查

- 构造 $\arccos$ 前将方向余弦裁剪到 $[-1,1]$；
- 检查所有目标到 UE 和接收站的距离分母不为零；
- 保证 $\gamma_n>0$，并统一 Gamma 分布的 rate/scale 约定；
- 保证 $\boldsymbol\Sigma_{h,n}$ 和几何先验协方差为 Hermitian/实对称正定矩阵；
- 用线性方程求解代替显式矩阵求逆；
- 对 $\mathbf C_{h,n}$、$\mathbf C_{\xi,n}$ 及采用高斯投影时的几何消息协方差使用 Cholesky 分解或对称化处理；
- 以当前 Taylor 中心为参考对周期相位增量解缠；
- 每次线谱相位中心更新后重新计算 $\boldsymbol\Phi_n$、$\boldsymbol\Phi_{\theta,n}$ 和 $\boldsymbol\Phi_{\tau,n}$；
- 检查 $\Delta\boldsymbol\theta_n$ 和 $\Delta\boldsymbol\tau_n$ 是否位于 Taylor 信赖域，并记录阻尼系数和局部 ELBO 变化；
- 几何消息传递必须使用空腔消息或外信息，避免在两个模块之间重复累计同一信息；
- 对非线性 Dirac 因子的推前消息和反向消息检查归一化、有限矩以及多峰性；
- 对解析矩阵矩与数值采样结果做小规模交叉验证，确认协方差项实现正确。

## 8. 实验整理

### 8.1 固定配置

- [ ] UE、两个接收站和 $K$ 个目标的二维坐标；
- [ ] 两个 ULA 的朝向、阵元数、阵元间距和载频；
- [ ] 导频子载波索引、子载波间隔、均匀抽取间隔和导频符号；
- [ ] 散射系数协方差 $\boldsymbol\Sigma_{h,n}$；
- [ ] Gamma 先验参数 $a_{\gamma,0},b_{\gamma,0}$；
- [ ] 2D-FFT 过采样倍数、峰值抑制半径和目标数设定；
- [ ] 路径关联门限与冲突消解规则；
- [ ] LM 初始阻尼、停止门限和最大迭代次数；
- [ ] 三类几何先验协方差 $\mathbf C_{p,k}^{(0)}$、$\mathbf C_{\mathrm{UE}}^{(0)}$、$\sigma_{\Delta t,0}^2$；
- [ ] Taylor 展开的信赖域、阻尼规则和线谱模块停止条件；
- [ ] 非线性几何消息采用粒子、sigma-point 还是局部 Laplace 投影；
- [ ] 两模块仅执行一次还是采用外信息反馈迭代，以及相应的消息停止条件。

### 8.2 建议对比方法

1. 仅 2D-FFT + 几何 LM，不进行贝叶斯精化；
2. Taylor 线谱估计后直接将角度和时延后验均值代入几何方程；
3. Taylor 线谱估计 + 单次非线性几何消息传递；
4. 带几何外信息反馈的两模块迭代消息传递；
5. 固定噪声精度、固定 UE 位置或固定时钟偏差的消融方法；
6. 原位置网格稀疏方法可作为历史基线，但不属于当前概率模型。

### 8.3 评价指标

- 目标位置 RMSE 或多目标 Chamfer distance；
- UE 位置 RMSE；
- 时钟偏差 RMSE；
- 归一化空间相位和时延相位的粗估计误差，以及反归一化后的物理角度和测量时延误差；
- 散射系数 NMSE；
- 噪声精度估计误差；
- 负对数似然、ELBO 和收敛迭代次数；
- 后验可信区间覆盖率；
- 运行时间与峰值内存。

### 8.4 正确性测试

- [ ] 无噪声单径情况下，2D-FFT 峰值映射到正确的归一化空间相位和时延相位，并能反归一化为正确的方向余弦和秒单位测量时延；
- [ ] 使用真实角度和时延时，LM 能恢复目标、UE 与时钟偏差；
- [ ] 数值差分验证 $\boldsymbol\Phi_{\theta,n}$ 和 $\boldsymbol\Phi_{\tau,n}$；
- [ ] 验证 Taylor 近似误差随增量范数呈二阶衰减，并测试信赖域边界；
- [ ] 固定几何时，$q(\boldsymbol h_n)$ 与标准线性复高斯后验一致；
- [ ] Gamma 更新中的期望残差包含均值和协方差两部分；
- [ ] 线谱模块每轮变分更新后局部 ELBO 不降低，或满足规定的接受准则；
- [ ] Dirac 因子的前向消息满足 $m_{f_n\rightarrow z}(\boldsymbol z)=m_{\mathrm{LS},n}(\boldsymbol g_n(\boldsymbol z))$；
- [ ] 关闭外信息反馈时，两阶段结果与直接顺序实现一致；开启反馈时不存在信息自循环；
- [ ] 两个接收站交换顺序不改变对称场景的最终结果；
- [ ] 路径排列改变后，经数据关联得到的几何结果保持一致；
- [ ] 几何秩亏时能够给出诊断，而不是输出虚假的高置信度结果。

## 9. 当前代码接口与后续实现边界

当前仓库已经具备：

- `CoarseEstimation2DFFT.m`：过采样 2D-FFT、二维峰值抑制和亚栅格插值；
- `SpecEstimation2D.m`：兼容旧实验调用的 FFT 包装接口；
- `estimate_bias_by_coherence.m`：双站路径关联和自适应 LM 粗估计。

后续需要新增或重构的核心模块为：

```text
build_parametric_path_matrix
build_path_derivative_matrices
build_linespectral_taylor_model
initialize_bayesian_model
update_scattering_posterior
update_noise_precision_posterior
update_angle_delay_posterior
form_linespectral_extrinsic_message
evaluate_linespectral_elbo
propagate_dirac_geometry_message
project_nonlinear_geometry_message
form_geometry_extrinsic_message
check_message_convergence
summarize_geometry_posterior
```

现有 `JSLES*.m`、`OGVAMP.m` 及位置网格感知矩阵函数实现的是原稀疏网格模型，不能直接视为本概率模型的变分实现。若保留，应在实验中标记为历史基线；当前技术路线的新推断模块应围绕 $K$ 列参数化路径矩阵 $\boldsymbol\Phi_n$、其角度/时延列导数矩阵和局部 Taylor 观测模型重新实现。

## 10. 实现前仍需固定的模型约定

1. $K$ 是已知场景参数，还是需要由 2D-FFT/模型选择联合估计；
2. $\boldsymbol\Sigma_{h,n}$ 是固定先验协方差，还是需要增加超先验并学习；
3. 几何先验协方差如何由 FFT 峰宽、LM Hessian、GPS 误差或经验参数确定；
4. Taylor 线性化允许的最大归一化角度增量、最大归一化时延增量以及阻尼/信赖域参数；
5. 非线性 Dirac 消息采用粒子、sigma-point 还是局部 Laplace 投影，以及是否保留多峰结构；
6. 是否在线谱模块与几何模块之间迭代，以及外信息的计算和阻尼规则；
7. 是否对目标位置后验采用目标间独立的块对角协方差，还是保留数据关联造成的相关性；
8. 时钟偏差使用秒 $\Delta t$ 还是距离偏差 $b=c\Delta t$ 作为内部变量；
9. 当 2D-FFT 漏检、虚警或两条路径不可分辨时，路径数和标签如何处理。

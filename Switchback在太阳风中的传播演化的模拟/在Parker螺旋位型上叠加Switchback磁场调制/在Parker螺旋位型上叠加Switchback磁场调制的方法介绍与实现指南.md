
# 在Parker螺旋位型上叠加Switchback磁场调制的方法介绍与实现指南

从 Parker 螺旋出发，构造单个/簇状 switchbacks（多种矢势方案），做并行补偿保证局点保模，再在球壳体积上做散度清除（Helmholtz 投影／球谐-径向 Poisson 求解）。每一步给出公式、实现要点、数值注意事项与可直接用的伪码/参数建议，便于立刻编码或移植到你的 MHD 流程中。

## 目录（快速导航）

1. 目标与总体策略
2. 基础：Parker 螺旋基场 $\mathbf B_0$（球坐标）
3. 构造 switchback 的思路与两条主路线
4. 矢势（$A_\theta$）构造方法（多种可选）
5. 计算横向扰动 $\delta\mathbf B_\perp=\nabla\times\mathbf A$ 的显式式子
6. 并行补偿：点点保证 $|\mathbf B|=|\mathbf B_0|$（根号解与分支）
7. 若需严格无散度：Helmholtz 投影（球壳域的谱-径向 Poisson 解法）
8. 完整算法流程与伪码（逐步）
9. 验证与调参建议（诊断量、常见问题）
10. 附：常用参数表与 Python 实现要点

---

## 1. 目标与总体策略（一句话）

在球壳体积上以 Parker 螺旋为背景，通过构造局域矢势 $A_\theta$ 生成无散度的横向扭转，再做并行分量修正以精确保模，最后用 Helmholtz 投影（求解 Poisson）消除残余散度，从而得到既**近乎保模**又\*\*（经投影）严格无散度\*\*的 switchback 磁场用于合成/数值初始化。

---

## 2. 基础：Parker 螺旋基场（球坐标）

球坐标 $(r,\theta,\phi)$，太阳自转角速度 $\Omega$，太阳风速度 $V_{\rm sw}$（近常数）：

$$
\boxed{%
B_r^{(0)}(r)=B_\ast\Big(\frac{r_\ast}{r}\Big)^2,\qquad
B_\phi^{(0)}(r,\theta)=-B_r^{(0)}\,\frac{\Omega r\sin\theta}{V_{\rm sw}},\qquad
B_\theta^{(0)}=0.}
$$

写成向量：

$$
\mathbf B_0 = B_r^{(0)}\;\hat e_r + B_\phi^{(0)}\;\hat e_\phi,\qquad
|\mathbf B_0|=\sqrt{(B_r^{(0)})^2+(B_\phi^{(0)})^2}.
$$

---

## 3. 构造 switchback 的两条主路线（选择其一或组合）

A. **直接旋转法（SO(3) 旋转）** — 最简单、严格保模：
对每点将 $\mathbf B_0$ 绕局部轴 $\hat n(\mathbf x)$ 旋转角 $\psi(\mathbf x,t)$：

$$
\mathbf B(\mathbf x)=\mathcal R(\hat n,\psi)\,\mathbf B_0
$$

优点：模长精确不变；易用于合成时序。缺点：局域旋转随空间变化会产生散度（需清散）并且不显式给出矢势（难直接保证无散度）。

B. **矢势叠加法（推荐用于与 MHD 初值配合）** — 先构造严格无散度的横向扰动 $\delta\mathbf B_\perp=\nabla\times(A_\theta\hat e_\theta)$，再做并行补偿保证模长，最后投影清散。优点：可控、易数值实现、天然无散度横向部分；适合多簇叠加与球壳投影。缺点：并行补偿与投影步骤必须谨慎处理以保证稳定。

本流程采用 B（矢势叠加 + 并行补偿 + 投影）。

---

## 4. 矢势 $A_\theta$ 的几种实用构造（可叠加多簇）

通式（球坐标）：

$$
\mathbf A = A_\theta(r,\theta,\phi,t)\,\hat e_\theta.
$$

常用组件（单簇）：

$$
A_\theta(r,\theta,\phi,t)=a_0(r,\theta)\;R\!\Big(\frac{r-r_c-V_{\rm sb}t}{L_r}\Big)\;W\!\Big(\frac{\theta-\theta_c}{L_\theta}\Big)\;W_\phi\!\Big(\frac{\phi-\phi_c}{L_\phi}\Big).
$$

* $R,W,W_\phi$：平滑窗（高斯 $e^{-x^2/2}$、tanh、或 compact bump）。
* 幅度定标建议（使横向扰动与角度 $\psi$ 一致）：

  $$
  \boxed{\,a_0(r,\theta)=|\mathbf B_0(r,\theta)|\,L_r\,\sin\psi(r,\theta)\,}
  $$

  这里的 $\psi$ 是目标局地旋转角，必须在 $[-π,π]$ 内来保证 $|\delta B_\perp|\le|B_0|$.
* 多簇：把 $A_\theta=\sum_i A_{\theta,i}$（每簇独立参数）叠加。**注意**：若要严格对应旋转合成，优选把每簇视作单独旋转再合成，但矢势叠加在工程上已足够。

**常见配方**（黄道面 $\theta=\pi/2$ 简化）：

* Gaussian blob: $R(\xi)=\exp(-\xi^2/2)$, $W_\phi(\eta)=\exp(-\eta^2/2)$.
* Erf-based integral shape (用于反号保证)：见之前“逆向构造”处方，用积分函数 $\mathcal S(\Delta\phi)=L_\phi\,\mathrm{erf}(\Delta\phi/(\sqrt2 L_\phi))$ 生成直接控制 $\delta B_r$ 的 A。

---

## 5. 由 $A_\theta$ 得到 $\delta\mathbf B_\perp=\nabla\times\mathbf A$（严格无散度）

球坐标公式（完全形式）：

$$
\boxed{\begin{aligned}
\delta B_r &= \frac{1}{r\sin\theta}\frac{\partial}{\partial\phi}\big(A_\theta\sin\theta\big)
\approx \frac{1}{r}\partial_\phi A_\theta \quad(\text{黄道面 }\sin\theta\approx1),\\[4pt]
\delta B_\theta &= 0,\\[4pt]
\delta B_\phi &= \frac{1}{r}\frac{\partial}{\partial r}\big(rA_\theta\big).
\end{aligned}}
$$

将 $\delta B$ 从球坐标投影到笛卡尔以便并行补偿与投影步骤。

---

## 6. 并行补偿（点点精确保模）：根号公式与分支选择

目标：得到在每点满足 $|\mathbf B|=|\mathbf B_0|$ 的场（先不考虑散度）：
设 $\delta\mathbf B_\perp$ 已知（笛卡尔），令 $b_0=\hat{\mathbf b}_0=\mathbf B_0/|B_0|$，$s^2=|\delta\mathbf B_\perp|^2$，则并行分量 $\delta B_\parallel$ 应满足

$$
2|B_0|\delta B_\parallel + \delta B_\parallel^2 + s^2 =0.
$$

解得

$$
\boxed{\,\delta B_\parallel = -|B_0| \pm \sqrt{|B_0|^2 - s^2}\, }.
$$

**分支选择（关键）**：若已知目标角 $\psi$（由构造），取 sign = $\operatorname{sgn}(\cos\psi)$ 从而恢复
$\delta B_\parallel=|B_0|(\cos\psi-1)$. 若不知道 $\psi$，可用连续性启发：选 $\sigma=\operatorname{sgn}(\mathbf B_0\cdot(\mathbf B_0+\delta\mathbf B_\perp))$。
**约束**：必须保证 $s^2 \le |B_0|^2$。若 $s^2 > |B_0|^2$（说明 A 设计过强），需缩小 $a_0$ 或把 $s^2$ 截断为 $ |B_0|^2-\varepsilon$。

并行分量向量：

$$
\delta\mathbf B_\parallel = \delta B_\parallel \; b_0.
$$

构造未投影的保模场：

$$
\mathbf B^{(c)} = \mathbf B_0 + \delta\mathbf B_\perp + \delta\mathbf B_\parallel,
\qquad |\mathbf B^{(c)}|=|B_0| \ (\text{点点}).
$$

---

## 7. Helmholtz 投影：在球壳上做 $\Phi$ 的球谐-径向 Poisson 求解（将 $\nabla\cdot\mathbf B^{(c)}$ 清除）

目标 PDE：

$$
\nabla^2\Phi(r,\theta,\phi)=\rho(r,\theta,\phi),\qquad \rho=\nabla\cdot\mathbf B^{(c)}.
$$

投影后 $\mathbf B=\mathbf B^{(c)}-\nabla\Phi$ 且 $\nabla\cdot\mathbf B=0$.

**高效方法（推荐）**：对角向做球谐展开，对径向做有限差分或谱法（混合谱-径向法）：

1. 球谐展开：

$$
\Phi(r,\theta,\phi)=\sum_{\ell=0}^{L_{\max}}\sum_{m=-\ell}^\ell \phi_{\ell m}(r) Y_{\ell m}(\theta,\phi),
\qquad
\rho(r,\theta,\phi)=\sum_{\ell m}\rho_{\ell m}(r) Y_{\ell m}(\theta,\phi).
$$

2. 对每对 $(\ell,m)$ 得到径向 ODE：

$$
\frac{1}{r^2}\frac{d}{dr}\Big(r^2\frac{d\phi_{\ell m}}{dr}\Big)-\frac{\ell(\ell+1)}{r^2}\phi_{\ell m}(r)=\rho_{\ell m}(r).
$$

3. 离散化（径向网格 $r_i$）：用二阶中心差分得到三对角线系统（Thomas 算法可解）。边界条件常选：

   * Neumann–Neumann（$\partial_r\Phi=0$ 在内外边界）——维持边界法向场不变，但需体积兼容（$\int_V \rho\,dV=0$），若不满足可均匀修正 $\rho$；或
   * Dirichlet–Neumann 混合（更稳定，无兼容性约束）。
4. 特别处理 $\ell=0$ 模：需保证总体源兼容或用 Dirichlet 条件。
5. 反变换：求得 $\phi_{\ell m}(r_i)$ 后在每径向层做逆球谐变换得到 $\Phi(r_i,\theta_j,\phi_k)$。
6. 计算 $\nabla\Phi$：径向导数由差分（或直接由 $\phi_{\ell m}'$ 得到）；角向导数可在系数空间用球谐导数表达（更精确）或在物理格点差分。最后 $\mathbf B=\mathbf B^{(c)}-\nabla\Phi$。

**数值细节/兼容性**：

* 计算初期检查 $I=\int_V \rho\,dV$。若 $|I|>$ 容忍度，做微小修正 $\rho\leftarrow\rho-I/V$（均匀减偏），或改用 Dirichlet 在一端。
* 使用高效球谐库 (`shtns`、`pyshtools`) 提速；无库时可用 `scipy.special.sph_harm` 做参考实现（角向慢）。
* 投影后 $|\mathbf B|$ 会微小偏离 $|\mathbf B_0|$（幅度与 $\nabla\Phi$ 尺度有关）；若偏差太大，可迭代：用投影后 $\mathbf B$ 作为新的背景，重新做并行补偿→投影，迭代 1–3 次通常足够。

---

## 8. 完整算法流程（逐步伪码）

下面给出端到端流程伪码，适用于球壳体积 $r\in[r_{\rm in},r_{\rm out}]$ 的离散网格（$N_r\times N_\theta\times N_\phi$）：

**输入**：网格、Parker 参数、switchback 簇参数（若干簇），$L_{\max}$、径向网格、投影 BC 等。
**输出**：投影后的 $\mathbf B$（无散度、近乎保模）与诊断量。

```text
# 0. 预计算与设置
compute B0(r,theta,phi) on grid  # Parker spiral
compute |B0| and unit vectors b0hat

# 1. 构造 A_theta(r,theta,phi): sum of clusters
for each grid point:
    A_theta = sum_i a0_i(r,theta)*R((r-rci-Vi*t)/Lr_i)*Wtheta((theta-thci)/Lth_i)*Wphi((phi-phci)/Lph_i)

# 2. 计算 deltaB_perp (球坐标) via:
#    deltaBr = (1/(r sin theta)) d_phi (A_theta sin theta)
#    deltaBphi = (1/r) d_r (r A_theta)
compute derivatives (spectral or FD). Convert deltaB_perp to Cartesian.

# 3. 点点并行补偿 to enforce |B|=|B0|
for each grid point:
    s2 = dot(deltaB_perp_cart, deltaB_perp_cart)
    if s2 > |B0|^2 - eps: s2 = |B0|^2 - eps   # regularize
    choose sigma = sign(cos(psi)) if psi known else heuristic
    deltaB_par = -|B0| + sigma*sqrt(|B0|^2 - s2)
    deltaB_par_cart = deltaB_par * b0hat
    Bc = B0_cart + deltaB_perp_cart + deltaB_par_cart
store Bc

# 4. 计算散度 rho = div(Bc) on grid (same discrete operator as later)
compute total integral I = sum(rho*dV)
if abs(I) > tol:
    rho -= I/V_total   # enforce compatibility (simple fix)

# 5. 球谐投影 (for each radial shell r_i): rho_lm(r_i) = SHT[rho(r_i,.,.)]
#    choose Lmax <= N_theta-1 and <= N_phi/2
for l in 0..Lmax:
  for m in -l..l:
    assemble tridiagonal system for phi_lm(r) from radial ODE:
       d/dr(r^2 d phi/dr) - l(l+1) phi = r^2 rho_lm(r)
    apply BC and solve (Thomas)

# 6. 逆变换：for each r_i reconstruct Phi(r_i,theta,phi) from phi_lm(r_i)
# 7. compute grad Phi (r,theta,phi): radial derivative by FD; angular by coeffs or FD
# 8. final B = Bc - grad Phi
# 9. diagnostics:
compute divB_final (should be ~0)
compute pointwise |B|-|B0| (measure deviation)
optionally iterate: recompute b0hat from B, redo step 3-8 for refinement
```

---

## 9. 验证与调参建议（要检查的指标）

* $\max|\nabla\cdot\mathbf B|$（投影后应接近求解器容差，例如 $10^{-12}$ 至 $10^{-8}$ 的相对值，取决于浮点/求解器）。
* 点对点与统计的 $|\mathbf B|-|\mathbf B_0|$（投影前为 0，投影后应是小量；若太大，说明投影引入强梯度）。
* 在簇中心检查 $B_r$ 是否按设计翻转（若设计目标是翻转）——用几何判断 $\alpha+\psi$ 或直接读取数值。
* 能量/通量守恒检查（若把此场作为 MHD 初值，监控边界法向与全域磁通）。
* 若并行补偿频繁触发截断 $s^2>|B0|^2$，减小 $a_0$ 或调整 $\psi$ 分布。

---

## 10. 附：常用参数与实现要点

* Parker：$\Omega\approx2.7\times10^{-6}\,\mathrm{rad/s}$，$V_{\rm sw}=400$–$700$ km/s。
* 簇尺度：$L_r$ 取 $10^4\!-\!10^6$ km（近太阳小，远日际大）；角向 $L_\theta,L_\phi$ 取 $0.01\!-\!0.5$ rad。
* 角幅：$\psi_0$ 取 $0.5\pi$（强横向）或接近 $\pi$（欲实现翻转并用并行补偿）。
* 投影：若能使用 `shtns` 或 `pyshtools` 性能最好；径向方程用 Thomas 三对角解法。
* 网格：角向采样 $N_\theta$（Gauss-Legendre）与 $N_\phi$（uniform）须能解析到 $L_{\max}$（$L_{\max}\lesssim N_\theta-1$）。
* 数值稳定：对 $s^2$ 做安全截断 $\min(s^2, |B0|^2-\varepsilon)$；为避免 NaN，设 $\varepsilon\sim 10^{-12}|B0|^2$。

---

## 附：关键代码片段（伪 Python）

（仅示意核心步骤；完整脚本可基于此扩展）

```python
# 1. compute A_theta on spherical grid (sum of Gaussians)
A = compute_A_theta(...)

# 2. derivatives -> deltaBr, deltaBphi
deltaBr = (1.0/(r*sinth)) * dphi( A*sinth )
d_rA = dr( r*A )
deltaBphi = (1.0/r) * d_rA

# 3. convert to cart: deltaB_perp_cart = deltaBr*er + deltaBphi*ephi
deltaB_perp_cart = to_cartesian(deltaBr, 0.0, deltaBphi)

# 4. parallel correction
s2 = dot(deltaB_perp_cart, deltaB_perp_cart)
s2 = min(s2, B0norm**2 - eps)
sigma = np.sign(np.cos(psi))  # if psi known, else heuristic
deltaB_par = -B0norm + sigma*np.sqrt(B0norm**2 - s2)
deltaB_par_cart = deltaB_par[...,None] * b0hat
Bc = B0cart + deltaB_perp_cart + deltaB_par_cart

# 5. divergence
rho = divergence(Bc)

# 6. SHT per r: rho_lm = sht_forward(rho(r_i,:,:))
# 7. for each (l,m) solve radial tridiag -> phi_lm(r)
# 8. inverse SHT to get Phi(r_i,:,:)
# 9. compute grad Phi and B = Bc - gradPhi
```


<!---
---

如果你愿意，我可以接下来把这套流程**写成可以运行的参考实现**，两种实现供选：

* A：**参考实现（纯 Python + SciPy）** — 可直接运行但角向变换较慢（使用 `scipy.special.sph_harm`）。
* B：**高性能实现（Python + shtns / pyshtools）** — 高效且适合较大网格/实际工程。

或者我可以把上面流程直接转换为**C/Fortran 风格伪码**以方便移植进你现有的 MHD 代码库（BATS-R-US/WAM-I等）。你想要哪种下一步？



---
# 在Parker Spiral磁场中引入磁场kink/switchback，构建扰动磁场的磁矢势 + 并行补偿 + 投影校正，写一下完整的表达式

把**完整的数学处方**与\*\*数值流程（伪码）\*\*写清楚，覆盖你要的三步：

1. 用矢势 $A_\theta$ 构造严格无散度的横向扰动 $\delta\mathbf B_\perp$；
2. 用并行分量 $\delta B_\parallel$ 精确保模（点对点）得到 $\mathbf B^{(c)}$；
3. 求解 Poisson 投影 $\nabla^2\Phi=\nabla\!\cdot\mathbf B^{(c)}$ 得到 $\mathbf B=\mathbf B^{(c)}-\nabla\Phi$（严格无散度，且“近乎保模”）。

我同时给出黄道面（2D）简化公式、数值注意事项与伪代码，方便直接实现。

---

## 一、符号与前置定义（连续形式）

使用球坐标 $(r,\theta,\phi)$ 与直角笛卡尔 $(x,y,z)$。
基场（Parker）：

$$
\mathbf B_0(r,\theta)=B_r^{(0)}(r)\,\hat e_r + B_\phi^{(0)}(r,\theta)\,\hat e_\phi,\qquad
|\mathbf B_0|=\sqrt{(B_r^{(0)})^2+(B_\phi^{(0)})^2}.
$$

矢势只保留 $\theta$ 分量：

$$
\boxed{\ \mathbf A = A_\theta(r,\theta,\phi,t)\,\hat e_\theta\ }.
$$

由旋度公式（球坐标）得到横向扰动（严格无散度）：

$$
\boxed{\begin{aligned}
\delta B_r &= \frac{1}{r\sin\theta}\,\partial_\phi\big(A_\theta\sin\theta\big)
\;\approx\;\frac{1}{r}\,\partial_\phi A_\theta\quad(\text{在黄道面 }\sin\theta\approx1),\\[4pt]
\delta B_\theta &= 0,\\[4pt]
\delta B_\phi &= \frac{1}{r}\,\partial_r\big(rA_\theta\big).
\end{aligned}}
$$

构造 $A_\theta$ 的示例（单簇中心 $r_c,\phi_c$）：

$$
A_\theta(r,\theta,\phi,t)=a_0(r,\theta)\,R\!\Big(\frac{r-r_c-V_{\rm sb}t}{L_r}\Big)\,
W\!\Big(\frac{\phi-\phi_c}{L_\phi}\Big)\,g(\theta),
$$

其中 $R,W$ 为平滑窗（如高斯 $e^{-x^2/2}$），且幅度定标取

$$
a_0(r,\theta)=|\mathbf B_0(r,\theta)|\,L_r\,\sin\psi(r,\theta)
$$

以保证横向扰动量与角度 $\psi$ 一致（见后）。

---

## 二、并行补偿（点点精确保模）

在每一点（$r,\theta,\phi$ 或在笛卡尔网格上的点）计算：

1. 先求出横向扰动（球坐标分量） $\delta\mathbf B_\perp=(\delta B_r,0,\delta B_\phi)$ 并转换到笛卡尔 $\delta\mathbf B_\perp^{\text{cart}}$。

2. 求背景与单位向量

$$
\mathbf B_0^{\text{cart}},\qquad b_0=\frac{\mathbf B_0^{\text{cart}}}{|\mathbf B_0|}.
$$

3. 求横向模长（笛卡尔）：

$$
s^2 \equiv |\delta\mathbf B_\perp^{\text{cart}}|^2.
$$

**判别与正则化**（避免数值误差）：

$$
\varepsilon\ll1\ (\text{例如 }10^{-12}|B_0|^2).
$$

若 $s^2>|\mathbf B_0|^2-\varepsilon$，把 $s^2 \leftarrow \min(s^2,\,|\mathbf B_0|^2-\varepsilon)$。

4. 选分支并给并行分量（恢复精确旋转关系）：
   令 $\cos\psi$ 的符号 $ \sigma = \operatorname{sgn}(\cos\psi)$ —— 若你有目标角度 $\psi$，则 $\sigma=\operatorname{sgn}(\cos\psi)$。若没有，可用
   $\sigma=\operatorname{sgn}\big(\mathbf B_0\cdot(\mathbf B_0+\delta\mathbf B_\perp^{\text{cart}})\big)$ 作为经验判据（保连续性）。
   然后

$$
\boxed{\ \delta B_\parallel = -|\mathbf B_0| + \sigma\sqrt{\,|\mathbf B_0|^2 - s^2\,}\ }.
$$

并行分量向量（笛卡尔）：

$$
\boxed{\ \delta\mathbf B_\parallel^{\text{cart}} = \delta B_\parallel\, b_0. \ }
$$

5. 构造修正后的场（未投影）：

$$
\boxed{\ \mathbf B^{(c)} = \mathbf B_0^{\text{cart}} + \delta\mathbf B_\perp^{\text{cart}} + \delta\mathbf B_\parallel^{\text{cart}}. \ }
$$

按构造，逐点满足 $|\mathbf B^{(c)}| = |\mathbf B_0|$（模长严格守恒，若上一步有正则化则近似成立）。

---

## 三、Helmholtz 投影：恢复严格无散度

计算散度（在笛卡尔网格上）：

$$
\rho(\mathbf x) \equiv \nabla\!\cdot\mathbf B^{(c)}(\mathbf x).
$$

解 Poisson 方程（域 $D$）：

$$
\boxed{\ \nabla^2\Phi(\mathbf x) = \rho(\mathbf x),\qquad \text{(边界条件见下)}\ }.
$$

然后

$$
\boxed{\ \mathbf B(\mathbf x)=\mathbf B^{(c)}(\mathbf x)-\nabla\Phi(\mathbf x),\quad \nabla\!\cdot\mathbf B=0.}
$$

**边界条件建议**（实际实现关键）：

* 周期域：用 FFT 解 Poisson（最简单），取周期 BC。
* 非周期/外推：用 Neumann BC $\partial_n\Phi = 0$（保持外法向通量不变），或用 Dirichlet $\Phi=0$ 若你愿意固定边界值。
* 若你要把场嵌入 MHD 模式，通常把投影域选择为与 MHD 网格一致并采用与 MHD 相同的边界策略（外推/指定）。

数值求解方法：FFT（周期），多重网格/共轭梯度（一般边界），迭代求解器（高效且稳健）。

---

## 四、黄道面（2D）专用简化公式与实现细节

在黄道面 $\theta=\pi/2$ 上，球坐标公式简化为：

$$
\delta B_r = \frac{1}{r}\partial_\phi A_\theta(r,\phi),\qquad
\delta B_\phi = \frac{1}{r}\partial_r(rA_\theta).
$$

实现流程：

1. 在二维极坐标网格或笛卡尔网格上评估 $A_\theta(r,\phi)$（可直接用例子窗函数与幅度定标）。
2. 用有限差分/谱方法计算 $\partial_\phi A_\theta,\ \partial_r(rA_\theta)$。
3. 把 $\mathbf B_0$ 与 $\delta\mathbf B_\perp$ 转到笛卡尔 $(B_x,B_y)$，做并行补偿（与上文完全相同，但在 2D 中 $b_0$ 也在平面内）。
4. 求 $\rho=\partial_x B_x^{(c)}+\partial_y B_y^{(c)}$，解二维 Poisson $\nabla^2\Phi=\rho$，然后 $(B_x,B_y)=(B_x^{(c)}-\partial_x\Phi,\;B_y^{(c)}-\partial_y\Phi)$。

---

## 五、伪代码（实现顺序）

```text
# 先在网格上准备 B0 (cartesian) 与 A_theta (在球/极坐标上)
for each grid point:
    compute B0_cart, B0_norm, b0hat
    compute A_theta and its derivatives: dA_dphi, d_rA_dr
    compute deltaB_perp_sph = (deltaBr = (1/r)*dA_dphi, 0, deltaBphi = (1/r)*d_rA_dr)
    convert deltaB_perp_sph -> deltaB_perp_cart

    s2 = dot(deltaB_perp_cart, deltaB_perp_cart)
    s2 = min(s2, B0_norm**2 - eps)   # regularize to avoid negative inside sqrt

    # determine sigma (sign) 
    if have_target_angle psi:
        sigma = sign(cos(psi))
    else:
        sigma = sign( dot(B0_cart, B0_cart + deltaB_perp_cart) )  # heuristic

    deltaB_par = -B0_norm + sigma * sqrt(B0_norm**2 - s2)
    deltaB_par_cart = deltaB_par * b0hat

    Bc = B0_cart + deltaB_perp_cart + deltaB_par_cart
    store Bc on grid

# Now compute divergence of Bc (discrete)
rho = div(Bc)

# Solve Poisson: Laplacian(Phi) = rho  with chosen BC
Phi = solve_poisson(rho, BC)

# Final divergence-free field
B = Bc - grad(Phi)
```

---

## 六、数值/物理注意事项与调参建议

1. **保证可解性**：必须在并行补偿步保证 $s^2 \le B_0^2$。如果你的 $A_\theta$ 设计会产生更大的 $s$，那表明你在横向上构造的扰动超过了同点能通过纯旋转实现的范围 —— 需要缩小 $a_0$ 或改为直接用旋转方法。
2. **正则化 $\varepsilon$**：避免浮点负数问题。
3. **分支选择 $\sigma$**：若你已经用角度函数 $\psi$ 设计 A 的幅度与位相，直接取 $\sigma=\operatorname{sgn}(\cos\psi)$ 最稳；如果没有，经验选择（见伪码）通常能保持连续性。
4. **投影对模长的影响**：Poisson 投影会修改模长：$|\mathbf B|$ 不再严格等于 $|\mathbf B_0|$；但若 $\nabla\Phi$ 很小（即 $\nabla\!\cdot\mathbf B^{(c)}$ 本身已小），偏离很小。若需要更高精度，可迭代：再做一次并行补偿（用新 B 估计），再投影，如此交替数次，直到模长与无散度都满意。
5. **边界对结果的影响**：投影解的边界条件会影响 $\Phi$ 以及最终场；选取与 MHD 模拟相容的 BC（或大域盒、周期）能最小化人工效应。
6. **效率**：对大网格优先用多重网格/FFT 求解 Poisson；二维情况更快，三维需更多算力。

---

## 七、示例：黄道面上一个单簇的解析处方（总结）

取黄道面 $\theta=\pi/2$。指定：

* $A_\theta(r,\phi)=a_0(r)\,e^{-(r-r_c)^2/(2L_r^2)}\,\mathrm{erf}\!\big((\phi-\phi_c)/(\sqrt2 L_\phi)\big)$（或高斯窗），
* 计算 $\delta B_r=(1/r)\partial_\phi A_\theta$、$\delta B_\phi=(1/r)\partial_r(rA_\theta)$，
* 转换并做并行补偿与投影，如上。

---

需要的话我可以：

* 把上面的伪码翻成**可直接运行的 Python（NumPy）脚本**：包含 A\_theta 定义、差分求导、并行补偿、用 FFT 解 Poisson（周期 BC）并展示投影前后模长与散度的对比图；或
* 给出**三维实现模板**（含 multigrid Poisson 接口伪码），适合插入你的 MHD 预处理管线。

你想先看 2D（黄道面）可运行脚本，还是直接要 3D 版本？

--->

<!---
---
# 在球壳上做 Helmholtz 投影（Poisson 解）校正清除磁场散度的数值方案

下面把**在球壳体积上**对 $\mathbf B^{(c)}$ 做 Helmholtz 投影（解 Poisson 方程 $\nabla^2\Phi=\nabla\!\cdot\mathbf B^{(c)}$ 并取 $\mathbf B=\mathbf B^{(c)}-\nabla\Phi$）的**完整、可实现**数值方案写清楚。重点采用对球面角用**球谐分解（spherical harmonic）**、径向用**有限差分（或谱）** 的混合谱-径向方法——这种方法在球壳上既准确又高效，且易于并行化。下面给出连续方程、每一步的离散化、边界条件选择、l=0 的特殊处理、实现要点与伪码/实现提示。

---

## 1. 问题与目标（回顾）

要解的方程（在球壳域 $r\in[r_{\rm in},r_{\rm out}],\ \theta\in[0,\pi],\ \phi\in[0,2\pi)$）：

$$
\nabla^2\Phi(r,\theta,\phi)=\rho(r,\theta,\phi),\qquad
\rho\equiv\nabla\!\cdot\mathbf B^{(c)}.
$$

投影后场

$$
\mathbf B=\mathbf B^{(c)}-\nabla\Phi,\quad\nabla\!\cdot\mathbf B=0.
$$

我们选择对 $(\theta,\phi)$ 做球谐变换，写

$$
\Phi(r,\theta,\phi)=\sum_{\ell=0}^{L_{\max}}\sum_{m=-\ell}^{\ell}\phi_{\ell m}(r)\,Y_{\ell m}(\theta,\phi),
\qquad
\rho(r,\theta,\phi)=\sum_{\ell m}\rho_{\ell m}(r)\,Y_{\ell m}(\theta,\phi).
$$

利用球谐的正交性，原 PDE 对每对 $(\ell,m)$ 退化为一个一维径向常系数二阶 ODE：

$$
\boxed{\ 
\frac{1}{r^2}\frac{d}{dr}(r^2\frac{d\phi_{\ell m}}{dr})
-\frac{\ell(\ell+1)}{r^2}\,\phi_{\ell m}(r)
=\rho_{\ell m}(r).\ }
$$

等价写法：

$$
\frac{d}{dr}(r^2\phi_{\ell m}')-\ell(\ell+1)\phi_{\ell m}=r^2\rho_{\ell m}(r).
$$

所以：问题转为对每个 $(\ell,m)$ 在径向解一个二阶常微分方程。

---

## 2. 球谐系数 $\rho_{\ell m}(r)$ 的计算（离散网格）

* 在给定三维网格上（通常是 $N_r\times N_\theta\times N_\phi$，$\theta$ 用 Gauss–Legendre 节点、$\phi$ 用均匀网格最合适），在每个径向层 $r_i$：

  $$
  \rho_{\ell m}(r_i)=\int_0^{2\pi}\!\int_0^\pi \rho(r_i,\theta,\phi)\,Y_{\ell m}^\ast(\theta,\phi)\,\sin\theta\,d\theta d\phi.
  $$

* 离散化用 Gauss-Legendre + trapezoidal in $\phi$（或现成的 SHT 库，如 `shtools`, `pyshtools`，或基于 FFT + associated Legendre 机制），计算得到 $\rho_{\ell m}(r_i)$ 对每个 $i$ 和 $(\ell,m)$。

实现要点：

* 选 $L_{\max}$ 不要超过角向网格可解析的最大模式（奈奎斯特）：若 $\theta$ 使用 $N_\theta$ 高斯点，通常 $L_{\max}\lesssim N_\theta-1$；$\phi$ 网格 $N_\phi$ 决定能表示的 $m$ 最大值。
* 用库（例如 `pyshtools`、`shtns`）可大幅提升速度；若无库也可用 `scipy.special.sph_harm` 做直接投影（慢）。

---

## 3. 径向 ODE 的离散化与边界条件

### 3.1 连续 ODE：

$$
\frac{d}{dr}\!\big(r^2\phi'\big)-\ell(\ell+1)\phi = r^2\rho_{\ell m}(r).
$$

### 3.2 边界条件的选择（实用建议）

常用两类 BC：

* **Neumann–Neumann（保边界法）**：在内外边界保持法向场不变，即

  $$
  \frac{\partial\Phi}{\partial r}\bigg|_{r=r_{\rm{in}}}=0,
  $$

  $$
  \frac{\partial\Phi}{\partial r}\bigg|_{r=r_{\rm{out}}}=0,
  $$

  因为 $B_n = B^{(c)}_n - \partial_n\Phi$，若要保持 $B_n$ 在边界不变，需要 $\partial_n\Phi=0$.
  但 Neumann–Neumann 对应的 Poisson 问题需要兼容性条件：总体源积分为零（见下）。
* **Dirichlet–Neumann（常用且稳定）**：例如在内边界设 $\Phi(r_{\rm in})=0$（固定势，等价于允许边界法向 $B$ 变化），外边界 Neumann（$\partial_r\Phi=0$）或相反。Dirichlet–Neumann 消除了 Neumann 兼容性问题，但会改变边界法向分布（需要评估是否可接受）。

**建议**（若希望尽量不改动边界法向场）：优先尝试 Neumann–Neumann **且**先检查兼容性（下文）。若兼容性差异很小（$<10^{-12}$），可以继续；若不满足，则采用 Dirichlet 在其中一侧或通过微小修正使兼容性成立（见 3.4）。

### 3.3 分析兼容性（Gauss 定理）

积分两边：

$$
\int_V \rho\,dV = \int_V \nabla\!\cdot\mathbf B^{(c)}\,dV = \oint_{\partial V}\mathbf B^{(c)}\!\cdot\!d\mathbf S .
$$

若选 Neumann–Neumann（$\partial_n\Phi=0$ 两边），则左边必须为 0 才有解。实际中若 $\mathbf B^{(c)}$ 是通过 $\nabla\times A$ + 并行补偿得到，理论上应满足
$\int_V \nabla\cdot\mathbf B^{(c)} dV = \oint B^{(c)}_n dS - \oint B^0_n dS$ — 但数值截断、映射误差会产生小非零残差。

**处理办法**：

* 直接计算 $I=\int_V \rho\,dV$. 若 $|I|$ 小于容忍度（比如机器精度乘以体积或目标误差），直接继续 Neumann。
* 若 $I$ 不小，给 $\rho$ 做微小调整：定义 $\tilde\rho=\rho - I/ V$（均匀减后使体积分为 0）。这是最常见的数值修正（会把残差自由度平均分配到域内，改变极微小）。
* 更理想的方法是修正边界法向：在边界上把 $B^{(c)}_n$ 做微小修正使通量匹配原始通量——这在工程上更复杂，一般不用，除非你必须严格保持边界 $B_n$。

### 3.4 径向差分离散（二阶中心差分，均匀网格示例）

用径向网格 $r_i,\ i=0\ldots N_r-1$，间距 $\Delta r$（可用非均匀，但推导更长）。在点 $r_i$ 用二阶差分离散方程：

定义 $q_i\equiv r_i^2$. 离散近似：

$$
\frac{1}{\Delta r}\Big( q_{i+\tfrac12}\frac{\phi_{i+1}-\phi_i}{\Delta r} - q_{i-\tfrac12}\frac{\phi_i-\phi_{i-1}}{\Delta r}\Big)- \ell(\ell+1)\phi_i = r_i^2\rho_{\ell m}(r_i),
$$

其中 $q_{i\pm1/2}=\tfrac12(q_{i}\!+\!q_{i\pm1})$（或更精确用 $r_{i\pm1/2}^2$）。化为三对角线形式：

$$
a_i\phi_{i-1} + b_i\phi_i + c_i\phi_{i+1} = f_i,
$$

系数（均匀 $\Delta r$ 的简洁形式）：

$$
\begin{aligned}
a_i &= \frac{q_{i-\tfrac12}}{\Delta r^2},\\
c_i &= \frac{q_{i+\tfrac12}}{\Delta r^2},\\
b_i &= -\big(a_i + c_i + \ell(\ell+1)\big),\\
f_i &= r_i^2\,\rho_{\ell m}(r_i).
\end{aligned}
$$

（注意单位/因子一致性：上面 $b_i$ 的 - 看你的左/右移项符号习惯调整；确保实现时检验维数。）

**边界差分**（Neumann）：

* Neumann 在 $i=0$（内边界）写 $\phi'(r_0)=0$ 近似 $(\phi_1-\phi_{-1})/(2\Delta r)=0$；用虚点消除得到 $\phi_{-1}=\phi_1$。代入得到 $i=0$ 的离散式（变系数）。
* 类似在 $i=N_r-1$（外边界）做处理 $\phi_{N_r}=\phi_{N_r-2}$。

得到对每 $(\ell,m)$ 的一个 $N_r\times N_r$ 三对角线系统，解法用 Thomas 算法（tridiag）或可向量化 BLAS/LAPACK。

---

## 4. l=0 模式（球谐平均）的特殊处理

$\ell=0,m=0$ 情况：
ODE 变为

$$
\frac{1}{r^2}\frac{d}{dr}\!\big(r^2\phi_{00}'\big)=\rho_{00}(r).
$$

积分得到

$$
r^2\phi_{00}'(r)=\int_{r_{\rm in}}^r r'^2\rho_{00}(r')\,dr' + C.
$$

若选 Neumann BC 两端为 0，则需 $\int_{r_{\rm in}}^{r_{\rm out}} r^2\rho_{00}(r)\,dr = 0$ 才能兼容。实现上：

* 若已用上面“均匀减法”保证总体体积分为 0，则该方程可解。
* 否则可以对 $\ell=0$ 用 Dirichlet（例如设 $\phi_{00}(r_{\rm in})=0$）再向外积分，这会允许调整边界通量。

总结：先检验 $\int_V \rho\,dV$，对不为 0 的做微调，是最实用的做法。

---

## 5. 反变换并构造 $\nabla\Phi$

* 解出每个 $\phi_{\ell m}(r_i)$ 后，在每个径向层重构 $\Phi(r_i,\theta_j,\phi_k)$ 用逆球谐变换：

  $$
  \Phi(r_i,\theta,\phi)=\sum_{\ell m}\phi_{\ell m}(r_i) Y_{\ell m}(\theta,\phi).
  $$
* 计算梯度 $\nabla\Phi$（球坐标形式）：

  $$
  \begin{aligned}
  \partial_r\Phi(r_i,\theta_j,\phi_k) &\approx \frac{\phi(r_{i+1})-\phi(r_{i-1})}{2\Delta r},\\
  \frac{1}{r}\partial_\theta\Phi &\text{ and } \frac{1}{r\sin\theta}\partial_\phi\Phi
  \end{aligned}
  $$

  其中 $\partial_\theta,\partial_\phi$ 可由球谐系数解析求导（更精确且高效）：

  * 利用球谐的角导数表示：可直接在系数空间得到 $\partial_\theta Y_{\ell m}$ 与 $\partial_\phi Y_{\ell m}=im Y_{\ell m}$ 的乘法关系，从而在系数层面直接算出角导数的重构系数，避免在实空间做差分。
* 把 $\nabla\Phi$ 转为笛卡尔或直接以球坐标形式从 $\mathbf B^{(c)}$ 中减去（两种都可）。

---

## 6. 算法总体流程（伪码）

```text
Given Bc on grid (r_i, theta_j, phi_k):

1. Compute divergence rho = div(Bc) on grid (use same finite-volume/discrete operator as Bc construction).

2. Check total integral I = ∑ rho * dV.
   if abs(I) > tol:
       # enforce compatibility (simple fix)
       rho = rho - I / V_total

3. For each radial layer r_i:
      compute spherical harmonic transform of rho(r_i, theta, phi) -> { rho_lm(r_i) } for all (l,m)

4. For each (l=0..Lmax, m=-l..l):
      assemble tridiagonal system for phi_lm(r) from discrete radial ODE
      apply BCs (Neumann/Dirichlet as chosen)
      solve for phi_lm(r) (Thomas or LAPACK)

5. For each radial layer r_i:
      reconstruct Phi(r_i,theta,phi) by inverse spherical harmonic transform from {phi_lm(r_i)}

6. Compute gradient ∇Phi on grid:
      - radial derivative by FD in r
      - angular derivatives via spherical harmonic coefficient relations or central differences
   Convert ∇Phi to Cartesian if needed.

7. Compute B = Bc - ∇Phi

8. Optionally: compute residual divergence of B for diagnostics (should be ~0 up to solver tolerance)
```

---

## 7. 实现细节、性能与数值建议

* **球谐库**：强烈建议用 `shtns`（C + Python 接口）或 `pyshtools`，因为直接用 `scipy.special.sph_harm` 在大网格上会很慢。若用 FFT+Legendre 配对法，角向变换复杂度约 $O(N_\theta N_\phi \log N_\phi)$。
* **并行化**：对不同 $(\ell,m)$ 的径向方程可独立求解，或对不同 $m$ 并行；角向 SHT 可用 MPI 并行分工。
* **选择 $L_{\max}$**：设 $L_{\max}\leq \min(N_\theta-1,\;N_\phi/2)$。实验上 $L_{\max}$ 取到 $N_\theta-1$ 会更精确但费时。
* **径向网格**：建议非均匀（如靠近内边界更细）以解析近日冕尺度；对应差分公式要用非均匀格点版本（上面给出均匀示例，但可扩展）。
* **边界 BC 一致性**：若你要把该场作为 MHD 初始场，一般希望保留内边界 $B_r$（从观测给定），因此内边界应坚持 $\partial_r\Phi=0$（Neumann）或至少保证投影不改动已指定的法向值。外边界若是开边可选 Neumann。若 Neumann 兼容性难满足，用体积分布修正策略。
* **验证**：运行后检查两项关键量：

  1. $\max|\nabla\!\cdot\mathbf B|$ （应该接近 0, 例如 $<10^{-10}$ 相对量），
  2. 点对点 $|\mathbf B|-|\mathbf B_0|$ 的误差（投影会引入小幅偏差，记录并观察是否可接受）。
* **迭代精化**：若投影后模长偏差太大，可做：并行补偿（用投影后 B 更新 b0hat）→再求解一次投影，迭代 1–3 次通常可以同时降低散度与保模误差。

---

## 8. 简短 Python/伪实现提示（工具）

* 角向 SHT：`pyshtools` 或 `shtns`（推荐 `shtns` 性能最佳）；若没有，角向用 `scipy.special.lpmv` + FFT 可实现但慢。
* 径向 tridiag：`scipy.linalg.solve_banded` 或自写 Thomas 算法。
* Poisson 检查：在周期角向下，`l=0,m=0` 的兼容性很关键。
* 若你需要，我可以把上面流程**写成可运行的 3D Python 脚本**（示例使用 `shtns` 或在纯 NumPy 下用 `scipy.special` 的简化版本），并演示在你现有 Bc 上做投影的完整数值结果和诊断图（散度/模长偏差）。

---

如果你希望继续，我可以做两件事中的一件（任选）：

1. 给出**完整可运行的 Python 代码**（基于 `numpy` + `scipy`，角向用 `scipy.special.sph_harm` 实现的较慢但可运行的参考实现；或基于 `shtns` 的高效版本），包含：从 Bc 计算 ρ、球谐投影、径向三对角求解、逆变换、构造 ∇Φ、并给出投影前后诊断。
2. 给出**离散化的具体系数表达式**（用于非均匀径向网格）和伪码，便于你在 Fortran/C++ 或 MHD 代码中直接移植。

你想我把哪种实现直接写出来？ (若选 1，请说明你偏好是否安装 `shtns`，或仅用 SciPy/NumPy。）

--->
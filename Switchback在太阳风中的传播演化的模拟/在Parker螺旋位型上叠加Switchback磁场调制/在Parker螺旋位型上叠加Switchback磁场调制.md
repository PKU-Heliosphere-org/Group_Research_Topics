

---
# 球壳上在 Helmholtz 投影（Poisson 解）校正清楚磁场散度的数值方案

好——下面把**在球壳体积上**对 $\mathbf B^{(c)}$ 做 Helmholtz 投影（解 Poisson 方程 $\nabla^2\Phi=\nabla\!\cdot\mathbf B^{(c)}$ 并取 $\mathbf B=\mathbf B^{(c)}-\nabla\Phi$）的**完整、可实现**数值方案写清楚。重点采用对球面角用**球谐分解（spherical harmonic）**、径向用**有限差分（或谱）** 的混合谱-径向方法——这种方法在球壳上既准确又高效，且易于并行化。下面给出连续方程、每一步的离散化、边界条件选择、l=0 的特殊处理、实现要点与伪码/实现提示。

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
\frac{1}{r^2}\frac{d}{dr}\!\Big(r^2\frac{d\phi_{\ell m}}{dr}\Big)
-\frac{\ell(\ell+1)}{r^2}\,\phi_{\ell m}(r)
=\rho_{\ell m}(r).\ }
$$

等价写法：

$$
\frac{d}{dr}\!\big(r^2\phi_{\ell m}'\big)-\ell(\ell+1)\phi_{\ell m}=r^2\rho_{\ell m}(r).
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
  \left.\frac{\partial\Phi}{\partial r}\right|_{r=r_{\rm in}}=0,\qquad
  \left.\frac{\partial\Phi}{\partial r}\right|_{r=r_{\rm out}}=0,
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
\frac{1}{\Delta r}\Big( q_{i+\tfrac12}\frac{\phi_{i+1}-\phi_i}{\Delta r} - q_{i-\tfrac12}\frac{\phi_i-\phi_{i-1}}{\Delta r}\Big)
- \ell(\ell+1)\phi_i = r_i^2\rho_{\ell m}(r_i),
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

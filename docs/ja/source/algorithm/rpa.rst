.. highlight:: none

.. _algorithm_sec:

乱雑位相近似法
==========================

乱雑位相近似(RPA)では相互作用のない状態を出発点に、電子相関効果による一体の演算子の揺らぎの応答を検出します。
UHF近似ではあらかじめ初期配置を予想しておく必要があるのに対して、RPA法では2次転移により生じる秩序相を推定することが可能です。
H-waveでは松原振動数を利用したRPA近似法を実装しており、解析接続によって実験で観測される動的な物理量との比較も行うことが可能です。

以下、アルゴリズムについて掲載します。
H-waveのRPAモードでは以下のHamiltonianを取り扱います。

.. math::
    \begin{aligned}
     {\cal H}&={\cal H}_0+{\cal H}_{\rm int},\\
     {\cal H}_0&=\sum_{\langle i\alpha;j\beta \rangle}
      (t_{ij}^{\alpha \beta}c_{i\alpha}^{\dagger}
      c_{j\beta}^{\mathstrut}+\mbox{H.c.}),\\
     {\cal H}_{\rm int}&=\sum_{ij}\sum_{\alpha, \alpha', \beta, \beta'}W_{ij}^{\beta\beta',\alpha\alpha'}\left(
      c_{i\alpha}^{\dagger}c_{i\alpha'}c_{j\beta'}^{\dagger}c_{j\beta}+\mbox{H.c.}\right)
    \end{aligned}

ここで、以下のフーリエ変換

.. math::
    \begin{aligned}
    c_{i\alpha}
    =\frac{1}{\sqrt{N_L}}\sum_{\bf{k}}
    e^{i \bf{k}\cdot \bf{r}_{i}}c_{\bf{k},\alpha}^{\mathstrut},
    \end{aligned}

を行うと、Hamiltonianは以下のように書き換えられます。

.. math::
    \begin{aligned}
     {\cal H}&=\sum_{\bf{k}\alpha\beta}
     (\varepsilon_{\alpha\beta}(\bf{k})c_{\bf{k}\alpha}^{\dagger}
     c_{\bf{k}\beta}^{\mathstrut}+\mbox{H.c.}) \nonumber\\
    &+\frac{1}{2N_L}\sum_{\bf{k}\bf{k}'\bf{q}}\sum_{\alpha\beta\alpha'\beta'}
     W^{\beta\beta',\alpha\alpha'}_{\bf{q}}
     c_{\bf{k}+\bf{q},\alpha}^{\dagger}
      c_{\bf{k},\alpha'}^{\mathstrut}
      c_{\bf{k}'-\bf{q},\beta'}^{\dagger}
      c_{\bf{k}',\beta}^{\mathstrut}
    \end{aligned}

として求められます。さて、RPA近似では :math:`{\cal H}_0` に対して、電子相関効果による密度揺らぎを検出します。
そのため、相互作用による散乱は :math:`{\cal H}_0` が対角化された空間で行う必要があるので、相互作用の項は以下のように近似されます。

.. math::
    \begin{aligned}
    &W^{\beta\beta',\alpha\alpha'}_{\bf{q}}c_{\bf{k}+\bf{q},\alpha}^{\dagger}c_{\bf{k},\alpha'}^{\mathstrut}
    c_{\bf{k}'-\bf{q},\beta'}^{\dagger} c_{\bf{k}',\beta}^{\mathstrut}\nonumber\\
    &\sim W^{\beta\beta',\alpha\alpha'}_{\bf{q}} \sum_{\gamma, \gamma'}
    u_{\alpha \gamma, \bf{k}+\bf{q}}^* d_{\bf{k}+\bf{q},\gamma}^{\dagger}
    u_{\alpha' \gamma, \bf{k}} d_{\bf{k},\gamma}^{\mathstrut}
    u_{\beta' \gamma', \bf{k}'-\bf{q}}^* d_{\bf{k}'-\bf{q},\gamma'}^{\dagger}
    u_{\beta  \gamma', \bf{k}'}d_{\bf{k}',\gamma'}^{\mathstrut}.
    \end{aligned}

ここで、

.. math::
    \begin{aligned}
    c_{\bf{k},\alpha} = \sum_{\gamma} u_{\alpha \gamma, \bf{k}} d_{\bf{k}, \gamma}
    \end{aligned}

であり、 :math:`d_{\bf{k}, \gamma}` は :math:`{\cal H}_0` を対角化する消滅演算子を表します( :math:`\gamma` は固有値のindex)。
このとき、一体の既約グリーン関数は以下のように与えられます。

.. math::
    \begin{aligned}
     G^{(0)\alpha\beta}_{\gamma}({\bf k}, i\omega_{n})=
      \frac{u^{\alpha\gamma}({\bf k})u^{*\beta\gamma}({\bf k})}{i\epsilon_{n}-\xi^{\gamma}({\bf k})+\mu}.
    \end{aligned}

既約感受率は対角化された成分で閉じる必要があるため、以下のように与えられます。

.. math::
    \begin{aligned}
     X^{(0)\alpha\alpha', \beta\beta'}({\bf q},i\omega_n)=
      -\frac{T}{N_L}
      \sum_{\gamma=1}^{n_{\rm orb}}\sum_{{\bf k},n}
      G^{(0)\alpha\beta}_{\gamma}({\bf k}+{\bf q}, i\omega_m+ i\epsilon_{n})
      G^{(0)\beta'\alpha'}_{\gamma}({\bf k}, i\epsilon_{n}),
    \end{aligned}

この既約感受率を用いることで、RPAで得られる感受率が以下のように得られます。

.. math::
    \begin{aligned}
    X^{\alpha\alpha', \beta\beta'}(q)&=
    X^{(0)\alpha\alpha', \beta\beta'}(q) - \sum_{\alpha_1'\beta_1'}
    X^{(0)\alpha\alpha', \beta_1\beta_1'}(q) W^{\beta_1\beta_1', \alpha_1\alpha_1'}_{\bf q}X^{\alpha_1 \alpha_1' , \beta \beta'}(q),
    \end{aligned}

ここで、 :math:`\alpha \alpha'` などをまとめて一つのindexにすると行列形式で表すことができ、
最終的には以下のような形式で感受率を得ることができます。

.. math::
    \begin{aligned}
     \hat{X}(q)&=\hat{X}^{(0)}(q)+\hat{X}^{(0)}(q)\hat{W}(q)\hat{X}(q)\nonumber\\
     &=\left[\hat{I}+\hat{X}^{(0)}(q)\hat{W}(q)\right]^{-1}\hat{X}^{(0)}(q).
    \end{aligned}

なお、より高次な相関効果を考慮する手法としてvertex補正の考慮などがあります。詳細については、例えばこちらの文献 [1]_ を参考にしてください。

.. [1] `K. Yoshimi, T. Kato, H. Maebashi, J. Phys. Soc. Jpn. {\bf 78}, 104002 (2009). <https://journals.jps.jp/doi/10.1143/JPSJ.78.104002>`_

.. highlight:: none

.. _algorithm_sec:

乱雑位相近似法
==========================

乱雑位相近似(RPA)では相互作用のない状態を出発点に、電子相関効果による一体の演算子の揺らぎの応答を検出します。
UHF近似ではあらかじめ初期配置を予想しておく必要があるのに対して、RPA法では2次転移により生じる秩序相を推定することが可能です。
H-waveでは松原振動数を利用したRPA法を実装しており、解析接続によって実験で観測される動的な物理量との比較も行うことが可能です。

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
     {\cal H}&=\sum_{{\bf k}\alpha\beta}
     (\varepsilon_{\alpha\beta}({\bf k})c_{{\bf k}\alpha}^{\dagger}
     c_{{\bf k}\beta}^{\mathstrut}+\mbox{H.c.}) \nonumber\\
    &+\frac{1}{2N_L}\sum_{{\bf k} {\bf k}'{\bf q}}\sum_{\alpha\beta\alpha'\beta'}
     W^{\beta\beta',\alpha\alpha'}_{{\bf q}}
     c_{{\bf k}+{\bf q},\alpha}^{\dagger}
      c_{{\bf k},\alpha'}^{\mathstrut}
      c_{{\bf k}'-{\bf q},\beta'}^{\dagger}
      c_{{\bf k}',\beta}^{\mathstrut}
    \end{aligned}

RPAでは :math:`{\cal H}_0` に対して、電子間相互作用を介した密度揺らぎの効果を考慮します。
具体的には、 :math:`{\cal H}_0` が対角化されるような軌道・スピンの混成基底を用いて、相互作用の項を以下のように近似します。

.. math::
    \begin{aligned}
    &W^{\beta\beta',\alpha\alpha'}_{\bf{q}}c_{\bf{k}+\bf{q},\alpha}^{\dagger}c_{\bf{k},\alpha'}^{\mathstrut}
    c_{\bf{k}'-\bf{q},\beta'}^{\dagger} c_{\bf{k}',\beta}^{\mathstrut}\nonumber\\
    &\sim W^{\beta\beta',\alpha\alpha'}_{\bf{q}} \sum_{\gamma, \gamma'}
    (u_{\alpha \gamma, \bf{k}+\bf{q}}^* d_{\bf{k}+\bf{q},\gamma}^{\dagger}
    u_{\alpha' \gamma, \bf{k}} d_{\bf{k},\gamma}^{\mathstrut})
    (u_{\beta' \gamma', \bf{k}'-\bf{q}}^* d_{\bf{k}'-\bf{q},\gamma'}^{\dagger}
    u_{\beta  \gamma', \bf{k}'}d_{\bf{k}',\gamma'}^{\mathstrut}) .
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
    X^{(0)\alpha\alpha', \beta\beta'}(q) - \sum_{\alpha_1,\alpha_1', \beta_1,\beta_1'}
    X^{(0)\alpha\alpha', \beta_1\beta_1'}(q) W^{\beta_1\beta_1', \alpha_1\alpha_1'}_{\bf q}X^{\alpha_1 \alpha_1' , \beta \beta'}(q),
    \end{aligned}

ここで、 :math:`\alpha \alpha'` などをまとめて一つのindexにすると行列形式で表すことができ、
最終的には以下のような形式で感受率を得ることができます。

.. math::
    \begin{aligned}
     \hat{X}(q)&=\hat{X}^{(0)}(q)-\hat{X}^{(0)}(q)\hat{W}(q)\hat{X}(q)\nonumber\\
     &=\left[\hat{I}+\hat{X}^{(0)}(q)\hat{W}(q)\right]^{-1}\hat{X}^{(0)}(q).
    \end{aligned}

    
上記の実装では、軌道とスピンを統一した一般化軌道として取り扱いました。計算の実行に必要な配列のうち、 感受率( :math:`X^{(0)\alpha\alpha', \beta\beta'}({\bf q},i\omega_n), X^{\alpha\alpha', \beta\beta'}({\bf q},i\omega_n)` )が一番大きなサイズの多次元配列となり、そのサイズは :math:`N_{\rm orb}^4 N_{\rm spin}^4 N_k N_{\omega}` で与えられ、サイズが大きくなるとメモリコスト、計算量が増大します。以下で説明するように、軌道とスピンを分離することで感受率の多次元配列のサイズを減らすことができます。
H-waveのRPAモードで取り扱う二体相互作用では、軌道とスピンを分離することで、

.. math::
    \begin{aligned}
    & W^{\beta\sigma_1\sigma_1',\alpha\sigma\sigma'}_{\bf{q}}c_{\bf{k}+\bf{q},\alpha \sigma}^{\dagger}c_{\bf{k},\alpha \sigma'}^{\mathstrut}
    c_{\bf{k}'-\bf{q},\beta\sigma_1'}^{\dagger} c_{\bf{k}',\beta\sigma_1}^{\mathstrut}    \end{aligned}

と書けます。軌道に対しては同一の軌道での散乱となるため、既約感受率は

.. math::
    \begin{aligned}
     X^{(0)\alpha, \beta}_{\sigma\sigma'\sigma_1\sigma_1'}({\bf q},i\omega_n)=
      -\frac{T}{N_L}
      \sum_{\gamma=1}^{n_{\rm orb}}\sum_{{\bf k},n}
      G^{(0)\alpha\beta}_{\sigma\sigma_1', \gamma}({\bf k}+{\bf q}, i\omega_m+ i\epsilon_{n})
      G^{(0)\beta\alpha}_{\sigma_1\sigma', \gamma}({\bf k}, i\epsilon_{n}),
    \end{aligned}

となり、 :math:`N_{\rm orb}^2 N_{\rm spin}^4 N_k N_{\omega}` にサイズを抑えることができます。このとき、RPAで得られる感受率は

.. math::
    \begin{aligned}
    X^{\alpha, \beta}_{\sigma\sigma'\sigma_1\sigma_1'}(q)&=
    X^{(0)\alpha, \beta}_{\sigma\sigma'\sigma_1\sigma_1'}(q) - \sum_{\alpha_1'\beta_1'}
    X^{(0)\alpha, \alpha_2}_{\sigma\sigma'\sigma_2\sigma_2'}(q) W^{\alpha_2, \alpha_3}_{\sigma_2\sigma_2', \sigma_3\sigma_3'}({\bf q})X^{\alpha_3, \beta}_{\sigma_3\sigma_3',\sigma_1\sigma_1'}(q),
    \end{aligned}

となります。 :math:`\alpha\sigma\sigma'` を一つのindexとみなせば、行列形式にすることができ、一般化軌道の場合と同様に、

.. math::
    \begin{aligned}
     \hat{X}(q)&=\hat{X}^{(0)}(q)-\hat{X}^{(0)}(q)\hat{W}(q)\hat{X}(q)\nonumber\\
     &=\left[\hat{I}+\hat{X}^{(0)}(q)\hat{W}(q)\right]^{-1}\hat{X}^{(0)}(q).
    \end{aligned}

と書けることがわかります。以上が一般的なRPAの定式化になります。

上述の近似では既約感受率の計算を

.. math::
    \begin{aligned}
     X^{(0)\alpha, \beta}_{\sigma\sigma'\sigma_1\sigma_1'}({\bf q},i\omega_n)=
      -\frac{T}{N_L}
      \sum_{\gamma=1}^{n_{\rm orb}}\sum_{{\bf k},n}
      G^{(0)\alpha\beta}_{\sigma\sigma_1', \gamma}({\bf k}+{\bf q}, i\omega_m+ i\epsilon_{n})
      G^{(0)\beta\alpha}_{\sigma_1\sigma', \gamma}({\bf k}, i\epsilon_{n})\nonumber
    \end{aligned}

として行っています。この場合、対角化した成分の和が必要となり、計算コストが多くかかってしまいます。
そのため、先行研究の多くは一体グリーン関数を

.. math::
    \begin{aligned}
     G^{(0)\alpha\beta}_{\sigma\sigma'}({\bf k}, i\omega_{n}) = \sum_{\gamma=1}^{n_{\rm orb}} G^{(0)\alpha\beta}_{\sigma\sigma', \gamma}({\bf k}, i\omega_{n})
    \end{aligned}

のように近似し、既約感受率を

.. math::
    \begin{aligned}
     X^{(0)\alpha, \beta}_{\sigma\sigma'\sigma_1\sigma_1'}({\bf q},i\omega_n)=
      -\frac{T}{N_L}
      \sum_{{\bf k},n}
      G^{(0)\alpha\beta}_{\sigma\sigma_1'}({\bf k}+{\bf q}, i\omega_m+ i\epsilon_{n})
      G^{(0)\beta\alpha}_{\sigma_1\sigma'}({\bf k}, i\epsilon_{n})\nonumber
    \end{aligned}

として計算して高速化する場合が多いです。
この既約感受率を用いた計算では、対角化成分が混在してしまう状況で近似精度が悪くなりますが、
バンド交差による :math:`\gamma` への技術的な対応を行う必要がないというメリットもあります。
先行研究との比較をするためにも、H-Waveではこの手法を採用しています(グリーン関数と既約感受率を正しく取り扱うモードについても実装する予定です)。
なお、より高次な相関効果を考慮する手法としてvertex補正の考慮などがあります。詳細については、例えばこちらの文献 [1]_ を参考にしてください。

.. [1] `K. Yoshimi, T. Kato, H. Maebashi, J. Phys. Soc. Jpn. 78, 104002 (2009). <https://journals.jps.jp/doi/10.1143/JPSJ.78.104002>`_

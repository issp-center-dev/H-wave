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

.. note::

   この :math:`e^{+i\bf{k}\cdot\bf{r}}` 規約（wannier90 形式の符号、
   UHFk と共通）は、issue #133 以降 RPA モジュールのすべての実空間係数
   構築（:math:`R \to k/q` 変換）が従う規約です：
   :math:`\varepsilon({\bf k}) = \sum_{\bf R} t({\bf R})
   e^{+i {\bf k}\cdot{\bf R}}` および :math:`W({\bf q}) =
   \sum_{\bf R} W({\bf R}) e^{+i {\bf q}\cdot{\bf R}}`（感受率計算
   内部の畳み込み変換は自己逆な対であり影響を受けません）。この修正以前
   は、非スピン軌道経路が全体で逆符号を使用しており---自己整合的な
   :math:`{\bf k} \to -{\bf k}` の再ラベルであり、保存された
   ``chi0q`` / ``chiq`` は、テンソルが FFT グリッド上で
   :math:`{\bf q} \to -{\bf q}` に対して要素毎に偶でない限り、運動量
   ラベルが反転していました---、スピン軌道経路は transfer と相互作用で
   2 つの符号が混在し、``chiq`` が :math:`\chi_0({\bf q})` と
   :math:`W(-{\bf q})` から解かれていました：
   :math:`W({\bf q}) \neq W(-{\bf q})` となる相互作用（方向性ボンド）
   では、旧スピン軌道 ``chiq`` は再ラベルではなく誤りです。この修正以降
   に書かれる運動量空間ファイルは ``momentum_convention = "e_plus_ikR"``
   を持ち、ローダは無印の旧ファイルを、内容が :math:`{\bf q}` 偶
   （両規約が一致する場合）でない限り拒否します。

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

    
.. note::

   **相互作用と感受率の添字対の規約について。**
   上の2つの量は、各ペアの内部で互いに逆向きの添字順序を用いており、
   その変換は RPA 方程式の一部です。

   既約感受率の2つのペアスロットには以下の双一次形式が対応します。
   左のペアは生成演算子の添字が **後**、右のペアは生成演算子の添字が
   **前** です。

   .. math::
      X^{(0)\alpha\alpha',\beta\beta'}(q) \;\sim\;
      \Big\langle \big(c^{\dagger}_{\alpha'}c^{\mathstrut}_{\alpha}\big)(q)\;;\;
      \big(c^{\dagger}_{\beta}c^{\mathstrut}_{\beta'}\big)(-q)\Big\rangle

   これは上の Green 関数の積 :math:`G^{\alpha\beta}(k+q)
   G^{\beta'\alpha'}(k)` が表すものです。一方、相互作用
   :math:`W^{\beta\beta',\alpha\alpha'}_{\bf q}` が掛かる演算子は
   :math:`c^{\dagger}_{\alpha}c^{\mathstrut}_{\alpha'}
   c^{\dagger}_{\beta'}c^{\mathstrut}_{\beta}` であり、**どちらのペアでも**
   感受率とは逆の順序になっています。したがって RPA
   方程式に入る行列 :math:`\hat W(q)` は、相互作用テンソルの
   **各ペア内で添字を入れ替えた** ものになります。

   .. math::
      \hat{W}(q)_{(\beta\beta'),(\alpha\alpha')}
      = W^{\beta'\beta,\,\alpha'\alpha}_{\bf q}

   これは H-wave の論文 (`arXiv:2308.00324 [cond-mat.str-el]
   <https://arxiv.org/abs/2308.00324>`_) の式 (16) と (20) の間に現れる
   並べ替えと同じもので、2つのペアの両方が対象です。

   この変換は密度型（ペア対角）成分に対しては恒等変換であり、より一般に
   実数のエルミート閉じた宣言に対しても値が変わりません（転置先のスロットが
   同じ値を持つため）。効いてくるのは **複素の** ペア交差型相互作用、
   すなわち複素のエルミート閉じた ``PairHop`` の場合だけで、変換を省くと
   複素共役なハミルトニアンに対する感受率が得られてしまいます。
   H-wave は ring および ladder の求解でバーテックスを組み立てる際に
   この変換を適用します。保存される ``chiq``/``chi0q`` と相互作用ファイルの
   規約は、それぞれの節に記載のとおりで変わりません。

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

ブロック対角化最適化
*****************************

相互作用ハミルトニアンがブロック対角構造を持つ場合
（例：スピン保存や軌道間結合がない場合）、
RPA方程式を各ブロックで独立に解くことができ、
計算コストを大幅に削減できます。

ブロック構造は相互作用行列の接続性を解析して自動検出されます:

1. 全k点にわたって相互作用ハミルトニアンの絶対値を合計し、接続パターン行列を得る。
2. 非ゼロの非対角要素（閾値: :math:`10^{-12}` ）から隣接グラフを構築する。
3. ラベル伝播（union-findアルゴリズム）により連結成分を求める。

行列が :math:`m` 個のブロック（サイズ :math:`n_1, n_2, \ldots, n_m` ）に
分解される場合、RPA方程式の計算コストは
:math:`O(N^3)` から
:math:`O(n_1^3 + n_2^3 + \cdots + n_m^3)` に削減されます。
ここで :math:`N = n_1 + n_2 + \cdots + n_m` です。

この最適化は自動的に適用され、ユーザーに対して透過的です。


横感受率（はしごダイアグラム）
*******************************************

標準的な（リングダイアグラム）RPA感受率に加えて、
H-waveは横感受率 :math:`\chi_{+-}(\mathbf{q})` を計算できます。
これはスピン反転相関
:math:`\langle S^+(\mathbf{q}) S^-(-\mathbf{q}) \rangle` を記述します。

横方向の裸感受率は

.. math::

   X^{(0)}_{+-,\alpha\gamma;\beta\delta}(\mathbf{q}, i\omega_n)
   = -\frac{T}{N_L} \sum_{\mathbf{k},n}
     G_{\alpha\beta,\uparrow}(\mathbf{k}+\mathbf{q}, i\omega_m + i\varepsilon_n)\,
     G_{\delta\gamma,\downarrow}(\mathbf{k}, i\varepsilon_n)

横方向の頂点 :math:`W_{+-}` は、
縦チャネルのハートリー（フォック交換）頂点の交差から得られます:

.. math::

   W_{+-} = W_{\uparrow\uparrow\uparrow\uparrow} - W_{\downarrow\downarrow\uparrow\uparrow}^{\rm crossed}

横方向の頂点は、相互作用テンソルの異スピンブロックとスピン反転ブロックのみ
から構成されます。同スピンブロックは寄与しません。同スピン相互作用は横方向
ループの上向き伝播関数と下向き伝播関数を結べないため、自己エネルギーは生じ
ますが頂点は生じません。

軌道ペアは2つの宣言の平均で対称化されます。相互作用ファイルでは同一の演算子
を2通りに書けるためです （ :math:`n_a n_b = n_b n_a` 、Exchange では
:math:`X_{ab} = X_{ba}` ）。これは UHFk が用いている規約と同じです。平均を取る
相手は相互作用型ごと（同値には、その型が占めるスロット族ごと）に異なります。密度-密度型と Exchange は単純な転置との平均、
PairHop は共役転置との平均です。PairHop の2つの宣言は同一係数ではなくエルミート
共役の組 （ :math:`P_{ba} = P_{ab}^{*}` ）だからです。したがって複素エルミート
閉じの Exchange の物理的結合 :math:`(J_{01} + J_{10})/2` は実数になり、複素
エルミート閉じの PairHop は完全な複素値を保ちます。

オンサイト相互作用に対する頂点は以下のとおりです:

- ``CoulombIntra`` :math:`U` : :math:`W_{+-} = -U`
- ``CoulombInter`` :math:`V` : :math:`W_{+-} = -V`
- ``Hund`` :math:`J` : :math:`W_{+-} = 0`
- ``Exchange`` :math:`J` : :math:`W_{+-} = -(J + J^{\rm T})/2`
- ``Ising`` :math:`I` : :math:`W_{+-} = +I` （wannier90 形式の k 空間ソルバーはすべてこの規格化で Ising ファイル
  を読みます。UHFk の因子 1/4 の不一致は issue #106 で解消済みです。
  別系統の実空間 UHFr リーダーは S^z 規約を保持します）
- ``PairLift`` :math:`J` : :math:`W_{+-} = 0`
- ``PairHop`` :math:`J` : :math:`W_{+-} = -J`

.. note::

   これらの値は、横方向チャネルを厳密対角化と照合する前に公開されていた値とは
   異なります。従来の記載は7種のうち4種が誤りで、1種が欠落していました。
   ``CoulombInter`` 、 ``Hund`` 、 ``Ising`` 、 ``Exchange`` を含む計算で得られた
   横方向帯磁率は再計算してください。影響を受けるのは ``chiq_pm`` のみで、
   ``chiq`` 、自己エネルギー、Eliashberg 頂点には波及しません。

``calc_type = "ring+ladder"`` では、縦方向の計算より前に、 **組み立てられた
横方向頂点が** （副格子で折り畳んだ格子上で） :math:`q` **に依存しないこと**
を相対許容値 :math:`10^{-10}` で検証し、満たさない入力を拒否します。横方向の
ペア :math:`c^\dagger_{i a \uparrow} c_{j b \downarrow}` はオフサイト項に
対して非局所となり、その頂点は :math:`q` のみの関数では表現できないためです。
実際にはオフサイトの ``CoulombInter`` 、 ``Ising`` 、 ``Exchange`` が拒否され、
オフサイトの ``Hund`` と ``PairLift`` は横方向頂点が消えるため受理されます。
なお、互いに打ち消し合う・値が食い違う宣言の組は、より早い読み込み時に
拒否されます（issue #93: 宣言ファイルはエルミート共役で閉じている必要が
あります）。また ``SubShape`` による折り畳みでスーパーセル内に収まる
サイト間ペアはセル内軌道ペアとなり、表現可能として受理されます。縦方向（ ``ring`` ）チャネルは影響を受けません。

.. warning::

   オフサイトの ``PairHop`` は相互作用の読み込み時にこの検査より前に暗黙に
   破棄されるため、拒否も反映もされません。RPA 計算でオフサイトの
   ``PairHop`` に依存しないでください。また、対角（同一軌道）の PairHop は
   密度項 :math:`2P\, n_\uparrow n_\downarrow` を意味しますが、読み込み時の
   係数は :math:`2P` ではなく :math:`P` として（縦・横チャネルで一貫して）
   扱われます。この縮退エントリの検証は別途追跡されています。

横方向のRPA感受率は

.. math::

   \hat{X}_{+-}(\mathbf{q})
   = \left[\hat{I} + \hat{X}^{(0)}_{+-}(\mathbf{q})\, \hat{W}_{+-}\right]^{-1}
     \hat{X}^{(0)}_{+-}(\mathbf{q})

横チャネルの計算を有効にするには、入力TOMLファイルで
``calc_type = "ring+ladder"`` と設定します。
これには ``general`` 計算スキームが必要です（自動選択されます）。

.. note::

   スピン軌道モードでハミルトニアンが本当にスピンを混合する場合
   （スピン軌道相互作用など）、横チャネルはバブルの :math:`S_z`
   保存ブロック :math:`G_\uparrow G_\downarrow` のみを抽出します。
   スピン混合のクロス項は含まれず、警告が出力されます。したがって
   スピン混合系の横感受率は、現在の実装では近似となります。


スピン軌道モード
*****************************

H-waveはスピンと軌道のインデックスが
ブロック分離ではなくインターリーブされるスピン軌道モードをサポートしています。

通常モードでは、複合インデックスは :math:`i = s \cdot n_{\rm orb} + a`
（スピンブロック優先）であり、
:math:`s = 0, 1` はスピンインデックス、
:math:`a = 0, \ldots, n_{\rm orb}-1` は軌道インデックスです。
スピン軌道モードでは、インデックスは :math:`i = 2a + s`
（インターリーブ）であり、スピン軌道相互作用を自然に扱えます。

スピン軌道モードは入力TOMLファイルで ``enable_spin_orbital = true``
と設定して有効化します。このモードでは:

- ハミルトニアンはスピン保存を仮定せず、
  完全な :math:`2n_{\rm orb} \times 2n_{\rm orb}` 空間で構成されます。
- 全ての相互作用型（ ``CoulombIntra`` 、 ``CoulombInter`` 、 ``Hund`` 、 ``Exchange`` 、
  ``Ising`` 、 ``PairLift`` 、 ``PairHop`` ）がサポートされます。
- 可能な場合、ブロック対角化最適化が自動的に適用されます。
- ``squashed`` 計算スキームもスピン軌道系で利用可能です。
- 幾何情報ファイル（``geom.dat``）の ``Norbit`` はスピン軌道の総数
  （= 物理軌道数 × 2 = Wannier90 の ``num_wann``）を表し、UHFk と同じ規約です。

.. note::

   **移行上の注意（RPA）：** スピン軌道入力の幾何 ``Norbit`` はスピン軌道数になりました。
   既存の RPA スピン軌道計算の ``geom.dat`` の ``Norbit`` は2倍にしてください。


.. [1] `K. Yoshimi, T. Kato, H. Maebashi, J. Phys. Soc. Jpn. 78, 104002 (2009). <https://journals.jps.jp/doi/10.1143/JPSJ.78.104002>`_

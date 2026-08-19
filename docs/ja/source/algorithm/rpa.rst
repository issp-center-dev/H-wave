.. highlight:: none

.. _algorithm_sec:

乱雑位相近似法
==========================

乱雑位相近似(RPA)では相互作用のない状態を出発点に、電子相関効果による一体の演算子の揺らぎの応答を検出します。
UHF近似ではあらかじめ初期配置を予想しておく必要があるのに対して、RPA法では2次転移により生じる秩序相を推定することが可能です。
H-waveでは松原振動数を利用したRPA法を実装しており、解析接続によって実験で観測される動的な物理量との比較も行うことが可能です。

.. note::

   RPAソルバーと1回反復のFLEXソルバーについて、どの設定がそれぞれで
   実行可能か、また両者が実行できる場合にどの程度一致するかを、小規模な
   参照入力に対してまとめた対応表があります。この表は日本語版には未収録
   のため、英語版マニュアルの
   `RPA and FLEX: support and equivalence <../../../en/html/algorithm/rpa_flex_equivalence.html>`_
   を参照してください。

以下、アルゴリズムについて掲載します。
H-waveのRPAモードでは以下のHamiltonianを取り扱います。

.. math::
    \begin{aligned}
     {\cal H}&={\cal H}_0+{\cal H}_{\rm int},\\
     {\cal H}_0&=\sum_{i\alpha;j\beta}
      t_{ij}^{\alpha \beta}c_{i\alpha}^{\dagger}
      c_{j\beta}^{\mathstrut},\\
     {\cal H}_{\rm int}&=\frac{1}{2}\sum_{ij}\sum_{\alpha, \alpha', \beta, \beta'}W_{ij}^{\beta\beta',\alpha\alpha'}
      c_{i\alpha}^{\dagger}c_{i\alpha'}c_{j\beta'}^{\dagger}c_{j\beta}
    \end{aligned}

入力の表には、1つの相互作用が2通りの並び順の両方で現れます。オンサイトの
クーロン項であれば、up-down のスロットと down-up のスロットの両方に入ります。
上の和は4つの添字を制限なく走るので、同じ項を2度数えることになり、
:math:`1/2`\ がそれを打ち消します。トランスファーの表も同様に
:math:`t_{\alpha\beta}(R)=t_{\beta\alpha}(-R)^{*}`\ を満たすように与えます。
これも両向きが揃っているということなので、 :math:`{\cal H}_0`\ に
エルミート共役項を別途書き足す必要はありません。

ここで、以下のフーリエ変換

.. math::
    \begin{aligned}
    c_{i\alpha}
    =\frac{1}{\sqrt{N_L}}\sum_{\bf{k}}
    e^{i \bf{k}\cdot \bf{r}_{i}}c_{\bf{k},\alpha}^{\mathstrut},
    \end{aligned}

を行うと、Hamiltonianは以下のように書き換えられます。

.. note::

   RPA モジュールの実空間係数構築（:math:`R \to k/q`\ 変換）は、すべて
   この :math:`e^{+i\bf{k}\cdot\bf{r}}`\ 規約に従います（wannier90 形式の
   符号、UHFk と共通）：
   :math:`\varepsilon({\bf k}) = \sum_{\bf R} t({\bf R})
   e^{+i {\bf k}\cdot{\bf R}}` および :math:`W({\bf q}) =
   \sum_{\bf R} W({\bf R}) e^{+i {\bf q}\cdot{\bf R}}`\ （感受率計算
   内部の畳み込み変換は自己逆な対であり影響を受けません）。1.0.x
   では、非スピン軌道経路が全体で逆符号を使用しており---自己整合的な
   :math:`{\bf k} \to -{\bf k}`\ の再ラベルであり、保存された
   ``chi0q`` / ``chiq``\ は、テンソルが FFT グリッド上で
   :math:`{\bf q} \to -{\bf q}`\ に対して要素毎に偶でない限り、運動量
   ラベルが反転していました---、スピン軌道経路は transfer と相互作用で
   2 つの符号が混在し、``chiq`` が :math:`\chi_0({\bf q})`\ と
   :math:`W(-{\bf q})`\ から解かれていました：
   :math:`W({\bf q}) \neq W(-{\bf q})`\ となる相互作用（方向性ボンド）
   では、旧スピン軌道 ``chiq``\ は再ラベルではなく誤りです。現行版
   で書かれる運動量空間ファイルは ``momentum_convention = "e_plus_ikR"``
   を持ち、ローダは無印の旧ファイルを、内容が :math:`{\bf q}`\ 偶
   （両規約が一致する場合）でない限り拒否します。

.. math::
    \begin{aligned}
     {\cal H}&=\sum_{{\bf k}\alpha\beta}
     \varepsilon_{\alpha\beta}({\bf k})c_{{\bf k}\alpha}^{\dagger}
     c_{{\bf k}\beta}^{\mathstrut} \nonumber\\
    &+\frac{1}{2N_L}\sum_{{\bf k} {\bf k}'{\bf q}}\sum_{\alpha\beta\alpha'\beta'}
     W^{\beta\beta',\alpha\alpha'}_{{\bf q}}
     c_{{\bf k}+{\bf q},\alpha}^{\dagger}
      c_{{\bf k},\alpha'}^{\mathstrut}
      c_{{\bf k}'-{\bf q},\beta'}^{\dagger}
      c_{{\bf k}',\beta}^{\mathstrut}
    \end{aligned}

RPAでは :math:`{\cal H}_0`\ に対して、電子間相互作用を介した密度揺らぎの効果を考慮します。
具体的には、 :math:`{\cal H}_0`\ が対角化されるような軌道・スピンの混成基底を用いると、相互作用の項は厳密な基底変換により以下のように書けます。

.. math::
    \begin{aligned}
    &W^{\beta\beta',\alpha\alpha'}_{\bf{q}}c_{\bf{k}+\bf{q},\alpha}^{\dagger}c_{\bf{k},\alpha'}^{\mathstrut}
    c_{\bf{k}'-\bf{q},\beta'}^{\dagger} c_{\bf{k}',\beta}^{\mathstrut}\nonumber\\
    &= W^{\beta\beta',\alpha\alpha'}_{\bf{q}}
    \sum_{\gamma_1 \gamma_2 \gamma_1' \gamma_2'}
    (u_{\alpha \gamma_1, \bf{k}+\bf{q}}^* d_{\bf{k}+\bf{q},\gamma_1}^{\dagger}
    u_{\alpha' \gamma_2, \bf{k}} d_{\bf{k},\gamma_2}^{\mathstrut})
    (u_{\beta' \gamma_1', \bf{k}'-\bf{q}}^* d_{\bf{k}'-\bf{q},\gamma_1'}^{\dagger}
    u_{\beta  \gamma_2', \bf{k}'}d_{\bf{k}',\gamma_2'}^{\mathstrut}) .
    \end{aligned}

各双一次形式は2つの独立なバンド添字を持ちます。

ここで、

.. math::
    \begin{aligned}
    c_{\bf{k},\alpha} = \sum_{\gamma} u_{\alpha \gamma, \bf{k}} d_{\bf{k}, \gamma}
    \end{aligned}

であり、 :math:`d_{\bf{k}, \gamma}` は :math:`{\cal H}_0`\ を対角化する消滅演算子を表します( :math:`\gamma`\ は固有値のindex)。
このとき、一体の既約グリーン関数は以下のように与えられます。

.. math::
    \begin{aligned}
     G^{(0)\alpha\beta}_{\gamma}({\bf k}, i\epsilon_{n})=
      \frac{u^{\alpha\gamma}({\bf k})u^{*\beta\gamma}({\bf k})}{i\epsilon_{n}-\xi^{\gamma}({\bf k})+\mu}.
    \end{aligned}

完全な一体 Green 関数はこれらのバンド寄与の和

.. math::
    \begin{aligned}
     G^{(0)\alpha\beta}({\bf k}, i\epsilon_{n})
      = \sum_{\gamma=1}^{n_{\rm orb}} G^{(0)\alpha\beta}_{\gamma}({\bf k}, i\epsilon_{n}),
    \end{aligned}

であり、既約感受率はこれから作られる粒子・正孔バブル

.. math::
    \begin{aligned}
     X^{(0)\alpha\alpha', \beta\beta'}({\bf q},i\omega_m)=
      -\frac{T}{N_L}
      \sum_{{\bf k},n}
      G^{(0)\alpha\beta}({\bf k}+{\bf q}, i\omega_m+ i\epsilon_{n})
      G^{(0)\beta'\alpha'}({\bf k}, i\epsilon_{n}),
    \end{aligned}

として与えられます。ここに含まれる2つのバンド和は独立で、バブルの粒子と正孔は
異なるバンドに乗ることができます。

H-wave は両方の和を保持します。バンド添字は現れた直後に和を取って消えるため、
相互作用は一貫して軌道基底で扱われ、対角化は Green 関数の構成にのみ用いられます。

この既約感受率を用いることで、RPAで得られる感受率が以下のように得られます。

.. math::
    \begin{aligned}
    X^{\alpha\alpha', \beta\beta'}(q)&=
    X^{(0)\alpha\alpha', \beta\beta'}(q) - \sum_{\alpha_1,\alpha_1', \beta_1,\beta_1'}
    X^{(0)\alpha\alpha', \beta_1\beta_1'}(q)
    W^{\beta_1'\beta_1, \alpha_1'\alpha_1}_{\bf q}
    X^{\alpha_1 \alpha_1' , \beta \beta'}(q),
    \end{aligned}

ここで、 :math:`\alpha \alpha'`\ などをまとめて一つのindexにすると行列形式で表すことができ、
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
      \Big\langle \big(c^{\dagger}_{\alpha'}c^{\mathstrut}_{\alpha}\big)(-q)\;;\;
      \big(c^{\dagger}_{\beta}c^{\mathstrut}_{\beta'}\big)(q)\Big\rangle

   これは上の Green 関数の積 :math:`G^{\alpha\beta}(k+q)
   G^{\beta'\alpha'}(k)` が表すものです。一方、相互作用
   :math:`W^{\beta\beta',\alpha\alpha'}_{\bf q}`\ が掛かる演算子は
   :math:`c^{\dagger}_{\alpha}c^{\mathstrut}_{\alpha'}
   c^{\dagger}_{\beta'}c^{\mathstrut}_{\beta}` であり、**どちらのペアでも**
   感受率とは逆の順序になっています。したがって上の RPA 方程式に入る
   バーテックスは、相互作用テンソルの **各ペア内で添字を入れ替えた**
   もの（式中の :math:`W^{\beta_1'\beta_1,\alpha_1'\alpha_1}_{\bf q}`\ ）
   であり、行列としては次のように書けます。

   .. math::
      \hat{W}(q)_{(\beta\beta'),(\alpha\alpha')}
      = W^{\beta'\beta,\,\alpha'\alpha}_{\bf q}

   この変換は密度型（ペア対角）成分に対しては恒等変換であり、より一般に
   実数のエルミート閉じた宣言に対しても値が変わりません（転置先のスロットが
   同じ値を持つため）。効いてくるのは **複素の** ペア交差型相互作用、
   すなわち複素のエルミート閉じた ``PairHop``\ の場合だけで、変換を省くと
   複素共役なハミルトニアンに対する感受率が得られてしまいます。
   H-wave は ring の求解でバーテックスを組み立てる際にこの変換を適用します。
   transverse (ladder) の組み立てはこの変換を通しません。こちらは相互作用
   テンソル自身を組み替えるため、ハミルトニアン規約のまま受け取ります。
   保存される ``chiq``/``chi0q``\ と相互作用ファイルの規約は、それぞれの節に
   記載のとおりで変わりません。

上記の実装では、軌道とスピンを統一した一般化軌道として取り扱いました。計算の実行に必要な配列のうち、 感受率( :math:`X^{(0)\alpha\alpha', \beta\beta'}({\bf q},i\omega_m), X^{\alpha\alpha', \beta\beta'}({\bf q},i\omega_m)` )が一番大きなサイズの多次元配列となり、そのサイズは :math:`N_{\rm orb}^4 N_{\rm spin}^4 N_k N_{\omega}`\ で与えられ、サイズが大きくなるとメモリコスト、計算量が増大します。以下で説明するように、軌道とスピンを分離することで感受率の多次元配列のサイズを減らすことができます。
H-waveのRPAモードで取り扱う二体相互作用では、軌道とスピンを分離することで、

.. math::
    \begin{aligned}
    & W^{\beta\sigma_1\sigma_1',\alpha\sigma\sigma'}_{\bf{q}}c_{\bf{k}+\bf{q},\alpha \sigma}^{\dagger}c_{\bf{k},\alpha \sigma'}^{\mathstrut}
    c_{\bf{k}'-\bf{q},\beta\sigma_1'}^{\dagger} c_{\bf{k}',\beta\sigma_1}^{\mathstrut}    \end{aligned}

と書けます。軌道に対しては同一の軌道での散乱となるため、既約感受率は

.. math::
    \begin{aligned}
     X^{(0)\alpha, \beta}_{\sigma\sigma'\sigma_1\sigma_1'}({\bf q},i\omega_m)=
      -\frac{T}{N_L}
      \sum_{{\bf k},n}
      G^{(0)\alpha\beta}_{\sigma\sigma_1'}({\bf k}+{\bf q}, i\omega_m+ i\epsilon_{n})
      G^{(0)\beta\alpha}_{\sigma_1\sigma'}({\bf k}, i\epsilon_{n}),
    \end{aligned}

となり、 :math:`N_{\rm orb}^2 N_{\rm spin}^4 N_k N_{\omega}`\ にサイズを抑えることができます。このとき、RPAで得られる感受率は

.. math::
    \begin{aligned}
    X^{\alpha, \beta}_{\sigma\sigma'\sigma_1\sigma_1'}(q)&=
    X^{(0)\alpha, \beta}_{\sigma\sigma'\sigma_1\sigma_1'}(q)
    - \sum_{\alpha_2,\alpha_3}
      \sum_{\sigma_2\sigma_2'\sigma_3\sigma_3'}
    X^{(0)\alpha, \alpha_2}_{\sigma\sigma'\sigma_2\sigma_2'}(q)
    W^{\alpha_2, \alpha_3}_{\sigma_2'\sigma_2, \sigma_3'\sigma_3}({\bf q})
    X^{\alpha_3, \beta}_{\sigma_3\sigma_3',\sigma_1\sigma_1'}(q),
    \end{aligned}

となります。 :math:`\alpha\sigma\sigma'`\ を一つのindexとみなせば、行列形式にすることができ、一般化軌道の場合と同様に、

.. math::
    \begin{aligned}
     \hat{X}(q)&=\hat{X}^{(0)}(q)-\hat{X}^{(0)}(q)\hat{W}(q)\hat{X}(q)\nonumber\\
     &=\left[\hat{I}+\hat{X}^{(0)}(q)\hat{W}(q)\right]^{-1}\hat{X}^{(0)}(q).
    \end{aligned}

と書けることがわかります。以上が一般的なRPAの定式化になります。

なお、より高次な相関効果を考慮する手法としてvertex補正の考慮などがあります。詳細については、例えばこちらの文献 [1]_ を参考にしてください。

セクター構造とブロック分解
*****************************

感受率が複数のブロックに分解できる場合、RPA 方程式は各ブロックで独立に解くことが
でき、計算コストを大幅に削減できます。

:math:`{\cal H}_0`\ が保存量 :math:`Q`\ を持つ場合、Green 関数はそれについて
ブロック対角になり、感受率も同じラベルでブロックに分かれます。ここで相関関数の
各辺のラベルとは、その双一次形式が与える :math:`Q`\ の変化 :math:`\Delta Q`\ の
ことで、左辺は :math:`A^{\dagger}`\ と読み :math:`A`\ のラベルを付けます。この規約
のもとで選択則は :math:`\Delta Q_{\rm L} = \Delta Q_{\rm R}`\ となります。
スピンが典型例で、
:math:`[{\cal H}_0, S_z]=0`\ なら密度（ :math:`\Delta S_z=0`\ ）とスピン反転
（ :math:`\Delta S_z=\pm 1`\ ）のチャネルが厳密に分離されます。これが縦型（ring）と
横型（ladder）を別々に解ける根拠です。ブロックごとに解くことは近似ではなく、
全体を解いた結果と一致します。

ただしこの分解が成り立つには、 **相互作用も同じ量を保存している** 必要があります。
:math:`X^{(0)}`\ がブロック対角でも :math:`\hat{W}`\ がそうでなければ、積
:math:`\hat{X}^{(0)}\hat{W}`\ もその逆行列もブロック対角ではなく、相互作用の
2次以降でセクターが混ざります。スピン軌道相互作用は :math:`S_z`\ を保存させないため
:math:`X^{(0)}`\ 自体に混合を持ち込み、 ``PairLift`` は :math:`{\cal H}_0`\ が
:math:`S_z`\ を保存していても相互作用の側で保存を破ります。

このため H-wave はブロック構造を仮定せず、相互作用行列と裸の感受率の **両方** の
接続性から数値的に検出します。

1. 全k点にわたって相互作用ハミルトニアンの絶対値を、また運動量と振動数にわたって
   裸の感受率の絶対値を合計し、接続パターン行列を得る。
2. 非ゼロの非対角要素（閾値: :math:`10^{-12}`\ ）から隣接グラフを構築する。
3. ラベル伝播（union-findアルゴリズム）により連結成分を求める。

行列が :math:`m`\ 個のブロック（サイズ :math:`n_1, n_2, \ldots, n_m`\ ）に
分解される場合、RPA方程式の計算コストは
:math:`O(N^3)`\ から
:math:`O(n_1^3 + n_2^3 + \cdots + n_m^3)`\ に削減されます。
ここで :math:`N = n_1 + n_2 + \cdots + n_m`\ です。

この最適化は自動的に適用され、ユーザーに対して透過的です。


横感受率（はしごダイアグラム）
*******************************************

標準的な（リングダイアグラム）RPA感受率に加えて、
H-waveは横感受率 :math:`\chi_{+-}(\mathbf{q})`\ を計算できます。
これはスピン反転相関を記述します。 :math:`S^{\pm}`\ の記法は曖昧さを伴うため、
双一次形式そのもので書くと、格納される配列は

.. math::

   \chi_{+-,\alpha\gamma;\beta\delta}(\mathbf{q}) =
   \Big\langle \big(c^{\dagger}_{\gamma\downarrow}
   c^{\mathstrut}_{\alpha\uparrow}\big)(-\mathbf{q}) \;;\;
   \big(c^{\dagger}_{\beta\uparrow}
   c^{\mathstrut}_{\delta\downarrow}\big)(\mathbf{q}) \Big\rangle

であり、添字対の規約は縦感受率と同じです。

.. _rpa_which_array:

どの配列がどのチャネルを保持するか
------------------------------------------

:math:`\chi_{+-}`\ は専用の配列 ``chiq_pm``\ に出力され、これが作られるのは
``calc_type = "ring+ladder"``\ （``calc_scheme = "general"``\ が必要）の
ときだけです。配列 ``chiq``\ は常に縦方向（ring）の結果です。

そもそも ``chiq``\ がそのようなスロットを持つかどうかはスキームによります。
``calc_scheme = "general"``\ では各々がスピンを担う対を添字とするため、対が
スピン非対角となるスロットを常に持ちます。 ``"reduced"``\ は密度対の成分のみを
格納するため、そのようなスロットはそもそも存在しません。

**これらのスロットは横感受率ではなく、バブルがスピンフリーあるいはスピン
対角なバブルの膨張で得られる限り計算されません。恒等的にゼロです。** 膨張は
同スピンのスロットしか作らないためです。これは ``enable_spin_orbital``\ を
使わないすべての計算に加えて、 ``enable_spin_orbital``\ を使っていても一体
ハミルトニアンがスピンフリーあるいはスピン対角と判定される計算にも当てはまり
ます。決めているのは設定ではなく、検出されたスピン構造です。感受率配列の
ゼロは計算結果として読めてしまいますが、ここでは不在を意味します。2軌道の
オンサイト模型では、縦方向のスロットが :math:`\max|\chi| = 1.53`\ に達する
一方、対がスピン非対角なスロットはすべて厳密に :math:`0`\ です。

例外は ``calc_scheme = "general"``\ での真にスピンフルな計算です。
``enable_spin_orbital``\ を使い、かつハミルトニアンが実際にスピンを混ぜる場合、
バブルは一般化軌道添字の上で直接構築されるため、これらのスロットは計算され、
一般にゼロではありません。その場合でも、横感受率は ``chiq_pm``\ が保持する
ものであって、これらのスロットが保持するものではありません。

ラダーを加えても縦方向の答えは変わりません。同一の入力に対して ``"ring"``\ と
``"ring+ladder"`` の ``chiq``\ はビット単位で一致し、``"ring+ladder"``\ が
加えるのは別配列の ``chiq_pm``\ です。したがって ``calc_type``\ の選択は、
横方向チャネルを計算するかどうかの選択であって、縦方向についての近似の
選択ではありません。

.. note::

   SU(2) 対称な常磁性状態では、横方向の応答は対称性から縦方向の結果に
   帰着します（ :math:`\chi_{+-} = 2\chi_{zz}`\ ）。ここで :math:`\chi_{zz}`
   は縦方向の結果から

   .. math::

      \chi_{zz} = \frac{1}{4}\sum_{\sigma\sigma'}\sigma\sigma'\,
      \chi^{\sigma\sigma'} ,\qquad \sigma,\sigma' = \pm 1

   として作ります（ :math:`\chi^{\sigma\sigma'}`\ は対がスピン対角な
   スロット）。空であるスピン非対角スロットから読むのではありません。
   対称性が成り立つ範囲では、これは ``chiq_pm``\ を浮動小数点の丸め誤差の範囲で
   再現するため（下記の入力での実測：要素ごとの **絶対** 差の最大が
   :math:`2.2\times 10^{-16}`\ 、これに対し :math:`\chi_{+-}`\ のピークは
   :math:`1.4`\ ）、この目的にラダーは不要です。

   SU(2) が破れた途端にこの関係は成り立たず、しかも小さなパラメータで
   制御されません。ずれは見積もって許容できる補正ではありません。
   1軌道正方格子（ :math:`U = 4` 、 :math:`4\times 4`\ 、
   :math:`N_{\rm mat} = 64`\ ）に ``coeff_extern = 0.35``\ ---2つのスピンの
   間の分裂は :math:`0.7`\ ---を加えると、静的な推定値は ``chiq_pm``\ から
   ピーク波数において :math:`T = 0.5` で :math:`+2\%` 、 :math:`T = 1.0`\ で
   :math:`-4\%` 、 :math:`T = 0.3` で :math:`-25\%` 、 :math:`T = 0.2`\ で
   :math:`-37\%`\ ずれます。符号は両方向に振れ、温度を下げるほど大きくなり、
   :math:`{\bf q} = 0`\ かつ :math:`T = 0.2`\ では 2 倍以上外れます。さらに
   :math:`T = 0.3`\ では推定値の最大値が真の最大値と異なる波数に現れるため、
   ピーク位置さえ信用できません。磁場、磁気秩序変数、スピン軌道相互作用が
   ある場合、横方向チャネルは計算する必要があります。

.. admonition:: 上記の数値の再現方法
   :class: note

   いずれの計算も ``calc_scheme = "general"`` 、 ``mu = 0.0``\ 、
   ``CellShape = [4, 4, 1]`` 、 ``SubShape = [1, 1, 1]`` 、 ``Nmat = 64``
   を用い、他は同一の入力に対して ``calc_type = "ring"``\ と
   ``calc_type = "ring+ladder"``\ を比較します。

   2軌道の数値は ``tests/rpa/input_2orb`` （ ``geom.dat``\ 、
   ``transfer.dat`` 、 ``coulombintra.dat``\ ）を :math:`T = 0.5`\ で用います。
   ``CoulombIntra``\ のみとしてください。オフサイトの ``CoulombInter``\ は
   ``calc_type = "ring+ladder"``\ では拒否されます。磁場下の数値は
   ``tests/rpa/input`` （ ``geom.dat`` 、 ``transfer.dat``\ 、
   ``coulombintra.dat`` 、 ``extern.dat``\ ）を ``coeff_extern = 0.35``
   で用います。

   静的な応答は松原添字 :math:`N_{\rm mat}/2 = 32`\ で読みます。添字
   :math:`0`\ ではありません。周波数軸は
   :math:`\omega_n = (2n - N_{\rm mat})\pi/\beta`\ であり、添字 :math:`0`\ は
   最も負の周波数だからです。百分率は :math:`\chi_{+-}`\ 自身がピークを取る
   波数で評価した
   :math:`(\chi^{\rm inferred} - \chi_{+-})/\chi_{+-}`\ です。

横方向の裸感受率は

.. math::

   X^{(0)}_{+-,\alpha\gamma;\beta\delta}(\mathbf{q}, i\omega_m)
   = -\frac{T}{N_L} \sum_{\mathbf{k},n}
     G^{(0)}_{\alpha\beta,\uparrow}(\mathbf{k}+\mathbf{q}, i\omega_m + i\epsilon_{n})\,
     G^{(0)}_{\delta\gamma,\downarrow}(\mathbf{k}, i\epsilon_{n})

横方向の頂点 :math:`W_{+-}`\ は、
縦チャネルのハートリー（フォック交換）頂点の交差から得られます:

.. math::

   W_{+-} = W^{\rm spin-flip} - \left[W^{\rm cross-spin}\right]^{\rm crossed}

このバーテックスの組み立てに入る相互作用テンソルは **ハミルトニアン規約** の
ままです。ring の求解が適用する添字対の変換（前節の注記を参照）をここに適用しては
いけません。上の crossing がテンソル自身を組み替えるためです。厳密な置換は
``_assemble_transverse_vertex``\ の実装が定義します。

横方向の頂点は、相互作用テンソルの異スピンブロックとスピン反転ブロックのみ
から構成されます。同スピンブロックは寄与しません。同スピン相互作用は横方向
ループの上向き伝播関数と下向き伝播関数を結べないため、自己エネルギーは生じ
ますが頂点は生じません。

軌道ペアは2つの宣言の平均で対称化されます。相互作用ファイルでは同一の演算子
を2通りに書けるためです （ :math:`n_a n_b = n_b n_a`\ 、Exchange では
:math:`X_{ab} = X_{ba}`\ ）。これは UHFk が用いている規約と同じです。平均を取る
相手は相互作用型ごと（同値には、その型が占めるスロット族ごと）に異なります。密度-密度型と Exchange は単純な転置との平均、
PairHop は共役転置との平均です。PairHop の2つの宣言は同一係数ではなくエルミート
共役の組 （ :math:`P_{ba} = P_{ab}^{*}`\ ）だからです。したがって複素エルミート
閉じの Exchange の物理的結合 :math:`(J_{01} + J_{10})/2`\ は実数になり、複素
エルミート閉じの PairHop は完全な複素値を保ちます。

:math:`W_{+-}`\ は軌道ペアを添字とする行列であり、各相互作用は自分が占める
スロットにのみ寄与します。頂点が使うのは相互作用テンソルの **異スピンブロック**
（:math:`\uparrow\uparrow\downarrow\downarrow`\ と
:math:`\downarrow\downarrow\uparrow\uparrow`\ ）と **スピン反転ブロック**
（:math:`\uparrow\downarrow\uparrow\downarrow`\ と
:math:`\downarrow\uparrow\downarrow\uparrow`\ ）の2つだけです。
同スピンブロックは寄与しません。同スピン相互作用は横型ループの上向き伝播関数と
下向き伝播関数を結べないため、自己エネルギーは生じても頂点は生じないからです。

オンサイト相互作用が与える成分は以下のとおりです。

.. list-table::
   :header-rows: 1
   :widths: 22 30 26

   * - 相互作用
     - 占めるブロック
     - :math:`W_{+-}`\ の成分
   * - ``CoulombIntra`` :math:`U`
     - 異スピン
     - :math:`-U`
   * - ``CoulombInter`` :math:`V`
     - 異スピンと同スピン（後者は不使用）
     - :math:`-V`
   * - ``Ising`` :math:`I`
     - 同上
     - :math:`+I`
   * - ``PairHop`` :math:`J`
     - 異スピン
     - :math:`-J`
   * - ``Exchange`` :math:`J`
     - スピン反転
     - :math:`-(J + J^{\rm T})/2`
   * - ``Hund`` :math:`J`
     - 同スピンのみ
     - :math:`0`
   * - ``PairLift`` :math:`J`
     - 二重スピン反転（ :math:`W_{+-}`\ では不使用 ）
     - :math:`0`

``Hund`` と ``PairLift``\ がゼロになるのは値が小さいからではなく、
頂点が使う2つのブロックの外に成分があるためです。

``Ising``\ の規格化について：wannier90 形式の k 空間ソルバーはすべてこの規格化で
Ising ファイルを読みます。UHFk の因子 1/4 の不一致はバージョン 2.0 で
解消済みです。別系統の実空間 UHFr リーダーは :math:`S^z`\ 規約を保持します。

.. note::

   これらの値は、横方向チャネルを厳密対角化と照合する前に公開されていた値とは
   異なります。従来の記載は7種のうち4種が誤りで、1種が欠落していました。
   ``CoulombInter`` 、 ``Hund`` 、 ``Ising`` 、 ``Exchange``\ を含む計算で得られた
   横方向帯磁率は再計算してください。影響を受けるのは ``chiq_pm``\ のみで、
   ``chiq``\ 、自己エネルギー、Eliashberg 頂点には波及しません。

``calc_type = "ring+ladder"``\ では、縦方向の計算より前に、 **組み立てられた
横方向頂点が** （副格子で折り畳んだ格子上で） :math:`q` **に依存しないこと**
を相対許容値 :math:`10^{-10}`\ で検証し、満たさない入力を拒否します。横方向の
ペア :math:`c^\dagger_{i a \uparrow} c_{j b \downarrow}`\ はオフサイト項に
対して非局所となり、その頂点は :math:`q`\ のみの関数では表現できないためです。
実際にはオフサイトの ``CoulombInter`` 、 ``Ising`` 、 ``Exchange``\ が拒否され、
オフサイトの ``Hund`` と ``PairLift``\ は横方向頂点が消えるため受理されます。
なお、互いに打ち消し合う・値が食い違う宣言の組は、より早い読み込み時に
拒否されます（バージョン 2.0 以降: 宣言ファイルはエルミート共役で閉じている必要が
あります）。また ``SubShape``\ による折り畳みでスーパーセル内に収まる
サイト間ペアはセル内軌道ペアとなり、表現可能として受理されます。縦方向（ ``ring``\ ）チャネルは影響を受けません。

.. warning::

   オフサイトの ``PairHop``\ は相互作用の読み込み時にこの検査より前に暗黙に
   破棄されるため、拒否も反映もされません。RPA 計算でオフサイトの
   ``PairHop``\ に依存しないでください。また、対角（同一軌道）の PairHop は
   密度項 :math:`2P\, n_\uparrow n_\downarrow`\ を意味しますが、読み込み時の
   係数は :math:`2P`\ ではなく :math:`P`\ として（縦・横チャネルで一貫して）
   扱われます。この縮退エントリの検証は別途追跡されています。

スピンを本当に混合する（genuinely spin-mixing）スピン軌道ハミルトニアンを
除く全てのモード（スピンなし、スピン対角、および :math:`S_z`\ を保存する
スピン軌道の場合）では、横方向のRPA感受率は上記の自己完結した頂点
:math:`W_{+-}`\ から

.. math::

   \hat{X}_{+-}(\mathbf{q})
   = \left[\hat{I} + \hat{X}^{(0)}_{+-}(\mathbf{q})\, \hat{W}_{+-}\right]^{-1}
     \hat{X}^{(0)}_{+-}(\mathbf{q})

として得られます。本当にスピンを混合するスピン軌道の場合は扱いが異なります。
下の注を参照してください。

横チャネルの計算を有効にするには、入力TOMLファイルで
``calc_type = "ring+ladder"``\ と設定します。
これには ``general``\ 計算スキームが必要です（自動選択されます）。

.. note::

   スピン軌道モードでハミルトニアンが本当にスピンを混合する場合
   （スピン軌道相互作用など）、上記の式は適用されません。代わりに、
   縦チャネルのスピンフル求解と同じ反対称化された頂点（direct 相互作用
   ＋オンサイトの exchange クロス項。設定リファレンスの
   ``spinful_vertex_exchange``\ を参照）で全スピン軌道空間をドレスし、
   そのうえで :math:`\chi_{+-}(\mathbf{q})`\ をドレスされたテンソルの
   スピン反転ブロックとして抽出します。したがってスピン混合のクロス項が
   含まれるのは、この反対称化された頂点のレベルまでであり、厳密では
   ありません。*オフサイト* 相互作用の exchange クロス項は独立な2つの
   フェルミオン運動量に依存するため単一の頂点 :math:`W(\mathbf{q})`\ では
   表現できず、実装されていません。相互作用の種類のうち、オフサイト宣言が
   このドレッシングに到達しうるのは ``Hund``\ と ``PairLift``\ のみです
   （他の種類は上述の通りオフサイトでは拒否されるか破棄されます）。その
   ため、これらの横方向への寄与は部分的にしか含まれず、該当する宣言が
   ある場合には警告が出力されます。この抽出は反対称化された縦チャネルの
   求解を再利用するため、``spinful_vertex_exchange = false``
   （ring のみとの互換性のためのスイッチ。設定リファレンスを参照）と
   ``calc_type = "ring+ladder"``\ 、かつ本当にスピンを混合する
   ハミルトニアンという組み合わせは拒否されます。


スピン軌道モード
*****************************

H-wave はスピン軌道モードをサポートしています。このモードでは入力ファイルの
スピン・軌道インデックスがブロック分離ではなくインターリーブされた形で与えられます。
ソルバーは読み込み時に内部のスピンブロック順へ並べ替えており、保存される感受率も
その順序です（ ``index_convention = "spin_block"``\ として記録されます）。

通常モードでは、複合インデックスは :math:`i = s \cdot n_{\rm orb} + a`
（スピンブロック優先）であり、
:math:`s = 0, 1`\ はスピンインデックス、
:math:`a = 0, \ldots, n_{\rm orb}-1`\ は軌道インデックスです。
スピン軌道モードでは、インデックスは :math:`i = 2a + s`
（インターリーブ）であり、スピン軌道相互作用を自然に扱えます。

スピン軌道モードは入力TOMLファイルで ``enable_spin_orbital = true``
と設定して有効化します。このモードでは:

- ハミルトニアンはスピン保存を仮定せず、
  完全な :math:`2n_{\rm orb} \times 2n_{\rm orb}`\ 空間で構成されます。
- 全ての相互作用型（ ``CoulombIntra`` 、 ``CoulombInter`` 、 ``Hund`` 、 ``Exchange``\ 、
  ``Ising`` 、 ``PairLift`` 、 ``PairHop``\ ）がサポートされます。
- 可能な場合、ブロック分解が自動的に適用されます。
- ``reduced`` と ``general``\ の両計算スキームがスピン軌道系で利用可能です。
- 幾何情報ファイル（``geom.dat``\ ）の ``Norbit``\ はスピン軌道の総数
  （= 物理軌道数 × 2 = Wannier90 の ``num_wann``\ ）を表し、UHFk と同じ規約です。

.. note::

   **移行上の注意（RPA）：** スピン軌道入力の幾何 ``Norbit``\ はスピン軌道数になりました。
   既存の RPA スピン軌道計算の ``geom.dat`` の ``Norbit``\ は2倍にしてください。


.. [1] `K. Yoshimi, T. Kato, H. Maebashi, J. Phys. Soc. Jpn. 78, 104002 (2009). <https://journals.jps.jp/doi/10.1143/JPSJ.78.104002>`_

.. highlight:: none

uhfk_to_mvmc.py — UHFk → mVMC PairProduct ブリッジ
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

``tools/uhfk_to_mvmc.py`` は H-wave UHFk の SCF 結果を、mVMC の
``InOrbital`` / ``InOrbitalAntiParallel`` / ``InOrbitalGeneral`` 初期
波動関数ファイル (``zqp_orbital_uhfk.dat``) に変換する。これにより
mVMC の PairProduct 状態を H-wave UHF の Slater 行列で初期化できる。

変換の規約と構成法は :doc:`../../algorithm/uhfk_to_mvmc` を参照。
本節では利用方法を述べる。

対応範囲
^^^^^^^^

軌道と格子
""""""""""

- 単軌道 (``norb_orig = 1``)。
- 副格子フォールドは ``SubShape`` が各方向で ``CellShape`` を割り切る
  任意の値に対応する。``SubShape`` を省略した場合は ``CellShape`` に
  フォールバックする (``uhfk.py`` の ``_init_lattice`` と一致)。

変換経路
""""""""

SCF 解の性質に応じて 3 つの経路がある。経路は CLI が自動的に選択する
(:ref:`後述のとおり <uhfk_to_mvmc_dispatch>`)。

AntiParallel 経路
   Sz 固定 (``2Sz = 0``) の解を扱う。``orbitalidx.def`` の 3 列または
   4 列形式を消費する。

General 経路
   Sz 固定で ``2Sz != 0`` の解、および Sz 非固定で mixed block を
   含まない解を扱う。``orbitalidx_general.def`` の 6 列形式を消費する。
   利用者は StdFace で ``orbitalidx_general.def`` を生成する必要がある
   (``2Sz != 0`` の場合は ``stan.in`` の ``2Sz`` を非 0 に設定する)。

General-SOC 経路
   Sz 非固定で mixed block を含む解、すなわちスピン軌道結合を有効に
   した解 (``enable_spin_orbital = true``) を扱う。Zeeman、Rashba、
   Dresselhaus、および一般の :math:`\sigma_x` / :math:`\sigma_y` 型
   1 体結合をサポートする。Sz を保存しない ``F[up_i, down_j]`` /
   ``F[down_i, up_j]`` を出力する。

スピン不均衡な Slater 状態については、canonical な
``(k, partner(k))`` ブロックで same-spin excess pair を発行することに
より、同スピン pair 成分 (``F[up, up]``、``F[down, down]``) をサポート
する。

境界条件
""""""""

周期境界条件と反周期境界条件のいずれも決定論的にサポートする。

スピン軌道結合と反周期境界条件と ``SubShape > [1, 1, 1]`` が同時に
指定された場合に限り、対応する格子形状と方向の組み合わせが以下に
限定される。

==============================  ==========================================
格子形状                        対応する反周期方向
==============================  ==========================================
``CellShape = [6, 4, 1]``       単方向 (x / y / z のいずれか 1 方向)
``SubShape  = [2, 2, 1]``
``CellShape = [4, 4, 4]``       多方向 (xy / xz / yz / xyz)
``SubShape  = [2, 2, 2]``
==============================  ==========================================

次のいずれかに該当する場合はこの制限を受けず、無条件にサポートされる。

- スピン軌道結合を用いない
- ``SubShape = [1, 1, 1]``
- スピン軌道結合を用いるが全方向が周期境界

上記のいずれにも当てはまらず、かつ表の組み合わせでもない指定は、
dispatch 前に拒否される。エラーメッセージは対応する組み合わせを
列挙する。新しい組み合わせを追加するには、fixture と gate 検証を
先に用意する必要がある。

その他
""""""

- T = 0 への Slater 投影を行う。Fermi 準位付近に分数占有が残る有限
  温度 SCF は拒否し、より小さい ``T`` での再計算を促す。
- ``params[idx]`` に微小な一様乱数を加算し、rank-deficient な
  :math:`F` に対する mVMC の Pfaffian Slater 評価の特異性を回避する。
  既定の振幅は ``1e-8`` で、``--epsilon-noise`` で変更できる。
  ``ComplexType 1`` では実部と虚部の双方に、``ComplexType 0`` では
  実部のみに加算する。mVMC 本体の ComplexUHF と同じ手法である
  (``mVMC-1.4.0/src/ComplexUHF/output.c:274``)。

制限
^^^^

- 多軌道 (``norb_orig > 1``) は対象外である。
- 2 体の Sz 非保存相互作用 (spin-flip Coulomb、Hund coupling、
  pair hopping) は対象外である。対応する 2 体項は CoulombIntra
  (オンサイト :math:`U`) のみである。
- 反周期境界の twist は :math:`\theta` の各成分が :math:`\{0, \pi\}`
  の場合のみ検証されている。一般の複素 twist は対象外である。
- 上表にない、スピン軌道結合と反周期境界と ``SubShape > [1, 1, 1]``
  の組み合わせは拒否される。

ワークフロー
^^^^^^^^^^^^

1. StdFace で mVMC 入力一式 (``orbitalidx.def`` 等) を生成する。
   ``CalcMode = 2`` と、副格子並進対称性に応じた ``Lsub`` を設定する。
   反周期境界条件は ``phase0 = 180.0`` を指定する。
2. H-wave UHFk SCF を実行し、``eigen.npz`` に加えて ``occupation.npz``
   を ``[file.output]`` で要求する::

       [file.output]
         path_to_output = "output"
         eigen          = "eigen.npz"
         green          = "green.npz"
         occupation     = "occupation.npz"
         onebodyg       = "greenone.dat"

3. ブリッジを実行する::

       python tools/uhfk_to_mvmc.py \
           --input        input.toml \
           --eigen        output/eigen.npz \
           --occupation   output/occupation.npz \
           --geometry     geometry_uhf.dat \
           --orbitalidx   mvmc_inputs/orbitalidx.def \
           --output       mvmc_inputs/zqp_orbital_uhfk.dat \
           --check-density \
           --onebodyg-uhf output/greenone.dat

4. mVMC の ``namelist.def`` に出力ファイルを登録する。経路に応じて
   ``InOrbital``、``InOrbitalAntiParallel``、または
   ``InOrbitalGeneral`` のいずれかのキーを用いる::

       InOrbitalGeneral zqp_orbital_uhfk.dat

   これにより mVMC が PairProduct パラメータを本ファイルで初期化する。

.. _uhfk_to_mvmc_dispatch:

経路の選択
^^^^^^^^^^

CLI は ``--orbitalidx`` と ``input.toml`` を parse し、次の 3 つの値の
組で経路を選択する。

``is_antiparallel_metadata``
   ``orbitalidx.def`` のメタデータが AntiParallel の要件を満たすか。
``orbitalidx_format``
   ``orbitalidx.def`` の列形式 (``antiparallel`` または ``general``)。
``is_soc_mode``
   ``input.toml`` の ``enable_spin_orbital``。

対応は次のとおりである。

==========  ==================  ==========  ====================
metadata    format              SOC         経路
==========  ==================  ==========  ====================
``True``    ``antiparallel``    ``False``   AntiParallel
``True``    ``general``         ``False``   forced-General
``False``   ``general``         ``False``   General
``False``   ``antiparallel``    ``False``   拒否
任意        ``general``         ``True``    General-SOC
任意        ``antiparallel``    ``True``    拒否 (6 列が必須)
==========  ==================  ==========  ====================

列は上から ``is_antiparallel_metadata``、``orbitalidx_format``、
``is_soc_mode`` に対応する。

``is_antiparallel_metadata`` は次の条件を **すべて** 満たしたときに
のみ ``True`` となる。

- ``input.toml`` の ``2Sz`` が明示的に 0 である
- ``N_up == N_down`` である
- ``column_spin`` の値が ``{0, 1}`` に含まれる
- ``column_mu_group`` の unique 数が 2 である
- ``column_spin`` と ``column_mu_group`` が全単射である

forced-General 経路は、AntiParallel のメタデータを持つ入力に対して
6 列形式が与えられた場合の分岐である。占有集合が
``(k_row, local_band)`` の pair-closure を満たしていれば、
``F[up, down]`` は AntiParallel 経路の :math:`F` を :math:`10^{-12}`
で再現する。満たしていない場合 (spin canting に由来する up-up excess
など) は警告を出力するが、出力自体は正しい ``InOrbitalGeneral`` 状態
となる。ただしその状態は AntiParallel 経路では表現できない。

``(False, antiparallel)`` が拒否されるのは、AntiParallel の要件を
満たさない解を 3 列または 4 列形式で表現できないためである。StdFace を
``2Sz`` 非 0 または Zeeman 駆動の Sz 非固定設定で実行し、
``orbitalidx_general.def`` を生成し直すこと。

BoundaryCondition の契約
^^^^^^^^^^^^^^^^^^^^^^^^

ブリッジは ``BoundaryCondition`` を H-wave の共有 helper
``normalize_boundary_condition`` に委譲して正規化する。受理する形式は
次のとおりで、大文字小文字を区別せず、前後の空白は除去される。

- 周期境界: ``"p"``、``"periodic"``
- 反周期境界: ``"ap"``、``"antiperiodic"``

これ以外の文字列は dispatch 前に ``ValueError`` となる。未知の値を
暗黙に解釈する fallback は存在しない。キーを省略した場合は全方向
周期境界を既定とする。

また ``eigen.npz`` の ``twist_offset`` を ``input.toml`` の
``BoundaryCondition`` と照合し、入力と eigen の組み合わせが古い場合を
排除する。

trans.def の出力
^^^^^^^^^^^^^^^^

スピン軌道結合を用いる場合、mVMC の ``vmcdry.out`` が生成する
``trans.def`` はスピン対角成分のみを保持し、:math:`s \neq t` の項を
黙って落とす。この不足を埋めるため、ブリッジが H-wave の
``Transfer.dat`` を読み、``(i, s, j, t, re, im)`` 形式の ``trans.def``
を出力する。写像規約は :doc:`../../algorithm/uhfk_to_mvmc` を参照。

必要な CLI フラグは ``--transfer`` (``Transfer.dat`` の読み込み) と
``--emit-trans`` (``trans.def`` の出力) である。

スピン軌道結合と ``SubShape > [1, 1, 1]`` を併用する場合は、さらに
``--emit-orbitalidx`` が必要である。StdFace が生成する
``orbitalidxgen.def`` は、折り畳み格子の下では Sz を保存しない pair
クラスを表現できず、クラスを過剰に統合してしまうためである。この
フラグを指定すると、ブリッジがクラスを統合しない
``orbitalidx_general.def`` を出力する。

クラス一致性チェック
^^^^^^^^^^^^^^^^^^^^

General 経路では ``aggregate_general_orbital_params`` が、平均化する
**前** に、``orbitalidx_general.def`` の各クラスに割り当てられた符号付き
:math:`F` 成分が ``class_consistency_tol`` (既定 :math:`10^{-8}`)
以内で一致しているかを検査する。一致しない場合は、該当する索引と
観測された最大残差とともに ``ClassInconsistencyError`` を送出する。

この検査は、StdFace が生成したクラスが仮定する対称性を Slater 状態が
守っていない場合 (対称なハミルトニアンに対して自発的対称性の破れた
UHF 基底状態が得られた場合など) に、異なる値が黙って平均化される
ことを防ぐ。

密度行列チェック
^^^^^^^^^^^^^^^^

``--check-density`` を指定すると、変換した波動関数から 1 体密度行列を
構築し、H-wave の結果と element-wise に比較する。許容差は
:math:`10^{-10}` で、超過は致命的エラーとする。不一致は H-wave の
境界条件処理、ブリッジ、または geometry の前提のいずれかに誤りが
あることを示す。

比較の基準は組み合わせによって異なる。

- 通常は :math:`(k, -k)` pair 構築から得た密度行列を、H-wave の
  物理基底の ``greenone.dat`` と比較する。
- スピン軌道結合と ``SubShape > [1, 1, 1]`` を併用する場合は、
  ``compare_against_green_sublattice`` に切り替える。H-wave の
  ``green_sublattice`` を ``gauge_lift`` で物理基底に持ち上げ、
  実際に出力する Slater 行列から作った :math:`\overline{A} A^{T}` と
  element-wise に照合する。この組み合わせでは ``greenone.dat`` の
  折り畳み経路を基準として用いない。

後者は出力そのものを検査対象とするため、出力経路の退行が検出されずに
通過することはない。

検証
^^^^

変換結果は 7 つの gate で検証する。各 gate は独立した照合であり、
すべてが通ったときにのみ変換が妥当と判断する。

G0-writer-check
   rank-lift ノイズを無効化した経路の :math:`F` が、集約済みの
   ``(mapping, params)`` と一致すること。書き出し処理そのものの健全性を
   確認する。許容差 :math:`10^{-10}`。
G1
   出力する Slater 行列から作った密度が、``gauge_lift`` で物理基底に
   持ち上げた ``green_sublattice`` と一致すること。
   許容差 :math:`10^{-10}`。
G2a-emitted-F
   出力した :math:`F` を skew-SVD で射影した密度が、ComplexUHF の
   1 体 Green 関数と一致すること。許容差 :math:`10^{-6}`。
G2a-in-memory-A
   メモリ上の Slater 行列から作った密度が ComplexUHF と一致すること。
   許容差 :math:`10^{-6}`。
G2b
   ``gauge_lift`` で持ち上げた ``green_sublattice`` が ComplexUHF と
   一致すること。許容差 :math:`10^{-6}`。
G3
   mVMC の :math:`\langle H \rangle` と H-wave の ``Energy_Total`` の
   相対差が 1 % 以下であること。
G4
   合成要素が現在の SCF 上で維持され、mutation を加えた場合に
   閾値 :math:`T_M = \max(10^{-5},\, 0.10 |G_\mathrm{base}|)` を
   超えて差が現れること。変換が格子のトポロジーに依存していることを
   確認する。

G2 系の 3 つは、一致するだけでなく **収縮** することを要求する。
ComplexUHF は H-wave の収束密度から seed するため、seed が既に許容差を
満たしていると、seed をそのまま返すソルバでも通過してしまう。そこで
各 G2 は、seed 時点で基準から少なくとも許容差の 10 倍離れていることと、
収束後に許容差の内側に入ることの両方を要求し、双方の値と収縮率を記録
する。また比較の前に非有限値を拒否する。NaN との比較は両方の境界判定を
偽にするため、そのままでは素通りするためである。

.. note::

   検証に用いる H-wave の実行は ``flag_fock = true`` である必要がある。
   比較対象の ``ComplexUHF`` はオンサイト交換項をコンパイル時に固定して
   おり (``src/ComplexUHF/include/Def.h`` の ``#define Fock 1``)、
   実行時スイッチを持たない。``flag_fock = false`` で得た解は異なる
   平均場汎関数の停留点であり、ComplexUHF のいかなる不動点とも一致
   しない。

gate は次のコマンドで実行する::

    bash tests/validation/uhfk_mvmc_pairproduct/run.sh <ケース名>

各 gate は ``<GATE 名> PASS mode=...`` の形式で始まる行を出力する。

終了コード
^^^^^^^^^^

- ``0`` — 成功
- ``2`` — fail-fast ガードが入力を拒否した。対象外のモード、Sz 非固定
  の SCF、有限温度における分数占有の残留、``orbitalidx.def`` や
  geometry や境界条件の不整合などが該当する。失敗した検査の名称を
  stderr に出力する。
- ``3`` — ``--check-density`` が許容差を超える不一致を検出した

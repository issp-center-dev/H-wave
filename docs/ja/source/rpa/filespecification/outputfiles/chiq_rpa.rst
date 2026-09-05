.. highlight:: none

.. _Subsec:chiq_rpa:

chiq, chi0q
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

感受率行列および既約感受率行列の計算結果を Numpy zip 形式で出力します。
ファイル名は、環境設定ファイルの中の\ ``file.output``\ セクションでキーワード\ ``chiq``\ および\ ``chi0q``\ を用いて指定された文字列(以下\ ``chiq_str``\ と記述)を用いて、\ ``chiq_str.npz``\ という名前で出力されます。

ファイルの内容は、以下のキーにバインドされる複数の配列データからなります。

- ``chiq``\ または\ ``chi0q``:

  感受率行列または既約感受率行列のデータ。データ形式は以下の節で説明します。

- ``momentum_convention``:

  フーリエ符号の来歴マーカー\ ``"e_plus_ikR"``\ （バージョン 2.0 以降）。運動量
  ラベルは\ :math:`M(k) = \sum_R M(R) e^{+ikR}`\ に従います。異なる値を
  記録したファイルはローダが拒否します。このフィールド導入前のファイル
  は、内容が\ :math:`q \to -q`\ に対して要素毎に偶（両規約が一致する
  場合）のときに限り受理されます。

- ``freq_index``:

  出力する松原振動数の値または範囲はパラメータ\ ``matsubara_frequency``\ で指定されます。出力データの配列インデックスと、実際の松原振動数のラベルの対応付けを\ ``freq_index``\ に格納します。

- ``wavevector_unit``\ および\ ``wavevector_index``:

  波数ベクトルの情報を格納します。詳細は\ :ref:`UHFk の出力ファイル<Subsec:eigen_uhfk.dat>`\ を参照してください。

- ``index_convention``:

  スピン軌道の軸が spin*norb+orb（スピンブロック）順で並んでいることを示す文字列マーカー(常に\ ``"spin_block"``)です。これは UHFk のインターリーブ(2*orb+spin)出力とは異なります。このマーカーを持たない\ ``chi0q``\ ファイル(規約修正前に出力されたもの)は、スピン軌道モードで再読み込みする際に拒否されます。

- ``coeff_tail``:

  感受率の計算時に使用された\ ``coeff_tail``\ パラメータ(Matsubara 高周波数裾補正)の値です。裾補正は中程度の\ ``Nmat``\ において\ :math:`\chi_0(\mathbf{q})`\ を\ :math:`O(1)`\ で変えるため、このファイルを読み込む計算(``hwave_sc``\ など)は、自身の\ ``coeff_tail``\ 設定が記録値と異なる場合に警告を出します。事前計算した\ ``chi0q_init``\ ファイルをそのまま出力する場合には入力ファイルの値が引き継がれます。生成元がこの値を記録していない場合(旧バージョン)、および一様格子の裾補正が適用されない IR 松原基底(``matsubara_basis = "ir"``\ 、densified 出力を含む)での FLEX 出力の場合には、キー自体が省略されます。

副格子を指定している場合は、出力されるデータは副格子を単位とした感受率の値です。波数ベクトルおよび軌道のインデックスは副格子に読み替えます。
   
``chi0q``\ の出力ファイルは、計算済み既約感受率データとして\ ``file.input``\ セクションの\ ``chi0q_init``\ に指定して使用できます。


スキームの来歴
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

バージョン 2.0 以降、RPA および FLEX が出力するすべての\ ``chiq``/``chi0q``\ ファイルには、以下の3つの文字列フィールドが 0 次元配列として格納されます（読み込みには\ ``.item()``\ を使用してください）。

- ``calc_scheme``: そのファイルの計算に実際に使用されたスキーム（``reduced`` | ``general``）。
- ``calc_scheme_requested``: 入力で指定された値（``auto`` | ``reduced`` | ``general``）。
- ``scheme_resolution``: どのように決定されたかを示します。閉じた語彙です（変更する場合は出力形式の変更として扱います）: ``explicit``, ``auto:ring_ladder``, ``auto:general_only``, ``auto:no_discarded_content``, ``auto:exact:diagonal_transfer``, ``auto:exact:folded_diagonal``, ``auto:mixed:transfer``, ``auto:mixed:extern``, ``auto:mixed:trans_mod``, ``auto:mixed:green_init``, ``auto:flex_forcing``\ 。

1.0.x で出力されたファイルにはこれらのフィールドはありません。読み込み側でもこれらを必須とはしていません。


chi0q のデータ形式
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

``chi0q``\ のデータ形式は、スピン軌道相互作用および外場の有無と、\ ``mode.calc_scheme``\ パラメータの値によって以下の形をとります。

- スピン非依存 (spin-free):

  ``enable_spin_orbital``\ パラメータが\ ``false``\ の場合、もしくは\ ``true``\ であっても\ :math:`T_{\tilde\alpha\tilde\beta}(k)`\ がスピンについて対角かつ対称な場合で、外場がないとき、既約感受率行列のスピン非依存部分を出力します。

  - ``calc_scheme = general``\ の場合、データ配列は\ ``ndarray(l,q,a,ap,b,bp)``\ で、インデックスは以下のとおりです

    - ``l``: 松原振動数のラベル。インデックスとラベルの対応は前述の\ ``freq_index``\ で与えられます。

    - ``q`` : 波数ベクトルのインデックス\ :math:`[q_x\ q_y\ q_z]`\ を1次元化したインデックスで、\ :math:`q = q_z + N_z \cdot (q_y + N_y \cdot q_x)`\ となります。

    - ``a``, ``ap``, ``b``, ``bp``\ はスピンを含まない軌道インデックス\ :math:`\alpha`, :math:`\alpha^\prime`, :math:`\beta`, :math:`\beta^\prime`\ に対応します。

  - ``calc_scheme = reduced``\ の場合、データ配列は\ ``ndarray(l,q,a,b)``\ です。インデックスの意味は上記と同じです。

- スピン対角 (spin-diagonal):

  ``enable_spin_orbital``\ パラメータが\ ``false``\ で外場がある場合、もしくは\ ``enable_spin_orbital``\ パラメータが\ ``true``\ で\ :math:`T_{\tilde\alpha\tilde\beta}(k)`\ がスピンについて対角な場合に、既約感受率行列の spin up/down 成分を出力します。

  - ``calc_scheme = general``\ の場合、データ配列は\ ``ndarray(s,l,q,a,ap,b,bp)``\ となります。\ ``s = 0 (1)``\ はそれぞれ spin up (down) 成分を表し、それ以外のインデックスは上記と同じです。
    
  - ``calc_scheme = reduced``\ の場合、データ配列は\ ``ndarray(s,l,q,a,b)``\ となります。

- スピン依存 (spinful):

  ``enable_spin_orbital``\ パラメータが\ ``true``\ で Transfer項が一般的な形の場合、一般化軌道をインデックスとする既約感受率行列を出力します。

  - ``calc_scheme = general``\ の場合、データ配列は\ ``ndarray(l,q,a,ap,b,bp)``\ となります。\ ``a``, ``ap``, ``b``, ``bp``\ はスピンを含む一般化軌道インデックス\ :math:`\tilde\alpha`, :math:`\tilde\alpha^\prime`, :math:`\tilde\beta`, :math:`\tilde\beta^\prime`\ に対応します。

  - ``calc_scheme = reduced``\ の場合、データ配列は\ ``ndarray(l,q,a,b)``\ となります。\ ``a``, ``b``\ はスピンを含む一般化軌道インデックス\ :math:`\tilde\alpha`, :math:`\tilde\beta`\ に対応します。

  

chiq のデータ形式
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

``chiq``\ のデータ形式は\ ``calc_scheme``\ パラメータの値によって以下の形をとります。

- ``calc_scheme = general``\ の場合、データ配列は\ ``ndarray(l,q,a,ap,b,bp)``\ となります。\ ``a``, ``ap``, ``b``, ``bp``\ はスピンを含む一般化軌道インデックス\ :math:`\tilde\alpha`, :math:`\tilde\alpha^\prime`, :math:`\tilde\beta`, :math:`\tilde\beta^\prime`\ に対応します。

- ``calc_scheme = reduced``\ の場合、データ配列は\ ``ndarray(l,q,a,b)``\ となります。\ ``a``, ``b``\ はスピンを含む一般化軌道インデックス\ :math:`\tilde\alpha`, :math:`\tilde\beta`\ に対応します。

``chiq``\ は縦方向（ring）の感受率です。\ ``calc_scheme = general``\ では添字の対がスピン非対角となるスロットを持ちます。\ ``reduced``\ は密度対の成分のみを格納するためそのようなスロットを持ちません。これらのスロットは横感受率では **ありません** 。またバブルがスピンフリーあるいはスピン対角なバブルの膨張で得られる限り **計算されません** 。膨張は同スピンのスロットしか作らないため、残りは恒等的にゼロのままです。これは\ ``enable_spin_orbital``\ を使わないすべての計算に加えて、\ ``enable_spin_orbital``\ を使っていても一体ハミルトニアンがスピンフリーあるいはスピン対角と判定される計算にも当てはまります。このゼロを横方向応答の計算結果として読まないでください。（真にスピンフルな\ ``general``\ の計算では、バブルは一般化軌道添字の上で直接構築されるため、これらのスロットは計算され、一般にゼロではありません。その場合でも横感受率は\ ``chiq_pm``\ です。）横方向チャネルの求め方は\ :ref:`アルゴリズム<rpa_which_array>`\ を参照してください。

.. admonition:: ``calc_scheme = "squashed"``\ からの移行（2.0 で削除）
   :class: note

   ``squashed``\ は\ ``reduced``\ と同じ感受率を数倍のコストで計算しており、
   その8軸出力のうちスピン添字の対が異なるスロットは構造的にゼロでした。
   この設定は現在、起動時にエラーで失敗します。\ ``calc_scheme = "reduced"``
   を使用してください。

   旧形式の8軸ファイルを読み込む解析スクリプトは、以下のように変換できます::

       chiq8 = data["chiq"]           # (l, q, s1, s2, a, s3, s4, b)
       l, q, _, _, norb, _, _, _ = chiq8.shape
       nd = 2 * norb
       chiq4 = np.zeros(chiq8.shape[:2] + (nd, nd), dtype=chiq8.dtype)
       for s1 in (0, 1):
           for s3 in (0, 1):
               chiq4[:, :, s1*norb:(s1+1)*norb, s3*norb:(s3+1)*norb] = \
                   chiq8[:, :, s1, s1, :, s3, s3, :]

   対の対角成分（``s1 == s2``\ 、``s3 == s4``\ ）のスロットのみがデータを持ち、
   8軸配列の他のスロットはすべて厳密にゼロです。旧\ ``squashed``\ 実行で
   生成された感受率ファイルは、\ ``reduced``\ の下で\ ``chi0q_init``\ として
   引き続き読み込めます。ただし他の旧ファイルと同様の来歴検証（例えば
   ``momentum_convention``\ マーカー）を経る必要があります。両スキームは
   常に同一のバブル表現を共有していたためです。

``calc_type = ring+ladder``\ の場合、\ ``chiq``\ ファイルには配列\ ``chiq_pm``\ も追加で格納されます。これは横感受率\ :math:`\chi_{+-}(q)`\ を保持します。配列形式は\ ``ndarray(l,q,a,ap,b,bp)``\ で、\ ``a``, ``ap``, ``b``, ``bp``\ はスピン自由度を含 **まない** 軌道インデックス\ :math:`\alpha`, :math:`\gamma`, :math:`\beta`, :math:`\delta`\ です（スピン構造は\ :math:`+-`\ のラベルで既に定まっているため）。したがって、対応する軸が一般化（スピン軌道）インデックスを走る\ ``chiq``\ よりも小さい配列になります。縦方向の\ ``chiq``\ はラダーの有無に影響されません。同一の入力に対して\ ``calc_type = ring``\ の結果とビット単位で一致します。

``longitudinal_bond_channels = true``\ （実験的機能、\ ``calc_type = "ring"``\ 、
:ref:`rpa_longitudinal_bond`\ を参照）の場合、\ ``chiq``\ ファイルにはボンド分解した
縦方向のオブジェクトが\ ``longitudinal_bond_``\ 接頭辞の下に追加で格納されます。
ファイル内の配列\ ``chiq``\ は上書き **されません** 。標準（Hartree のみ）の ring の
結果のままで、交換交差を含めた静的感受率は新しいキーの下にのみ格納されます。
ボンドチャネル数を\ :math:`B`\ （チャネル0はオンサイト\ :math:`R = 0`\ ）、
:math:`n_d = n_{\rm orb}^2`\ 、\ :math:`N_D = B\, n_d`\ として、

- ``longitudinal_bond_chi_s``\ 、\ ``longitudinal_bond_chi_c``\ : ドレスされた静的
  スピン・電荷感受率。複素\ ``ndarray(q, I, J)``\ で、ボンド優先の添字は
  :math:`I = m\, n_d + l_1 n_{\rm orb} + l_2`\ 、
  :math:`J = m'\, n_d + l_3 n_{\rm orb} + l_4`\ です
  （\ ``longitudinal_bond_index_order``\ にこの規則が文字列として記録されます）。
- ``longitudinal_bond_chiq_s_static``\ 、\ ``longitudinal_bond_chiq_c_static``\ :
  その\ :math:`(m = 0, m' = 0)`\ ブロック。軌道添字の複素\ ``ndarray(q, a, b, c, d)``\ で、
  標準の ring の静的スピン・電荷チャネルにオフサイト交換交差を含めたものに相当します。
- ``longitudinal_bond_delta_r``\ （整数\ ``(B, 3)``\ の変位。チャネル0が先頭）、
  ``longitudinal_bond_reverse``\ （整数\ ``(B,)``\ 、\ :math:`-R`\ のチャネル）、
  ``longitudinal_bond_spatial_shape``\ （整数\ ``(3,)``\ 、\ ``CellShape``\ ）、
  ``longitudinal_bond_q_convention``\ （文字列:
  ``q = 2*pi*(n_x/N_x, n_y/N_y, n_z/N_z), C-order flattened``\ ）、
  ``longitudinal_bond_spin_mode``\ （文字列、\ ``spin-free``\ ）、
  ``longitudinal_bond_normalization``\ （文字列:
  ``chi_bar = -(T/N) sum_k G G, per site``\ ）、
  ``longitudinal_bond_types``\ （文字列配列。常に\ ``CoulombInter, Hund, Ising``\ の
  3種で、宣言の有無によらずボンドブロックを構成する種類の組。未宣言の種類の
  ブロックはゼロ）、\ ``longitudinal_bond_max_shells``\ （整数。未指定なら\ ``-1``\ ）、
  ``longitudinal_bond_cond_min_s`` / ``longitudinal_bond_cond_min_c``\ （浮動小数。
  スピン／電荷の RPA 分母の\ :math:`q`\ にわたる最小の条件数スコア。不安定性の
  下限に達すると実行は拒否されます）、\ ``longitudinal_bond_schema``\ （整数、\ ``1``\ ）。


データ読み込みの例
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

以下にデータの読み込み方の例を示します。

.. code-block:: python

    import numpy as np
    data = np.load('output/chiq.npz')

    chiq = data['chiq']
    freq_index = data['freq_index']

.. raw:: latex

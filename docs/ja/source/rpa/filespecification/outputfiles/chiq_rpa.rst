.. highlight:: none

.. _Subsec:chiq_rpa:

chiq, chi0q
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

感受率行列および既約感受率行列の計算結果を Numpy zip 形式で出力します。
ファイル名は、環境設定ファイルの中の ``file.output`` セクションでキーワード ``chiq`` および ``chi0q`` を用いて指定された文字列(以下 ``chiq_str`` と記述)を用いて、 ``chiq_str.npz`` という名前で出力されます。

ファイルの内容は、以下のキーにバインドされる複数の配列データからなります。

- ``chiq`` または ``chi0q``:

  感受率行列または既約感受率行列のデータ。データ形式は以下の節で説明します。

- ``momentum_convention``:

  フーリエ符号の来歴マーカー ``"e_plus_ikR"``\ （バージョン 2.0 以降）。運動量
  ラベルは :math:`M(k) = \sum_R M(R) e^{+ikR}` に従います。異なる値を
  記録したファイルはローダが拒否します。このフィールド導入前のファイル
  は、内容が :math:`q \to -q` に対して要素毎に偶（両規約が一致する
  場合）のときに限り受理されます。

- ``freq_index``:

  出力する松原振動数の値または範囲はパラメータ ``matsubara_frequency`` で指定されます。出力データの配列インデックスと、実際の松原振動数のラベルの対応付けを ``freq_index`` に格納します。

- ``wavevector_unit`` および ``wavevector_index``:

  波数ベクトルの情報を格納します。詳細は :ref:`UHFk の出力ファイル<Subsec:eigen_uhfk.dat>` を参照してください。

- ``index_convention``:

  スピン軌道の軸が spin*norb+orb（スピンブロック）順で並んでいることを示す文字列マーカー(常に ``"spin_block"``)です。これは UHFk のインターリーブ(2*orb+spin)出力とは異なります。このマーカーを持たない ``chi0q`` ファイル(規約修正前に出力されたもの)は、スピン軌道モードで再読み込みする際に拒否されます。

- ``coeff_tail``:

  感受率の計算時に使用された ``coeff_tail`` パラメータ(Matsubara 高周波数裾補正)の値です。裾補正は中程度の ``Nmat`` において :math:`\chi_0(\mathbf{q})` を :math:`O(1)` で変えるため、このファイルを読み込む計算(``hwave_sc`` など)は、自身の ``coeff_tail`` 設定が記録値と異なる場合に警告を出します。事前計算した ``chi0q_init`` ファイルをそのまま出力する場合には入力ファイルの値が引き継がれます。生成元がこの値を記録していない場合(旧バージョン)、および一様格子の裾補正が適用されない IR 松原基底(``matsubara_basis = "ir"``、densified 出力を含む)での FLEX 出力の場合には、キー自体が省略されます。

副格子を指定している場合は、出力されるデータは副格子を単位とした感受率の値です。波数ベクトルおよび軌道のインデックスは副格子に読み替えます。
   
``chi0q`` の出力ファイルは、計算済み既約感受率データとして ``file.input`` セクションの ``chi0q_init`` に指定して使用できます。


chi0q のデータ形式
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

``chi0q`` のデータ形式は、スピン軌道相互作用および外場の有無と、 ``mode.calc_scheme`` パラメータの値によって以下の形をとります。

- スピン非依存 (spin-free):

  ``enable_spin_orbital`` パラメータが ``false`` の場合、もしくは ``true`` であっても :math:`T_{\tilde\alpha\tilde\beta}(k)` がスピンについて対角かつ対称な場合で、外場がないとき、既約感受率行列のスピン非依存部分を出力します。

  - ``calc_scheme = general`` の場合、データ配列は ``ndarray(l,q,a,ap,b,bp)`` で、インデックスは以下のとおりです

    - ``l``: 松原振動数のラベル。インデックスとラベルの対応は前述の ``freq_index`` で与えられます。

    - ``q`` : 波数ベクトルのインデックス :math:`[q_x\ q_y\ q_z]` を1次元化したインデックスで、 :math:`q = q_z + N_z \cdot (q_y + N_y \cdot q_x)` となります。

    - ``a``, ``ap``, ``b``, ``bp`` はスピンを含まない軌道インデックス :math:`\alpha`, :math:`\alpha^\prime`, :math:`\beta`, :math:`\beta^\prime` に対応します。

  - ``calc_scheme = reduced`` または ``squashed`` の場合、データ配列は ``ndarray(l,q,a,b)`` です。インデックスの意味は上記と同じです。   

- スピン対角 (spin-diagonal):

  ``enable_spin_orbital`` パラメータが ``false`` で外場がある場合、もしくは ``enable_spin_orbital`` パラメータが ``true`` で :math:`T_{\tilde\alpha\tilde\beta}(k)` がスピンについて対角な場合に、既約感受率行列の spin up/down 成分を出力します。

  - ``calc_scheme = general`` の場合、データ配列は ``ndarray(s,l,q,a,ap,b,bp)`` となります。 ``s = 0 (1)`` はそれぞれ spin up (down) 成分を表し、それ以外のインデックスは上記と同じです。
    
  - ``calc_scheme = reduced`` または ``squashed`` の場合、データ配列は ``ndarray(s,l,q,a,b)`` となります。

- スピン依存 (spinful):

  ``enable_spin_orbital`` パラメータが ``true`` で Transfer項が一般的な形の場合、一般化軌道をインデックスとする既約感受率行列を出力します。

  - ``calc_scheme = general`` の場合、データ配列は ``ndarray(l,q,a,ap,b,bp)`` となります。 ``a``, ``ap``, ``b``, ``bp`` はスピンを含む一般化軌道インデックス :math:`\tilde\alpha`, :math:`\tilde\alpha^\prime`, :math:`\tilde\beta`, :math:`\tilde\beta^\prime` に対応します。

  - ``calc_scheme = reduced`` の場合、データ配列は ``ndarray(l,q,a,b)`` となります。 ``a``, ``b`` はスピンを含む一般化軌道インデックス :math:`\tilde\alpha`, :math:`\tilde\beta` に対応します。

  

chiq のデータ形式
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

``chiq`` のデータ形式は ``calc_scheme`` パラメータの値によって以下の形をとります。

- ``calc_scheme = general`` の場合、データ配列は ``ndarray(l,q,a,ap,b,bp)`` となります。 ``a``, ``ap``, ``b``, ``bp`` はスピンを含む一般化軌道インデックス :math:`\tilde\alpha`, :math:`\tilde\alpha^\prime`, :math:`\tilde\beta`, :math:`\tilde\beta^\prime` に対応します。

- ``calc_scheme = reduced`` の場合、データ配列は ``ndarray(l,q,a,b)`` となります。 ``a``, ``b`` はスピンを含む一般化軌道インデックス :math:`\tilde\alpha`, :math:`\tilde\beta` に対応します。

- ``calc_scheme = squashed`` の場合、データ配列は ``ndarray(l,q,s1,s2,a,s3,s4,b)`` となります。 ``a``, ``b`` は軌道インデックス :math:`\alpha`, :math:`\beta` に対応し、 ``s1`` 〜 ``s4`` はスピンインデックス :math:`\sigma`, :math:`\sigma^\prime`, :math:`\sigma_1`, :math:`\sigma_1^\prime` にそれぞれ対応します。

``chiq`` は縦方向（ring）の感受率です。 ``calc_scheme = general`` および ``squashed`` では添字の対がスピン非対角となるスロットを持ちます（ ``calc_scheme = reduced`` は密度対の成分のみを格納するため、そのようなスロットは存在しません）。これらのスロットは横感受率では **ありません** 。またバブルがスピンフリーあるいはスピン対角なバブルの膨張で得られる限り **計算されません** 。膨張は同スピンのスロットしか作らないため、残りは恒等的にゼロのままです。これは ``enable_spin_orbital`` を使わないすべての計算に加えて、 ``enable_spin_orbital`` を使っていても一体ハミルトニアンがスピンフリーあるいはスピン対角と判定される計算にも当てはまります。このゼロを横方向応答の計算結果として読まないでください。（真にスピンフルな計算では、バブルは一般化軌道添字の上で直接構築されるため、これらのスロットは計算され、一般にゼロではありません。その場合でも横感受率は ``chiq_pm`` です。）横方向チャネルの求め方は :ref:`アルゴリズム<rpa_which_array>` を参照してください。

``calc_type = ring+ladder`` の場合、 ``chiq`` ファイルには配列 ``chiq_pm`` も追加で格納されます。これは横感受率 :math:`\chi_{+-}(q)` を保持します。配列形式は ``ndarray(l,q,a,ap,b,bp)`` で、 ``a``, ``ap``, ``b``, ``bp`` はスピン自由度を含 **まない** 軌道インデックス :math:`\alpha`, :math:`\gamma`, :math:`\beta`, :math:`\delta` です（スピン構造は :math:`+-` のラベルで既に定まっているため）。したがって、対応する軸が一般化（スピン軌道）インデックスを走る ``chiq`` よりも小さい配列になります。縦方向の ``chiq`` はラダーの有無に影響されません。同一の入力に対して ``calc_type = ring`` の結果とビット単位で一致します。


データ読み込みの例
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

以下にデータの読み込み方の例を示します。

.. code-block:: python

    import numpy as np
    data = np.load('output/chiq.npz')

    chiq = data['chiq']
    freq_index = data['freq_index']

.. raw:: latex

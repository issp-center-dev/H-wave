.. _initial_states:

***********************************
初期状態の候補（プレビュー）
***********************************

UHF 計算が収束する状態は初期グリーン関数に依存します。
``python -m hwave.initial`` は UHFk の入力に対して初期状態の候補を用意し、
候補ごとに実行ディレクトリを作ります。各ディレクトリで ``hwave`` を実行すると、
このツールが結果を集めて自由エネルギーの順に並べます。

.. note::

   これはプレビュー版です。コマンド、オプション、ファイル形式は今後変わる
   可能性があります。このリリースは二つの tier を提供します。large tier は
   既定の候補の集合で、セルで表せるすべての q でのスピン・電荷のパターンと、
   上限で削られない予約行 ``para``\ 、 ``ferro``\ 、 ``random`` からなります。
   small tier はこれに、常磁性状態の静的な解析から得た軌道のパターンと
   不安定性の行を加えます（ :ref:`initial_states_small_tier` を参照）。
   より大きな模型のための medium tier と、較正したコストモデルによる tier の
   選択は、今後のリリースで加わります。

The selected state is the lowest free energy among the candidates tried. It is
not a proof of the ground state.
（選ばれる状態は、試した候補の中で自由エネルギーが最も低いものです。基底状態で
あることの証明ではありません。）

使い方の流れ
============

.. code-block:: bash

   $ python -m hwave.initial prepare INPUT_DIR --out WORK
   $ bash WORK/commands.txt          # または: cd WORK && parallel < commands.txt
   $ python -m hwave.initial collect WORK

すべてのオプションと終了コードは ``python -m hwave.initial --help`` および
``python -m hwave.initial <command> --help`` で表示されます。

``prepare``
-----------

``INPUT_DIR`` は UHFk の入力 TOML ファイルをちょうど一つ含むディレクトリ、
またはその TOML ファイル自身です。オプション:

- ``--out WORK``\ （必須）: 作業ディレクトリ。
- ``--subshape Bx,By,Bz``\ （複数指定可）: すべての候補をそれぞれの
  ``SubShape`` で評価します。省略時は入力の ``SubShape`` です。
- ``--extra-candidates FILE.json``: 追加の候補。
  ``{"family": "user" | "literature", "pattern": ..., "label": ...,
  "citation": ...}`` の JSON リストです。 ``pattern`` はパッケージと一緒に
  インストールされる JSON スキーマ ``hwave/initial/patternspec.schema.json``
  （Python では
  ``importlib.resources.files("hwave.initial") / "patternspec.schema.json"``\ ）
  に従います。軌道の番号は物理軌道の番号です（spin-orbital モードでも同じ
  です）。
- ``--t-hint T``: ヒント温度。単位は入力のエネルギーの単位です。small tier
  はこの温度でヒントを評価します（ :ref:`initial_states_small_tier` を参照）。
  変わるのは不安定性の行とその順序だけで、他の候補は変わらず、候補の実行は
  入力の ``T`` のままです。large tier では記録されるだけです。省略時は、
  ``T > 0`` なら入力の ``T``\ 、 ``T = 0`` ならバンドと k メッシュから決まる
  ``max(W/(2 L), 1e-4 W)`` です。 ``W`` はフェルミ準位付近のバンドの
  エネルギー幅（1e-3 以上）、 ``L`` は ``CellShape`` の 2 より長い軸のうち
  最も短いものです。使った規則は要約に表示されます。
- ``--amplitude A``\ （初期状態の振幅の要求値、(0, 1) の範囲、既定値 0.5）、
  ``--max-candidates N``\ （ ``SubShape`` ごとの生成行数、既定値 32）、
  ``--orbital-site FILE``\ （物理軌道ごとにサイトのラベルを一つ並べた JSON
  リスト）。
  ``T = 0`` では、二つの準位がフェルミ準位で交差するちょうどその場で、
  初期状態を作る場が要求した振幅に達することがあります。このとき初期状態は
  その交差点で作ります。交差する準位を一様に占有して要求を満たせばその占有を、
  満たさなければ交差を越えた直後の占有を使います。そのため ``coverage.json``
  の ``amplitude_used`` が要求値を上回ることがあります。スピンについて対称な
  模型では、スピン成分を持たないパターン（para、電荷、軌道）から、上向きと
  下向きのスピン成分が一致する初期状態ができます。
- ``--tier auto|large|small``: 候補の集合（ ``medium`` は今後のリリースの
  ために予約されています）。既定値の ``auto`` は、模型が十分に小さければ
  small tier を、そうでなければ large tier を選びます
  （ :ref:`initial_states_auto_tier` を参照）。
- ``--para-check``\ （small tier）: 事前の常磁性計算を一様な軌道占有からも
  始め、自由エネルギーの低い方を使います（ :ref:`initial_states_small_tier`
  を参照）。 ``--tier large`` と一緒には指定できません。 ``auto`` が large
  tier を選んだ場合は、警告を出して無視します。
- ``--stage0``\ （small tier）: 事前の常磁性計算を ``prepare`` の中で実行
  せず、バッチシステム向けに独立した実行ディレクトリとして書き出します
  （ :ref:`initial_states_preliminary_run` を参照）。large tier を指定した
  場合や large tier が選ばれた場合は拒否します。

``prepare`` は WORK に次を書きます。

- ``coverage.json``: 実行前の網羅性レポート（後述）。
- ``manifest.json``, ``candidates.json``: 用意した内容。各実行の
  ``input.toml``\ 、 ``initial.npz`` とコピーした入力の SHA-256 を含みます。
- ``commands.txt``: WORK に移動する最初の一行と ``#`` で始まる数行の
  コメントの後に、実行ごとに一行、WORK からの相対パスで書かれます（行番号で
  行を選ぶジョブ配列では最初の数行を飛ばしてください）。
- 実行ごとのディレクトリ。 ``input.toml``, ``initial.npz``, ``inputs/`` を
  含みます。

各実行の ``input.toml`` は入力を次のように書き換えたものです。
``[file.input] path_to_input = ""`` と ``initial = "initial.npz"``\ 、
``[mode.param] SubShape`` の設定（既存のキーはその綴りのまま置き換え）、
``initial_uhf`` と ``initial_mode`` の削除、入力ファイル（相互作用の
ファイル、 ``geometry_uhf``, ``onebodyg_uhf``\ ）を ``inputs/`` に
**コピー**\ （ ``[file.input.interaction] path_to_input = "inputs"``\ ）、
``[file.output] energy`` がなければ ``"energy.dat"`` を追加。元の TOML の
コメントは残りません。出力の絶対パスや ``..`` を含むパス、
``green = ""`` や ``energy = ""``\ 、実行ディレクトリ自身のファイルを上書き
する出力名は拒否されます。入力はリンクではなくコピーなので、WORK を移動
したり元のファイルを削除したりしても、実行と collect ができます。

キーはソルバーと同じ大文字・小文字の規則で読みます。 ``[mode]``,
``[file]``, ``[file.input]``, ``[file.output]`` の各テーブルとそのキーは
大文字・小文字を区別し、 ``[mode.param]`` と ``[file.input.interaction]``
のキーは区別しません。区別するテーブルで既知のキーの綴りが違う場合（例:
``[file.input] Initial``\ ）は、ソルバーがそのキーを無視するため、正しい
綴りを示して拒否します。

同じ WORK に同じコマンド・同じ入力で ``prepare`` を再実行しても何も
しません。コマンドや入力が違う場合、または ``prepare`` が書いたものでない
空でないディレクトリは、違いを示して拒否します。

複素数の係数を受け付けるのは ``PairHop`` だけです。 ``CoulombInter``\ 、
``Hund``\ 、 ``Ising``\ 、 ``Exchange``\ 、 ``PairLift``\ 、まとめた
``Coulomb`` の表に複素数の係数があると、 ``prepare`` は項、ボンド、軌道を
示して拒否します。これらの項では、係数が複素数のときソルバーの平均場
ポテンシャルが説明どおりのハミルトニアンのものにならず、その上に作った
候補は意味を持たないためです。

``2Sz`` を固定した場合、ソルバーは相互作用項のうちスピンを保存する部分
だけを使います。z 軸まわりのスピン保存を破る相互作用項（例: ``PairLift``\ ）
があると、 ``prepare`` はその項を名指しした警告を記録し、要約に表示します。

.. _initial_states_small_tier:

small tier
----------

small tier は large tier の集合に二つの family を加えます。

- ``orbital``: サイト上で隣り合う相互作用軌道 ``a``, ``b`` の軌道偏極
  ``n_a - n_b``\ 。一様なものとサイトの符号パターンを付けたものを、セルで
  表せるすべての q で作ります（baseline の行と同じく ``--max-candidates``
  に数えます）。
- ``instability``: 常磁性状態の静的な線形応答の解析から得たパターン
  （ヒント）。ヒントは不安定性の行の順序を決めるだけで、他の family の
  候補を削ることはありません。

ヒントは四つの段階で計算します。

1. **事前の常磁性計算。** UHFk を基本セル（ ``SubShape = [1, 1, 1]``\ ）で
   ヒント温度 ``T_hint`` において、入力の ``Ncond``\ 、 ``2Sz``\ 、
   ``flag_fock`` と相互作用の表を使い、裸の一体ハミルトニアンの自由な密度
   から実行します。 ``EPS`` は 12 以上に引き上げます（ ``para_eps`` として
   記録）。計算が収束し、模型の対称性を 1e-8 以内で保つときだけ、その結果を
   ヒントの基準にします。対称性とは、一体ハミルトニアンと相互作用の両方が
   保つスピン回転、 ``2Sz`` を固定した場合の各セクターの電子数、模型が持つ
   ときの spin-orbital モードでの時間反転です。そうでなければヒントは理由
   （ ``hint.skip_reason``\ ）とともに省略され、すべての初期状態は裸の一体
   ハミルトニアンの上に作られます。 ``--para-check`` を付けると、一様な軌道
   占有から二つ目の計算を始め、自由エネルギーの低い方を使います。差が
   セルあたり 1e-8 を超えると警告します（もう一方の出発点は準安定な軌道
   秩序状態に達したということです）。基準は ``WORK/para_hf/``\ （
   ``green.npz``\ 、 ``para_hf.json`` と計算）に保存されるので、 ``collect``
   で計算し直す必要はありません。すべての family の初期状態は、その
   ハミルトニアン ``H0 + Sigma'`` を各セルに折り畳んだものの上に作られ
   （ ``seed_base = "para_hf"``\ ）、 ``prepare`` と ``collect`` のラベルも
   それを基準にします。
2. **バーテックス** ``Gamma(q)``\ : ソルバーの Hartree-Fock ポテンシャルの
   密度による微分です。ソルバーが読み書きする一体のスロット（ペア空間
   ``P``\ 。相互作用軌道に限る）の上で作ります。入力の ``flag_fock`` に
   従い、Hartree 項は常に、Fock 項は実行が Fock 項を含むときだけ入ります。
   そのため ``CoulombInter`` のボンドからのボンド秩序のパターンは Fock 項が
   あるときだけ現れます。
3. **裸の感受率** ``chi0(q)``\ : ``T_hint`` での常磁性バンドの静的な
   Lindhard 応答を、同じペア空間の上で求めます。 ``q = 0`` では、ソルバーが
   反復ごとに化学ポテンシャルを合わせ直すのと同じく、各セクターの電子数
   （ ``2Sz`` を固定したときはスピンごと、それ以外は系全体で、相互作用の
   ない軌道も含みます）を固定した応答を使います。その電子数を一様に変える
   モードは不安定性になりません（ ``lambda = 0`` ）。 ``2Sz`` を固定したとき
   は一様な ``S^z`` もこれに含まれます。
4. **固有値** ``lambda``\ : ``CellShape`` のメッシュのすべての q で
   ``K(q) = -chi0(q) Gamma(q)`` の正の固有値のうち大きい方から
   ``top_k = 3`` 個を求めます（スピン回転の確認に通れば、電荷と ``S^z`` の
   セクターを別々に解き、打ち切りの前にまとめます）。各固有ベクトル ``x``
   から、共役な場 ``O = Gamma x`` を初期状態のパターンとする不安定性の行が
   一つできます。

``lambda`` はそのパターンの静的な Stoner 型の因子です。 ``lambda >= 1``
なら ``T_hint`` での常磁性状態は Hartree-Fock の範囲でそのパターンに対して
不安定で、 ``lambda`` が大きいほどその傾向が強いことを表します。 ``T_hint``
で ``lambda < 1`` となるのは普通のことで、入力の ``T`` でそのパターンが秩序
化しないという意味ではありません。不安定性の各行には ``rank``\ （順序での
位置）、 ``lambda``\ 、 ``lambda_imag``\ （バーテックスがエルミートでない
場合を除き 0。その場合は報告されます）が記録されます。行は ``lambda`` の
大きい順に並びます。互いに相対 1e-8 以内の固有値は一つの縮退した組とし、
決まった規則（決まった順序の q、次にその q での固有値の順位）で並べ、一つの q
の中の縮退した組には標準基底を与えます。そのため行とその id は固有値
ソルバーの丸め誤差に依存しません。ヒントの行はすべて一覧に載ります。
``SubShape`` ごとに実行されるのは最大 ``quota_instability = 12`` 行で、
残りの行にはその理由（ ``quota: instability``\ 、または q を表せるセルを
伴う ``unrepresentable``\ ）が書かれます。

**ヒント温度はヒントをぼかします。** 各遷移は ``(f(e) - f(e'))/(e' - e) <=
1/(4 T_hint)`` の重みで ``chi0`` に入るため、フェルミ準位から ``T_hint``
程度以内の構造や、 ``T_hint`` より小さいエネルギー差から来る構造（小さな
ギャップ、不完全なネスティング、分裂した van Hove ピーク）は抑えられます。
``lambda`` は下がり、ピークは広がり、最大の q がずれ、近い行の順序が入れ
替わることがあります。要約には、使った値とともにこの説明が書かれます。

``coverage.json`` のヒントの記録（ ``hint``\ ）には、ヒント温度とその規則、
常磁性計算の状態（ ``para_status``\ 、 ``para_free_energy``\ 、
``para_residual``\ 、 ``para_polarisation``\ 、 ``para_eps``\ ）、
バーテックスがエルミートかどうか、方法、q の集合と各 q の状態（ ``ok``\ 、
``fewer than top_k positive eigenvalues``\ 、 ``no positive eigenvalue``\ ）
が含まれます。

.. _initial_states_auto_tier:

tier の選択
-----------

このリリースでは ``--tier auto`` は模型だけで決まるので、同じ入力はどの
ホストでも同じ tier になります。ペア空間のスロット数が ``|P| <= 400``\ 、
メッシュの点の数が ``N_q <= 4096``\ 、かつ ``N_q |P|^3 <= 1e11`` なら small
tier、そうでなければ警告を出して large tier です。 ``--tier small`` を明示
すると、これらの上限を超えても警告を出して実行します。ホストのメモリは tier
を選びません。small tier のメモリの見積もり
``16 (|P|^2 + N_q top_k |P| + N_k nb^2 + N_k nb |P|)`` バイト（ ``N_k`` は
k 点の数、 ``nb`` はバンドの数）が利用できるメモリを超える場合に拒否する
だけです。 ``coverage.json`` の ``tier`` の記録に、サイズ、見積もり、規則が
書かれます。

.. _initial_states_preliminary_run:

バッチシステムでの事前の常磁性計算
----------------------------------

small tier の常磁性計算は通常 ``prepare`` の中で実行されます。大きな模型
では、独立したジョブとして実行できます。

.. code-block:: bash

   $ python -m hwave.initial prepare INPUT_DIR --out WORK --tier small --stage0
   $ bash WORK/stage0_commands.txt      # 上で終了コード 3: 先にこれを実行
   $ python -m hwave.initial prepare INPUT_DIR --out WORK --tier small --stage0

最初の ``prepare`` は ``WORK/stage0_para/``\ （ ``--para-check`` を付けると
``WORK/stage0_para_check/`` も）を候補の実行と同じ方法で用意し
（ ``SubShape = [1, 1, 1]``\ 、 ``T = T_hint``\ 、 ``EPS`` は 12 以上）、
``WORK/stage0_commands.txt`` を書いて終了コード 3 で終わります。同じ
``prepare`` コマンドを再実行すると、それらの入力を ``manifest.json`` と
照合し（変わったファイルは名前を示して拒否）、計算の結果を読みます。収束
し、かつ模型の対称性を保っていれば（ :ref:`initial_states_small_tier` の
手順 1）それが基準になり、候補の実行が書かれます。実行されたが失敗した、
収束しなかった、またはそれらの対称性を破った場合は、初期状態は裸の一体
ハミルトニアンで作られ、ヒントはその理由とともに省略されます。まだ実行されていない、または終わっていない
場合は、 ``prepare`` はコマンドをもう一度表示し、何も変えずに終了コード 3 で
終わります。 ``--stage0`` はコマンドの一部です。これを付けて用意した WORK
には、これを付けたときだけ再び入ります。

初期状態ファイル
----------------

各 ``initial.npz`` には UHFk が読むグリーン関数と来歴情報が入っています。
来歴情報は、候補と実行の id、パターンのラベル・q・位相・符号、
``SubShape``\ 、 ``CellShape``\ 、 ``Ncond``\ 、 ``2Sz``\ 、使った場の強さ
（ ``h_used``\ ）、入力の SHA-256、初期状態の基準です。 ``seed_base`` は
``"para_hf"``\ （small tier の常磁性の基準）または ``"bare"`` で、基準を
使った場合はその温度（ ``seed_base_para_T``\ ）と入力のダイジェスト
（ ``seed_base_para_digest``\ ）も記録されます。ライブラリの書き出し関数は、
別の入力や別のヒント温度に属する基準を拒否します。

候補の実行
----------

``bash WORK/commands.txt`` はどのディレクトリからでも実行できます。各行は
実行ディレクトリに入り、初期状態ファイルを確認してから
``hwave input.toml > hwave.log 2>&1`` を実行します。GNU parallel やジョブ
配列では、WORK から各行を実行してください。

``check``
---------

.. code-block:: bash

   $ python -m hwave.initial check CONFIG.toml [NPZ] [--require-provenance]

``CONFIG.toml`` が読む初期状態ファイル（または ``NPZ``\ ）を ``hwave`` の
実行前に確認します。確認項目は来歴情報、 ``SubShape``\ 、形状、エルミート
性、固有値の範囲、電子数、ソルバーのブロック構造です。来歴情報のない
ファイルは未検証として受け付けますが、 ``--require-provenance`` を付ける
と終了コード 4 で拒否します。

``collect``
-----------

``collect WORK`` はすべての実行を読み、 ``WORK/coverage.json``\ （実行後）と
``WORK/summary.txt`` を書きます。繰り返し実行できます。各実行の状態は次の
とおりです。

- ``seed_unverified``: ``input.toml`` が用意した ``initial.npz`` を指して
  いない、コピーした入力や npz が ``prepare`` の書いたものと違う、または
  確認で拒否された場合。
- それ以外で ``scf_status.dat`` があれば、それで決まります:
  ``converged``, ``not_converged``\ 、SCF を一度も行っていなければ
  ``failed``\ 。状態ファイルがあるのに読めない場合は ``failed``
  （ ``"scf_status.dat unreadable: ..."``\ ）となり、順位付けされません。
- 状態ファイルがない場合: ``not_run``\ （ ``hwave.log`` も出力もない）、
  ``convergence_unverified``\ （ ``hwave.log`` がなく出力がある。別の
  リストに印付きで並べます）、 ``failed``\ （トレースバックやソルバーの
  エラー）、 ``hwave.log`` から判定した ``not_converged`` または
  ``converged``\ 。途中で止められた実行（状態ファイルもトレースバックも
  ない）は ``not_converged`` となり、理由は
  ``hwave.log has no convergence message`` です。

収束した実行の出力は読めなければならず、読めない場合は ``failed`` です。
失敗した実行や収束しなかった実行があっても ``collect`` の終了コードは 0
です。各実行の状態はレポートで確認してください。

``summary.txt`` には、モデル（軌道、電子数、 ``T``\ 、 ``2Sz``\ 、
``CellShape``\ ）、tier とヒント温度とその規則、ラベルの閾値、行の数と
試さなかった行の理由、順位の表（id、ラベル、セルあたりの自由エネルギー、
最も低い行とのサイトあたりの差、状態、最終ラベル）、印の付いた行と警告、
後述の制限、このページの冒頭に引用した文が書かれます。

終了コード
==========

======  ====================================================================
コード  意味
======  ====================================================================
0       完了（ ``check``: 検証済み、または ``--require-provenance`` なしでの
        未検証）
1       エラー、または確認でファイルが拒否された
2       引数のエラー
3       ``prepare --stage0``: 事前の常磁性計算の実行待ち。
        ``bash WORK/stage0_commands.txt`` の後、同じ ``prepare`` を再実行する
4       ``check --require-provenance``: ファイルに来歴情報がない
======  ====================================================================

網羅性レポート
==============

``coverage.json`` には候補と ``SubShape`` の組み合わせごとに一行があります。行は
すべて残り、削られることはありません。試さなかった行にはその理由が書かれ
ます（対称性で他の状態と等しい場合は ``duplicate of <run_id> (<map>)``\ 、
``quota: <family>``\ 、または理由と提案するセルを伴う ``unrepresentable``\ 。
例: ``q_x denominator 3 does not divide B_x = 2`` と 3 × 3 の提案）。

- family は ``para`` と ``ferro``\ （予約行、 ``SubShape`` ごとに一回実行）、
  ``random``\ 、 ``baseline``\ （セルで表せる q でのスピン・電荷パターン）、
  ``orbital`` と ``instability``\ （small tier）、 ``user`` と
  ``literature``\ （ ``--extra-candidates`` から。上限で削られません）です。
- 候補は初期状態のパターンで名付けられます（例: ``spin-z ferro-orbital
  q=(1/2,1/2,0) φ=0``\ ）。実行後の状態には ``final_label`` が付き、重みの
  10% 以上を持つ成分を名指しします。逆のパリティの重みが 10% 未満なら時間
  反転（TR）について even または odd、そうでなければ mixed と呼びます。
  秩序ベクトルはそれが表すサイト上の成分の直後に置かれ、ボンドの成分は
  その後に続きます（例: ホッピングの一様な変化を伴うネール状態は
  ``spin-z ferro-orbital q=(1/2,1/2,0) + bond-charge``\ ）。閾値はレポートに
  記録されます。
- ``collect`` は収束した実行を、 ``T``\ 、セルあたりの ``Ncond``\ 、
  ``CellShape``\ 、 ``2Sz`` が等しい比較グループの中で、サイトあたり（サイト
  数が定まらない場合はセルあたり）の自由エネルギーの順に並べます。互いに
  サイトあたり 1e-4（入力のエネルギーの単位）以内の行には
  ``near-degenerate`` の印が付きます。収束
  しなかった行や失敗した行は理由とともに並べ、順位は付けません。

制限
====

UHF は秩序と秩序化温度を過大評価し、k メッシュとセルは有限です。次の状態は
どの候補にも含まれません: 複数の q の重ね合わせ（double-Q, triple-Q,
渦の結晶）、残した二つの位相以外のスライド位相、縮退した固有空間での標準
基底以外の組み合わせ、√3 × √3 のような回転したセルと与えた ``SubShape``
以外のセル、large tier での軌道・ボンドのパターンとヒント以外のボンドの
パターン、候補あたり複数の振幅、スピン回転の確認に通らない場合のスパイラル、
分母が 16 を超える q。要約にもこれらの制限が書かれます。

ヒントが表すのは常磁性状態だけで、静的に、 ``T_hint`` で評価したものです。
一次転移や強結合で到達する状態には、ヒントがまったくないこともあります。
そのため予約行と baseline の行は常に実行され、ヒントが行を削ることは
ありません。小さな ``lambda`` はその状態を否定する根拠ではなく、最大の
``lambda`` が結果を選ぶのでもありません。結果を選ぶのは、試した候補の中での
収束した実行の自由エネルギーだけです。

コマンドラインが対応するのは UHFk の入力のみです。ライブラリ関数
``hwave.initial.initial_state`` は、同じ格子の UHFk 入力とサイトの対応表から
UHFr の ``Initial`` ファイルも書けます。spin-orbital モードは UHFk のみ対応
します。

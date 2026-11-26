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
   可能性があります。このリリースが提供するのは large tier のみです。large
   tier は既定の候補の集合で、セルで表せるすべての q でのスピン・電荷の
   パターンと、上限で削られない予約行 ``para``\ 、 ``ferro``\ 、 ``random``
   からなります。常磁性状態の不安定性解析から選んだより小さな候補の集合は、
   今後のリリースで加わります。

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
- ``--t-hint T``: ヒント温度。単位は入力のエネルギーの単位です。レポートに
  記録されるだけで、候補の集合や初期状態を変えることはありません。省略時は、
  ``T > 0`` なら入力の ``T``\ 、 ``T = 0`` ならバンドと k メッシュから決まる
  ``max(W/(2 L), 1e-4 W)`` です。 ``W`` はフェルミ準位付近のバンドの
  エネルギー幅（1e-3 以上）、 ``L`` は ``CellShape`` の 2 より長い軸のうち
  最も短いものです。使った規則は要約に表示されます。
- ``--amplitude A``\ （初期状態の振幅の要求値、(0, 1) の範囲、既定値 0.5）、
  ``--max-candidates N``\ （ ``SubShape`` ごとの生成行数、既定値 32）、
  ``--orbital-site FILE``\ （物理軌道ごとにサイトのラベルを一つ並べた JSON
  リスト）、 ``--tier``\ （このリリースでは ``auto`` は ``large`` に
  なります）。
  ``T = 0`` では、二つの準位がフェルミ準位で交差するちょうどその場で、
  初期状態を作る場が要求した振幅に達することがあります。このとき初期状態は
  その交差点で作ります。交差する準位を一様に占有して要求を満たせばその占有を、
  満たさなければ交差を越えた直後の占有を使います。そのため ``coverage.json``
  の ``amplitude_used`` が要求値を上回ることがあります。スピンについて対称な
  模型では、スピン成分を持たないパターン（para、電荷、軌道）から、上向きと
  下向きのスピン成分が一致する初期状態ができます。

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

``2Sz`` を固定した場合、ソルバーは相互作用項のうちスピンを保存する部分
だけを使います。z 軸まわりのスピン保存を破る相互作用項（例: ``PairLift``\ ）
があると、 ``prepare`` はその項を名指しした警告を記録し、要約に表示します。

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
3       今後のリリースのために予約（事前の常磁性計算の実行待ち）。現在は
        返しません
4       ``check --require-provenance``: ファイルに来歴情報がない
======  ====================================================================

網羅性レポート
==============

``coverage.json`` には候補と ``SubShape`` の組ごとに一行があります。行は
すべて残り、削られることはありません。試さなかった行にはその理由が書かれ
ます（対称性で他の状態と等しい場合は ``duplicate of <run_id> (<map>)``\ 、
``quota: <family>``\ 、または理由と提案するセルを伴う ``unrepresentable``\ 。
例: ``q_x denominator 3 does not divide B_x = 2`` と 3 × 3 の提案）。

- family は ``para`` と ``ferro``\ （予約行、 ``SubShape`` ごとに一回実行）、
  ``random``\ 、 ``baseline``\ （セルで表せる q でのスピン・電荷パターン）、
  ``user`` と ``literature``\ （ ``--extra-candidates`` から。上限で削られ
  ません）です。
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
以外のセル、large tier での軌道・ボンドのパターン、候補あたり複数の振幅、
スピン回転の確認に通らない場合のスパイラル、分母が 16 を超える q。要約にも
これらの制限が書かれます。

コマンドラインが対応するのは UHFk の入力のみです。ライブラリ関数
``hwave.initial.initial_state`` は、同じ格子の UHFk 入力とサイトの対応表から
UHFr の ``Initial`` ファイルも書けます。spin-orbital モードは UHFk のみ対応
します。

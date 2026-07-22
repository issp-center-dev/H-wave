.. highlight:: none

occupation.npz
~~~~~~~~~~~~~~

収束した SCF の占有数、化学ポテンシャル、温度、列ごとのメタデータを
``occupation.npz`` として出力する。キー名は
``[file.output].occupation`` で設定可能(既定 ``occupation.npz``)。
``tools/uhfk_to_mvmc.py`` ブリッジ等の消費者は ``eigen.npz`` と
合わせてこの情報を読み、SCF を再計算せずに有限温度占有を T=0 の
Slater 行列に射影する。

キー
^^^^

- ``occupation`` (shape ``(nvol, nd)``, ``float64``)
    SCF 最終反復で得られた (k 点, バンド) ごとの Fermi-Dirac 重み
    :math:`f(\epsilon_{k,n})`。T=0 は 0 か 1、T>0 では分数。

- ``mu`` (shape ``(n_mu_groups,)``, ``float64``)
    mu-group ごとの化学ポテンシャル。Sz-fixed (``2Sz = 0``) では 2 グループ
    (up と down)、Sz-free では 1 グループ。

- ``T`` (scalar ``float64``)
    SCF で用いた温度。

- ``column_spin`` (shape ``(nd,)``, ``int64``)
    各固有ベクトル列のスピン性格。``0`` = up のみ block、``1`` = down のみ block、
    ``-1`` = mixed (Sz-free)。

- ``column_mu_group`` (shape ``(nd,)``, ``int64``)
    各列の mu-group index。``mu[column_mu_group[n]]`` で列 ``n`` の
    化学ポテンシャルが引ける。

規約
^^^^

列レイアウトは ``eigen.npz`` の ``eigenvector`` と一致する: 列は
``UHFk._init_block_structure`` の block 順に並び、行は元の ``nd``
indexspace に置かれる。

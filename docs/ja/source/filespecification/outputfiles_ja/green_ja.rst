.. highlight:: none

.. _Subsec:cgcisajs:

green
~~~~~~~~~~

OneBodyGで指定された一体グリーン関数\ :math:`\langle c_{i\sigma_1}^{\dagger}c_{j\sigma_2}\rangle`\ の計算結果を出力します。
ファイル名は環境設定ファイルの中の ``file.output`` セクションでキーワード ``green`` を用いて指定することができます。
以下にファイル例を記載します。

::

    0 0 0 0 0.9517526553947047 0.0
    0 0 1 0 -0.03971951040016314 0.0
    0 0 2 0 0.09202884754223833 0.0
    0 0 3 0 -0.039719448981075135 0.0
    0 0 4 0 0.09202884754219534 0.0
    0 0 5 0 -0.03971947216145664 0.0
    0 0 6 0 0.09202884753253462 0.0
    0 0 7 0 0.09202884754259735 0.0
    0 1 0 1 0.04824734460529617 0.0
    0 1 1 1 0.03971951040016307 0.0
    …

ファイル形式
^^^^^^^^^^^^

-  :math:`[`\ int01\ :math:`]`  :math:`[`\ int02\ :math:`]`  :math:`[`\ int03\ :math:`]`  :math:`[`\ int04\ :math:`]`  :math:`[`\ double01\ :math:`]`  :math:`[`\ double02\ :math:`]`

パラメータ
^^^^^^^^^^

-  :math:`[`\ int01\ :math:`]`, :math:`[`\ int03\ :math:`]`

   **形式 :** int型

   **説明 :**
   サイト番号を指定する整数。\ :math:`[`\ int01\ :math:`]`\ が\ :math:`i`\ サイト、\ :math:`[`\ int03\ :math:`]`\ が\ :math:`j`\ サイトを表します。

-  :math:`[`\ int02\ :math:`]`, :math:`[`\ int04\ :math:`]`

   **形式 :** int型

   | **説明 :**
     スピンを指定する整数。\ :math:`[`\ int02\ :math:`]`\ が\ :math:`\sigma_1`\ 、\ :math:`[`\ int03\ :math:`]`\ が\ :math:`\sigma_2`\ に対応します。
   | 0: アップスピン
   | 1: ダウンスピン
   | を表します。

-  :math:`[`\ double01\ :math:`]`, :math:`[`\ double02\ :math:`]`

   **形式 :** double型

   | **説明 :**
     :math:`\langle c_{i\sigma_1}^{\dagger}c_{j\sigma_2}\rangle`\ の値を表します。
   | :math:`[`\ double01\ :math:`]`\ が実部、\ :math:`[`\ double02\ :math:`]`\ が虚部を表します。

.. raw:: latex

   \newpage

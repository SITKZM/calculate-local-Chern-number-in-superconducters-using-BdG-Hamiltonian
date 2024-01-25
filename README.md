# calculate-local-Chern-number-in-superconducters-using-BdG-Hamiltonian
## 日本語
超伝導体における局所チャーン数[1]を計算するためのコードです。
ローカライザーLの定義がclassAのものになっていますが、classDのもの[2]でも同じ結果になることを数値的には確かめています。

### ファイルの説明
ファイル"coordinates.txt"には、1列目に各サイトのx座標、2列目にy座標が入っています。

ファイル"hop.txt"には、1列目に飛び移り先サイトのインデックス、2列目に飛び移り元サイトのインデックス、3列目にホッピングの単位方向ベクトルのx成分、4列目にホッピングの単位方向ベクトルのy成分が入っています。

ファイル"pair_potential.txt"には、1列目に各サイトのペアポテンシャルの絶対値、2列目にその位相が入っています。

ファイル"particle_number.txt"には、1列目に各サイトのアップスピン電子数の期待値、2列目にダウンスピン電子数の期待値が入っています。

## English
This code is for calculating the local Chern number [1] in superconductors.
The definition of localizer L is for classA, but we have numerically confirmed that the same result is obtained for classD [2].

### File Description
The file "coordinates.txt" contains the x-coordinates of each site in the first column and the y-coordinates in the second column.

The file "hop.txt" contains the index of the destination site in the first column, the index of the source site in the second column, the x component of the unit direction vector of hopping in the third column, and the y component of the unit direction vector of hopping in the fourth column.

The file "pair_potential.txt" contains the absolute value of the pair potential for each site in the first column and its phase in the second column.

The file "particle_number.txt" contains the expected number of up-spin electrons for each site in the first column and the expected number of down-spin electrons in the second column.

## 参考文献（References）
[1] A. Cerjan and T. A. Loring, Physical Review B 106, 064109 (2022).

[2] T. A. Loring, Annals of Physics 356, 383 (2015).

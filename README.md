# 3dPanner

HRIRのデータを保持するSOFAファイルをもとに3Dパンニングを行うアプリケーションの開発です。
開発のフレームワークにはJUCEを使用しており(https://juce.com/)
、外部スクリプトは既存の3Dパンナ―のOrbiter(https://github.com/superkittens/Orbiter?tab=readme-ov-file)
を参考にしました。




## 外部ライブラリ

[libBasicSOFA](https://github.com/superkittens/libBasicSOFA) を使用しSOFAファイルのデータを取得するようにしました。[Orbiter](https://github.com/superkittens/Orbiter?tab=readme-ov-file)のIssuesの項目にある通り、
Theta/Phi/Radiusが必ずとれる状態でなかったので少しスクリプトを改変しています。




## 使用可能なSOFAファイル

https://sofacoustics.org/data/database/fhk/
このリンクのsofaファイルは読み込めることが確認されている。読み込めないものの場合はplayボタンを押しても無音になってしまっている。

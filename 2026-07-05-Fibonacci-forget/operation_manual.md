# B13地震解析 実行手順書 v1.0

**作成日**: 2026-04-20  
**対象マシン**: moroc-HP-EliteDesk-800-G5-SFF（Ubuntu 24.04）  
**IPアドレス**: 192.168.11.16

---

## 環境構成

```
~/b13/
  ├── scripts/
  │   ├── b13phase.py        # B13ライブラリ
  │   ├── b13_pipeline.py    # 波形→CSV変換
  │   ├── b13_to_db.py       # CSV→PostgreSQL登録
  │   └── knet_download.py   # K-NETダウンロード
  ├── data/
  │   ├── *.zip              # K-NETダウンロードZIP
  │   └── results/           # b13_pipeline出力CSV
  └── work/                  # 作業用
```

**PostgreSQL**
- DB名: b13db
- ユーザー: moroc
- ホスト: localhost:5432

---

## 標準作業フロー

### 1. K-NETからデータ取得
K-NETサイト（要ログイン）で条件検索してZIPをダウンロード:
```
https://www.kyoshin.bosai.go.jp/ja/eqdownload/
```
ダウンロードしたZIPを `~/b13/data/` に置く。

### 2. B13解析（ZIP→CSV）
```bash
python3 ~/b13/scripts/b13_pipeline.py \
  ~/b13/data/20260326231857_csv.zip \
  ~/b13/data/results/20260326231857.csv
```

出力例:
```
結果: 506点  有意=129(25.5%)  正=334(66.0%)
```

### 3. DBに登録（CSV→PostgreSQL）
```bash
python3 ~/b13/scripts/b13_to_db.py \
  --csv ~/b13/data/results/20260420165300.csv \
  --event 20260420165300
```

**イベントIDの形式**: `YYYYMMDDHHmmSS`（例: `20110407233241`）

### 4. DB確認
```bash
python3 -c "
import psycopg2
conn=psycopg2.connect(host='localhost',dbname='b13db',user='moroc',password='b13research2026')
cur=conn.cursor()
cur.execute('SELECT event_id, count(*) FROM b13_results GROUP BY event_id ORDER BY event_id')
for row in cur.fetchall(): print(row)
conn.close()
"
```

### 5. CSVをClaudeに持ち込む
`~/b13/data/results/` のCSVファイルをClaudeにアップロード。
Moran's I・統計検定・可視化・解釈はClaudeが担当。

---

## 複数イベントの一括処理

```bash
# ZIPを全部処理
for zip in ~/b13/data/*_csv.zip; do
  name=$(basename $zip _csv.zip)
  out=~/b13/data/results/${name}.csv
  if [ ! -f $out ]; then
    python3 ~/b13/scripts/b13_pipeline.py $zip $out
  fi
done

# DB一括登録
for csv in ~/b13/data/results/*.csv; do
  name=$(basename $csv .csv)
  python3 ~/b13/scripts/b13_to_db.py --csv $csv --event $name
done
```

---

## トラブルシューティング

**PostgreSQLが起動していない場合:**
```bash
sudo systemctl start postgresql
```

**psqlに接続できない場合（Pythonは通る）:**
```bash
# パスワード入力時にタイポに注意
psql -U moroc -d b13db -h localhost
# パスワード: b13research2026
```

**パイプラインが0点を返す場合:**
- 記録が短すぎる（15秒未満）→ 正常除外
- ASCII形式の場合は観測点数が少ない（3〜5点）→ 正常

---

## プロトコルv1.0 固定パラメータ（変更禁止）

```python
FS=100, W=500, STR=50, EXP_H=4/13
N_SHUF=200, SEED=42
BP_LOW=0.5, BP_HIGH=10.0, BP_ORDER=4
```

これを変えると過去の結果との比較ができなくなる。

---

## WindowsからSSH接続

```
ssh moroc@192.168.11.16
```

ファイル転送（WinSCP推奨）:
- ホスト: 192.168.11.16
- ユーザー: moroc
- ポート: 22

---

**END OF MANUAL v1.0**

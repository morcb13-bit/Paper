## ローカルAIとローカルSQLサーバーにフィボナッチ忘却を実装　　ー　公開ノート　

## 苦痛：iclaudeは賢いのに、昨日のことを全く覚えていない

それは、前に解決したでしょう！？といったところでclaudeには記憶がないし、私の記憶も曖昧である。

Claudeは賢い。しかしセッションが終われば記憶がリセットされる。

これはプライバシーと安全性の観点から正しい設計だ。しかし継続的な研究作業には致命的な弱点になる。
私はB13という物理理論の研究を続けながら、セッションをまたぐたびに「引継書」を書いもらっていた。v1から始まったそれは今やv111まで積み上がっている。
しかし、引継書だけでB13の思想を伝えるには知識が全然足りない。
それを少しでも補う手段はないかと悩んだ末、思いついたのがローカルサーバーでの引継書の管理である。

## 解決策の設計

3つの要素を組み合わせた：
### クラウドAI（Claude）    ← 高精度な推論・実装・執筆
      ↕ テキスト
### ローカルAI（Ollama）    ← 記憶の補完・引継書の検索
      ↕ SQL
### PostgreSQL              ← 引継書の永続保存・構造化

※ローカルモデルは精度がクラウドAIより劣るため、最終判断はClaude等で行う。これが現実的な分業です。

## 環境構築

### ハードウェア

HP EliteDesk 800 G5 SFF（中古）
Ubuntu 24.04 LTS
RAM 32GB（Ollama用）

※ GPUなしでは会話には遅すぎるが、引継書を自動的に要約して登録したり、文書の検索といった用途ならこれで十分だ。

### インストール

bash# PostgreSQL 16
sudo apt install postgresql-16

### Ollama

curl -fsSL https://ollama.ai/install.sh | sh

### 日本語モデル

ollama pull lucas2024/llama-3-elyza-jp-8b:q5_k_m

### DB作成

sqlCREATE DATABASE ＜データベース名＞;
CREATE USER ＜userID＞ WITH PASSWORD ＜パスワード＞;
GRANT ALL ON DATABASE b13db TO ＜データベース名＞;

### 引継書テーブルの設計

sqlCREATE TABLE handover (
    id SERIAL PRIMARY KEY,
    version VARCHAR(20),
    filename TEXT,
    content TEXT,
    created_at TIMESTAMP DEFAULT NOW()
);

### セッション終了時に引継書を登録：

pythonimport psycopg2

conn = psycopg2.connect(
    host='localhost', dbname='****',
    user='****', password='****'
)
cur = conn.cursor()
cur.execute(
    "INSERT INTO handover (version, content) VALUES (%s, %s)",
    ('v111', open('handover_v111.md').read())
)
conn.commit()

## 核心：フィボナッチ忘却という設計

人間の記憶は「均等に忘れる」のではなく、「急速に忘れ、ゆっくり残る」。それを離散化したのがフィボナッチ重みです。
新しい記憶 → 重み377
1つ前     → 重み233
2つ前     → 重み144
3つ前     → 重み 89
...（フィボナッチ数列に従って減衰）
容量を超えたら古いものをまとめてバンドル化し、圧縮して次の層へ：
Layer 0: [w=377] v111  [w=233] v110  [w=144] v109 ...
Layer 1: [Bundle] v90-v100 の要約
Layer 2: [Bundle] v70-v90 の要約

### 実装：

pythonFIB_WEIGHTS = [377, 233, 144, 89, 55, 34, 21, 13, 8, 5, 3, 2, 1, 1]

def recalculate_weights():
    """file_mtime降順（新しい順）でFIB_WEIGHTSを割り当て"""
    conn = get_conn()
    cur = conn.cursor()
    cur.execute("""
        SELECT id FROM fib_memory
        WHERE is_bundle = FALSE
        ORDER BY file_mtime DESC NULLS LAST
    """)
    rows = cur.fetchall()
    for i, (id,) in enumerate(rows):
        w = FIB_WEIGHTS[min(i, len(FIB_WEIGHTS)-1)]
        cur.execute("UPDATE fib_memory SET weight=%s WHERE id=%s", (w, id))
    conn.commit()

### ローカルAIに引継書を読ませる

pythonimport requests

def ask_ollama(question):
    cur.execute("""
        SELECT title, summary, content
        FROM fib_memory
        ORDER BY weight DESC LIMIT 3
    """)
    rows = cur.fetchall()
    context = "\n---\n".join([
        f"{r[0]}\n{r[1] or r[2][:300]}" for r in rows
    ])
    prompt = f"以下はB13研究の引継書です。\n\n{context}\n\n質問: {question}\n回答:"
    resp = requests.post('http://localhost:11434/api/generate', json={
        "model": "lucas2024/llama-3-elyza-jp-8b:q5_k_m",
        "prompt": prompt, "stream": False
    }, timeout=300)
    return resp.json()['response']

### 実際の動作例：

質問: 次にやることを教えてください
回答: 20110407233241（最大余震M7.1）と20050816135518（2005年宮城沖M7.2）
      のB13解析とDB登録が次のステップです。

### 実証された分業

あるとき、K-NETのダウンロードスクリプトがどこにあるか分からなくなった。ローカルAIに引継書を検索させた。
返ってきたコードの断片はURLが壊れていてロジックも不完全だった。しかし曖昧な記憶から断片を拾ってきた。
Claudeがその断片と過去チャットを照合して正確なスクリプトの在処を特定した。

### まとめ

課題解決策クラウドAIの記憶なし引継書（テキスト）で補う引継書が増えすぎるフィボナッチ重み付きDB管理何を読ませるか重要度順に上位3件を選択断片的な記憶ローカルAIが検索・Claudeが検証
クラウドAIの弱点を補うために、ローカルAIに引継書を読ませている。その引継書を書いているのも、クラウドAIと一緒だ。そして別のクラウドAI（ChatGPT）が記事を査読している。
記憶のないAIたちが、テキストファイルを介して協力している。

### 今後の拡張

意味検索（embedding）× weightの組み合わせで、関連性と新しさの両方で引けるようになります。これは次のステップです。

構築に使ったスクリプト類は順次公開予定です。

現状は、圧縮やデータの削除は実装していないので危険はない。ubuntuのローカルサーバーを構築できる環境にある人は、自分の作ったドキュメントをまとめて放り込んでいけばいいだけだ。

そうすれば、ローカルAIで自在に検索できる。

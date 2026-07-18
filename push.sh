#!/usr/bin/env bash
# 便捷推送脚本：SSH 加好 key 后直接跑此脚本推送
# 用法: ./push.sh [remote] [branch]
set -euo pipefail
cd "$(dirname "$0")"
REMOTE="${1:-origin}"
BRANCH="${2:-feature/grpc-server-impl}"

echo "[push] 测试 SSH 连通性..."
if ! ssh -T -o StrictHostKeyChecking=no -o ConnectTimeout=10 git@github.com 2>&1 | grep -q 'successfully authenticated'; then
  echo "[push] ✗ SSH 认证失败。请先把公钥加到 GitHub:"
  echo "      cat ~/.ssh/id_ed25519.pub"
  echo "      粘贴到 https://github.com/settings/ssh/new"
  exit 1
fi

echo "[push] SSH 认证成功，推送 ${BRANCH} 到 ${REMOTE}..."
git push "$REMOTE" "$BRANCH"
echo "[push] ✓ 推送完成"

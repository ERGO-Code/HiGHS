# ⚠️ 待推送提交

GitHub HTTPS (443) 在当前网络环境下被阻断，已切换到 SSH (22) 协议。

## 待推送的提交
- `e4d8a2c` refactor: 重要示例迁到 examples/grpc/

## 推送步骤

### 1. 添加 SSH 公钥到 GitHub
```bash
cat ~/.ssh/id_ed25519.pub
```
复制输出，粘贴到 https://github.com/settings/ssh/new

### 2. 推送
```bash
./push.sh
# 或直接
git push origin feature/grpc-server-impl
```

## 验证 SSH 认证
```bash
ssh -T git@github.com
# 成功提示: Hi awakeningofthetrailblazer! You've successfully authenticated...
```

## 备注
- SSH key 已生成: `~/.ssh/id_ed25519` (ed25519, 无密码)
- SSH config 已配置: `~/.ssh/config` (github.com 走 22 端口)
- remote 已切换为 SSH: `git@github.com:awakeningofthetrailblazer/HiGHServer.git`
- 若 SSH 也不通，可改用其他网络环境后执行 `git push`

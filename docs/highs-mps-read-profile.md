# HiGHS 大型 MPS 读取性能：可复现的改进线索

2026-10-04，Intel Core i7-14700，Release 构建；本地 HiGHS 1.15.1 / hipdlp2
分支 `d3e4aca596ba733593453cc45f5f0744ac8472c7`。这是本地版本的测量，尚未
在上游 HiGHS 最新版重复，不能直接视为上游所有版本的性能结论。

## 结论与范围

`Highs::readModel()` 的大型 MPS 读取有明显 CPU 开销，发生在 HiPDLP 或其他
求解器运行之前。这里的“读取”包含文本解析、稀疏矩阵构建、输入检查和模型注册，
**不包含 presolve、HiPDLP scaling 或 GPU 初始化**。

三者均关闭 presolve，并读取同一个已经解压的 canonical `.mps` 文件。此处
没有 gzip 解压成本。输入归一化和 HiGHS 导入结果已有逐字段一致性审计。

| 模型 | 非零元 | HiGHS 读取，秒 | HPR-LP-C 读取，秒 | cuPDLPx 读取，秒 |
|---|---:|---:|---:|---:|
| a2864 | 20,078,717 | 29.564 | 21.752 | 7.385 |
| neso-2005 | 36,876,017 | 48.395 | 37.996 | 16.483 |

数据来自[八例比较](hipdlp-expanded-gpu.md)。各入口都返回可供求解的内存模型，
但内部构建过程不同；这不是三个解析函数逐条操作相同的微基准。测量为单次运行。

对约 426 MiB 的 `a2864` 文件，单独拆分 HiGHS 路径得到：

| 阶段 | Wall time，秒 | CPU time，秒 |
|---|---:|---:|
| `readModelFromFile`：解析与矩阵构建，含 reader 清理 | 28.630 | 28.622 |
| `passModel`：模型注册与验证 | 1.411 | 1.409 |

CPU 与 wall time 接近，说明此样本主要耗在 CPU 处理，磁盘等待不是主要部分。
这是独立复测，所以两项之和不要求与前表的 29.564 秒完全相等。

## 解析器来源与改进方向

cuPDLPx 使用仓库自带的
[`src/mps_parser.c`](https://github.com/MIT-Lu-Lab/cuPDLPx/blob/083b90e482f9f0e8229eac965d189a1ba386ca20/src/mps_parser.c)，
核心为 `FastLineReader`：4 MiB 缓冲区、普通文件 `fread`、gzip 文件 `gzread`、
原地分词、名称哈希表和稀疏矩阵构建。它与 PSLP presolver 是分开的模块。

当前 HiGHS 路径是 `Highs::readModel` → `FilereaderMps` → `HMpsFF`，随后
`passModel` 调用 `assessLp` 等输入检查。源码中可以看到逐行字符串处理、trim、
名称查找、数值字符串复制以及中间矩阵条目构建。**这些只是待采样的候选热点，
目前没有把 28.63 秒精确归因到某一个函数，也未证明单纯放大缓冲区即可解决。**

建议先在上游版本复现，再采样 `HMpsFF::parseCols`、`getMpsLine`、`getValue`、
矩阵构建和输入检查；评估减少字符串分配、复用解析缓冲区、直接构造稀疏数组等
方案。改进需继续保留重复项、名称、空列、边界、D 指数、异常输入及压缩文件语义。

## 复现

工具为 `hipdlp-scaling-study/tools/highs_read_profile.cpp`，只读取和注册模型，
没有调用 `run()` 或任何 presolver。传入本地已审计的 `.mps` 路径：

```sh
c++ -O3 -std=c++11 -Ihighs -Ibuild-scaling-cuda \
  docs/hipdlp-scaling-study/tools/highs_read_profile.cpp \
  -Lbuild-scaling-cuda/lib -Wl,-rpath,"$PWD/build-scaling-cuda/lib" \
  -lhighs -o /tmp/highs-read-profile
/tmp/highs-read-profile /tmp/hipdlp-expanded-gpu/models/a2864-canonical.mps
```

[分阶段原始数据](hipdlp-scaling-study/data/2026-10-04-read-profile/a2864.json)；
[输入哈希与审计](hipdlp-scaling-study/data/2026-10-04-expanded-gpu/manifest.json)。

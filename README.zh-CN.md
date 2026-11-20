# FastGA Joint Index

从 FASTA 直接构建的多基因组（2–64 个）无损联合索引。完整保留参考基因组的种子位置，
用 minimizer 锚点发现的精确拷贝块表示另一基因组的共享位置，只单独存储剩余位置。
读取时恢复完整位置及出现次数，继续调用 FastGA 的播种、链化和局部比对。

核心 C/C++ 代码与已验证版本逐字节一致；本仓库将其独立构建、测试和文档化，
不再依赖原 SXGC 大仓库。当前为实验性方法，联合模式固定 k=40，支持 2–64 个输入。

## 主要来源与归属

**本项目的主要来源是 PFP（Prefix-Free Parsing）和
[SXGC](https://github.com/pangenome/sxgc)。** PFP 提供可选的 Big-BWT 后缀排序后端；共享精确
拷贝由 minimizer 锚点查找器发现；共享索引代码最初在 SXGC 工作区中作为 FastGA 集成开发，再整理为
本独立仓库。[CORE_SOURCES.json](CORE_SOURCES.json) 记录了提取的 C/C++ 源码校验和。

| 来源 | 在本项目中的作用 |
|---|---|
| [PFP：Boucher 等，2019](https://doi.org/10.1186/s13015-019-0148-5)；实现为 [Big-BWT](https://gitlab.com/manzai/Big-BWT) | 可选的 PFP 后缀排序后端 |
| [SXGC](https://github.com/pangenome/sxgc) | 共享索引集成、构建和验证流程的原始研发工作区 |
| [FastGA：Myers 与 Zhou](https://github.com/thegenemyers/FASTGA) | 原有采样、播种、链化、局部比对、GDB 和输出格式 |
| [libsais](https://github.com/IlyaGrebnov/libsais) | PFP 后端中解析与字典的后缀数组和 LCP 构建 |

本仓库整理了 FastGA 兼容的共享种子表示、精确拷贝与残余位置编码、读取适配及
单命令入口。SXGC 的 suffixient-array 形式化证明适用于其自身结构；本适配器以
全量记录、种子对和比对结果测试提供正确性证据。使用时请标注 PFP、SXGC 两项
主要来源，并同时引用 FastGA 及实际使用的上游实现。完整引用见文末及
[第三方说明](THIRD_PARTY_NOTICES.md)。

## 构建与使用

构建依赖 Linux/macOS、CMake ≥3.20、C11/C++17 编译器、zlib、pthread、patch；
Linux 另需 OpenSSL ≥1.1 开发文件，
macOS 使用系统 CommonCrypto 计算 SHA-256。JSON 库及原有上游源码已固定版本保存，
构建不下载源码。**独立索引命令 `GIXmake-shared` 与 `FastGA-shared` 共用 C++
构建器，原生排序和 PFP 排序均不调用 Python；原 Python 索引构建器已删除。**
**构建仅使用 CMake，不依赖 Python**；源码校验由 CMake 完成，构建记录由 C++ 工具生成。
Python ≥3.11 仅用于可选测试、基准分析和 FASTA 准备工具；保留 `scripts/merge_fasta.py`。

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --parallel 4
export PATH="$PWD/build:$PATH"

FastGA-shared -T24 -pafx A.fa B.fa > A_B.paf
FastGA-shared -T24 -pafx A.fa > A.paf
```

完整测试需安装 `requirements-benchmark.txt` 并配置 `-DSJ_PYTHON_TESTS=ON`，然后重新构建并运行
`ctest --test-dir build --output-on-failure`；默认仅启用原生 codec 测试。
基准脚本现在同时输出索引大小与完整构建时长：默认四轮独立从 FASTA 构建，报告原始秒数、
中位数及联合索引相对原生/v3 的耗时比例，区分双样本、合并自比对和全部视图。
合并视图计时包含新建 GDB，复用原生索引时未实测的构建时间标为空；旧实验不补造计时数据。
仅比较构建可使用 `--construction-only`：完成全部构建正确性检查后报告 `MEASURED`，
如实保留联合索引更大的结果，不宣称已通过空间或速度验收。
CMake 迁移验证见 [构建与计时检查](benchmarks/cmake-build-verification.json)。
当前读入器的块缓存复用、批量恢复和位图查找见 [实现说明](docs/reader-optimization.md)，
完整模拟数据及 sanitizer 结果见 [读入器正确性记录](benchmarks/reader-implementation-verification.json)。

以上直接读取你的 FASTA：两个输入自动构建 v4 联合索引并比对；一个输入构建
v3 独立压缩索引并自比对，也支持 `.fna`、`.fasta` 及其 gzip 文件。
建索引日志写入 stderr，stdout 只有比对结果。默认结束后清理临时索引；加 `-k`
会在当前目录的 `.fastga-shared/` 保留索引，并按完整输入内容、掩码、参数和
构建工具校验和复用，多个进程同时调用时只构建一份。`-P<目录>` 指定临时空间，
`-T1` 到 `-T32` 均可使用；低线程时顺序写入完整种子流以遵守预算。
首次运行的总耗时包含建索引。`-1:` 输出会自动保留所需 GDB/索引。

原有 `.gix` 参数直接交给原比对内核。移动工具时请一并保留原生程序及 `pfp-tools/`；
完整说明见[命令行文档](docs/cli.md)。

也可以继续显式构建、保存索引后比对：

```sh
mkdir -p indexes
GIXmake-shared \
  --storage joint --reference A --threads 24 \
  --input A=genome-A.fa --input B=genome-B.fa \
  --output indexes/joint
FastGA-shared -T24 -pafx indexes/joint/A.gix indexes/joint/B.gix > A-to-B.paf
```

`genome-A.fa`、`genome-B.fa` 指你自己的输入，也支持 `.fa.gz`；输出目录必须不存在。
`--threads` 是构建总预算，范围为 1–32；低线程时原生排序采用顺序种子流，
PFP 排序至少需要 2 个线程。默认使用原生排序；添加
`--backend pfp` 可替换排序阶段。两种模式均用 minimizer 锚点查找共享拷贝块，仅 pfp 后端使用 PFP 解析。
`build/` 内记录真实编译命令和源文件/可执行文件校验和。

`GIXmake-shared --help` 列出完整参数；默认 `--storage auto` 在两个及以上输入时选择
联合索引（pfp 后端只支持两个），其余选择 source v3，联合参考默认为第一个输入。
三个及以上输入时，每个样本从它前面的所有样本中查找拷贝，建议把最完整的组装放在第一个；
`--copy-index-mib N` 把拷贝查找的 minimizer 表分块建立，降低峰值内存，索引内容不变。
任意两个样本的视图可直接比对，例如
`FastGA-shared -T24 -pafx indexes/joint3/A.gix indexes/joint3/C.gix`。四种存储格式均可选择
原生/PFP 排序。单命令建索引加比对也可用 `--index-backend pfp` 切换排序后端。

单命令建索引加比对时，`--index native` 改为建立上游 FastGA 的 GDB/GIX 文件，
默认的 `--index joint` 建立 FastGA-shared 索引。两者用同一个比对程序，PAF 逐字节
相同；上游索引占用的磁盘约为 3–5 倍，构建更快。

合并自比对复用同一物理索引：

```sh
python3 scripts/merge_fasta.py \
  --input A=genome-A.fa --input B=genome-B.fa --output merged.fa
FAtoGDB merged.fa merged.1gdb
gixshare view --index indexes/joint --samples A,B \
  --gdb merged.1gdb --output merged.gix
FastGA-shared -T24 -pafx merged.gix > self.paf
```

保留索引伴随文件以及对应 GDB/BPS 序列；使用本版本 `FastGA-shared`。
基因组或掩码变化后重新建索引。更详细的参数和运行说明见 [英文主页](README.md)。

## 当前版本与原版 FastGA 对比

下表为真实测量，数据是完整的 *Macaca nemestrina* 两个单倍型（PR00030）：
- hap1 作参考，hap2 作查询，参数 `-T24 -pafx`。
- 运行环境是一台多人共用的 Linux 主机，256 逻辑核。
- 原版是未改动的上游构建 `FastGA-baseline`；当前版本是使用最小化子锚点 copy 查找器的 `FastGA-shared`。

| 指标 | 原版 FastGA | FastGA-shared | 差异 |
|---|---:|---:|---:|
| PAF 输出 | 540,544 行，14,770,377,539 B | 两次均与原版逐字节一致（`cmp`） | 无 |
| 索引文件（不含 GDB 序列） | 64,450,238,138 B（两份 `.gix` + `.ktab.*`） | 13,034,351,452 B（全部联合索引文件） | 为原版的 20.2%，缩小 4.94 倍 |
| 端到端墙钟（建索引 + 比对） | 2,149 s | 2,286 s；2,332 s | +6.3%；+8.5% |
| 端到端 CPU（用户态 + 内核态） | 9,029 核秒（CPU 利用率 420%） | 18,186；19,165 核秒 | 2.01 倍；2.12 倍 |
| 端到端峰值内存 | 19,210,152 KiB（18.32 GiB） | 19,739,836；19,702,008 KiB（18.83；18.79 GiB） | +2.8%；+2.6% |

两种方式的 GDB 序列文件（`.bps`，两基因组合计 1,543,856,064 B）完全相同，不计入上表的索引大小。

当前版本本身的数据：
- 建索引用 584 s，其中 copy 查找 94 s、合并 445 s。
- 单独建索引时峰值内存为 15,515,908 KiB（14.80 GiB），所以端到端的峰值出现在比对阶段。
- 精确 copy 覆盖查询序列 97.60% 的碱基，剩余 3.951 亿行 residual 索引记录。

边界：
- 每行都只跑了一次，而且主机是共用的：原版运行期间负载在 85–200，当前版本两次运行期间负载在 65–161。因此墙钟差距只是近似值。
- 原版只用到 420% 的 CPU，当前版本多出的 CPU 大多落在原本空闲的核上。
- 原版的建索引时间没有和比对分开测量。
- 人类 CHM13 对 HG002 的对比（建索引和比对分开计时）还在运行，这里不报告。

## 已有真实结果与边界

这些数据来自独立整理前的原实现；公共 GDB 序列空间单独统计。

| 完整输入 | 原生两份索引 | 独立 v3 两份索引 | 联合索引，含三个视图 |
|---|---:|---:|---:|
| FastGA 上游 HAP 示例 | 1.955 GB | 0.632 GB | 0.500 GB |
| 人类 HG002 + HG005 | 65.704 GB | 23.121 GB | 11.779 GB |

人类联合索引完成了全量记录核验，合并视图包含 **4,362,338,887** 条逻辑记录。
两个方向的首轮完整种子与 PAF/CIGAR 相同。人类比对已按用户要求停止：
**尚未完成四轮速度比较，也没有完整人类合并自比对结果。**

完整 HAP 示例的双向和自比对验收均通过，四轮耗时中位数分别为 v3 的
1.178、1.145、1.053 倍。该测量使用 12 逻辑 CPU 的 macOS 主机、`-T24` 和
未控制的 OS 缓存；HAP 示例结果不能当作人类基因组速度结果。

2026-10-29 又对完整 HAP 做了**同一份联合 v4 索引的旧、新读入器对照**：四轮总耗时中位数
在 A→B、B→A、合并自比对上分别缩短 **7.4%、4.5%、6.5%**，全量记录、种子与排序后
PAF/CIGAR 一致，索引文件哈希未改变。环境为 Apple M2 Pro、12 逻辑核、16 GiB、`-T24`，
OS 缓存未控制；同轮配对比值的中位数对应缩短 7.1%、4.3%、4.6%，自比对存在单轮变慢，
所有轮次均保留。详见 [阶段耗时、内存及敏感性记录](docs/reader-optimization.md#complete-hap-reader-comparison-2026-10-29)
和 [原始计时摘要](benchmarks/reader-hap.json)。这组结果比较两个联合索引读入器版本。

[结果摘要](benchmarks/) 保留真实数字及来源摘要校验和，不包含个人机器路径、
SSH 信息、基因组 FASTA、索引大文件或二进制程序。所有生成的小数据都明确标记
为 **SIMULATED TEST DATA**，仅用于正确性测试。

方法和文件格式见 [方法说明](docs/method.md)、[格式说明](docs/format.md)，
复现与检查流程见 [reproducibility](docs/reproducibility.md)。上游及第三方许可证
按原文保留；新增代码的发行许可证尚未由仓库所有者指定。

## 引用

1. Boucher, C., Gagie, T., Kuhnle, A., Langmead, B., Manzini, G., and Mun, T.
   (2019). *Prefix-free parsing for building big BWTs*. Algorithms for Molecular
   Biology **14**, 13。[DOI：10.1186/s13015-019-0148-5](https://doi.org/10.1186/s13015-019-0148-5)。
   [Big-BWT 实现](https://gitlab.com/manzai/Big-BWT) 在本项目固定为
   `fa9a1663af10661d92a47d0e1f07218577bfd8bd`。
2. *SXGC*：[pangenome/sxgc 软件与研究仓库](https://github.com/pangenome/sxgc)。
   引用 SXGC 相关工作时记录实际使用的提交；本次提取的集成源码快照由
   [CORE_SOURCES.json](CORE_SOURCES.json) 标识。
3. Myers, G., and Zhou, C. *FastGA: A Fast Genome Aligner*。
   [源码](https://github.com/thegenemyers/FASTGA) 在本项目固定为
   `178afcbba2f79abc1686e3bd2b8a916969467e83`，另附可审阅的集成补丁。
4. *libsais*：[源码](https://github.com/IlyaGrebnov/libsais)，版本 2.12.1，
   固定为 `c284ba843dd16c159acb976f17bc19382793ea75`。

可直接复用的 BibTeX 见 [REFERENCES.bib](REFERENCES.bib)。发表实验时同时记录
本仓库提交与依赖版本；本集成尚无单独发表的论文或 DOI。

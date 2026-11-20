# FastGA Joint Index

从 FASTA 直接构建的多基因组（2–64 个）无损联合索引。完整保留第一个基因组的种子位置，
后面每个基因组用 minimizer 锚点发现的、来自前面各基因组的精确拷贝块表示共享位置，只单独存储剩余位置。
读取时恢复完整位置及出现次数，继续调用 FastGA 的播种、链化和局部比对。
比对结果与上游 FastGA 逐字节一致，索引是上游索引的 0.18–0.36 倍；人类基因组对和合并基因组自比对
比上游快 4–5 倍，小型跨物种基因组对为 0.7–1.35 倍（见[实测结果](#实测结果)）。

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

## 为什么索引小、比对快

**索引小是因为共享序列只存一次；比对快主要来自重写的比对引擎。** 联合索引完整保存第一个基因组的种子位置，
后面每个基因组都从前面所有基因组拷贝：minimizer 锚点查找器（k=21、w=20）在由前面各基因组序列组成的
伪参考上找出至少 40 bp 的精确共享片段，片段内的每个 40-mer 种子出现位置都直接继承、不再存储，
继承前逐个核对两个基因组对该种子的选取一致。只有拷贝之外的种子才写成新行，每行带一个权重，
记录它代表多少次出现；六个人类单倍型的 133.5 亿次出现只需存 29.1 亿行。索引仍是无损的：
读取时通过拷贝表把每一行还原成全部出现位置，所以下面所有比对都与上游 FastGA 逐字节一致，
每个样本解码出的行也与上游 GIXmake 完全相同。

速度主要来自 FastGA-shared 的比对引擎：共享的计算令牌池让同时计算的线程不超过 `-T`，
加载下一批种子与搜索并行而不争抢核心，种子合并本身也做了优化。即使用上游索引，FastGA-shared
也比上游 FastGA 快，最高 5.4 倍（小型基因组对为 1.0–1.6 倍）。读取联合索引有额外解码开销：
人类基因组对和合并基因组为 0–16%，比对只需 25–40 s 的小型基因组对为 2–46%。建索引时所有基因组
流经同一个写入器，不产生中间临时文件，但仍比建上游索引慢，所以代价是一次性的、更长的建索引时间，
换来小得多的索引和更快的比对。

## 实测结果

数据均为完整真实组装（登录号见 [`benchmarks/current/samples.tsv`](benchmarks/current/samples.tsv)），
在一台 256 逻辑核的共用 Linux 主机上测量，每个程序 `-T24` 并绑定 24 个核，比对参数
`-pafx -f10 -c85 -s1000 -l100 -i0.7`。上游 FastGA 为未改动的构建，FastGA-shared 为提交 `d4efe6a`。
每次比对排序后的 PAF 与上游 FastGA 的 SHA-256 和记录数都相同。每种配置、每个方向只跑一次（第 1 轮）；
同样配置的第 2 轮在 15–70 s 的小型基因组对上最多相差 61%，在人类基因组对上相差 10–23%，
所以两列 FastGA-shared 之间这种量级的差异没有意义。索引大小包括索引目录的全部文件，
其中 GDB 序列文件两种方案完全相同。上游各列测于 2026-11-16/18，FastGA-shared 各列测于 2026-11-20。
数据见 [`benchmarks/current/`](benchmarks/current/)，图见 [基因组对](docs/figures/pairs.pdf)、[合并基因组](docs/figures/merge.pdf)。

**基因组对。** 比对时间为 A→B 与 B→A 的平均值。

| 基因组对（A / B） | 上游索引 A+B | 联合索引 | 上游 FastGA | FastGA-shared，上游索引 | FastGA-shared，联合索引 | 建索引，上游 A+B | 建索引，联合 | 峰值内存，上游 / 联合 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| *W. vidua* / *P. parallelum* | 25.84 GiB | 8.41 GiB（0.33×） | 38.6 s | 25.6 s | 28.7 s | 44.2 s | 93.7 s | 3.0 / 5.9 GiB |
| *P. vaucheri* / *A. gentilis* | 31.40 GiB | 10.24 GiB（0.33×） | 30.5 s | 26.5 s | 29.9 s | 58.8 s | 153.7 s | 4.7 / 7.3 GiB |
| *G. mutica* / *A. ochrocephala* | 37.70 GiB | 12.14 GiB（0.32×） | 41.2 s | 25.6 s | 37.4 s | 69.2 s | 134.6 s | 11.0 / 9.2 GiB |
| *A. centrarchus* / *P. auritus* | 31.02 GiB | 10.66 GiB（0.34×） | 23.5 s | 23.7 s | 32.9 s | 58.2 s | 121.9 s | 2.4 / 6.7 GiB |
| *C. striatus* / *E. lucius* | 22.78 GiB | 7.15 GiB（0.31×） | 30.3 s | 24.6 s | 25.0 s | 52.0 s | 91.7 s | 6.1 / 5.3 GiB |
| HG002 / HG005（人类，父系） | 62.56 GiB | 11.47 GiB（0.18×） | 5,426 s | 997 s | 1,125 s | 114.6 s | 274.7 s | 11.6 / 14.3 GiB |

六个联合索引是同时构建的，构建时间包含相互争用内存带宽和 `/dev/shm` 的影响。

**合并基因组自比对。** N 个跨物种基因组 M1–MN（2.19、2.16、1.48、1.27、1.79、1.46 Gbp）合并成一个 FASTA
后自比对；联合索引是一个 N 样本索引加一个合并视图。

| N | N 个独立上游索引 | 合并基因组的上游索引 | 联合索引 | 上游 FastGA | FastGA-shared，上游索引 | FastGA-shared，联合索引 | 建索引：N 个上游 / 合并上游 / 联合 | 峰值内存，上游 / 联合 |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 3 | 61.68 GiB | 62.78 GiB | 20.62 GiB（0.33×） | 9,304 s | 1,903 s | 1,884 s（4.9×） | 136 / 85 / 264 s | 18.2 / 18.0 GiB |
| 4 | 75.35 GiB | 76.40 GiB | 26.23 GiB（0.35×） | 10,921 s | 2,110 s | 2,446 s（4.5×） | 161 / 174 / 552 s | 18.7 / 22.5 GiB |
| 5 | 93.33 GiB | 95.60 GiB | 33.11 GiB（0.35×） | 11,113 s | 2,071 s | 2,372 s（4.7×） | 196 / 162 / 791 s | 23.6 / 27.6 GiB |
| 6 | 108.27 GiB | 111.52 GiB | 39.37 GiB（0.36×） | 10,449 s | 2,365 s | 2,594 s（4.0×） | 224 / 288 / 802 s | 23.5 / 33.9 GiB |

每个 N 的两次 FastGA-shared 自比对同时启动。联合索引的峰值内存在 N=4–6 时来自建索引，N=3 时来自比对；
上游的峰值来自比对。

**六个人类单倍型合建一个索引。** HPRC Year 1 的 HG002、HG005、HG00438 父系与母系单倍型
（[`benchmarks/current/hprc6.json`](benchmarks/current/hprc6.json)）：

| 指标 | 数值 |
|---|---:|
| 建索引（`GIXmake-shared -T24`，六个一起） | 683 s，峰值内存 24.3 GiB |
| 索引大小 | 15.66 GB（14.59 GiB） |
| 存储行数 / 出现次数 | 2,908,024,513 / 13,345,046,896（21.8%） |
| 每个样本用 `gixshare verify-equivalent` 与上游 GIXmake 对照 | 全部通过，共 13,345,046,896 行 |
| 从该索引比对 HG002P → HG005P | 1,157 s，峰值内存 9.4 GiB，492,284 条记录与上游一致 |

后面的单倍型有 91–95% 的出现位置是继承来的。作为对照，一个人类单倍型的上游索引约 31 GiB（见上面的 HG002/HG005）。

[`benchmarks/current/`](benchmarks/current/) 保存上面各表的数据，不包含个人机器路径、
基因组 FASTA、索引大文件或二进制程序。所有生成的小数据都明确标记
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

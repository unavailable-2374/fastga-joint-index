# FastGA Joint Index

从 FASTA 直接构建的双基因组无损联合索引。完整保留参考基因组的种子位置，
用 PFP 发现的精确拷贝块表示另一基因组的共享位置，只单独存储剩余位置。
读取时恢复完整位置及出现次数，继续调用 FastGA 的播种、链化和局部比对。

核心 C/C++ 代码与已验证版本逐字节一致；本仓库将其独立构建、测试和文档化，
不再依赖原 SXGC 大仓库。当前为实验性方法，联合模式固定 k=40、两个输入。

## 主要来源与归属

**本项目的主要来源是 PFP（Prefix-Free Parsing）和
[SXGC](https://github.com/pangenome/sxgc)。** PFP 提供短语解析，用于发现共享精确
拷贝的锚点；共享索引代码最初在 SXGC 工作区中作为 FastGA 集成开发，再整理为
本独立仓库。[CORE_SOURCES.json](CORE_SOURCES.json) 记录了提取的 C/C++ 源码校验和。

| 来源 | 在本项目中的作用 |
|---|---|
| [PFP：Boucher 等，2019](https://doi.org/10.1186/s13015-019-0148-5)；实现为 [Big-BWT](https://gitlab.com/manzai/Big-BWT) | 短语解析和共享拷贝锚点；可选的 PFP 后缀排序后端 |
| [SXGC](https://github.com/pangenome/sxgc) | 共享索引集成、构建和验证流程的原始研发工作区 |
| [FastGA：Myers 与 Zhou](https://github.com/thegenemyers/FASTGA) | 原有采样、播种、链化、局部比对、GDB 和输出格式 |
| [libsais](https://github.com/IlyaGrebnov/libsais) | 短语 token 的后缀数组和 LCP 构建 |

本仓库整理了 FastGA 兼容的共享种子表示、精确拷贝与残余位置编码、读取适配及
单命令入口。SXGC 的 suffixient-array 形式化证明适用于其自身结构；本适配器以
全量记录、种子对和比对结果测试提供正确性证据。使用时请标注 PFP、SXGC 两项
主要来源，并同时引用 FastGA 及实际使用的上游实现。完整引用见文末及
[第三方说明](THIRD_PARTY_NOTICES.md)。

## 构建与使用

依赖 Linux/macOS、CMake ≥3.20、Python ≥3.11、C11/C++17 编译器、zlib、pthread、patch。
上游源代码已随仓库固定版本保存，构建不需要下载依赖。

```sh
python3 scripts/build.py --jobs 4
export PATH="$PWD/build:$PATH"
ctest --test-dir build --output-on-failure

FastGA-shared -T24 -pafx A.fa B.fa > A_B.paf
FastGA-shared -T24 -pafx A.fa > A.paf
```

以上直接读取你的 FASTA：两个输入自动构建 v4 联合索引并比对；一个输入构建
v3 独立压缩索引并自比对，也支持 `.fna`、`.fasta` 及其 gzip 文件。
建索引日志写入 stderr，stdout 只有比对结果。默认结束后清理临时索引；加 `-k`
会在当前目录的 `.fastga-shared/` 保留索引，并按完整输入内容、掩码、参数和
构建工具校验和复用，多个进程同时调用时只构建一份。`-P<目录>` 指定临时空间，
`-T1` 到 `-T32` 均可使用；低线程时顺序写入完整种子流以遵守预算。
首次运行的总耗时包含建索引。`-1:` 输出会自动保留所需 GDB/索引。

原有 `.gix` 参数直接交给原比对内核。移动工具时请保留整个 `build/`，包括
`fastga-runtime/`；完整说明见[命令行文档](docs/cli.md)。

也可以继续显式构建、保存索引后比对：

```sh
mkdir -p indexes
python3 scripts/fasta_index.py \
  --storage joint --reference A --threads 24 \
  --input A=genome-A.fa --input B=genome-B.fa \
  --output indexes/joint
FastGA-shared -T24 -pafx indexes/joint/A.gix indexes/joint/B.gix > A-to-B.paf
```

`genome-A.fa`、`genome-B.fa` 指你自己的输入，也支持 `.fa.gz`；输出目录必须不存在。
`--threads` 是构建总预算，范围为 1–32；低线程时原生排序采用顺序种子流，
PFP 排序至少需要 2 个线程。默认使用原生排序；添加
`--backend pfp` 可替换排序阶段。两种模式均使用 PFP 查找共享拷贝块。
`build/` 内记录真实编译命令和源文件/可执行文件校验和。

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

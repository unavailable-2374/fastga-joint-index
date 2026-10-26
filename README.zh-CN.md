# FastGA Joint Index

从 FASTA 直接构建的双基因组无损联合索引。完整保留参考基因组的种子位置，
用 PFP 发现的精确拷贝块表示另一基因组的共享位置，只单独存储剩余位置。
读取时恢复完整位置及出现次数，继续调用 FastGA 的播种、链化和局部比对。

核心 C/C++ 代码与已验证版本逐字节一致；本仓库将其独立构建、测试和文档化，
不再依赖原 SXGC 大仓库。当前为实验性方法，固定 k=40、两个输入。

## 构建与使用

依赖 Linux/macOS、CMake ≥3.20、Python ≥3.11、C11/C++17 编译器、zlib、pthread、patch。
上游源代码已随仓库固定版本保存，构建不需要下载依赖。

```sh
python3 scripts/build.py --jobs 4
export PATH="$PWD/build:$PATH"
ctest --test-dir build --output-on-failure

mkdir -p indexes
python3 scripts/fasta_index.py \
  --storage joint --reference A --threads 24 \
  --input A=genome-A.fa --input B=genome-B.fa \
  --output indexes/joint
FastGA-shared -T24 -pafx indexes/joint/A.gix indexes/joint/B.gix > A-to-B.paf
```

`genome-A.fa`、`genome-B.fa` 指你自己的输入，也支持 `.fa.gz`；输出目录必须不存在。
`--threads` 是构建总预算，双输入时范围为 3–32。默认使用原生排序；添加
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

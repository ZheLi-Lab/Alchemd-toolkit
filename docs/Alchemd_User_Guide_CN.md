# Alchemd User Guide v1.2

---
**许可证**: 本文档属于 Alchemd 的一部分，采用 CC BY-NC-SA 4.0 许可证用于学术用途。商业用途需要单独授权。完整条款请参见 [LICENSE.md](../LICENSE.md)。
---

## 目录
1. [简介](#简介)
2. [快速开始](#快速开始)
3. [核心模块详解](#核心模块详解)
4. [完整工作流程](#完整工作流程)
5. [常见问题](#常见问题)

## 简介

Alchemd 是一个基于 Python 的 CADD（计算辅助药物设计）应用程序，专门用于炼金术变换的自由能微扰（FEP）计算。该工具包提供了从输入准备到分子动力学模拟再到结果分析的完整工作流程。

### 项目结构

```
Alchemd/
├── core/                    # 核心模块
│   ├── prepare_file.py     # 输入文件准备
│   ├── run_md.py          # 分子动力学模拟
│   └── analyze_result.py  # 结果分析
├── gui/                    # 图形界面
├── dependencies/           # 依赖工具
└── docs/                   # 文档
```

### 计算方法
- **Amber**: 分子动力学引擎
- **CAR (Convergence-Adaptive Roundtrip Method)**: 自适应收敛方法，默认启用

## 快速开始

### 标准工作流程

```bash
# 1. 准备输入文件
python core/prepare_file.py -auto -p protein.pdb --cs-align -d -l ligands.sdf -o ./output

# 2. 运行MD计算（使用CAR方法）
python core/run_md.py -ln 50 -p ./output/prepare -auto

# 3. 分析结果
python core/analyze_result.py -auto -w ./output/run
```

## 核心模块详解

## 1. prepare_file.py - 输入文件准备

### 基本语法
```bash
python core/prepare_file.py [参数]
```

### 主要参数

#### 必需参数

**`-p, --protein PROTEIN`**
- **作用**: 指定蛋白质PDB文件路径
- **类型**: 字符串
- **格式**: .pdb
- **示例**: `-p protein.pdb`
- **说明**: 蛋白质结构文件，将用作受体分子

**`-o, --output-dir OUTPUT_DIR`**  
- **作用**: 指定输出目录
- **类型**: 字符串
- **示例**: `-o ./output`
- **说明**: 所有处理文件的保存位置，必须指定

#### 配体输入方式（二选一）

**方式一：从SDF文件输入配体**

**`-l, --ligands LIGANDS`**
- **作用**: 指定配体SDF文件路径
- **类型**: 字符串
- **格式**: .sdf
- **示例**: `-l ligands.sdf`
- **要求**: 
  - 必须包含键合信息
  - 必须包含显式和隐式氢原子
  - 每个分子需有唯一的名称

**方式二：从现有文件夹加载**

**`-af, --advance-folder ADVANCE_FOLDER`**
- **作用**: 指定已组织好的prepare文件夹路径
- **类型**: 字符串
- **示例**: `-af ./existing_prepare`
- **说明**: 适用于已有pair分子PDB文件的情况

#### 配体筛选与配对

**`-r, --ref-ligand-list REF_LIGAND_LIST`**
- **作用**: 指定参考配体列表文件
- **类型**: 字符串
- **格式**: .lst
- **示例**: `-r ref_ligands.lst`
- **说明**: 
  - 每行一个配体名称
  - 名称需与SDF文件中的分子名称对应
  - 不指定时，所有配体均作为参考配体

**`-pl, --pair-list PAIR_LIST`**
- **作用**: 指定配体对列表文件
- **类型**: 字符串
- **格式**: .lst
- **示例**: `-pl pairs.lst`
- **说明**: 
  - 每行一个配体对，格式：`ligand1-ligand2`
  - 不指定时，使用LOMAP2自动生成配体对

#### 分子处理选项

**`--cs-align`**
- **作用**: 在分子生成时启用cs-align叠合分子
- **类型**: 布尔值
- **默认**: False
- **建议**: **强烈推荐启用**
- **说明**: 从参考分子生成配对分子的3D构象

**`-d, --dock`**
- **作用**: 启用分子对接
- **类型**: 布尔值  
- **默认**: False
- **建议**: **推荐启用**
- **说明**: 在分子生成后进行蛋白-配体对接

**`-npp, --not-prepare-protein`**
- **作用**: 跳过蛋白质预处理
- **类型**: 布尔值
- **默认**: False
- **建议**: 通常保持默认（即进行蛋白质预处理）

#### 配体对处理选项

**`-sp, --strict-pair-list`**
- **作用**: 严格遵循用户定义的配体对列表
- **类型**: 布尔值
- **默认**: False
- **条件**: 仅在同时使用 `--cs-align` 和 `-p` 时有效
- **说明**: 
  - 启用时：严格按照给定顺序处理配体对
  - 禁用时：自动优化配体对的处理顺序

**`-f, --flip`**
- **作用**: 设置FEP炼金转化方向，从大分子转化到小分子
- **类型**: 布尔值
- **默认**: False
- **说明**: 根据范德华体积自动调整配体对方向，确保从大分子到小分子的炼金转化

**`-as, --allow-sep`**
- **作用**: 允许多个不连通的配体对组
- **类型**: 布尔值
- **默认**: False
- **说明**: 处理存在多个独立配体对网络的情况

#### 特殊氨基酸支持

**`-npd, --nsaa-params-dir NSAA_PARAMS_DIR`**
- **作用**: 指定包含预生成的非标准氨基酸力场参数的目录
- **类型**: 字符串
- **示例**: `-npd ./nsaa_parameters`
- **说明**: 包含预生成的.frcmod和.prepi文件的目录路径。参数必须事先使用nsaa_preparer.py生成
- **重要**: 该参数使用现有参数，不执行参数生成。参数生成是独立的预处理步骤
- **详细配置**: 请参见[NSAA准备模块](#nsaa准备模块)章节了解完整的配置步骤和Gaussian依赖设置

#### 自动化选项

**`-auto`**
- **作用**: 自动准备拓扑和坐标文件
- **类型**: 布尔值
- **默认**: False
- **建议**: **推荐启用**
- **说明**: 在准备配体对文件后自动执行准备步骤（而非生成一个任务队列）

### 输出文件

执行成功后将生成以下目录和文件：

**`preprocess/`** - 预处理目录
- 存储蛋白质和配体的预处理文件
- 包含分子对接和坐标生成的中间文件

**`prepare/`** - 准备工作目录  
- 每个配体对一个子目录（如 `A-B/`）
- 每个子目录包含：
  - `A.pdb` - 参考配体PDB文件
  - `B.pdb` - 目标配体PDB文件
  - `submit.sh` - 准备脚本

**`prepare_queue.lst`** - 准备任务队列
- 每行包含工作目录和执行脚本路径
- 格式：`/path/to/pair/dir /path/to/submit.sh`

**`info.pickle`** - 元信息文件
- 包含配体对信息和分子数据

### 使用示例

**基本使用**：
```bash
python core/prepare_file.py -p protein.pdb -l ligands.sdf -o ./output --cs-align -d -auto
```

**指定参考配体**：
```bash
python core/prepare_file.py -p protein.pdb -l ligands.sdf -r ref_list.lst -o ./output --cs-align -d -auto
```

**使用自定义配体对**：
```bash
python core/prepare_file.py -p protein.pdb -l ligands.sdf -pl pairs.lst -sp -o ./output --cs-align -d -auto
```

### 输入文件格式详解

本节详细说明 `prepare_file.py` 中三个**重要可选**参数的文件格式要求和使用方法。

#### pairs.lst (--pair-list) 格式规范

**文件作用**：手动指定配体对列表，替代自动的LOMAP2配对生成

**文件格式**：
- 纯文本文件，推荐使用 `.txt` 或 `.lst` 扩展名
- 每行包含一个配体对，格式为：`ligandA-ligandB`
- 使用连字符 `-` 分隔两个配体名称
- 支持空行和注释（以 `#` 开头的行将被忽略）
- 文件编码：推荐使用UTF-8编码

**示例文件内容**：
```
# 配体对列表示例
mol1-mol2
mol2-mol3
mol3-mol4
mol1-mol4
mol5-mol6
```

**重要要求**：
- 配体名称必须与SDF文件中分子的 `_Name` 属性完全匹配
- 配体对具有双向对称性：`mol1-mol2` 与 `mol2-mol1` 被视为相同配对
- 与 `--strict-pair-list` 参数配合使用时，严格按照文件中的顺序处理配体对
- 每个配体对仅需指定一次，程序会自动处理双向关系

#### ref.lst (--ref-ligand-list) 格式规范

**文件作用**：指定作为参考配体的分子列表，用于坐标系对齐和结构生成

**文件格式**：
- 纯文本文件，推荐使用 `.txt` 或 `.lst` 扩展名
- 每行包含一个配体名称
- 支持空行和注释（以 `#` 开头的行将被忽略）
- 文件编码：推荐使用UTF-8编码

**示例文件内容**：
```
# 参考配体列表示例
mol1
mol3
mol5
```

**重要要求**：
- 配体名称必须与SDF文件中分子的 `_Name` 属性完全匹配
- 参考配体用于提供可靠的3D构象，建议选择具有实验结构或高质量构象的分子
- 未列在此文件中的分子将被标记为活性配体，需要进行构象生成和优化
- 如果不提供此文件，所有分子都被视为参考配体

#### 文件夹准备模式 (--advance-folder) 格式规范

**功能作用**：使用预处理好的配体对PDB文件，跳过自动的分子对齐和构象生成步骤

**目录结构要求**：
```
advance_folder/
├── mol1-mol2/
│   ├── mol1.pdb
│   └── mol2.pdb
├── mol2-mol3/
│   ├── mol2.pdb
│   └── mol3.pdb
├── mol3-mol4/
│   ├── mol3.pdb
│   └── mol4.pdb
└── mol5-mol6/
    ├── mol5.pdb
    └── mol6.pdb
```

**重要要求**：
- 每个子文件夹名称必须遵循 `ligandA-ligandB` 格式（必须包含连字符 `-`）
- 每个子文件夹内必须包含对应的两个PDB文件
- PDB文件名必须与配对名称中的分子名称完全匹配
- PDB文件必须包含完整的原子坐标信息和正确的残基信息
- 系统会自动从PDB文件中提取分子信息和电荷信息

**使用场景**：
- 已经通过外部工具（如分子对接、量子化学计算）生成了高质量的分子构象
- 需要使用特定的分子取向或构象进行FEP计算
- 希望跳过自动构象生成步骤以节省计算时间

**注意事项**：
- 该模式跳过了所有分子预处理步骤，用户需确保输入的PDB文件质量
- 不会执行分子对接或坐标系对齐，直接使用提供的坐标

#### 格式验证建议

**文件编码**：建议使用UTF-8编码保存所有文本文件

**路径使用**：
- 建议使用绝对路径避免路径解析问题
- 确保所有引用的文件和目录都存在且可访问
- 避免在文件路径和分子名称中使用特殊字符（如空格、中文字符等）

**命名规范**：
- 分子名称应使用英文字母、数字和下划线，避免使用连字符以外的特殊符号
- 文件夹名称遵循系统命名规范，避免过长路径

**命名一致性检查**：
```bash
# 验证SDF文件中的分子名称
grep "_Name" ligands.sdf

# 检查配体对文件格式
cat pairs.lst | grep -v "^#" | grep -v "^$"

# 验证参考配体列表
cat ref.lst | grep -v "^#" | grep -v "^$"

# 检查文件是否存在特殊字符
file -bi pairs.lst  # 检查文件编码
```

**常见错误及解决方案**：

1. **分子名称不匹配**：检查SDF文件中的`_Name`属性是否与列表文件中的名称完全一致
2. **文件格式错误**：确保使用UTF-8编码，检查是否存在隐藏字符
3. **路径问题**：使用绝对路径并确保文件存在
4. **权限问题**：确保文件具有可读权限

---

## 2. run_md.py - 分子动力学模拟

### 基本语法
```bash
python core/run_md.py [参数]
```

### 主要参数

#### 必需参数

**`-p, --preparation_dir PREPARATION_DIR`**
- **作用**: 指定prepare目录路径
- **类型**: 字符串
- **示例**: `-p ./output/prepare`
- **说明**: 
  - 必须是完成prepare_file.py处理的目录
  - 确保已执行完prepare_queue.lst中的所有任务

#### Lambda参数配置

**`-ln, --lambda_nums LAMBDA_NUMS`**
- **作用**: 设置Lambda窗口数量
- **类型**: 整数
- **默认**: 50
- **示例**: `-ln 50`
- **说明**: 
  - 实际Lambda点数为 `lambda_nums + 1`
  - 例如：`-ln 50` 产生51个Lambda值 (1.0, 0.98, ..., 0.02, 0.0)
  - 推荐值：50（对应Lambda间隔0.02）

#### 计算引擎配置

**`-m, --mode MODE`**
- **作用**: 指定分子动力学引擎
- **类型**: 字符串
- **默认**: `Amber`
- **示例**: `-m Amber`
- **说明**: 使用Amber引擎配合CAR方法（默认设置）

#### CAR自适应收敛方法配置

**`-ci, --CAR_input CAR_INPUT`**
- **作用**: 指定自定义CAR输入文件
- **类型**: 字符串
- **格式**: .txt
- **示例**: `-ci custom_car_settings.txt`
- **说明**: 不指定时使用默认CAR配置模板

#### 自动化选项

**`-auto`**
- **作用**: 准备完成后自动运行MD
- **类型**: 布尔值
- **默认**: False
- **建议**: **推荐启用**
- **说明**: 在本地环境自动**非并行**执行所有MD任务

### 输出文件

执行成功后将在preparation_dir的父目录生成：

**`run/`** - MD工作目录
- 每个配体对一个子目录
- 每个配体对包含4个边：
  - `lM2A/` - 溶液中M→A炼金转化
  - `lM2B/` - 溶液中M→B炼金转化  
  - `cM2A/` - 复合物中M→A炼金转化
  - `cM2B/` - 复合物中M→B炼金转化

**每个边目录包含**：
- 拓扑文件（`.prmtop`）
- 坐标文件（`.prmcrd`）
- Lambda配置文件（`.json`）
- 执行脚本（`submit.sh`）
- 输入配置文件

**`run_all_queue.lst`** - MD任务队列
- 包含所有MD计算任务
- 格式：`/path/to/edge/dir /path/to/submit.sh`

### 使用示例

**标准模式**：
```bash
python core/run_md.py -ln 50 -p ./output/prepare -auto
```

**自定义CAR配置**：
```bash
python core/run_md.py -ln 50 -p ./output/prepare -ci my_car_config.txt -auto
```

---

## 3. analyze_result.py - 结果分析

### 基本语法
```bash
python core/analyze_result.py [参数]
```

### 主要参数

#### 必需参数

**`-w, --work_dir WORK_DIR`**
- **作用**: 指定MD结果目录路径
- **类型**: 字符串
- **默认**: `./work`
- **示例**: `-w ./output/run`
- **说明**: run_md.py生成的run目录路径

#### 计算方法配置

**`-m, --mode MODE`**
- **作用**: 指定MD引擎类型
- **类型**: 字符串
- **默认**: `Amber`
- **说明**: 必须与run_md.py中的设置保持一致

#### 自动化选项

**`-auto`**
- **作用**: 自动计算能量并进行统计分析
- **类型**: 布尔值
- **默认**: False
- **建议**: **推荐启用**
- **说明**: 
  - 自动准备计算文件
  - 并行执行能量计算
  - 生成统计分析报告

### 并行处理配置

程序内部使用以下配置进行并行处理：
- **最大工作进程数**: `min(cpu_count, 4)`
- **最大重试次数**: 3次
- **任务超时时间**: 600秒

**自定义并行配置**：用户可根据计算资源修改脚本中的 `MAX_WORKER_COUNT` 参数以优化性能。该参数统一控制文件准备和计算任务的最大并行进程数。

### 输出文件

执行成功后将在work_dir的父目录生成：

**`calculate/`** - 计算工作目录
- 包含所有配体对的能量计算文件
- 每个计算单元的输入输出文件

**`results/`** - 分析结果目录
- **`results.csv`**: 所有配体对的能量分析结果
  - 包含ΔΔG值、不确定度、收敛性指标
- **`ene_detail.csv`**: 详细的能量分解分析
  - 包含各个Lambda窗口的详细能量数据
- **`summary.pdf`**: 完整分析报告
  - 包含能量图表、收敛性分析、统计信息
- **`summary_abstract.pdf`**: 简化分析报告
  - 关键结果和图表的摘要版本

**可选输出文件**：
- **`sus_pairs.lst`**: 可疑配体对列表
  - 包含收敛性差或异常结果的配体对
- **`err_pairs.lst`**: 错误配体对列表
  - 计算失败的配体对列表
- **`calculate_queue.lst`**: 计算任务队列

### 使用示例

**标准分析流程**：
```bash
python core/analyze_result.py -auto -w ./output/run
```

**手动控制分析流程**：
```bash
# 仅准备计算文件
python core/analyze_result.py -w ./output/run

# 之后可以手动执行计算和统计
```

---

## NSAA准备模块

### 模块概述

非标准氨基酸（NSAA）准备模块是Alchemd工具包中的独立预处理工具，通过量子力学计算为非标准氨基酸生成力场参数。该模块作为独立的预处理步骤运行，完全独立于主要的FEP计算工作流程。

**模块位置**: `core/prepare_input_tools/nsaa_preparer.py`

**关键特征**:
- 独立的预处理工具
- 需要在主工作流程之前单独执行
- 生成基于量子力学的力场参数
- 输出结构化参数目录供prepare_file.py使用

### 计算依赖

NSAA模块需要多个外部计算工具：

**必需依赖**：
- **Gaussian** (g03/g09/g16): 电子结构量子化学计算
- **pyautomd**: 自动化分子动力学参数生成接口
- **AmberTools**: antechamber、parmchk2、tleap用于力场参数处理
  - **重要**: 与只需要PMEMD的主安装指南不同，NSAA功能需要完整的AmberTools编译
  - 用户必须单独编译AmberTools，独立于主要的PMEMD安装
  - AmberTools编译包含必要的参数生成工具
- **GAUSS_SCRDIR**: 指向Gaussian临时文件目录的环境变量

### 两步NSAA工作流程

NSAA处理遵循强制性的两步程序，无法绕过或集成：

#### 第一步：独立的NSAA参数生成

```bash
python core/prepare_input_tools/nsaa_preparer.py protein_with_nsaa.pdb -o ./nsaa_parameters
```

**处理描述**：
1. **结构分析**: 扫描PDB文件识别非标准氨基酸残基
2. **残基提取**: 提取NSAA及邻近残基以保持上下文
3. **电荷计算**: 使用RDKit分子描述符计算形式电荷
4. **量子计算**: 调用pyautomd → Gaussian进行电子结构计算
5. **参数生成**: 产生.frcmod（力场修饰）和.prepi（残基拓扑）文件
6. **目录组织**: 输出参数到`nsaa_parameters/NSAA_NAME/`结构

**输出结构**：
```
nsaa_parameters/
├── NSAA_NAME_1/
│   ├── NSAA_NAME_1.frcmod
│   └── NSAA_NAME_1.prepi
└── NSAA_NAME_2/
    ├── NSAA_NAME_2.frcmod
    └── NSAA_NAME_2.prepi
```

#### 第二步：主工作流程中的参数使用

```bash
python core/prepare_file.py -p protein.pdb -l ligands.sdf -npd ./nsaa_parameters -o ./output
```

### 配置要求

#### 环境配置

**Gaussian环境设置**：
```bash
# 设置Gaussian暂存目录
export GAUSS_SCRDIR="/path/to/gaussian/scratch"

# 确保Gaussian可执行文件可访问
export g16root="/path/to/gaussian"
source $g16root/g16/bsd/g16.profile
```

#### 模板修改

NSAA准备需要模板配置以实现Gaussian集成：

**模板文件**: `core/templates/nsaa_prep.txt`

**必需修改**：
- 第11行：指定Gaussian版本（g03/g09/g16）
- 临时文件管理的暂存目录路径
- 计算资源分配参数

### 使用示例

#### 完整的两步过程

```bash
# 第一步：生成NSAA参数（独立预处理）
python core/prepare_input_tools/nsaa_preparer.py \
    protein_with_nsaa.pdb \
    -o ./nsaa_parameters \
    --no-delete  # 保留临时文件用于调试

# 第二步：在主工作流程中使用参数
python core/prepare_file.py \
    -auto \
    -p protein_with_nsaa.pdb \
    -l ligands.sdf \
    -npd ./nsaa_parameters \
    --cs-align -d \
    -o ./fep_project
```

#### 多NSAA系统

```bash
# 处理含有NSAAs的不同蛋白质
python core/prepare_input_tools/nsaa_preparer.py protein_A.pdb -o ./nsaa_params_A
python core/prepare_input_tools/nsaa_preparer.py protein_B.pdb -o ./nsaa_params_B

# 在各自的工作流程中使用参数
python core/prepare_file.py -p protein_A.pdb -l ligands_A.sdf -npd ./nsaa_params_A -o ./project_A
python core/prepare_file.py -p protein_B.pdb -l ligands_B.sdf -npd ./nsaa_params_B -o ./project_B
```

### 故障排除指南

#### 常见问题及解决方案

**问题**: "Environment variable GAUSS_SCRDIR is not set"
- **解决方案**: 配置Gaussian暂存目录：`export GAUSS_SCRDIR="/tmp/gaussian"`

**问题**: "pyautomd command not found"
- **解决方案**: 确保pyautomd已安装并在PATH中可访问

**问题**: "Gaussian calculation timeout"
- **解决方案**: 增加超时值或优化计算资源

**问题**: "No non-standard amino acids found"
- **解决方案**: 验证PDB文件包含不在标准列表中的实际非标准残基

#### 验证程序

**参数验证**：
1. 验证.frcmod文件包含合理的力常数
2. 检查.prepi文件的原子连接性
3. 使用小型tleap验证运行测试参数
4. 将生成的电荷与预期化学行为进行比较

### 技术限制

**当前限制**：
- 肽链连接的NSAAs需要手动封端
- 默认使用AM1-BCC电荷（RESP电荷需要额外的QM计算）
- 仅支持序列处理（无并行NSAA参数生成）
- Gaussian依赖限制了跨平台兼容性

**高级用户建议**：
- 考虑RESP电荷计算以提高精度
- 在可用时根据实验数据验证参数
- 参数生成前进行几何优化
- 电荷计算使用适当的溶剂化模型

---

## 完整工作流程

### 标准工作流程

**第一步：准备输入文件**
```bash
python core/prepare_file.py \
    -auto \
    -p protein.pdb \
    --cs-align \
    -d \
    -l ligands.sdf \
    -o ./my_project
```

**第二步：运行MD计算（使用CAR方法）**
```bash
python core/run_md.py \
    -ln 50 \
    -p ./my_project/prepare \
    -auto
```

**第三步：分析结果**
```bash
python core/analyze_result.py \
    -auto \
    -w ./my_project/run
```

### 高级工作流程选项

**使用自定义配体对**：
```bash
# 创建pairs.lst文件，每行一个配体对
echo "ligand1-ligand2" > pairs.lst
echo "ligand2-ligand3" >> pairs.lst

python core/prepare_file.py \
    -auto \
    -p protein.pdb \
    --cs-align \
    -d \
    -l ligands.sdf \
    -pl pairs.lst \
    -sp \
    -o ./custom_pairs_project
```

**处理非标准氨基酸**：
```bash
python core/prepare_file.py \
    -auto \
    -p protein_with_nsaa.pdb \
    --cs-align \
    -d \
    -l ligands.sdf \
    -npd ./nsaa_parameters \
    -o ./nsaa_project
```

**批处理模式（仅准备，不自动执行）**：
```bash
# 准备文件但不自动执行
python core/prepare_file.py \
    -p protein.pdb \
    --cs-align \
    -d \
    -l ligands.sdf \
    -o ./batch_project

# 准备MD但不自动执行
python core/run_md.py \
    -ln 50 \
    -p ./batch_project/prepare

# 分析但不自动计算
python core/analyze_result.py \
    -w ./batch_project/run
```

### 文件组织结构

完整项目的文件结构如下：
```
my_project/
├── preprocess/              # 预处理文件
│   ├── all_generate_mol.sdf # 所有生成分子的SDF文件
│   ├── dock.pdb            # 用于分子对接的蛋白质结构
│   ├── gen_file/           # 配体对生成详细信息文件
│   ├── gen_mol/            # 生成的分子构象文件
│   ├── given_mol/          # 参考分子mol文件
│   ├── pair/               # 配体对PDB文件
│   ├── protein.pdb         # 预处理后的蛋白质
│   └── rec_prepare/        # 蛋白质预处理过程文件
├── prepare/                # 准备工作目录
│   ├── ligand1-ligand2/   # 配体对目录
│   │   ├── ligand1.pdb
│   │   ├── ligand2.pdb
│   │   └── submit.sh
│   └── ...
├── run/                   # MD工作目录
│   ├── ligand1-ligand2/   # 配体对计算目录
│   │   ├── lM2A/         # 溶液中M→A炼金转化
│   │   ├── lM2B/         # 溶液中M→B炼金转化
│   │   ├── cM2A/         # 复合物中M→A炼金转化
│   │   └── cM2B/         # 复合物中M→B炼金转化
│   └── ...
├── calculate/             # 计算分析目录
├── results/               # 最终结果
│   ├── results.csv
│   ├── ene_detail.csv
│   ├── summary.pdf
│   └── summary_abstract.pdf
├── prepare_queue.lst      # 准备任务队列
├── run_all_queue.lst     # MD任务队列
├── calculate_queue.lst   # 计算任务队列
└── info.pickle          # 元信息文件
```

## 常见问题

### Q1: 如何选择合适的Lambda数量？
**A**: 推荐使用 `-ln 50`（51个Lambda点），对应Lambda间隔0.02。复杂变换可增至100，简单变换可减至25。

### Q2: 为什么必须使用CAR方法？
**A**: CAR（Convergence-Adaptive Roundtrip Method，自适应收敛方法）是Alchemd的核心技术，提供优异的采样效率和收敛性，确保FEP计算的可靠性。

### Q3: 如何处理计算失败的配体对？
**A**: 检查 `err_pairs.lst` 文件和日志（`run.log`, `calculate.log`）。常见问题包括分子结构缺陷、缺少氢原子或计算资源不足。

### Q4: 如何调整并行性能？
**A**: 在 `analyze_result.py` 脚本中修改 `MAX_WORKER_COUNT` 参数以适配计算资源。该参数控制并行进程池的最大工作进程数。

---

## 许可证和引用

本文档是 Alchemd 的一部分。许可证信息和引用要求请参考 [LICENSE.md](../LICENSE.md)。

### 引用格式

如果您在研究中使用了 Alchemd，请引用：

```bibtex
@software{alchemd2025,
  title={Alchemd: A toolkit for alchemical free energy perturbation calculations},
  year={2025},
  url={https://github.com/ZheLi-Lab/Alchemd-toolkit}
}
```

**方法学参考**：
- CAR 方法: https://pubs.acs.org/doi/10.1021/acs.jctc.4c00939
- cs-fep 实现: https://doi.org/10.1016/j.apsb.2024.06.021

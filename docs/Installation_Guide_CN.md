# Alchemd 安装指南

---
**许可证**: 本文档属于 Alchemd 的一部分，采用 CC BY-NC-SA 4.0 许可证用于学术用途。商业用途需要单独授权。完整条款请参见 [LICENSE.md](../LICENSE.md)。
---

[![License](https://img.shields.io/badge/license-CC%20BY--NC--SA--4.0-blue.svg)](../LICENSE.md)
[![Platform](https://img.shields.io/badge/platform-Linux-lightgrey.svg)](#system-requirements)
[![Python](https://img.shields.io/badge/python-3.9+-blue.svg)](#python-environment)
[![CUDA](https://img.shields.io/badge/CUDA-12.0-green.svg)](#cuda-installation)

> 本指南包含 Alchemd 的完整安装方法，包括自动化安装和详细的手动安装步骤

## 概述

本指南提供 Alchemd 的完整安装方法。推荐使用自动化脚本进行快速部署，也可根据特定需求选择手动安装进行定制化配置。

**注意**: 本指南已针对简化的 PMEMD 单组件架构进行更新，不再需要 AmberTools。

---

## 目录

- [第零章：自动安装指南](#第零章自动安装指南)
- [第一章：系统环境准备](#第一章系统环境准备)
- [第二章：Python 环境配置](#第二章python-环境配置)
- [第三章：外部工具配置](#第三章外部工具配置)
- [第四章：PMEMD 编译](#第四章pmemd-编译)
- [第五章：Alchemd 配置](#第五章alchemd-配置)
- [第六章：可选组件配置](#第六章可选组件配置)
- [第七章：常见问题解答](#第七章常见问题解答)
- [附录：故障排除](#附录故障排除)

---

## 第零章：自动安装指南

## 概述

Alchemd 提供了自动化安装脚本以简化复杂的环境配置过程。自动安装通过三个脚本按顺序执行系统环境准备、Python 环境配置和 PMEMD 编译，适用于标准 Linux 环境的快速部署。

**注意**：自动安装脚本适用于大多数标准配置，如遇特殊环境或依赖冲突，建议参考手动安装指南进行定制化配置。

---

## 0.1 前提条件

### 0.1.1 系统要求

- **操作系统**：Ubuntu 18.04+、CentOS 7+、Fedora、openSUSE 或 Arch Linux
- **CUDA 工具包**：CUDA 12.0（必需，用于 GPU 加速）
- **Amber 许可证**：需要向 https://ambermd.org/ 申请获取

### 0.1.2 必要软件

**Conda/Mamba 环境管理器**：
```bash
# 推荐安装 Miniconda
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh

# 或安装 Mamba（性能更优）
conda install mamba -n base -c conda-forge
```

**NVIDIA 驱动和 CUDA**：
- 确保安装了兼容的 NVIDIA 驱动程序
- CUDA 12.0 工具包必须可通过 `nvcc --version` 命令访问

### 0.1.3 网络连接优化

对于中国用户，可配置清华大学镜像源以提升下载速度：

```bash
# 配置 Conda 清华源
conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/pkgs/free/
conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/pkgs/main/
conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud/conda-forge/
conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud/bioconda/
```

---

## 0.2 安装步骤

### 0.2.1 系统环境准备

执行第一个脚本安装系统依赖和 GCC 工具链：

```bash
chmod +x install/setup_linux.sh
./install/setup_linux.sh
```

此脚本将自动：
- 检测 Linux 发行版和包管理器
- 安装基础开发包（build-essential、cmake、flex、bison 等）
- 配置 GCC 9+ 编译器环境
- 验证编译器版本和可用性

### 0.2.2 Python 环境配置

执行第二个脚本创建 Alchemd 专用 Python 环境：

```bash
chmod +x install/setup_conda.sh
./install/setup_conda.sh
```

此脚本将自动：
- 创建名为 `alchemd` 的 conda 环境（Python 3.9）
- 安装科学计算依赖（numpy、scipy、pandas、matplotlib 等）
- 安装另外的依赖（rdkit、openmm、parmed、lomap2 等）
- 下载并配置外部工具（WatVina等）
- 初始化 Alchemd 配置文件

### 0.2.3 PMEMD 编译

执行第三个脚本编译带 CUDA 支持的 PMEMD：

```bash
chmod +x install/build_amber.sh
./install/build_amber.sh /path/to/pmemd24
```

此脚本将自动：
- 应用 Alchemd 专用补丁到 PMEMD 源码
- 配置 CUDA 编译选项
- 执行并行编译和安装
- 设置环境变量集成

**参数说明**：
- `/path/to/pmemd24`：PMEMD 源码解压目录的路径
- 可选第二个参数指定编译并行度（默认 2 核）

---

## 0.3 验证安装

### 0.3.1 激活环境

```bash
conda activate alchemd
```

### 0.3.2 验证 Python 环境

```bash
# 测试核心模块
python core/prepare_file.py --help

# 验证关键依赖
python -c "import numpy, rdkit, openmm, parmed; print('Core dependencies OK')"
```

### 0.3.3 验证 PMEMD 安装

```bash
# 检查 PMEMD CUDA 版本
which pmemd.cuda_SPFP
pmemd.cuda_SPFP --help

# 运行环境测试脚本
./install/test_amber_env.sh
```

### 0.3.4 验证配置文件

```bash
# 检查配置文件生成
ls -la core/configs.toml

# 验证路径配置
python -c "from core.analyze_tools.SettingManager import SettingManager; s=SettingManager(); print('Config loaded successfully')"
```

---

## 0.4 故障排除

### 0.4.1 常见问题

**GCC 版本过低**：
脚本会自动尝试安装 GCC 11，如失败请手动安装或参考手动安装指南。

**CUDA 环境问题**：
确保 `nvcc --version` 可正常执行且显示 CUDA 12.0 版本。

**依赖包冲突**：
可尝试清理 conda 缓存：`conda clean --all`，然后重新运行脚本。

**网络下载失败**：
配置镜像源或使用代理网络环境。

### 0.4.2 重新安装

如需重新执行某个步骤：

```bash
# 重新配置 Python 环境（保留现有环境）
./install/setup_conda.sh

# 仅更新依赖（跳过环境创建）
./install/setup_conda.sh --skip-env

# 重新编译 PMEMD
./install/build_amber.sh /path/to/pmemd24
```

---

## 第一章：系统环境准备

### 1.1 Ubuntu 系统更新

首先更新系统软件包列表和已安装的软件包：

```bash
# 更新软件包列表
sudo apt update

# 升级已安装的软件包（可选，但推荐）
sudo apt upgrade -y

# 安装必要的系统工具
sudo apt install -y curl wget git vim
```

**验证**：
```bash
# 检查系统版本
lsb_release -a
# 应显示 Ubuntu 18.04+ 版本
```

### 1.2 基础开发工具安装

安装编译和开发所需的基础工具：

```bash
# 安装构建工具
sudo apt install -y build-essential

# 安装编译器和工具链
sudo apt install -y gcc g++ gfortran

# 安装构建系统
sudo apt install -y make cmake

# 安装其他必需工具
sudo apt install -y tcsh flex bison patch bc unzip

# 安装开发库
sudo apt install -y xorg-dev libz-dev libbz2-dev
```

**验证**：
```bash
# 检查工具是否安装成功
gcc --version      # 应显示 GCC 版本信息
g++ --version      # 应显示 G++ 版本信息
gfortran --version # 应显示 Gfortran 版本信息
make --version     # 应显示 Make 版本信息
cmake --version    # 应显示 CMake 版本信息
```

### 1.3 GCC 编译器安装（GCC 9+）

检查并确保 GCC 版本满足要求：

```bash
# 检查当前 GCC 版本
gcc --version

# 提取主版本号
gcc -dumpversion | cut -d. -f1
```

**如果版本 >= 9**：
```bash
echo "GCC 版本满足要求，继续下一步"
```

**如果版本 < 9**：
```bash
# 安装更新的 GCC 版本
sudo apt install -y gcc-11 g++-11 gfortran-11

# 设置默认版本（可选）
sudo update-alternatives --install /usr/bin/gcc gcc /usr/bin/gcc-11 60 \
    --slave /usr/bin/g++ g++ /usr/bin/g++-11 \
    --slave /usr/bin/gfortran gfortran /usr/bin/gfortran-11

# 验证版本
gcc --version
```

**依赖说明**：
- **GCC 9+**：编译 PMEMD 的最低要求
- **GCC 11+**：推荐版本，具有更好的优化和兼容性
- **Gfortran**：编译 Fortran 代码（PMEMD 核心）必需

### 1.4 CUDA 12.0 安装

**重要**：CUDA 12.0 是 Alchemd 的硬性要求，必须手动安装。

#### 1.4.1 检查 NVIDIA 驱动

```bash
# 检查 NVIDIA 驱动
nvidia-smi

# 如果命令不存在或报错，需要先安装 NVIDIA 驱动
```

#### 1.4.2 安装 NVIDIA 驱动（如果需要）

```bash
# 自动检测并安装推荐驱动
sudo ubuntu-drivers autoinstall

# 或手动安装特定版本
# sudo apt install -y nvidia-driver-535

# 重启系统
sudo reboot
```

#### 1.4.3 下载和安装 CUDA 12.0

**重要提示**：CUDA版本必须与GCC编译器版本兼容。CUDA 12.0建议使用GCC 9-11版本，避免使用过高版本的GCC（如GCC 12+）编译低版本的CUDA，可能导致兼容性问题。

```bash
# 访问 NVIDIA 官网下载 CUDA 12.0
# 官网链接：https://developer.nvidia.com/cuda-downloads
# 选择 Linux -> x86_64 -> Ubuntu -> 版本 -> runfile(local)
# 下载到 ~/Downloads 目录

cd ~/Downloads
# 运行安装程序（文件名可能因版本而异）
# 首先进行正常的CUDA安装
sudo sh cuda_12.0.0_*_linux.run

# 在安装界面中：
# - 取消选择 "Driver"（如果已安装驱动）
# - 保持选择 "CUDA Toolkit"
# - 确认安装路径为 /usr/local/cuda-12.0

# 如果需要额外安装NVIDIA驱动（在基本安装完成后）：
# sudo sh cuda_12.0.0_*_linux.run --silent --driver
```

#### 1.4.4 配置 CUDA 环境变量

```bash
# 编辑 ~/.bashrc
vim ~/.bashrc

# 添加以下行到文件末尾：
export CUDA_HOME=/usr/local/cuda-12.0
export PATH=$CUDA_HOME/bin:$PATH
export LD_LIBRARY_PATH=$CUDA_HOME/lib64:$LD_LIBRARY_PATH

# 保存并重新加载配置
source ~/.bashrc
```

**验证**：
```bash
# 检查 CUDA 安装
nvcc --version
nvidia-smi
# 两个命令都应该成功运行并显示版本信息
```

---

## 第二章：Python 环境配置

### 2.1 Conda/Mamba 安装

#### 2.1.1 选择安装方式

**推荐 Mamba**（更快的包管理器）：
```bash
# 访问官网下载最新版本的 Mambaforge
# 官网链接：https://github.com/conda-forge/miniforge
# 选择适合你系统的版本下载到 ~/Downloads 目录

# 安装 Mambaforge（以x86_64版本为例）
cd ~/Downloads
bash Mambaforge-Linux-x86_64.sh

# 按提示操作：
# - 同意许可协议
# - 确认安装路径（默认 ~/mambaforge）
# - 选择初始化 shell（推荐 yes）
```

**或者安装 Miniconda**：
```bash
# 访问官网下载最新版本的 Miniconda
# 官网链接：https://docs.conda.io/en/latest/miniconda.html
# 选择适合你系统的版本下载到 ~/Downloads 目录

cd ~/Downloads
# 安装 Miniconda（以x86_64版本为例）
bash Miniconda3-latest-Linux-x86_64.sh
```

#### 2.1.2 重启终端并验证

```bash
# 重启终端或重新加载配置
source ~/.bashrc

# 验证安装
mamba --version  # 如果安装了 Mamba
# 或
conda --version  # 如果安装了 Conda
```

**依赖说明**：
- **Conda/Mamba**：Python 包和环境管理
- **Python 3.9.21**：Alchemd 推荐的 Python 版本

### 2.2 虚拟环境创建

```bash
# 创建 Alchemd 专用环境
mamba create -n alchemd python=3.9.21
# 或使用 conda create -n alchemd python=3.9.21

# 激活环境
mamba activate alchemd
# 或 conda activate alchemd

# 验证环境
python --version  # 应显示 Python 3.9.21
which python      # 应显示 mambaforge/envs/alchemd/bin/python 路径
```

### 2.3 Python 依赖安装

#### 2.3.1 获取 Alchemd 源码

```bash
# 克隆项目（如果还没有）
cd ~
git clone <repository-url> Alchemd_release
cd Alchemd_release
```

#### 2.3.2 安装 Python 依赖

```bash
# 确保在 alchemd 环境中
mamba activate alchemd

# 安装依赖（使用 requirements.txt）
pip install -r requirements.txt

# 验证关键依赖
python -c "import numpy; print('NumPy:', numpy.__version__)"
python -c "import rdkit; print('RDKit:', rdkit.__version__)"
python -c "import openmm; print('OpenMM:', openmm.__version__)"
```

**依赖说明**：
- **NumPy, SciPy**：科学计算基础
- **RDKit**：分子处理和化学信息学
- **OpenMM**：分子动力学模拟引擎
- **Pandas**：数据处理
- **Matplotlib**：绘图和可视化

---

## 第三章：外部工具配置

### 3.1 WatVina 下载配置

```bash
# 确保在项目根目录
cd ~/Alchemd_release

# 创建 dependencies 目录
mkdir -p dependencies/watvina

# 下载 WatVina
wget -O dependencies/watvina/watvina \
    https://github.com/biocheming/watvina/releases/download/v20241125/watvina

# 设置执行权限
chmod +x dependencies/watvina/watvina

# 验证
./dependencies/watvina/watvina --help
```

### 3.2 相关项目依赖

```bash
# 进入 dependencies 目录
cd dependencies

# 下载 Weighted_cc 项目
wget -O weighted_cc.zip \
    https://github.com/zlisysu/Weighted_cc/archive/refs/heads/main.zip
unzip weighted_cc.zip

# 下载 Pyautomd 项目
wget -O pyautomd.zip \
    https://github.com/phloglucinol/Pyautomd/archive/refs/heads/main.zip
unzip pyautomd.zip

# 安装 Pyautomd
pip install ./Pyautomd-main/

# 下载并安装 Lomap（lomap2）
wget -O lomap.zip \
    https://github.com/OpenFreeEnergy/Lomap/archive/refs/heads/main.zip
unzip lomap.zip

# 安装 Lomap
pip install ./Lomap-main/

# 验证 Lomap 安装
python -c "import lomap; print('Lomap安装成功')"

# 清理下载文件
rm -f weighted_cc.zip pyautomd.zip lomap.zip
```

**依赖说明**：
- **Weighted_cc**: 分析工具
- **Pyautomd**: 自动准备分子工具
- **Lomap**: 生成perturbation map工具

---

## 第四章：PMEMD 编译

### 4.1 获取 Amber 许可证

1. **访问 Amber 官网**：https://ambermd.org/
2. **申请学术许可证**
3. **下载许可文件**：收到确认邮件后下载 PMEMD 包

### 4.2 PMEMD 源码准备

```bash
# 创建安装目录
mkdir -p ~/software
cd ~/software

# 解压 PMEMD 包（假设已下载到 ~/Downloads）
tar -xjf ~/Downloads/pmemd24.tar.bz2
# 这将创建 pmemd24 目录

# 验证目录结构
ls pmemd24/
# 应包含：build/, src/, 等目录
```

### 4.3 编译配置和执行

#### 4.3.1 预编译准备

```bash
# 进入 PMEMD 目录
cd ~/software/pmemd24

# 检查 build 目录
ls build/
# 应包含 run_cmake 脚本

# 应用 Alchemd 补丁
cp ~/Alchemd_release/dependencies/AmberFiles/patch.sh src/pmemd/src/
cp ~/Alchemd_release/dependencies/AmberFiles/ti.F90.patch src/pmemd/src/
cd src/pmemd/src/
chmod +x patch.sh
./patch.sh
cd ../../..
```

#### 4.3.2 配置编译选项

```bash
# 编辑编译脚本
cd build
cp run_cmake run_cmake.backup

# 修改 run_cmake 文件
sed -i 's/DCUDA=FALSE/DCUDA=TRUE/g' run_cmake
sed -i 's/DDOWNLOAD_MINICONDA=TRUE/DDOWNLOAD_MINICONDA=FALSE/g' run_cmake

# 验证修改
grep "DCUDA=TRUE" run_cmake
grep "DDOWNLOAD_MINICONDA=FALSE" run_cmake
```

#### 4.3.3 执行编译

```bash
# 确保在 alchemd 环境中
mamba activate alchemd

# 设置环境变量
export AMBERHOME=~/software/pmemd24
export CUDA_HOME=/usr/local/cuda-12.0

# 清理之前的构建（如果存在）
if [ -f "./clean_build" ]; then
    echo "y" | ./clean_build
else
    echo "y" | make clean || true
fi

# 运行 CMake 配置
./run_cmake

# 开始编译（使用 4 核）
make install -j4

# 编译过程可能需要 30-120 分钟，取决于系统性能
```

**故障排除**：
- 如果内存不足，减少并行度：`make install -j2`
- 如果 CUDA 错误，检查 CUDA_HOME 环境变量
- 如果 GCC 错误，确认版本 >= 9

#### 4.3.4 验证编译结果

```bash
# 检查生成的二进制文件
ls ~/software/pmemd24/bin/

# 应包含：
# - pmemd
# - pmemd.cuda
# - pmemd.cuda_SPFP

# 测试关键二进制文件
~/software/pmemd24/bin/pmemd.cuda_SPFP --help
```

### 4.4 环境变量设置

#### 4.4.1 自动环境设置

```bash
# 使用 Alchemd 提供的环境设置脚本
cd ~/Alchemd_release
./install/setup_amber_env.sh ~/software/pmemd24
```

#### 4.4.2 手动环境设置（可选）

```bash
# 编辑 conda 环境激活脚本
mkdir -p ~/mambaforge/envs/alchemd/etc/conda/activate.d/
cat > ~/mambaforge/envs/alchemd/etc/conda/activate.d/amber_init.sh << 'EOF'
#!/bin/bash
# PMEMD environment activation
export AMBERHOME="$HOME/software/pmemd24"
source "$AMBERHOME/amber.sh"
echo "PMEMD environment activated: $AMBERHOME"
EOF

chmod +x ~/mambaforge/envs/alchemd/etc/conda/activate.d/amber_init.sh

# 创建去激活脚本
mkdir -p ~/mambaforge/envs/alchemd/etc/conda/deactivate.d/
cat > ~/mambaforge/envs/alchemd/etc/conda/deactivate.d/amber_init.sh << 'EOF'
#!/bin/bash
# PMEMD environment deactivation
if [[ "$AMBERHOME" == "$HOME/software/pmemd24" ]]; then
    unset AMBERHOME
    echo "PMEMD environment deactivated"
fi
EOF

chmod +x ~/mambaforge/envs/alchemd/etc/conda/deactivate.d/amber_init.sh
```

#### 4.4.3 测试环境设置

```bash
# 重新激活环境
mamba deactivate
mamba activate alchemd

# 检查 AMBERHOME
echo $AMBERHOME
# 应显示：/home/用户名/software/pmemd24

# 测试 PMEMD 工具
which pmemd.cuda_SPFP
pmemd.cuda_SPFP --help
```

---

## 第五章：Alchemd 配置

### 5.1 配置文件设置

Alchemd 使用 TOML 格式的配置文件来管理各种工具和依赖的路径。我们提供了一个自动化配置工具来简化这个过程。

#### 5.1.1 使用自动配置工具

```bash
# 确保在 alchemd 环境中
mamba activate alchemd
cd ~/Alchemd_release

# 运行配置初始化工具
python core/init_setup.py
```

该工具将：
- 自动检测项目路径
- 验证依赖工具的可用性
- 生成 `core/configs.toml` 配置文件
- 提供配置完整性测试

#### 5.1.2 交互式配置（高级用户）

如果需要自定义路径设置：

```bash
# 交互式配置，可自定义各工具路径
python core/init_setup.py --interactive
```

#### 5.1.3 验证配置

配置完成后，可以检查生成的配置文件：

```bash
# 查看配置文件内容
cat core/configs.toml

# 应包含类似以下内容：
# [paths]
# AlchemdToolkit_Path = "/home/用户名/Alchemd_release"
# CAR_Path = "/home/用户名/Alchemd_release/dependencies/CAR-FEP"
# ...
```

### 5.2 测试验证

#### 5.2.1 基础功能测试

```bash
# 确保在 alchemd 环境中
mamba activate alchemd
cd ~/Alchemd_release

# 测试核心模块
python core/prepare_file.py --help
python core/run_md.py --help  
python core/analyze_result.py --help

# 测试 GUI（可选）
python gui/gui.py
```

#### 5.2.2 PMEMD 集成测试

**重要说明**：请严格按照以下顺序进行测试，确保环境设置正确：

```bash
# 1. 确保在 alchemd 环境中（编译完成后应该自动处于此环境）
echo $CONDA_DEFAULT_ENV  # 应显示 alchemd

# 2. 首先设置 PMEMD 环境变量集成
./install/setup_amber_env.sh /path/to/pmemd24

# 3. 重新激活环境以加载新的环境变量
conda deactivate
conda activate alchemd

# 4. 运行环境测试脚本
./install/test_amber_env.sh

# 5. 检查输出，应显示大部分测试通过
```

**故障排除**：
- 如果测试失败，确认 PMEMD 编译已完成且无错误
- 检查 `/path/to/pmemd24/bin/` 目录下是否存在 `pmemd.cuda_SPFP` 等文件
- 确认 CUDA 环境变量已正确设置

#### 5.2.3 示例数据测试

Alchemd 提供了示例数据用于测试完整工作流程：

```bash
# 查看示例数据
ls example/
# 应显示：
# ligands.sdf  pairs.txt  protein.pdb

# 使用示例数据测试准备步骤
python core/prepare_file.py -p example/protein.pdb -l example/ligands.sdf -o test_output --cs-align -d -auto

# 检查输出目录
ls test_output/
# 应生成准备好的输入文件
```

**示例数据说明**：
- `protein.pdb`：示例蛋白质结构文件
- `ligands.sdf`：包含多个配体分子的 SDF 文件
- `pairs.txt`：预定义的配体配对文件

**注意**：完整的 MD 模拟测试需要已编译的 PMEMD，通常计算时间较长，建议在生产环境中进行。

## 第六章：可选组件配置

### 6.1 NSAA模块依赖（可选）

NSAA（非标准氨基酸）准备模块是Alchemd的可选功能组件，专门用于处理包含非标准氨基酸的蛋白质系统。该模块需要Gaussian量子化学软件支持。

#### 6.1.1 Gaussian软件要求

**重要说明**：Gaussian软件需要用户自行获取和安装，Alchemd不提供Gaussian软件本身。

**支持的Gaussian版本**：
- Gaussian 03 (g03)
- Gaussian 09 (g09) 
- Gaussian 16 (g16) - 推荐版本

**许可证要求**：
- 学术用户：需要有效的学术许可证
- 商业用户：需要购买商业许可证
- 申请地址：https://gaussian.com/

#### 6.1.2 Gaussian安装指导

**注意**：以下仅为一般性指导，具体安装步骤请参考Gaussian官方文档。

```bash
# 创建Gaussian安装目录
sudo mkdir -p /opt/gaussian
cd /opt/gaussian

# 解压Gaussian软件包（假设已获得许可并下载）
# tar -xvf g16.tar  # 示例：实际文件名可能不同

# 设置环境变量
export g16root="/opt/gaussian"
export GAUSS_SCRDIR="/tmp"
export GAUSS_MEMDEF="2GB"

# 加载Gaussian环境
source $g16root/g16/bsd/g16.profile
```

#### 6.1.3 Alchemd中的NSAA配置

```bash
# 进入Alchemd项目目录
cd ~/Alchemd_release

# 编辑NSAA配置模板
vim core/templates/nsaa_prep.txt
```

在第11行找到Gaussian配置：
```
gaussian = g03
```

根据您安装的Gaussian版本修改：
```
# 配置示例
gaussian = g16           # 如果使用Gaussian 16
gaussian = g09           # 如果使用Gaussian 09
gaussian = /opt/gaussian/g16/g16  # 使用完整路径
```

#### 6.1.4 验证NSAA模块配置

```bash
# 激活alchemd环境
mamba activate alchemd

# 测试Gaussian可用性
which g16  # 或您配置的Gaussian命令
g16 < /dev/null  # 测试Gaussian是否能正常启动

# 测试NSAA模块
python -c "from core.prepare_input_tools.nsaa_preparer import NSAAProcessor; print('NSAA模块可用')"
```

#### 6.1.5 使用说明

NSAA模块通过prepare_file.py的`-npd`参数使用：

```bash
# 包含非标准氨基酸的蛋白质处理示例
python core/prepare_file.py \
    -p protein_with_nsaa.pdb \
    -l ligands.sdf \
    -npd ./nsaa_parameters \
    -o ./nsaa_project \
    --cs-align -d -auto
```

#### 6.1.6 故障排除

**Q: 找不到Gaussian命令**
```bash
# 检查环境变量
echo $g16root
echo $PATH | grep gaussian

# 重新加载Gaussian环境
source /opt/gaussian/g16/bsd/g16.profile
```

**Q: Gaussian许可证问题**
- 确认许可证文件位置正确
- 检查许可证是否在有效期内
- 联系Gaussian技术支持

**Q: 内存不足**
```bash
# 调整Gaussian内存设置
export GAUSS_MEMDEF="4GB"  # 根据系统配置调整
```

**注意事项**：
1. NSAA模块仅在处理包含非标准氨基酸的蛋白质时需要
2. 标准氨基酸的FEP计算不需要Gaussian
3. 量子化学计算较为耗时，建议在高性能环境中使用
4. 定期检查Gaussian许可证有效性

---

## 第七章：常见问题解答（FAQ）

### 7.1 Conda/Mamba 相关问题

#### Q: 遇到 CondaToSNonInteractiveError 错误

**问题描述**: "Terms of Service have not been accepted for the following channels. Please accept or remove them before proceeding"

**解决方案**:
```bash
# 方法一：接受指定频道的服务条款
conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/main
conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/r

# 方法二：查看详细帮助并手动处理
conda tos --help
```

**或者**: 查阅 Conda 官方文档了解如何接受服务条款。

### 7.2 编译相关问题

#### Q: WSL 环境下编译时出现段错误 (Segment Fault)

**问题描述**: 在 WSL 环境中编译时莫名出现段错误或编译失败

**解决方案**:
```bash
# 方法一：调整编译核数（默认为4）
./install/build_amber.sh /path/to/pmemd24_src 2  # 使用2核编译
./install/build_amber.sh /path/to/pmemd24_src 4  # 使用4核编译

# 方法二：重新执行编译
# WSL 有时会出现随机段错误，重新编译通常能解决
./install/build_amber.sh /path/to/pmemd24_src
```

**建议**: WSL 环境下推荐使用较少的编译核数（2-4核）以避免内存不足导致的问题。

### 7.3 安装验证问题

#### Q: build_amber.sh 编译完成后验证失败

**问题描述**: 编译完成但 pmemd 工具无法找到或使用

**正确的验证流程**:
```bash
# 1. 编译完成后，保持在 alchemd 环境中
# （build_amber.sh 脚本会自动保持环境激活状态）

# 2. 设置环境变量
./install/setup_amber_env.sh /path/to/pmemd24

# 3. 重新激活环境以加载新的环境变量
conda deactivate
conda activate alchemd

# 4. 验证安装
./install/test_amber_env.sh

# 5. 测试 PMEMD 工具
which pmemd.cuda_SPFP
pmemd.cuda_SPFP --help
```

**注意**: 不要在编译脚本执行完毕后立即退出 conda 环境，正确的顺序是先设置环境变量，再重新激活环境。

---

## 附录：故障排除

### A.1 编译问题

**Q: 编译时出现"CUDA not found"错误**
```bash
# 检查 CUDA 安装路径
ls /usr/local/cuda*

# 设置正确的 CUDA_HOME
export CUDA_HOME=/usr/local/cuda-12.0

# 确保路径在 PATH 中
export PATH=$CUDA_HOME/bin:$PATH
```

**Q: GCC 版本太旧**
```bash
# 安装新版本
sudo apt install gcc-11 g++-11 gfortran-11

# 更新默认版本
sudo update-alternatives --install /usr/bin/gcc gcc /usr/bin/gcc-11 60
```

### A.2 运行时问题

**Q: 找不到 pmemd.cuda_SPFP**
```bash
# 检查文件是否存在
ls $AMBERHOME/bin/pmemd*

# 如果不存在，重新编译 PMEMD
cd ~/software/pmemd24/build
make install -j4
```

**Q: GPU 无法使用**
```bash
# 检查 NVIDIA 驱动
nvidia-smi

# 测试 CUDA
cd $CUDA_HOME/samples/1_Utilities/deviceQuery
make
./deviceQuery
```

### A.3 环境问题

**Q: Conda 环境混乱**
```bash
# 重建环境
mamba deactivate
mamba env remove -n alchemd
mamba create -n alchemd python=3.9.21
mamba activate alchemd
pip install -r requirements.txt
```

**Q: 权限问题**
```bash
# 检查文件权限
ls -la ~/software/pmemd24/bin/

# 修复权限（如果需要）
chmod +x ~/software/pmemd24/bin/*
```

### A.4 获取帮助

如果遇到本指南未涵盖的问题：

1. **检查日志文件**：审查编译或运行时的详细错误信息
2. **查阅官方文档**：Amber 官方文档和 Alchemd 项目文档
3. **社区帮助**：通过 GitHub Issues 或相关论坛获取帮助
4. **回退到自动化脚本**：如果手动安装遇到困难，尝试使用项目提供的自动化安装脚本

### A.5 配置工具常见问题

**Q: init_setup.py 报告依赖缺失**
```bash
# 检查依赖目录
ls dependencies/

# 确保已解压必要的 tar 文件
cd dependencies
tar -xf CAR-FEP.tar
tar -xf genambRBFE.tar
tar -xf AlchemConvTools.tar
```

**Q: AMBERHOME 环境变量未设置**
```bash
# 确保已编译 PMEMD 并设置环境
source /path/to/pmemd24/amber.sh
echo $AMBERHOME
```

**Q: 配置文件生成失败**
```bash
# 检查权限
ls -la core/

# 手动创建目录（如果需要）
mkdir -p core

# 重新运行配置
python core/init_setup.py --auto
```

### A.6 配置验证失败

**Q: 依赖工具测试失败**
```bash
# 手动测试各个工具
./dependencies/watvina/watvina --help
python dependencies/Weighted_cc-main/wcc_main.py -h
python dependencies/CAR-FEP/segmented_converge_control.py -h
```

**Q: 重新配置系统**
```bash
# 删除旧配置并重新生成
rm -f core/configs.toml
python core/init_setup.py --interactive
```

---

## 总结

恭喜！您已经完成了 Alchemd 的手动安装。现在您应该拥有：

- ✅ 完整配置的 Ubuntu 系统环境
- ✅ GCC 9+ 编译器工具链
- ✅ CUDA 12.0 支持
- ✅ Python 3.9.21 和所有必需依赖
- ✅ 完全编译的 PMEMD（带 CUDA 支持）
- ✅ 正确配置的 Alchemd 环境

您现在可以开始使用 Alchemd 进行分子动力学模拟和自由能计算了！

### 后续步骤

- 阅读 Alchemd 用户指南了解具体使用说明
- 准备您的蛋白质和配体输入文件
- 开始您的第一个 FEP 计算项目

### 建议事项

- 定期备份您的配置和结果
- 保持 CUDA 驱动程序更新
- 在计算过程中监控系统资源使用情况

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
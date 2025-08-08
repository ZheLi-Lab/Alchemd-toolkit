# Alchemd

[![License](https://img.shields.io/badge/license-CC%20BY--NC--SA--4.0-blue.svg)](#license)
[![Platform](https://img.shields.io/badge/platform-Linux-lightgrey.svg)](#system-requirements)
[![Python](https://img.shields.io/badge/python-3.9+-blue.svg)](#installation)
[![CUDA](https://img.shields.io/badge/CUDA-12.0-green.svg)](#system-requirements)

> 用于自由能微扰计算的计算机辅助药物设计工具包

## 项目简介

Alchemd 是基于 Python 的相对结合自由能(RBFE)计算工具包。该软件提供分子准备、动力学模拟和热力学分析的集成平台，配备增强收敛算法。

## 系统要求

- **操作系统**：Linux (Ubuntu 18.04+) 或 Windows (需要 WSL2)
- **硬件**：支持 CUDA 的 NVIDIA 显卡 (推荐 RTX 4060 Ti+)
- **依赖**：CUDA 12.0、Conda/Mamba、GCC 9+
- **许可证**：需要 Amber 许可证 ([申请地址](https://ambermd.org/))
- **可选依赖**：Gaussian（用于NSAA非标准氨基酸模块）

## 安装

### 自动化安装

**先决条件**：在运行自动化安装脚本之前，请确保您已安装：
- **Conda 或 Mamba** 并正确配置
- **CUDA 12.0** 并可在系统 PATH 中访问

依次执行以下安装脚本：

```bash
# 1. 系统环境设置
chmod +x install/setup_linux.sh
./install/setup_linux.sh

# 2. Python 环境配置
chmod +x install/setup_conda.sh
./install/setup_conda.sh

# 3. PMEMD 编译 (需要 Amber 许可证)
chmod +x install/build_amber.sh
./install/build_amber.sh /path/to/pmemd24
```

### 安装指南

详细安装步骤请参考：
- [`docs/Installation_Guide_CN.md`](docs/Installation_Guide_CN.md)

## 快速开始

```bash
# 结构准备
python core/prepare_file.py -p protein.pdb -l ligands.sdf -o workspace --cs-align -d -auto

# 分子动力学模拟
python core/run_md.py -c -ln 50 -p ./workspace/prepare -m Amber -auto

# 结果分析
python core/analyze_result.py -auto -c -w ./workspace/run -m Amber

# GUI 界面
python gui/gui.py
```

## 文档

- **用户指南**：[`docs/Alchemd_User_Guide_CN.md`](docs/Alchemd_User_Guide_CN.md)
- **安装指南**：[`docs/Installation_Guide_CN.md`](docs/Installation_Guide_CN.md)
- **示例数据**：[`example/`](example/) 目录

## 引用

如果您在研究中使用 Alchemd，请引用：

```bibtex
@software{alchemd2025,
  title={Alchemd: A toolkit for alchemical free energy perturbation calculations},
  year={2025},
  url={https://github.com/ZheLi-Lab/Alchemd-toolkit}
}
```

**方法学参考文献：**
- CAR 方法：[DOI](https://pubs.acs.org/doi/10.1021/acs.jctc.4c00939)  
- cs-fep 实现：[DOI](https://doi.org/10.1016/j.apsb.2024.06.021)

## 许可证

Alchemd 采用双重许可：

- **非商业使用**：根据 [CC BY-NC-SA 4.0](LICENSE.md) 许可证 - 个人、学术和非营利组织可免费使用
- **商业使用**：需要单独的商业许可证 - 请联系我们获取许可

完整的条款和条件请参见 [LICENSE](LICENSE.md) 文件。

## 支持

如有问题和建议，请使用 [GitHub Issues](../../issues)。
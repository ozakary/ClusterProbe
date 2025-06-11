# ClusterProbe

[![Python](https://img.shields.io/badge/python-3.7+-blue.svg)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![ASE](https://img.shields.io/badge/ASE-3.0+-green.svg)](https://wiki.fysik.dtu.dk/ase/)

---

📄 Author: **Ouail Zakary**  
- 📧 Email: [Ouail.Zakary@oulu.fi](mailto:Ouail.Zakary@oulu.fi)  
- 🔗 ORCID: [0000-0002-7793-3306](https://orcid.org/0000-0002-7793-3306)  
- 🌐 Website: [Personal Webpage](https://cc.oulu.fi/~nmrwww/members/Ouail_Zakary.html)  
- 📁 Portfolio: [GitHub Portfolio](https://ozakary.github.io/)

---
A Python tool for analyzing the local environment around Xenon atoms in molecular clusters and identifying anomalous structures based on coordination number thresholds.

## 🔬 Overview

This tool is designed to perform sanity checks on molecular cluster structures containing Xenon atoms. It analyzes the local coordination environment of Xenon atoms by:

- Defining spherical cutoff regions around each Xenon atom
- Counting neighboring atoms within the specified cutoff radius
- Identifying clusters with anomalously low coordination numbers
- Optionally isolating problematic clusters for further investigation

## 🚀 Features

- **Automated Structure Analysis**: Processes multiple cluster files automatically
- **Flexible Parameters**: Customizable cutoff radius and coordination thresholds
- **Progress Tracking**: Real-time progress bars for all operations
- **Comprehensive Reporting**: Detailed statistics and summaries
- **Error Handling**: Robust error detection and reporting
- **File Management**: Optional sorting of anomalous clusters
- **ASE Integration**: Leverages the Atomic Simulation Environment for distance calculations

## 📋 Requirements

### Dependencies
- Python 3.7+
- ASE (Atomic Simulation Environment) 3.0+
- NumPy
- tqdm (for progress bars)

### File Structure Requirements
Your cluster files should be organized as follows:
```
base_directory/
├── cluster_1/
│   ├── coord_1ClusterAroundXeNew.xyz
│   └── [other cluster files...]
├── cluster_2/
│   ├── coord_2ClusterAroundXeNew.xyz
│   └── [other cluster files...]
├── cluster_N/
│   ├── coord_NClusterAroundXeNew.xyz
│   └── [other cluster files...]
└── [script files]
```

### XYZ File Format
```
570
HeaderEmpty
Xe     17.762500000   16.982800000   11.928100000
H      9.768510000   13.474700000   18.691100000
H      8.112340000   14.364500000   20.419700000
H      9.283990000   15.473600000   20.973700000
...
```

## 🛠️ Installation

### Method 1: Direct Installation (Recommended)
```bash
# Create a Python virtual environment where all packages will be installed
python -m venv cluster_probe_env
source cluster_probe_env/bin/activate

# Clone the repository
git clone https://github.com/ozakary/ClusterProbe.git
cd ClusterProbe

# Install dependencies
pip install -r requirements.txt
```

### Method 2: Using conda
```bash
# Create conda environment
conda create -n cluster-analysis python=3.8
conda activate cluster-analysis

# Install dependencies
conda install -c conda-forge ase numpy tqdm
```

### requirements.txt
```
ase>=3.20.0
numpy>=1.19.0
tqdm>=4.60.0
```

## 📖 Usage

### Basic Usage
```bash
python cluster_sanity_check.py
```

### Advanced Usage
```bash
# Custom parameters
python cluster_sanity_check.py --rcut 5.5 --num_atoms 8

# With anomaly sorting
python cluster_sanity_check.py --rcut 6.0 --num_atoms 10 --sort_out_anomalies yes

# Specify different directory
python cluster_sanity_check.py --base_dir /path/to/clusters --sort_out_anomalies yes
```

### Command Line Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `--rcut` | float | 6.0 | Cutoff radius for neighbor search (Å) |
| `--num_atoms` | int | 10 | Minimum number of neighbors for good clusters |
| `--sort_out_anomalies` | str | no | Move anomalous clusters (`yes`/`no`) |
| `--base_dir` | str | . | Base directory containing cluster folders |

### Help
```bash
python cluster_sanity_check.py --help
```

## 📊 Output Examples

### Console Output
```
Cluster Sanity Check Analysis
==================================================
Parameters:
  - Base directory: .
  - Cutoff radius: 6.0 Å
  - Minimum neighbors: 10
  - Sort out anomalies: no

Searching for cluster files...
Found 150 cluster files.

Processing clusters: 100%|████████████| 150/150 [00:45<00:00,  3.31it/s]

Detailed Results:
--------------------------------------------------------------------------------
cluster_1       |  570 atoms | 2 Xe | Neighbors: [12, 14       ] | GOOD
cluster_2       |  568 atoms | 2 Xe | Neighbors: [8, 11        ] | ANOMALOUS
cluster_3       |  572 atoms | 2 Xe | Neighbors: [13, 15       ] | GOOD
...

================================================================================
CLUSTER ANALYSIS SUMMARY
================================================================================
Analysis Parameters:
  - Cutoff radius: 6.00 Å
  - Minimum neighbors threshold: 10

File Statistics:
  - Total cluster files found: 150
  - Successfully analyzed: 150
  - Failed analyses: 0

Cluster Composition:
  - Total atoms per cluster: 565 - 575 (avg: 570.2)
  - Xe atoms per cluster: 2 - 2 (avg: 2.0)

Xenon Coordination Analysis:
  - Total Xe atoms analyzed: 300
  - Coordination numbers: 6 - 16 (avg: 11.8 ± 2.1)

Quality Assessment:
  - Good clusters: 142 (94.7%)
  - Anomalous clusters: 8 (5.3%)
```

### Generated Directory Structure (with `--sort_out_anomalies yes`)
```
base_directory/
├── cluster_1/
├── cluster_3/
├── cluster_4/
├── ...
└── bad_seeds/
    ├── cluster_2/
    ├── cluster_7/
    └── cluster_15/
```

## 🔧 Algorithm Details

### Coordination Analysis Process

1. **File Discovery**: Automatically locates all cluster files matching the pattern `cluster_*/coord_*ClusterAroundXeNew.xyz`

2. **Structure Reading**: Uses ASE to parse XYZ files and extract atomic coordinates

3. **Xenon Identification**: Locates all Xenon atoms in each cluster

4. **Distance Calculation**: For each Xenon atom:
   - Creates a spherical region with radius `rcut`
   - Calculates Euclidean distances to all other atoms
   - Counts atoms within the cutoff sphere

5. **Quality Assessment**: Flags clusters where any Xenon atom has fewer than `num_atoms` neighbors

6. **Optional Cleanup**: Moves entire cluster directories (preserving all files) to `bad_seeds` folder

### Technical Notes

- **No Periodic Boundary Conditions**: Calculations assume isolated clusters
- **Distance Metric**: Euclidean distance in 3D space
- **Self-Exclusion**: Xenon atoms don't count themselves as neighbors
- **Error Resilience**: Continues processing even if individual clusters fail

## 📈 Performance

- **Processing Speed**: ~3-5 clusters/second (depends on cluster size)
- **Memory Usage**: Minimal - processes one cluster at a time
- **Scalability**: Tested with 1000+ clusters without issues

## 🐛 Troubleshooting

### Common Issues

**Issue**: `No cluster files found`
```bash
# Solution: Check file pattern and directory structure
ls cluster_*/coord_*ClusterAroundXeNew.xyz
```

**Issue**: `ImportError: No module named 'ase'`
```bash
# Solution: Install ASE
pip install ase
# or
conda install -c conda-forge ase
```

**Issue**: `Permission denied when moving folders`
```bash
# Solution: Check write permissions
chmod +w cluster_*
```

### Debug Mode
For detailed debugging, modify the script to include:
```python
import logging
logging.basicConfig(level=logging.DEBUG)
```

## 📚 Examples and Tutorials

### Example 1: Default Quality Check
```bash
# Check all clusters with default parameters
python cluster_sanity_check.py
```

### Example 2: Advanced Quality Control
```bash
# Use smaller cutoff and higher coordination requirement
python cluster_sanity_check.py --rcut 5 --num_atoms 10 --base_dir ./path/to/cluster/folders/
```

### Example 3: Advanced Quanlity Control with Automated Cleanup
```bash
# Identify and isolate bad clusters automatically
python cluster_sanity_check.py --rcut 5 --num_atoms 10 --base_dir ./path/to/cluster/folders/ --sort_out_anomalies yes
```

## 📄 License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## 📞 Support

If you encounter any issues or have questions:

1. Check the [Issues](https://github.com/ozakary/ClusterProbe/issues) page
2. Create a new issue with:
   - Python version
   - Operating system
   - Complete error message
   - Sample input files (if possible)

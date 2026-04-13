# Barrier Certified Safety Learning Control  
### When Sum-of-Squares Programming Meets Reinforcement Learning  
📄 Paper: https://arxiv.org/pdf/2206.07915  

---

## Overview

This repository contains the official implementation of the paper:

**Barrier Certified Safety Learning Control: When Sum-of-Squares Programming Meets Reinforcement Learning**  
Published at **IEEE Conference on Control Technology and Applications (CCTA 2022)**.

This project develops a **safety-critical reinforcement learning framework** that integrates **sum-of-squares (SOS) optimization** and **barrier certificates** to ensure provable safety guarantees during learning.

The goal of this work is to enable **safe learning-based control systems** for applications where safety is essential, such as:

- Autonomous vehicles  
- Robotics systems  
- Safety-critical cyber-physical systems  
- Intelligent transportation systems  

---

## Motivation

Traditional reinforcement learning methods often lack **formal safety guarantees**, which limits their adoption in real-world safety-critical applications.

This project addresses that challenge by:

- Embedding **barrier certificates** into learning pipelines  
- Using **sum-of-squares optimization** to verify system safety  
- Providing **provable safe regions** during policy learning  
- Enabling scalable safety verification for nonlinear systems  

This work contributes toward building **reliable and trustworthy learning-enabled control systems**.

---

## Key Contributions

- Developed a **barrier-certified reinforcement learning framework**
- Combined **sum-of-squares optimization** with **learning-based control**
- Provided **formal safety guarantees** for nonlinear systems
- Demonstrated improved **safe region estimation**
- Enabled safer deployment of learning-enabled control policies

---

## Repository Structure

```
CCTA2022_SafeRL/
│
├── controllers/ # Safe control policy implementations
├── sos_solver/ # Sum-of-squares optimization modules
├── rl_training/ # Reinforcement learning components
├── experiments/ # Simulation and evaluation scripts
├── utils/ # Helper functions
│
├── main.py # Main execution entry
├── config.yaml # Experiment configurations
└── README.md
```

---

## Installation

#### Clone the repository:

```bash
git clone https://github.com/Wogwan/CCTA2022_SafeRL.git
cd CCTA2022_SafeRL
```

#### Install dependencies:

```
pip install -r requirements.txt

(Optional) Install additional solver dependencies if required:

pip install cvxpy
pip install sympy
```

#### Usage

Run the main training pipeline:

```
python main.py
```


#### Run experiments:

```
python experiments/run_experiment.py
```

#### Modify configuration parameters in `config.yaml`

The framework computes safe regions of attraction and learns safe control policies while ensuring safety constraints are satisfied.

Typical outputs include:

- Verified safe regions
- Policy learning trajectories
- Safety boundary visualization
- Performance metrics

This framework supports safety-critical applications including:

- Autonomous driving safety systems
- Robotics motion planning
- Aerospace control systems
- Intelligent transportation systems
- Cyber-physical infrastructure safety

#### Reproducibility

All experiments reported in the paper can be reproduced using:

```
python experiments/reproduce_results.py
```

Random seeds and configurations are provided for consistency.

## Citation

If you use this repository in your research, please cite:
```
@inproceedings{huang2022barrier,
  title={Barrier Certified Safety Learning Control: When Sum-of-Squares Programming Meets Reinforcement Learning},
  author={Huang, Hejun and others},
  booktitle={IEEE Conference on Control Technology and Applications (CCTA)},
  year={2022}
}
```

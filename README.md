# genomepaths
This repository contains an automated bioinformatics pipeline designed for the genome-scale metabolic reconstruction of microalgae, targeted screening of key pathways (carbon, nitrogen, and calcification), and interactive flux/presence network visualization.

---

## Overview

Reconstructing metabolic networks for eukaryotic microalgae poses unique challenges due to subcellular compartmentalization (chloroplasts, mitochondria) and complex evolutionary origins. This pipeline automates this process by combining structural homology alignment, constraint-based modeling, and dynamic functional annotation enrichment.

### Computational Workflow
1. **_De novo_ reconstruction:** Uses **CarveMe** to build a draft metabolic model from a protein FASTA file (`.faa`). Eukaryotic guidance is achieved by supplying orthology data via the `-i` (initialization) flag.
2. **Dynamic functional enrichment:** Parses **eggNOG-mapper** outputs dynamically to cross-reference and enrich missing Enzyme Commission (EC) numbers and KEGG Orthology (KO) terms, avoiding brittle hardcoded table indices.
3. **Targeted enzyme screening:** Performs multi-angle matching to identify and report the presence of critical enzymes:
   * **Carbonic Anhydrase** (EC 4.2.1.1) — Key for Carbon Concentrating Mechanisms (CCM) and Calcification.
   * **Urease** (EC 3.5.1.5) — Crucial for Nitrogen assimilation.
4. **Network visualization:** Maps the microalgae's metabolic repertoire onto the **RECON1 (Metabolic Core)** eukaryotic template using **Escher**, highlighting present pathways in green and absent ones in gray.

## Dependencies and prerequisites

The pipeline relies on both standalone bioinformatics binaries and Python libraries. 

### Core Requirements:
* **Python 3.11** (Recommended for stability across COBRApy and CarveMe packages).
* **Diamond** (Ultra-fast protein aligner used internally by CarveMe).
* **Linear Programming (LP) Solver:** **CPLEX** (via pip) or **Gurobi** is mandatory for CarveMe to optimize and prune the universal metabolic network.
* **Tkinter:** Python's standard GUI library for interactive file selection dialogs.

---

## Installation Guide

To avoid dependency conflicts, it is highly recommended to install this pipeline inside a dedicated **Conda** environment.

### Linux Setup (Terminal)

Open your terminal and run the following commands sequentially:

```bash
# 1. Create a clean dedicated environment
conda create -n environment_name python=3.11 -y

# 2. Activate the environment
conda activate environment_name

# 3. Install binary dependencies from specialized channels
conda install -c bioconda diamond -y
conda install -c conda-forge tk spyder-kernels -y

# 4. Install Python libraries and solvers
pip install cplex carveme cobra escher pandas

```

## Windows Setup (Anaconda Prompt)
Always open your Anaconda Prompt as an Administrator.
If you encounter an error stating 'powershell' is not recognized, your system PATH is misconfigured. You must add C:\Windows\System32\WindowsPowerShell\v1.0\ to your Windows System Environment Variables before running these commands.

```bash
# 1. Create a clean dedicated environment
conda create -n pipe_alga python=3.11 -y

# 2. Activate the environment
conda activate pipe_alga

# 3. Install binary dependencies
conda install -c bioconda diamond -y
conda install -c conda-forge tk spyder-kernels -y

# 4. Install Python libraries and solvers
pip install cplex carveme cobra escher pandas

```

## Running the Pipeline

### 1. Connect your IDE (e.g., Spyder):

Go to Tools > Preferences > Python interpreter.

Select Use the following Python interpreter.

Set the path to your new environment's Python binary:

**Linux:** /home/YOUR_USER/anaconda3/envs/environment_name/bin/python

**Windows:** C:\Users\YOUR_USER\anaconda3\envs\environment_name\python.exe

Restart your IPython console.

### 2. Execute the script:
Run the pipegenome.py file. Interactive native pop-up windows will appear mimicking R's choose.file() behavior:

Select your Microalgae Protein FASTA file (.faa).

Select your eggNOG-mapper Annotations file (.emapper.annotations) 

[Press Cancel if you don't have it].

---

## Cite us!

Paper not published yet :/

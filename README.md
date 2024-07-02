# Generalized Phylogenetic Diversity
This repository contains a rust crate to compute the PD statistics of a phylogenetic tree

## Installation

To install you must first have cargo and rustup installed:
```bash
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
```

After installing the above command you can run the following to install genpd:
```bash
cargo install --git=https://github.com/sriram98v/Generalized-PD
```

Alternatively, you can install seq_class by cloning this repository and building it locally:
```bash
git clone https://github.com/sriram98v/Generalized-PD
cd Generalized-PD
cargo install --path=./
```

## Usage
### Finding the Cophenetic distance between a pair of trees
To compute the cophenetic distance between a pair of trees, please create a single file with the extension ```.tre``` containing the two trees in Newick format (line-separated). The run the following command to compute the cophenetic distance with depth as the path function:
```bash
genpd PD min -f <PATH TO .TRE FILE> -n <NUM_TAXA>
```

Please refer the help page for details on how to use other path functions using:
```bash
genpd PD min -h
```
extern crate clap;

use phylo::tree::io::Newick;
use clap::{arg, Command};
use phylo::tree::simple_rtree::RootedMetaTree;
use PD::pd::phylogenetic_diversity::PhylogeneticDiversity;
use phylo::tree::SimpleRootedTree;
use itertools::Itertools;
use PD::pd::TreePDMap;
use std::fs::File;
use std::io::Read;

fn main(){
    let matches = Command::new("Phylogenetics Rust")
        .version("1.0")
        .author("Sriram Vijendran <vijendran.sriram@gmail.com>")
        .subcommand(Command::new("PD")
            .about("Compute Phylogenetic Diversity")
            .subcommand(Command::new("min")
                .about("Compute minPD")
                .arg(arg!(-f --file <TREE_FILE> "Input Tree File")
                    .required(true)
                    .value_parser(clap::value_parser!(String))
                )
                .arg(arg!(-n --num_taxa <NUM_TAXA> "Input Tree File")
                    .required(true)
                    .value_parser(clap::value_parser!(usize))
                )
            )
        )
        .about("CLI tool for quick tree operations")
        .get_matches();

        match matches.subcommand(){
            Some(("PD", sub_m)) => {
                match sub_m.subcommand(){
                    Some(("min", min_pd)) => {
                        let mut tree_file = File::open(min_pd.get_one::<String>("file").expect("required")).unwrap();
                        let num_taxa = min_pd.get_one::<usize>("num_taxa").expect("required");
                        let mut trees = String::new();

                        tree_file.read_to_string(&mut trees).unwrap();
                        let tree_string = trees.split("\n").collect_vec()[0];
                        let tree = SimpleRootedTree::from_newick(tree_string.as_bytes());
                        let mut tree_pd = TreePDMap::new(&tree);
                        tree_pd.precompute_minPDs();
                        println!("minPD: {}\nnormalized minPD: {}\nminPD set:{}\nnormalized minPD set:{}", tree_pd.get_minPD(num_taxa.clone()), tree_pd.get_norm_minPD(num_taxa.clone()), tree_pd.get_minPD_taxa_set(num_taxa.clone()).map(|x| tree.get_node_taxa(x).unwrap()).join(","), tree_pd.get_norm_minPD_taxa_set(num_taxa.clone()).map(|x| tree.get_node_taxa(x).unwrap()).join(","));
                        // dbg!("{}", tree);
                    }
                    _ => println!("No valid PD metric chosen! Refer help page (-h flag)")
                }
            },
            _ => {
                println!("No option selected! Refer help page (-h flag)");
            }
        }
}
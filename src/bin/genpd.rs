extern crate clap;

use clap::{arg, Command};
use itertools::Itertools;
use phylo::prelude::RootedTree;
use phylo::tree::io::Newick;
use phylo::tree::simple_rtree::RootedMetaTree;
use phylo::tree::SimpleRootedTree;
use std::cmp;
use std::fs::File;
use std::io::Read;
use PD::pd::phylogenetic_diversity::{binarize_tree, TreePDMap};
use PD::pd::TreePD;
use anyhow::Result;

fn main() -> Result<()>{
    let matches = Command::new("Phylogenetics Rust")
        .version("1.0")
        .author("Sriram Vijendran <vijendran.sriram@gmail.com>")
        .subcommand(
            Command::new("PD")
                .about("Compute Phylogenetic Diversity")
                .subcommand(
                    Command::new("min")
                        .about("Compute minPD")
                        .arg(
                            arg!(-f --file <TREE_FILE> "Input Tree File")
                                .required(true)
                                .value_parser(clap::value_parser!(String)),
                        )
                        .arg(
                            arg!(-n --num_taxa <NUM_TAXA> "Input Tree File")
                                .required(true)
                                .value_parser(clap::value_parser!(usize)),
                        ),
                )
                .subcommand(
                    Command::new("max")
                        .about("Compute maxPD")
                        .arg(
                            arg!(-f --file <TREE_FILE> "Input Tree File")
                                .required(true)
                                .value_parser(clap::value_parser!(String)),
                        )
                        .arg(
                            arg!(-n --num_taxa <NUM_TAXA> "Input Tree File")
                                .required(true)
                                .value_parser(clap::value_parser!(usize)),
                        ),
                )
                .subcommand(
                    Command::new("all_max")
                        .about("Compute maxPD")
                        .arg(
                            arg!(-f --file <TREE_FILE> "Input Tree File")
                                .required(true)
                                .value_parser(clap::value_parser!(String)),
                        )
                )
                .subcommand(
                    Command::new("all_min")
                        .about("Compute maxPD")
                        .arg(
                            arg!(-f --file <TREE_FILE> "Input Tree File")
                                .required(true)
                                .value_parser(clap::value_parser!(String)),
                        )
                ),

        )
        .subcommand(
            Command::new("gen")
                .about("Compute Phylogenetic Diversity")
                .subcommand(
                    Command::new("min").about("Compute minPD").arg(
                        arg!(-f --file <TREE_FILE> "Input Tree File")
                            .required(true)
                            .value_parser(clap::value_parser!(String)),
                    ),
                )
                .subcommand(
                    Command::new("max").about("Compute maxPD").arg(
                        arg!(-f --file <TREE_FILE> "Input Tree File")
                            .required(true)
                            .value_parser(clap::value_parser!(String)),
                    ),
                )
                .subcommand(
                    Command::new("all")
                        .about("Compute maxPD")
                        .arg(
                            arg!(-f --file <TREE_FILE> "Input Tree File")
                                .required(true)
                                .value_parser(clap::value_parser!(String)),
                        )
                ),
        )
        .about("CLI tool for quick tree operations")
        .get_matches();

    match matches.subcommand() {
        Some(("PD", sub_m)) => {
            match sub_m.subcommand() {
                Some(("min", min_pd)) => {
                    let mut tree_file =
                        File::open(min_pd.get_one::<String>("file").expect("required"))?;
                    let num_taxa = min_pd.get_one::<usize>("num_taxa").expect("required");
                    let mut trees = String::new();

                    tree_file.read_to_string(&mut trees).unwrap();
                    let tree_string = trees.split("\n").collect_vec()[0];
                    let mut tree = SimpleRootedTree::from_newick(tree_string.as_bytes())?;
                    if !tree.is_binary(){
                        binarize_tree(&mut tree);
                    }
                    let num_taxa = match *num_taxa==0{
                        true => {println!("setting k to n");tree.num_taxa()},
                        false => cmp::min(*num_taxa as usize, tree.num_taxa()),
                    };
                    let tree_pd = TreePD::new(&tree);
                    println!(
                        "minPD: {}\nnormalized minPD: {}\nminPD set:{}\nnormalized minPD set:{}",
                        tree_pd.get_minPD(num_taxa.clone()),
                        tree_pd.get_norm_minPD(num_taxa.clone()),
                        tree_pd
                            .get_minPD_taxa_set(num_taxa.clone())
                            .map(|x| tree.get_node_taxa(x).unwrap())
                            .join(","),
                        tree_pd
                            .get_norm_minPD_taxa_set(num_taxa.clone())
                            .map(|x| tree.get_node_taxa(x).unwrap())
                            .join(",")
                    );
                    // dbg!("{}", tree);
                },
                Some(("all_min", min_pd)) => {
                    let mut tree_file =
                        File::open(min_pd.get_one::<String>("file").expect("required"))?;
                    // let num_taxa = min_pd.get_one::<usize>("num_taxa").expect("required");
                    let mut trees = String::new();

                    tree_file.read_to_string(&mut trees).unwrap();
                    let tree_string = trees.split("\n").collect_vec()[0];
                    let mut tree = SimpleRootedTree::from_newick(tree_string.as_bytes())?;
                    if !tree.is_binary(){
                        binarize_tree(&mut tree);
                    }
                    // let num_taxa = match *num_taxa==0{
                    //     true => {println!("setting k to n");tree.num_taxa()},
                    //     false => *num_taxa as usize,
                    // };
                    let tree_taxa: usize = tree.num_taxa();

                    let tree_pd = TreePD::new(&tree);
                    for num_taxa in 3..tree_taxa+1{
                        println!(
                            "k: {}\nminPD: {}\nnormalized minPD: {}\nminPD set:{}\nnormalized minPD set:{}\n",
                            num_taxa,
                            tree_pd.get_minPD(num_taxa.clone()),
                            tree_pd.get_norm_minPD(num_taxa.clone()),
                            tree_pd
                                .get_minPD_taxa_set(num_taxa.clone())
                                .map(|x| tree.get_node_taxa(x).unwrap())
                                .join(","),
                            tree_pd
                                .get_norm_minPD_taxa_set(num_taxa.clone())
                                .map(|x| tree.get_node_taxa(x).unwrap())
                                .join(",")
                        );    
                    }
                    // dbg!("{}", tree);
                },
                Some(("max", max_pd)) => {
                    let mut tree_file =
                        File::open(max_pd.get_one::<String>("file").expect("required"))?;
                    let n_taxa = max_pd.get_one::<usize>("num_taxa").expect("required");
                    let mut trees = String::new();

                    tree_file.read_to_string(&mut trees).unwrap();
                    let tree_string = trees.split("\n").collect_vec()[0];
                    let mut tree = SimpleRootedTree::from_newick(tree_string.as_bytes())?;
                    if !tree.is_binary(){
                        binarize_tree(&mut tree);
                    }
                    let num_taxa = match *n_taxa==0{
                        true => {println!("setting k to n");tree.num_taxa()},
                        false => cmp::min(*n_taxa, tree.num_taxa()),
                    };

                    let tree_pd = TreePD::new(&tree);
                    println!(
                        "maxPD: {}\nnormalized maxPD: {}\nmaxPD set:{}\nnormalized maxPD set:{}",
                        tree_pd.get_maxPD(num_taxa.clone()),
                        tree_pd.get_norm_maxPD(num_taxa.clone()),
                        tree_pd
                            .get_maxPD_taxa_set(num_taxa.clone())
                            .map(|x| tree.get_node_taxa(x).unwrap())
                            .join(","),
                        tree_pd
                            .get_norm_maxPD_taxa_set(num_taxa.clone())
                            .map(|x| tree.get_node_taxa(x).unwrap())
                            .join(",")
                    );
                    // dbg!("{}", tree);
                },
                Some(("all_max", max_pd)) => {
                    let mut tree_file =
                        File::open(max_pd.get_one::<String>("file").expect("required"))?;
                    // let n_taxa = max_pd.get_one::<usize>("num_taxa").expect("required");
                    let mut trees = String::new();

                    tree_file.read_to_string(&mut trees).unwrap();
                    let tree_string = trees.split("\n").collect_vec()[0];
                    let mut tree = SimpleRootedTree::from_newick(tree_string.as_bytes())?;
                    if !tree.is_binary(){
                        binarize_tree(&mut tree);
                    }
                    // let num_taxa = match *n_taxa==0{
                    //     true => {println!("setting k to n");tree.num_taxa()},
                    //     false => *n_taxa,
                    // };

                    let tree_taxa: usize = tree.num_taxa();

                    let tree_pd = TreePD::new(&tree);
                    for num_taxa in 3..tree_taxa+1{
                        println!(
                            "k: {}\nmaxPD: {}\nnormalized maxPD: {}\nmaxPD set:{}\nnormalized maxPD set:{}\n",
                            num_taxa,
                            tree_pd.get_maxPD(num_taxa.clone()),
                            tree_pd.get_norm_maxPD(num_taxa.clone()),
                            tree_pd
                                .get_maxPD_taxa_set(num_taxa.clone())
                                .map(|x| tree.get_node_taxa(x).unwrap())
                                .join(","),
                            tree_pd
                                .get_norm_maxPD_taxa_set(num_taxa.clone())
                                .map(|x| tree.get_node_taxa(x).unwrap())
                                .join(",")
                        );
                    }
                    // dbg!("{}", tree);
                },

                _ => println!("No valid PD metric chosen! Refer help page (-h flag)"),
            }
        }
        Some(("gen", sub_m)) => {
            match sub_m.subcommand() {
                Some(("min", min_pd)) => {
                    let mut tree_file =
                        File::open(min_pd.get_one::<String>("file").expect("required"))?;
                    let mut trees = String::new();

                    tree_file.read_to_string(&mut trees).unwrap();
                    let tree_string = trees.split("\n").collect_vec()[0];
                    let mut tree = SimpleRootedTree::from_newick(tree_string.as_bytes())?;
                    if !tree.is_binary(){
                        binarize_tree(&mut tree);
                    }
                    let num_taxa = tree.num_taxa();
                    let tree_pd = TreePD::new(&tree);
                    match num_taxa >=3{
                        true => {
                            println!(
                                "minGenPD: {}\nminGenPD set: {}\nminGenPD set size: {}",
                                tree_pd.get_min_genPD(),
                                tree_pd
                                    .get_min_genPD_set()
                                    .map(|x| tree.get_node_taxa(x).unwrap())
                                    .join(","),
                                tree_pd.get_min_genPD_set().collect_vec().len()
                            );
        
                        },
                        false => {
                            println!(
                                "minGenPD: {}\nminGenPD set: {}\nminGenPD set size: {}",
                                0,
                                "",
                                0
                            );
        
                        }
                    }
                    // dbg!("{}", tree);
                },
                Some(("max", max_pd)) => {
                    let mut tree_file =
                        File::open(max_pd.get_one::<String>("file").expect("required"))?;
                    let mut trees = String::new();

                    tree_file.read_to_string(&mut trees).unwrap();
                    let tree_string = trees.split("\n").collect_vec()[0];
                    let mut tree = SimpleRootedTree::from_newick(tree_string.as_bytes())?;
                    if !tree.is_binary(){
                        binarize_tree(&mut tree);
                    }
                    let num_taxa = tree.num_taxa();
                    let tree_pd = TreePD::new(&tree);
                    match num_taxa >=3{
                        true => {
                            println!(
                                "maxGenPD: {}\nmaxGenPD set: {}\nmaxGenPD set size: {}",
                                tree_pd.get_max_genPD(),
                                tree_pd
                                    .get_max_genPD_set()
                                    .map(|x| tree.get_node_taxa(x).unwrap())
                                    .join(","),
                                tree_pd.get_max_genPD_set().collect_vec().len()
                            );
                
                        },
                        false => {
                            println!(
                                "maxGenPD: {}\nmaxGenPD set: {}\nmaxGenPD set size: {}",
                                0,
                                "",
                                0
                            );
        
                        }
                    }

                    // dbg!("{}", tree);
                },
                _ => println!("No valid PD metric chosen! Refer help page (-h flag)"),
            }
        }
        _ => {
            println!("No option selected! Refer help page (-h flag)");
        }
    };
    Ok(())
}

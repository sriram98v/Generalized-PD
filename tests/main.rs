use phylo::prelude::*;
use itertools::Itertools;
use PD::pd::{TreePD, phylogenetic_diversity::{PhylogeneticDiversity, TreePDMap}};

#[test]
fn pd() {
    let input_str: String = String::from("(((A:1,B:2):2,C:7):4,(D:1,E:2):5);");
    let tree = SimpleRootedTree::from_newick(input_str.as_bytes()).unwrap();
    let mut tree_pd = TreePD::new(&tree);
    let num_taxa = 2;

    let tree_taxa = tree.get_taxa_space().collect_vec();

    dbg!(
        tree_pd.get_minPD(num_taxa.clone()),
        tree_pd.get_norm_minPD(num_taxa.clone()),
        tree_pd
            .get_minPD_taxa_set(num_taxa.clone())
            .map(|x| tree.get_node_taxa(x).unwrap())
            .join(","),
        tree_pd
            .get_norm_minPD_taxa_set(num_taxa.clone())
            .map(|x| tree.get_node_taxa(x).unwrap())
            .join(","),
        tree_pd.get_min_genPD(),
        tree_pd.get_min_genPD_set()
            .map(|x| tree.get_node_taxa(x).unwrap())
            .join(","),
    );

    dbg!(
        tree_pd.get_maxPD(num_taxa.clone()),
        tree_pd.get_norm_maxPD(num_taxa.clone()),
        tree_pd
            .get_maxPD_taxa_set(num_taxa.clone())
            .map(|x| tree.get_node_taxa(x).unwrap())
            .join(","),
        tree_pd
            .get_norm_maxPD_taxa_set(num_taxa.clone())
            .map(|x| tree.get_node_taxa(x).unwrap())
            .join(","),
        tree_pd.get_max_genPD(),
        tree_pd.get_max_genPD_set()
            .map(|x| tree.get_node_taxa(x).unwrap())
            .join(","),
    );
}

use phylo::prelude::*;
use itertools::Itertools;
use PD::pd::{TreePDMap, phylogenetic_diversity::PhylogeneticDiversity};

#[test]
fn pd() {
    let input_str: String = String::from("(((A:1,B:1):2,C:1):3,D:1);");
    let tree = SimpleRootedTree::from_newick(input_str.as_bytes()).unwrap();
    let mut tree_pd = TreePDMap::new(&tree);
    let num_taxa = 2;

    let tree_taxa = tree.get_taxa_space().collect_vec();

    tree_pd.precompute_minPDs();

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

    tree_pd.precompute_maxPDs();
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

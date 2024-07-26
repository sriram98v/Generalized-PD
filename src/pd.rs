pub mod phylogenetic_diversity;
use crate::pd::phylogenetic_diversity::{PhylogeneticDiversity, TreePDMap};
use itertools::Itertools;
use phylo::prelude::*;
use phylo::tree::SimpleRootedTree;
use std::cmp::Ordering;

pub struct TreePD<'a> {
    tree: &'a SimpleRootedTree,
    precomputed_min: Vec<Vec<(f32, u32)>>,
    precomputed_min_set: Vec<Vec<Vec<usize>>>,
    precomputed_norm_min: Vec<Vec<(f32, u32)>>,
    precomputed_norm_min_set: Vec<Vec<Vec<usize>>>,
    precomputed_max: Vec<Vec<(f32, u32)>>,
    precomputed_max_set: Vec<Vec<Vec<usize>>>,
    precomputed_norm_max: Vec<Vec<(f32, u32)>>,
    precomputed_norm_max_set: Vec<Vec<Vec<usize>>>,
}

impl<'a> TreePD<'a> {
    pub fn new(tree: &'a SimpleRootedTree) -> Self {
        let (min, min_set, min_norm, min_norm_set) = tree.compute_norm_min();
        let (max, max_set, max_norm, max_norm_set) = tree.compute_norm_max();
        TreePD {
            tree,
            precomputed_min: min,
            precomputed_min_set: min_set,
            precomputed_norm_min: min_norm,
            precomputed_norm_min_set: min_norm_set,
            precomputed_max: max,
            precomputed_max_set: max_set,
            precomputed_norm_max: max_norm,
            precomputed_norm_max_set: max_norm_set,
        }
    }
}
impl<'a>  TreePDMap<'a> for TreePD<'a> {

    type Tree = SimpleRootedTree;

    fn reset(&mut self) {
        self.precomputed_min = Vec::new();
        self.precomputed_min_set = Vec::new();
        self.precomputed_norm_min = Vec::new();
        self.precomputed_norm_min_set = Vec::new();
        self.precomputed_max = Vec::new();
        self.precomputed_max_set = Vec::new();
        self.precomputed_norm_max = Vec::new();
        self.precomputed_norm_max_set = Vec::new();
    }

    fn get_tree(&self)->&Self::Tree {
        self.tree
    }

    fn get_minPD(
        &self,
        num_taxa: usize,
    ) -> TreeNodeWeight<'a, SimpleRootedTree>
    {
        self.precomputed_min[self.tree.get_root_id()][std::cmp::min(num_taxa, self.tree.num_taxa())]
            .0
    }

    fn get_minPD_node(
        &self,
        node_id: TreeNodeID<'a, SimpleRootedTree>,
        num_taxa: usize,
    ) -> TreeNodeWeight<'a, SimpleRootedTree>
    {
        self.precomputed_min[node_id][std::cmp::min(num_taxa, self.tree.num_taxa())].0
    }

    fn get_norm_minPD(
        &self,
        num_taxa: usize,
    ) -> TreeNodeWeight<'a, SimpleRootedTree>
    {
        self.precomputed_norm_min[self.tree.get_root_id()]
            [std::cmp::min(num_taxa, self.tree.num_taxa())]
        .0
    }

    fn get_minPD_taxa_set(
        &self,
        num_taxa: usize,
    ) -> impl Iterator<Item = TreeNodeID<'a, SimpleRootedTree>> {
        self.precomputed_min_set[self.tree.get_root_id()][num_taxa]
            .clone()
            .into_iter()
    }

    fn get_norm_minPD_taxa_set(
        &self,
        num_taxa: usize,
    ) -> impl Iterator<Item = TreeNodeID<'a, SimpleRootedTree>> {
        self.precomputed_norm_min_set[self.tree.get_root_id()][num_taxa]
            .clone()
            .into_iter()
    }

    fn get_min_genPD(
        &self,
    ) -> TreeNodeWeight<'a, SimpleRootedTree>
    {
        self.precomputed_norm_min[self.tree.get_root_id()]
            .iter()
            .filter(|x| x.0 != 0_f32 && x.0 != f32::INFINITY)
            .min_by(|x, y| x.0.total_cmp(&y.0))
            .unwrap()
            .0
    }

    fn get_min_genPD_set(
        &self,
    ) -> impl Iterator<Item = TreeNodeID<'a, SimpleRootedTree>> {
        let num_taxa = self.precomputed_norm_min[self.tree.get_root_id()]
            .iter()
            .enumerate()
            .filter(|x| x.1 .0 != 0_f32 && x.1 .0 != f32::INFINITY)
            .min_by(|x, y| x.1 .0.total_cmp(&y.1 .0))
            .unwrap()
            .0;
        self.precomputed_norm_min_set[self.tree.get_root_id()][num_taxa]
            .clone()
            .into_iter()
    }

    fn get_maxPD(
        &self,
        num_taxa: usize,
    ) -> TreeNodeWeight<'a, SimpleRootedTree>
    {
        self.precomputed_max[self.tree.get_root_id()][std::cmp::max(num_taxa, 1)]
            .0
    }

    fn get_maxPD_node(
        &self,
        node_id: TreeNodeID<'a, SimpleRootedTree>,
        num_taxa: usize,
    ) -> TreeNodeWeight<'a, SimpleRootedTree>
    {
        self.precomputed_max[node_id][std::cmp::max(num_taxa, self.tree.num_taxa()-1)].0
    }

    fn get_norm_maxPD(
        &self,
        num_taxa: usize,
    ) -> TreeNodeWeight<'a, SimpleRootedTree>
    {
        self.precomputed_norm_max[self.tree.get_root_id()]
            [std::cmp::max(num_taxa, self.tree.num_taxa()-1)]
        .0
    }

    fn get_maxPD_taxa_set(
        &self,
        num_taxa: usize,
    ) -> impl Iterator<Item = TreeNodeID<'a, SimpleRootedTree>> {
        self.precomputed_max_set[self.tree.get_root_id()][num_taxa]
            .clone()
            .into_iter()
    }

    fn get_norm_maxPD_taxa_set(
        &self,
        num_taxa: usize,
    ) -> impl Iterator<Item = TreeNodeID<'a, SimpleRootedTree>> {
        self.precomputed_norm_max_set[self.tree.get_root_id()][num_taxa]
            .clone()
            .into_iter()
    }

    fn get_max_genPD(
        &self,
    ) -> TreeNodeWeight<'a, SimpleRootedTree>
    {
        self.precomputed_norm_max[self.tree.get_root_id()]
            .iter()
            .filter(|x| x.0 != 0_f32 && x.0 != f32::INFINITY)
            .max_by(|x, y| x.0.total_cmp(&y.0))
            .unwrap()
            .0
    }

    fn get_max_genPD_set(
        &self,
    ) -> impl Iterator<Item = TreeNodeID<'a, SimpleRootedTree>> {
        let num_taxa = self.precomputed_norm_max[self.tree.get_root_id()]
            .iter()
            .enumerate()
            .filter(|x| x.1 .0 != 0_f32 && x.1 .0 != f32::INFINITY)
            .max_by(|x, y| x.1 .0.total_cmp(&y.1 .0))
            .unwrap()
            .0;
        self.precomputed_norm_max_set[self.tree.get_root_id()][num_taxa]
            .clone()
            .into_iter()
    }

}

impl PhylogeneticDiversity for SimpleRootedTree {

    fn compute_dp_table<'a>(
        &self,
        op: std::cmp::Ordering,
    ) -> (
        Vec<Vec<(f32, u32)>>,
        Vec<Vec<Vec<TreeNodeID<'a, Self>>>>,
        Vec<Vec<(f32, u32)>>,
        Vec<Vec<Vec<TreeNodeID<'a, Self>>>>,
    )
    {
        let op_fn = match op{
            Ordering::Greater => PartialOrd::ge,
            _ => PartialOrd::le,
        };

        let num_leaves = self.get_leaves().len();
        let mut delta_bar =
            vec![vec![(f32::INFINITY, 0_u32); num_leaves + 1]; self.get_nodes().len()];
        let mut delta_bar_sets: Vec<Vec<Vec<usize>>> =
            vec![vec![vec![]; num_leaves + 1]; self.get_nodes().len()];
        let mut delta_hat =
            vec![vec![(f32::INFINITY, 0_u32); num_leaves + 1]; self.get_nodes().len()];
        let mut delta_hat_sets: Vec<Vec<Vec<usize>>> =
            vec![vec![vec![]; num_leaves + 1]; self.get_nodes().len()];

        for node_id in self.get_node_ids() {
            delta_bar[node_id][0] = (0_f32, 0_u32);
            delta_hat[node_id][0] = (0_f32, 0_u32);
            if self.is_leaf(node_id){
                delta_bar[node_id][1] = (0_f32, 0_u32);
                delta_hat[node_id][1] = (0_f32, 0_u32);
                delta_bar_sets[node_id][1] = vec![node_id];
                delta_hat_sets[node_id][1] = vec![node_id];
            }
        }

        for node_id in self.postord_ids(self.get_root_id()) {
            if !self.is_leaf(node_id) {
                for i in 1..(std::cmp::min(num_leaves, self.get_cluster_size(node_id)) + 1)
                {
                    // Min pd
                    let mut min_bar = match op{
                        Ordering::Greater => f32::MIN,
                        _ => f32::MAX,
                    };
                    let mut min_hat = match op{
                        Ordering::Greater => f32::MIN,
                        _ => f32::MAX,
                    };
                    // min pd set
                    let mut min_bar_set: Vec<usize> = vec![];
                    // min norm pd set
                    let mut min_hat_set: Vec<usize> = vec![];
                    // min pd edge_count
                    let mut min_e_bar = 0_u32;
                    // min norm pd edge count
                    let mut min_e_hat = 0_u32;
                    let node_children = self.get_node_children_ids(node_id).collect_vec();
                    let x = node_children[0];
                    let y = node_children[1];
                    for r in 0..i + 1 {
                        let l = i - r;
                        if delta_bar[y][r].0 == f32::INFINITY || delta_bar[x][l].0 == f32::INFINITY
                        {
                            continue;
                        }
                        let val_bar = delta_bar[y][r].0
                            + (self.get_node(y).unwrap().get_weight().unwrap()
                                * (std::cmp::min(r, 1) as f32))
                            + delta_bar[x][l].0
                            + (self.get_node(x).unwrap().get_weight().unwrap()
                                * (std::cmp::min(l, 1) as f32));
                        let mut e = 0_u32;
                        let mut set = vec![];
                        if l == 0 {
                            e += delta_bar[y][r].1 + 1;
                            set = delta_bar_sets[y][r].clone();
                        } else if r == 0 {
                            e += delta_bar[x][l].1 + 1;
                            set = delta_bar_sets[x][l].clone();
                        } else {
                            e += delta_bar[x][l].1 + delta_bar[y][r].1 + 2;
                            set = delta_bar_sets[x][l].clone();
                            set.extend(delta_bar_sets[y][r].clone());
                        }
                        let val_hat = val_bar / (e as f32);
                        if op_fn(&val_bar, &min_bar) {
                            min_bar = val_bar;
                            min_e_bar = e;
                            min_bar_set = set.clone();
                        }
                        if op_fn(&val_hat, &min_hat) {
                            min_hat = val_hat;
                            min_e_hat = e;
                            min_hat_set = set;
                        }
                    }
                    delta_bar[node_id][i] = (min_bar, min_e_bar);
                    delta_hat[node_id][i] = (min_hat, min_e_hat);
                    delta_bar_sets[node_id][i] = min_bar_set;
                    delta_hat_sets[node_id][i] = min_hat_set;
                }
            }
        }

        (delta_bar, delta_bar_sets, delta_hat, delta_hat_sets)
    }

    // fn precompute_minPDs(&mut self) {
    //     let (delta_bar, delta_bar_sets, delta_hat, delta_hat_sets) = self.compute_norm_min();
    //     self.precomputed_min = delta_bar;
    //     self.precomputed_norm_min = delta_hat;
    //     self.precomputed_min_set = delta_bar_sets;
    //     self.precomputed_norm_min_set = delta_hat_sets;
    // }

    // fn compute_norm_min(
    //     &self,
    // ) -> (
    //     Vec<Vec<(f32, u32)>>,
    //     Vec<Vec<Vec<usize>>>,
    //     Vec<Vec<(f32, u32)>>,
    //     Vec<Vec<Vec<usize>>>,
    // ) {
    //     let num_leaves = self.tree.get_leaves().len();
    //     let mut delta_bar =
    //         vec![vec![(f32::INFINITY, 0_u32); num_leaves + 1]; self.tree.get_nodes().len()];
    //     let mut delta_bar_sets: Vec<Vec<Vec<usize>>> =
    //         vec![vec![vec![]; num_leaves + 1]; self.tree.get_nodes().len()];
    //     let mut delta_hat =
    //         vec![vec![(f32::INFINITY, 0_u32); num_leaves + 1]; self.tree.get_nodes().len()];
    //     let mut delta_hat_sets: Vec<Vec<Vec<usize>>> =
    //         vec![vec![vec![]; num_leaves + 1]; self.tree.get_nodes().len()];

    //     for node_id in self.tree.get_node_ids() {
    //         delta_bar[node_id][0] = (0_f32, 0_u32);
    //         delta_hat[node_id][0] = (0_f32, 0_u32);
    //         if self.tree.is_leaf(node_id){
    //             delta_bar[node_id][1] = (0_f32, 0_u32);
    //             delta_hat[node_id][1] = (0_f32, 0_u32);
    //             delta_bar_sets[node_id][1] = vec![node_id];
    //             delta_hat_sets[node_id][1] = vec![node_id];
    //         }
    //     }

    //     for node_id in self.tree.postord_ids(self.tree.get_root_id()) {
    //         if !self.tree.is_leaf(node_id) {
    //             for i in 1..(std::cmp::min(num_leaves, self.tree.get_cluster_size(node_id)) + 1)
    //             {
    //                 // Min pd
    //                 let mut min_bar = f32::MAX;
    //                 // Min norm pd
    //                 let mut min_hat = f32::MAX;
    //                 // min pd set
    //                 let mut min_bar_set: Vec<usize> = vec![];
    //                 // min norm pd set
    //                 let mut min_hat_set: Vec<usize> = vec![];
    //                 // min pd edge_count
    //                 let mut min_e_bar = 0_u32;
    //                 // min norm pd edge count
    //                 let mut min_e_hat = 0_u32;
    //                 let node_children = self.tree.get_node_children_ids(node_id).collect_vec();
    //                 let x = node_children[0];
    //                 let y = node_children[1];
    //                 for r in 0..i + 1 {
    //                     let l = i - r;
    //                     if delta_bar[y][r].0 == f32::INFINITY || delta_bar[x][l].0 == f32::INFINITY
    //                     {
    //                         continue;
    //                     }
    //                     let val_bar = delta_bar[y][r].0
    //                         + (self.tree.get_node(y).unwrap().get_weight().unwrap()
    //                             * (std::cmp::min(r, 1) as f32))
    //                         + delta_bar[x][l].0
    //                         + (self.tree.get_node(x).unwrap().get_weight().unwrap()
    //                             * (std::cmp::min(l, 1) as f32));
    //                     let mut e = 0_u32;
    //                     let mut set = vec![];
    //                     if l == 0 {
    //                         e += delta_bar[y][r].1 + 1;
    //                         set = delta_bar_sets[y][r].clone();
    //                     } else if r == 0 {
    //                         e += delta_bar[x][l].1 + 1;
    //                         set = delta_bar_sets[x][l].clone();
    //                     } else {
    //                         e += delta_bar[x][l].1 + delta_bar[y][r].1 + 2;
    //                         set = delta_bar_sets[x][l].clone();
    //                         set.extend(delta_bar_sets[y][r].clone());
    //                     }
    //                     let val_hat = val_bar / (e as f32);
    //                     if val_bar < min_bar {
    //                         min_bar = val_bar;
    //                         min_e_bar = e;
    //                         min_bar_set = set.clone();
    //                     }
    //                     if val_hat < min_hat {
    //                         min_hat = val_hat;
    //                         min_e_hat = e;
    //                         min_hat_set = set;
    //                     }
    //                 }
    //                 delta_bar[node_id][i] = (min_bar, min_e_bar);
    //                 delta_hat[node_id][i] = (min_hat, min_e_hat);
    //                 delta_bar_sets[node_id][i] = min_bar_set;
    //                 delta_hat_sets[node_id][i] = min_hat_set;
    //             }
    //         }
    //     }

    //     (delta_bar, delta_bar_sets, delta_hat, delta_hat_sets)
    // }

    // fn get_minPD(
    //     &self,
    //     num_taxa: usize,
    // ) -> TreeNodeWeight<'a, Self::Tree>
    // {
    //     self.precomputed_min[self.tree.get_root_id()][std::cmp::min(num_taxa, self.tree.num_taxa())]
    //         .0
    // }

    // fn get_minPD_node(
    //     &self,
    //     node_id: TreeNodeID<'a, Self::Tree>,
    //     num_taxa: usize,
    // ) -> TreeNodeWeight<'a, Self::Tree>
    // {
    //     self.precomputed_min[node_id][std::cmp::min(num_taxa, self.tree.num_taxa())].0
    // }

    // fn get_norm_minPD(
    //     &self,
    //     num_taxa: usize,
    // ) -> TreeNodeWeight<'a, Self::Tree>
    // {
    //     self.precomputed_norm_min[self.tree.get_root_id()]
    //         [std::cmp::min(num_taxa, self.tree.num_taxa())]
    //     .0
    // }

    // fn get_minPD_taxa_set(
    //     &self,
    //     num_taxa: usize,
    // ) -> impl Iterator<Item = TreeNodeID<'a, Self::Tree>> {
    //     self.precomputed_min_set[self.tree.get_root_id()][num_taxa]
    //         .clone()
    //         .into_iter()
    // }

    // fn get_norm_minPD_taxa_set(
    //     &self,
    //     num_taxa: usize,
    // ) -> impl Iterator<Item = TreeNodeID<'a, Self::Tree>> {
    //     self.precomputed_norm_min_set[self.tree.get_root_id()][num_taxa]
    //         .clone()
    //         .into_iter()
    // }

    // fn get_min_genPD(
    //     &self,
    // ) -> TreeNodeWeight<'a, Self::Tree>
    // {
    //     self.precomputed_norm_min[self.tree.get_root_id()]
    //         .iter()
    //         .filter(|x| x.0 != 0_f32 && x.0 != f32::INFINITY)
    //         .min_by(|x, y| x.0.total_cmp(&y.0))
    //         .unwrap()
    //         .0
    // }

    // fn get_min_genPD_set(
    //     &self,
    // ) -> impl Iterator<Item = TreeNodeID<'a, Self::Tree>> {
    //     let num_taxa = self.precomputed_norm_min[self.tree.get_root_id()]
    //         .iter()
    //         .enumerate()
    //         .filter(|x| x.1 .0 != 0_f32 && x.1 .0 != f32::INFINITY)
    //         .min_by(|x, y| x.1 .0.total_cmp(&y.1 .0))
    //         .unwrap()
    //         .0;
    //     self.precomputed_norm_min_set[self.tree.get_root_id()][num_taxa]
    //         .clone()
    //         .into_iter()
    // }

    // fn precompute_maxPDs(&mut self) {
    //     let (delta_bar, delta_bar_sets, delta_hat, delta_hat_sets) = self.compute_norm_max();
    //     self.precomputed_max = delta_bar;
    //     self.precomputed_norm_max = delta_hat;
    //     self.precomputed_max_set = delta_bar_sets;
    //     self.precomputed_norm_max_set = delta_hat_sets;
    // }

    // fn compute_norm_max(
    //     &self,
    // ) -> (
    //     Vec<Vec<(f32, u32)>>,
    //     Vec<Vec<Vec<usize>>>,
    //     Vec<Vec<(f32, u32)>>,
    //     Vec<Vec<Vec<usize>>>,
    // ) {
    //     let num_leaves = self.tree.get_leaves().len();
    //     let mut delta_bar =
    //         vec![vec![(f32::INFINITY, 0_u32); num_leaves + 1]; self.tree.get_nodes().len()];
    //     let mut delta_bar_sets: Vec<Vec<Vec<usize>>> =
    //         vec![vec![vec![]; num_leaves + 1]; self.tree.get_nodes().len()];
    //     let mut delta_hat =
    //         vec![vec![(f32::INFINITY, 0_u32); num_leaves + 1]; self.tree.get_nodes().len()];
    //     let mut delta_hat_sets: Vec<Vec<Vec<usize>>> =
    //         vec![vec![vec![]; num_leaves + 1]; self.tree.get_nodes().len()];

    //     for node_id in self.tree.get_node_ids() {
    //         delta_bar[node_id][0] = (0_f32, 0_u32);
    //         delta_hat[node_id][0] = (0_f32, 0_u32);
    //         if self.tree.is_leaf(node_id) {
    //             delta_bar[node_id][1] = (0_f32, 0_u32);
    //             delta_hat[node_id][1] = (0_f32, 0_u32);
    //             delta_bar_sets[node_id][1] = vec![node_id];
    //             delta_hat_sets[node_id][1] = vec![node_id];
    //         }
    //     }

    //     for node_id in self.tree.postord_ids(self.tree.get_root_id()) {
    //         if !self.tree.is_leaf(node_id) {
    //             for i in 1..std::cmp::max(num_leaves, self.tree.get_cluster_size(node_id)) + 1
    //             {
    //                 let mut max_bar = 0_f32;
    //                 let mut max_hat = 0_f32;
    //                 let mut max_bar_set: Vec<usize> = vec![];
    //                 let mut max_hat_set: Vec<usize> = vec![];
    //                 let mut max_e_bar = 0_u32;
    //                 let mut max_e_hat = 0_u32;
    //                 let node_children = self.tree.get_node_children_ids(node_id).collect_vec();
    //                 let x = node_children[0];
    //                 let y = node_children[1];
    //                 for r in 0..i + 1 {
    //                     let l = i - r;
    //                     if delta_bar[y][r].0 == f32::INFINITY || delta_bar[x][l].0 == f32::INFINITY
    //                     {
    //                         continue;
    //                     }
    //                     let val_bar = delta_bar[y][r].0
    //                         + (self.tree.get_node(y).unwrap().get_weight().unwrap()
    //                             * (std::cmp::min(r, 1) as f32))
    //                         + delta_bar[x][l].0
    //                         + (self.tree.get_node(x).unwrap().get_weight().unwrap()
    //                             * (std::cmp::min(l, 1) as f32));
    //                     let mut e = 0_u32;
    //                     let mut set: Vec<usize>;
    //                     if l == 0 {
    //                         e += delta_bar[y][r].1 + 1;
    //                         set = delta_bar_sets[y][r].clone();
    //                     } else if r == 0 {
    //                         e += delta_bar[x][l].1 + 1;
    //                         set = delta_bar_sets[x][l].clone();
    //                     } else {
    //                         e += delta_bar[x][l].1 + delta_bar[y][r].1 + 2;
    //                         set = delta_bar_sets[x][l].clone();
    //                         set.extend(delta_bar_sets[y][r].clone());
    //                     }
    //                     let val_hat = val_bar / (e as f32);
    //                     if val_bar > max_bar {
    //                         max_bar = val_bar;
    //                         max_e_bar = e;
    //                         max_bar_set = set.clone();
    //                     }
    //                     if val_hat > max_hat {
    //                         max_hat = val_hat;
    //                         max_e_hat = e;
    //                         max_hat_set = set;
    //                     }
    //                 }
    //                 delta_bar[node_id][i] = (max_bar, max_e_bar);
    //                 delta_hat[node_id][i] = (max_hat, max_e_hat);
    //                 delta_bar_sets[node_id][i] = max_bar_set;
    //                 delta_hat_sets[node_id][i] = max_hat_set;
    //             }
    //         }
    //     }

    //     (delta_bar, delta_bar_sets, delta_hat, delta_hat_sets)
    // }

    // fn get_maxPD(
    //     &self,
    //     num_taxa: usize,
    // ) -> TreeNodeWeight<'a, Self::Tree>
    // {
    //     self.precomputed_max[self.tree.get_root_id()][std::cmp::max(num_taxa, 1)]
    //         .0
    // }

    // fn get_maxPD_node(
    //     &self,
    //     node_id: TreeNodeID<'a, Self::Tree>,
    //     num_taxa: usize,
    // ) -> TreeNodeWeight<'a, Self::Tree>
    // {
    //     self.precomputed_max[node_id][std::cmp::max(num_taxa, self.tree.num_taxa()-1)].0
    // }

    // fn get_norm_maxPD(
    //     &self,
    //     num_taxa: usize,
    // ) -> TreeNodeWeight<'a, Self::Tree>
    // {
    //     self.precomputed_norm_max[self.tree.get_root_id()]
    //         [std::cmp::max(num_taxa, self.tree.num_taxa()-1)]
    //     .0
    // }

    // fn get_maxPD_taxa_set(
    //     &self,
    //     num_taxa: usize,
    // ) -> impl Iterator<Item = TreeNodeID<'a, Self::Tree>> {
    //     self.precomputed_max_set[self.tree.get_root_id()][num_taxa]
    //         .clone()
    //         .into_iter()
    // }

    // fn get_norm_maxPD_taxa_set(
    //     &self,
    //     num_taxa: usize,
    // ) -> impl Iterator<Item = TreeNodeID<'a, Self::Tree>> {
    //     self.precomputed_norm_max_set[self.tree.get_root_id()][num_taxa]
    //         .clone()
    //         .into_iter()
    // }

    // fn get_max_genPD(
    //     &self,
    // ) -> TreeNodeWeight<'a, Self::Tree>
    // {
    //     self.precomputed_norm_max[self.tree.get_root_id()]
    //         .iter()
    //         .filter(|x| x.0 != 0_f32 && x.0 != f32::INFINITY)
    //         .max_by(|x, y| x.0.total_cmp(&y.0))
    //         .unwrap()
    //         .0
    // }

    // fn get_max_genPD_set(
    //     &self,
    // ) -> impl Iterator<Item = TreeNodeID<'a, Self::Tree>> {
    //     let num_taxa = self.precomputed_norm_max[self.tree.get_root_id()]
    //         .iter()
    //         .enumerate()
    //         .filter(|x| x.1 .0 != 0_f32 && x.1 .0 != f32::INFINITY)
    //         .max_by(|x, y| x.1 .0.total_cmp(&y.1 .0))
    //         .unwrap()
    //         .0;
    //     self.precomputed_norm_max_set[self.tree.get_root_id()][num_taxa]
    //         .clone()
    //         .into_iter()
    // }
}

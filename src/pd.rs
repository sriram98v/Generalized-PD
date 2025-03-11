pub mod phylogenetic_diversity;
use crate::pd::phylogenetic_diversity::{RootedPhylogeneticDiversity, TreePDMap, binary_splits};
use itertools::Itertools;
use phylo::prelude::*;
use phylo::tree::SimpleRootedTree;
use phylogenetic_diversity::pascal_triangle;
use std::cmp::{Ordering, min, max};

pub struct TreePD<'a,T:NodeTaxa,W:EdgeWeight,Z:NodeWeight> {
    tree: &'a SimpleRootedTree<T,W,Z>,
    precomputed_min: Vec<Vec<(W, u32)>>,
    precomputed_min_set: Vec<Vec<Vec<usize>>>,
    precomputed_norm_min: Vec<Vec<(W, u32)>>,
    precomputed_norm_min_set: Vec<Vec<Vec<usize>>>,
    precomputed_max: Vec<Vec<(W, u32)>>,
    precomputed_max_set: Vec<Vec<Vec<usize>>>,
    precomputed_norm_max: Vec<Vec<(W, u32)>>,
    precomputed_norm_max_set: Vec<Vec<Vec<usize>>>,
    precomputed_avg: Vec<Vec<W>>,
}

impl<'a,T:NodeTaxa,W:EdgeWeight,Z:NodeWeight> TreePD<'a,T,W,Z> {
    pub fn new(tree: &'a SimpleRootedTree<T,W,Z>) -> Self {
        let (min, min_set, min_norm, min_norm_set) = tree.compute_norm_min();
        let (max, max_set, max_norm, max_norm_set) = tree.compute_norm_max();
        let avg = tree.compute_avg();
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
            precomputed_avg: avg,
        }
    }
}
impl<'a,T:NodeTaxa,W:EdgeWeight,Z:NodeWeight> TreePDMap for TreePD<'a,T,W,Z> {

    type Tree = SimpleRootedTree<T,W,Z>;

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
    ) -> TreeNodeWeight<Self::Tree>
    {
        self.precomputed_min[self.tree.get_root_id()][min(num_taxa, self.tree.num_taxa())]
            .0
    }

    fn get_minPD_node(
        &self,
        node_id: TreeNodeID<Self::Tree>,
        num_taxa: usize,
    ) -> TreeNodeWeight<Self::Tree>
    {
        self.precomputed_min[node_id][min(num_taxa, self.tree.num_taxa())].0
    }

    fn get_norm_minPD(
        &self,
        num_taxa: usize,
    ) -> TreeNodeWeight<Self::Tree>
    {
        self.precomputed_norm_min[self.tree.get_root_id()]
            [min(num_taxa, self.tree.num_taxa())]
        .0
    }

    fn get_minPD_taxa_set(
        &self,
        num_taxa: usize,
    ) -> impl Iterator<Item = TreeNodeID<Self::Tree>> {
        self.precomputed_min_set[self.tree.get_root_id()][num_taxa]
            .clone()
            .into_iter()
    }

    fn get_norm_minPD_taxa_set(
        &self,
        num_taxa: usize,
    ) -> impl Iterator<Item = TreeNodeID<Self::Tree>> {
        self.precomputed_norm_min_set[self.tree.get_root_id()][num_taxa]
            .clone()
            .into_iter()
    }

    fn get_min_genPD(
        &self,
    ) -> TreeNodeWeight<Self::Tree>
    {
        self.precomputed_norm_min[self.tree.get_root_id()]
            .iter()
            .enumerate()
            .filter(|(x,_)| *x > 2)
            .filter(|(_, x)| x.0 != W::zero() && x.0 != W::infinity())
            .min_by(|(_, x), (_, y)| x.0.partial_cmp(&y.0).unwrap())
            .unwrap()
            .1.0
    }

    fn get_min_genPD_set(
        &self,
    ) -> impl Iterator<Item = TreeNodeID<Self::Tree>> {
        let num_taxa = self.precomputed_norm_min[self.tree.get_root_id()]
            .iter()
            .enumerate()
            .filter(|x| x.0 > 2 && x.1.0 != W::zero() && x.1 .0 != W::infinity())
            .min_by(|x, y| x.1 .0.partial_cmp(&y.1 .0).unwrap())
            .unwrap()
            .0;
        self.precomputed_norm_min_set[self.tree.get_root_id()][num_taxa]
            .clone()
            .into_iter()
    }

    fn get_maxPD(
        &self,
        num_taxa: usize,
    ) -> TreeNodeWeight<Self::Tree>
    {
        self.precomputed_max[self.tree.get_root_id()][max(num_taxa, 1)]
            .0
    }

    fn get_maxPD_node(
        &self,
        node_id: TreeNodeID<Self::Tree>,
        num_taxa: usize,
    ) -> TreeNodeWeight<Self::Tree>
    {
        self.precomputed_max[node_id][max(num_taxa, self.tree.num_taxa()-1)].0
    }

    fn get_norm_maxPD(
        &self,
        num_taxa: usize,
    ) -> TreeNodeWeight<Self::Tree>
    {
        self.precomputed_norm_max[self.tree.get_root_id()]
            [max(num_taxa, self.tree.num_taxa()-1)]
        .0
    }

    fn get_maxPD_taxa_set(
        &self,
        num_taxa: usize,
    ) -> impl Iterator<Item = TreeNodeID<Self::Tree>> {
        self.precomputed_max_set[self.tree.get_root_id()][num_taxa]
            .clone()
            .into_iter()
    }

    fn get_norm_maxPD_taxa_set(
        &self,
        num_taxa: usize,
    ) -> impl Iterator<Item = TreeNodeID<Self::Tree>> {
        self.precomputed_norm_max_set[self.tree.get_root_id()][num_taxa]
            .clone()
            .into_iter()
    }

    fn get_max_genPD(
        &self,
    ) -> TreeNodeWeight<Self::Tree>
    {
        self.precomputed_norm_max[self.tree.get_root_id()]
            .iter()
            .enumerate()
            .filter(|(x,_)| *x > 2)
            .filter(|(_, x)| x.0 != W::zero() && x.0 != W::infinity())
            .max_by(|(_, x), (_, y)| x.0.partial_cmp(&y.0).unwrap())
            .unwrap()
            .1.0
    }

    fn get_max_genPD_set(
        &self,
    ) -> impl Iterator<Item = TreeNodeID<Self::Tree>> {
        let num_taxa = self.precomputed_norm_max[self.tree.get_root_id()]
            .iter()
            .enumerate()
            .filter(|x| x.0 > 2 && x.1.0 != W::zero() && x.1 .0 != W::infinity())
            .max_by(|x, y| x.1 .0.partial_cmp(&y.1 .0).unwrap())
            .unwrap()
            .0;
        self.precomputed_norm_max_set[self.tree.get_root_id()][num_taxa]
            .clone()
            .into_iter()
    }

    fn get_avgPD(
            &self,
            num_taxa: usize,
        ) -> TreeNodeWeight<Self::Tree> {
        self.precomputed_avg[self.tree.get_root_id()][num_taxa]
    }

}

impl<T:NodeTaxa,W:EdgeWeight,Z:NodeWeight> RootedPhylogeneticDiversity for SimpleRootedTree<T,W,Z> {

    fn compute_dp_table(
        &self,
        op: Ordering,
    ) -> (
        Vec<Vec<(W, u32)>>,
        Vec<Vec<Vec<TreeNodeID<Self>>>>,
        Vec<Vec<(W, u32)>>,
        Vec<Vec<Vec<TreeNodeID<Self>>>>,
    )
    {
        let op_fn = match op{
            Ordering::Greater => PartialOrd::ge,
            _ => PartialOrd::le,
        };

        let start_val = match op{
            Ordering::Greater => W::min_value(),
            _ => W::infinity(),
        };

        let num_leaves = self.get_leaves().len();
        let mut delta_bar: Vec<Vec<(W, u32)>> =
            vec![vec![(start_val, 0_u32); num_leaves + 1]; self.get_nodes().len()];
        let mut delta_bar_sets: Vec<Vec<Vec<usize>>> =
            vec![vec![vec![]; num_leaves + 1]; self.get_nodes().len()];
        let mut delta_hat: Vec<Vec<(W, u32)>> =
            vec![vec![(start_val, 0_u32); num_leaves + 1]; self.get_nodes().len()];
        let mut delta_hat_sets: Vec<Vec<Vec<usize>>> =
            vec![vec![vec![]; num_leaves + 1]; self.get_nodes().len()];

        for node_id in self.get_node_ids() {
            delta_bar[node_id][0] = (W::zero(), 0_u32);
            delta_hat[node_id][0] = (W::zero(), 0_u32);
            if self.is_leaf(node_id){
                delta_bar[node_id][1] = (W::zero(), 0_u32);
                delta_hat[node_id][1] = (W::zero(), 0_u32);
                delta_bar_sets[node_id][1] = vec![node_id];
                delta_hat_sets[node_id][1] = vec![node_id];
            }
        }

        for node_id in self.postord_ids(self.get_root_id()) {
            if !self.is_leaf(node_id) {
                for i in 1..(min(num_leaves, self.get_cluster_size(node_id)) + 1)
                {
                    // Min pd
                    let mut min_bar = match op{
                        Ordering::Greater => W::min_value(),
                        _ => W::max_value(),
                    };
                    let mut min_hat = match op{
                        Ordering::Greater => W::min_value(),
                        _ => W::max_value(),
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
                    let x_cluster_size = self.get_cluster_size(x);
                    let y_cluster_size = self.get_cluster_size(y);
                    for (l,r) in binary_splits(i, x_cluster_size, y_cluster_size){
                        let discount_l = self.get_node(y).unwrap().get_weight().is_some() as u32;
                        let discount_r = self.get_node(x).unwrap().get_weight().is_some() as u32;
                        let val_bar: W = delta_bar[y][r].0
                            + (self.get_node(y).unwrap().get_weight().unwrap_or(W::zero())*W::from(min(r, 1)).unwrap())
                            + delta_bar[x][l].0
                            + (self.get_node(x).unwrap().get_weight().unwrap_or(W::zero())
                                * W::from(min(l, 1)).unwrap());
                        let mut e = 0_u32;
                        let mut set: Vec<usize>;
                        if l == 0 {
                            e += delta_bar[y][r].1;
                            set = delta_bar_sets[y][r].clone();
                        } else if r == 0 {
                            e += delta_bar[x][l].1;
                            set = delta_bar_sets[x][l].clone();
                        } else {
                            e += delta_bar[x][l].1 + delta_bar[y][r].1 + discount_r +discount_l;
                            set = delta_bar_sets[x][l].clone();
                            set.extend(delta_bar_sets[y][r].clone());
                        }
                        let val_hat = val_bar / W::from(e).unwrap();
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

    fn compute_avg(
            &self,
        ) -> Vec<Vec<TreeNodeWeight<Self>>> {
        let num_leaves = self.num_taxa();
        let pascal = pascal_triangle(num_leaves as u32);
        let mut alpha = vec![vec![W::zero(); num_leaves + 1]; self.get_nodes().len()];
        let mut beta = vec![vec![W::zero(); num_leaves + 1]; self.get_nodes().len()];

        for node_id in self.postord_ids(self.get_root_id()){
            match self.is_leaf(node_id) {
                true => {
                    alpha[node_id][1] = W::zero();
                    beta[node_id][1] = W::zero();
                },
                _ => {
                    let node_children = self.get_node_children_ids(node_id).collect_vec();
                    let x = node_children[0];
                    let y = node_children[1];
                    let x_cluster_size = self.get_cluster_size(x);
                    let y_cluster_size = self.get_cluster_size(y);
                    let w_x = self.get_node(x).unwrap().get_weight().unwrap_or(W::zero());
                    let w_y = self.get_node(y).unwrap().get_weight().unwrap_or(W::zero());
                    for i in 1..(min(num_leaves, self.get_cluster_size(node_id)) + 1){
                        let s_x = beta[x][i] + W::from(pascal[x_cluster_size][i]).unwrap()*w_x;
                        let s_y = beta[y][i] + W::from(pascal[y_cluster_size][i]).unwrap()*w_y;
                        let mut s = s_x + s_y;
                        for r in max(1,i as isize-y_cluster_size as isize)..min(x_cluster_size as isize,i as isize-1)+1{
                            let l = (i as isize-r) as usize;
                            s = s + W::from(pascal[y_cluster_size][l]).unwrap()*beta[x][r as usize]
                                + W::from(pascal[x_cluster_size][r as usize]).unwrap()*beta[y][l]
                                + W::from(pascal[y_cluster_size][l]).unwrap()*W::from(pascal[x_cluster_size][r as usize]).unwrap()*(w_x+w_y);
                        }
                        beta[node_id][i] = s;
                        alpha[node_id][i] = s/W::from(pascal[x_cluster_size+y_cluster_size][i]).unwrap();
                    }
                },
            };
        }
        alpha
    }

}

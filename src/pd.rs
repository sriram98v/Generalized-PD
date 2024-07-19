pub mod phylogenetic_diversity;
use crate::pd::phylogenetic_diversity::PhylogeneticDiversity;
use itertools::Itertools;
use phylo::prelude::*;
use phylo::tree::SimpleRootedTree;

pub struct TreePDMap<'a> {
    tree: &'a SimpleRootedTree,
    precomputed_min: Vec<Vec<(f32, u32)>>,
    precomputed_min_set: Vec<Vec<Vec<usize>>>,
    precomputed_norm_min: Vec<Vec<(f32, u32)>>,
    precomputed_norm_min_set: Vec<Vec<Vec<usize>>>,
}

impl<'a> TreePDMap<'a> {
    pub fn new(tree: &'a SimpleRootedTree) -> Self {
        TreePDMap {
            tree,
            precomputed_min: vec![],
            precomputed_min_set: vec![],
            precomputed_norm_min: vec![],
            precomputed_norm_min_set: vec![],
        }
    }
}

impl<'a> PhylogeneticDiversity<'a> for TreePDMap<'a> {
    type Tree = SimpleRootedTree;

    fn get_tree(&self) -> &Self::Tree {
        self.tree
    }

    fn precompute_minPDs(&mut self) {
        let (delta_bar, delta_bar_sets, delta_hat, delta_hat_sets) = self.compute_norm_min();
        self.precomputed_min = delta_bar;
        self.precomputed_norm_min = delta_hat;
        self.precomputed_min_set = delta_bar_sets;
        self.precomputed_norm_min_set = delta_hat_sets;
    }

    fn compute_norm_min(
        &self,
    ) -> (
        Vec<Vec<(f32, u32)>>,
        Vec<Vec<Vec<usize>>>,
        Vec<Vec<(f32, u32)>>,
        Vec<Vec<Vec<usize>>>,
    ) {
        let num_leaves = self.tree.get_leaves().len();
        let mut delta_bar =
            vec![vec![(f32::INFINITY, 0_u32); num_leaves + 1]; self.tree.get_nodes().len()];
        let mut delta_bar_sets: Vec<Vec<Vec<usize>>> =
            vec![vec![vec![]; num_leaves + 1]; self.tree.get_nodes().len()];
        let mut delta_hat =
            vec![vec![(f32::INFINITY, 0_u32); num_leaves + 1]; self.tree.get_nodes().len()];
        let mut delta_hat_sets: Vec<Vec<Vec<usize>>> =
            vec![vec![vec![]; num_leaves + 1]; self.tree.get_nodes().len()];

        for node in self.tree.get_nodes() {
            delta_bar[node.get_id()][0] = (0_f32, 0_u32);
            delta_hat[node.get_id()][0] = (0_f32, 0_u32);
            if node.is_leaf() {
                delta_bar[node.get_id()][1] = (0_f32, 0_u32);
                delta_hat[node.get_id()][1] = (0_f32, 0_u32);
                delta_bar_sets[node.get_id()][1] = vec![node.get_id()];
                delta_hat_sets[node.get_id()][1] = vec![node.get_id()];
            }
        }

        for node_id in self.tree.postord_ids(self.tree.get_root_id()) {
            if !self.tree.is_leaf(node_id) {
                for i in 1..std::cmp::min(num_leaves, self.tree.get_cluster_size(node_id)) + 1
                {
                    let mut min_bar = f32::INFINITY;
                    let mut min_hat = f32::INFINITY;
                    let mut min_bar_set: Vec<usize> = vec![];
                    let mut min_hat_set: Vec<usize> = vec![];
                    let mut min_e_bar = 0_u32;
                    let mut min_e_hat = 0_u32;
                    let node_children = self.tree.get_node_children_ids(node_id).collect_vec();
                    let x = node_children[0];
                    let y = node_children[1];
                    for r in 0..i + 1 {
                        let l = i - r;
                        if delta_bar[y][r].0 == f32::INFINITY || delta_bar[x][l].0 == f32::INFINITY
                        {
                            continue;
                        }
                        let val_bar = delta_bar[y][r].0
                            + (self.tree.get_node(y).unwrap().get_weight().unwrap()
                                * (std::cmp::min(r, 1) as f32))
                            + delta_bar[x][l].0
                            + (self.tree.get_node(x).unwrap().get_weight().unwrap()
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
                        if val_bar < min_bar {
                            min_bar = val_bar;
                            min_e_bar = e;
                            min_bar_set = set.clone();
                        }
                        if val_hat < min_hat {
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

    fn get_minPD(
        &self,
        num_taxa: usize,
    ) -> TreeNodeWeight<'a, Self::Tree>
    {
        self.precomputed_min[self.tree.get_root_id()][std::cmp::min(num_taxa, self.tree.num_taxa())]
            .0
    }

    fn get_minPD_node(
        &self,
        node_id: TreeNodeID<'a, Self::Tree>,
        num_taxa: usize,
    ) -> TreeNodeWeight<'a, Self::Tree>
    {
        self.precomputed_min[node_id][std::cmp::min(num_taxa, self.tree.num_taxa())].0
    }

    fn get_norm_minPD(
        &self,
        num_taxa: usize,
    ) -> TreeNodeWeight<'a, Self::Tree>
    {
        self.precomputed_norm_min[self.tree.get_root_id()]
            [std::cmp::min(num_taxa, self.tree.num_taxa())]
        .0
    }

    fn backtrack_min(
        &self,
        node_id: TreeNodeID<'a, Self::Tree>,
        num_taxa: usize,
        taxaset: &mut Vec<TreeNodeID<'a, Self::Tree>>,
    ) {
        if self.tree.is_leaf(node_id) {
            taxaset.push(node_id);
        } else {
            let node_children = self.tree.get_node_children_ids(node_id).collect_vec();
            let left_child_id = node_children[0];
            let right_child_id = node_children[1];
            let left_edge_weight = self.tree.get_edge_weight(node_id, left_child_id).unwrap();
            let right_edge_weight = self.tree.get_edge_weight(node_id, right_child_id).unwrap();
            if self.get_minPD_node(node_id, num_taxa)
                == self.get_minPD_node(left_child_id, num_taxa) + left_edge_weight
            {
                if self.tree.is_leaf(left_child_id) {
                    taxaset.push(left_child_id);
                } else {
                    self.backtrack_min(left_child_id, num_taxa, taxaset)
                }
            } else if self.get_minPD_node(node_id, num_taxa)
                == self.get_minPD_node(right_child_id, num_taxa) + right_edge_weight
            {
                if self.tree.is_leaf(right_child_id) {
                    taxaset.push(right_child_id);
                } else {
                    self.backtrack_min(right_child_id, num_taxa, taxaset)
                }
            } else {
                // the following line needs to be checked!
                let ((left_pd, left_edges), (right_pd, right_edges)) = (1..num_taxa)
                    .map(|l| (l, num_taxa - l))
                    .map(|(l, r)| {
                        (
                            self.precomputed_min[left_child_id]
                                [std::cmp::min(l, self.tree.num_taxa())],
                            self.precomputed_min[right_child_id]
                                [std::cmp::min(r, self.tree.num_taxa())],
                        )
                    })
                    .filter(|(l, r)| l.0 != f32::INFINITY && r.0 != f32::INFINITY)
                    // .inspect(|x| println!("{:?}", &x))
                    .min_by(|x, y| (x.0 .0 + x.1 .0).total_cmp(&(y.0 .0 + y.1 .0)))
                    .unwrap();
                let num_left_taxa = (left_edges / 2) as usize + 1;
                let num_right_taxa = (right_edges / 2) as usize + 1;
                if num_left_taxa > 0 {
                    self.backtrack_min(left_child_id, num_left_taxa, taxaset);
                }
                if num_right_taxa > 0 {
                    self.backtrack_min(right_child_id, num_right_taxa, taxaset);
                }
            }
        }
    }

    fn get_minPD_taxa_set(
        &self,
        num_taxa: usize,
    ) -> impl Iterator<Item = TreeNodeID<'a, Self::Tree>> {
        self.precomputed_min_set[self.tree.get_root_id()][num_taxa]
            .clone()
            .into_iter()
    }

    fn get_norm_minPD_taxa_set(
        &self,
        num_taxa: usize,
    ) -> impl Iterator<Item = TreeNodeID<'a, Self::Tree>> {
        self.precomputed_norm_min_set[self.tree.get_root_id()][num_taxa]
            .clone()
            .into_iter()
    }

    fn get_min_genPD(
        &self,
    ) -> TreeNodeWeight<'a, Self::Tree>
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
    ) -> impl Iterator<Item = TreeNodeID<'a, Self::Tree>> {
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
}

use itertools::Itertools;
use phylo::prelude::*;

use std::{cmp::{max, min, Ordering}, collections::VecDeque};

pub trait TreePDMap
where 
    <Self::Tree as RootedTree>::Node: RootedWeightedNode + RootedMetaNode,
{
    type Tree: RootedWeightedTree + RootedMetaTree;
 
    fn get_tree(&self)->&Self::Tree;
    fn reset(&mut self);

    fn get_minPD(
        &self,
        num_taxa: usize,
    ) -> TreeNodeWeight<Self::Tree>;

    fn get_minPD_node(
        &self,
        node_id: TreeNodeID<Self::Tree>,
        num_taxa: usize,
    ) -> TreeNodeWeight<Self::Tree>;

    fn get_norm_minPD(
        &self,
        num_taxa: usize,
    ) -> TreeNodeWeight<Self::Tree>;


    fn get_minPD_taxa_set(
        &self,
        num_taxa: usize,
    ) -> impl Iterator<Item = TreeNodeID<Self::Tree>>;

    fn get_norm_minPD_taxa_set(
        &self,
        num_taxa: usize,
    ) -> impl Iterator<Item = TreeNodeID<Self::Tree>>;

    fn get_min_genPD(
        &self,
    ) -> TreeNodeWeight<Self::Tree>;

    fn get_min_genPD_set(
        &self,
    ) -> impl Iterator<Item = TreeNodeID<Self::Tree>>;

    fn get_maxPD(
        &self,
        num_taxa: usize,
    ) -> TreeNodeWeight<Self::Tree>;

    fn get_maxPD_node(
        &self,
        node_id: TreeNodeID<Self::Tree>,
        num_taxa: usize,
    ) -> TreeNodeWeight<Self::Tree>;

    fn get_norm_maxPD(
        &self,
        num_taxa: usize,
    ) -> TreeNodeWeight<Self::Tree>;

    fn get_maxPD_taxa_set(
        &self,
        num_taxa: usize,
    ) -> impl Iterator<Item = TreeNodeID<Self::Tree>>;

    fn get_norm_maxPD_taxa_set(
        &self,
        num_taxa: usize,
    ) -> impl Iterator<Item = TreeNodeID<Self::Tree>>;

    fn get_max_genPD(
        &self,
    ) -> TreeNodeWeight<Self::Tree>;

    fn get_max_genPD_set(
        &self,
    ) -> impl Iterator<Item = TreeNodeID<Self::Tree>>;

}

pub trait RootedPhylogeneticDiversity: RootedWeightedTree + Clusters
where 
    <Self as RootedTree>::Node: RootedWeightedNode
{
    fn compute_dp_table(
        &self,
        op: Ordering,
    ) -> (
        Vec<Vec<(TreeNodeWeight<Self>, u32)>>,
        Vec<Vec<Vec<TreeNodeID<Self>>>>,
        Vec<Vec<(TreeNodeWeight<Self>, u32)>>,
        Vec<Vec<Vec<TreeNodeID<Self>>>>,
    );

    fn compute_norm_min(
        &self,
    ) -> (
        Vec<Vec<(TreeNodeWeight<Self>, u32)>>,
        Vec<Vec<Vec<TreeNodeID<Self>>>>,
        Vec<Vec<(TreeNodeWeight<Self>, u32)>>,
        Vec<Vec<Vec<TreeNodeID<Self>>>>,
    ) {
        self.compute_dp_table(Ordering::Less)
    }

    fn compute_norm_max(
        &self,
    ) -> (
        Vec<Vec<(TreeNodeWeight<Self>, u32)>>,
        Vec<Vec<Vec<TreeNodeID<Self>>>>,
        Vec<Vec<(TreeNodeWeight<Self>, u32)>>,
        Vec<Vec<Vec<TreeNodeID<Self>>>>,
    ) {
        self.compute_dp_table(Ordering::Greater)
    }
}

pub trait UnrootedPhylogeneticDiversity: RootedWeightedTree + Clusters
where 
    <Self as RootedTree>::Node: RootedWeightedNode
{
    fn compute_dp_table(
        &self,
        op: Ordering,
    ) -> (
        Vec<Vec<(TreeNodeWeight<Self>, u32)>>,
        Vec<Vec<Vec<TreeNodeID<Self>>>>,
        Vec<Vec<(TreeNodeWeight<Self>, u32)>>,
        Vec<Vec<Vec<TreeNodeID<Self>>>>,
    );

    fn compute_norm_min(
        &self,
    ) -> (
        Vec<Vec<(TreeNodeWeight<Self>, u32)>>,
        Vec<Vec<Vec<TreeNodeID<Self>>>>,
        Vec<Vec<(TreeNodeWeight<Self>, u32)>>,
        Vec<Vec<Vec<TreeNodeID<Self>>>>,
    ) {
        self.compute_dp_table(Ordering::Less)
    }

    fn compute_norm_max(
        &self,
    ) -> (
        Vec<Vec<(TreeNodeWeight<Self>, u32)>>,
        Vec<Vec<Vec<TreeNodeID<Self>>>>,
        Vec<Vec<(TreeNodeWeight<Self>, u32)>>,
        Vec<Vec<Vec<TreeNodeID<Self>>>>,
    ) {
        self.compute_dp_table(Ordering::Greater)
    }
}

/// Enumerate all partitions of +ve integer k into two parts
pub fn binary_splits(k: usize, u1_max: usize, u2_max: usize)->impl Iterator<Item=(usize,usize)>{
    let mut out_vec = vec![];
    for l in max(0_i32, k as i32 - (u2_max+1) as i32) as usize..min(k, u1_max)+1{
        out_vec.push((l, k-l));
    }
    out_vec.into_iter()
}

/// Arbitrarily binarize a non-binary tree
pub fn binarize_tree(tree: &mut SimpleRootedTree){
    let mut stack = tree.postord_ids(tree.get_root_id()).collect::<VecDeque<_>>();
    loop{
        let n_id=stack.pop_front();
        match n_id{
            Some(node_id) => {
                let node_children = tree.get_node_children_ids(node_id).collect_vec();
                match node_children.len().cmp(&3){
                    Ordering::Less => {continue;},
                    _ => {
                        let node_children_ids = tree.get_node_children_ids(node_id).collect_vec();
                        let child1 = node_children_ids[0];
                        let child2 = node_children_ids[1];
                        let split_node_id = tree.next_id();
                        let split_node = Node::new(split_node_id);
                        // added split vertex to edge into first child
                        tree.split_edge((node_id, child1), split_node);
                        // deleted edge into second child
                        tree.delete_edge(node_id, child2);
                        // added edge from split vertex to second child
                        tree.set_child(split_node_id, child2);
                        // add node_id back to stack
                        stack.push_front(node_id);
                    }
                }                
            },
            None => {break}
        }
    }
}
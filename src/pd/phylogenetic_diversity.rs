use phylo::prelude::*;

use std::cmp::{Ordering, min, max};

pub trait TreePDMap<'a>
where 
    <Self::Tree as RootedTree<'a>>::Node: RootedWeightedNode + RootedMetaNode<'a>,
{
    type Tree: RootedWeightedTree<'a> + RootedMetaTree<'a>;
 
    fn get_tree(&self)->&Self::Tree;
    fn reset(&mut self);

    fn get_minPD(
        &self,
        num_taxa: usize,
    ) -> TreeNodeWeight<'a, Self::Tree>;

    fn get_minPD_node(
        &self,
        node_id: TreeNodeID<'a, SimpleRootedTree>,
        num_taxa: usize,
    ) -> TreeNodeWeight<'a, SimpleRootedTree>;

    fn get_norm_minPD(
        &self,
        num_taxa: usize,
    ) -> TreeNodeWeight<'a, Self::Tree>;


    fn get_minPD_taxa_set(
        &self,
        num_taxa: usize,
    ) -> impl Iterator<Item = TreeNodeID<'a, SimpleRootedTree>>;

    fn get_norm_minPD_taxa_set(
        &self,
        num_taxa: usize,
    ) -> impl Iterator<Item = TreeNodeID<'a, SimpleRootedTree>>;

    fn get_min_genPD(
        &self,
    ) -> TreeNodeWeight<'a, SimpleRootedTree>;

    fn get_min_genPD_set(
        &self,
    ) -> impl Iterator<Item = TreeNodeID<'a, SimpleRootedTree>>;

    fn get_maxPD(
        &self,
        num_taxa: usize,
    ) -> TreeNodeWeight<'a, Self::Tree>;

    fn get_maxPD_node(
        &self,
        node_id: TreeNodeID<'a, SimpleRootedTree>,
        num_taxa: usize,
    ) -> TreeNodeWeight<'a, SimpleRootedTree>;

    fn get_norm_maxPD(
        &self,
        num_taxa: usize,
    ) -> TreeNodeWeight<'a, Self::Tree>;

    fn get_maxPD_taxa_set(
        &self,
        num_taxa: usize,
    ) -> impl Iterator<Item = TreeNodeID<'a, SimpleRootedTree>>;

    fn get_norm_maxPD_taxa_set(
        &self,
        num_taxa: usize,
    ) -> impl Iterator<Item = TreeNodeID<'a, SimpleRootedTree>>;

    fn get_max_genPD(
        &self,
    ) -> TreeNodeWeight<'a, SimpleRootedTree>;

    fn get_max_genPD_set(
        &self,
    ) -> impl Iterator<Item = TreeNodeID<'a, SimpleRootedTree>>;

}

pub trait RootedPhylogeneticDiversity: for<'a> RootedWeightedTree<'a>+ for <'a> Clusters<'a>
where 
    for <'a> <Self as RootedTree<'a>>::Node: RootedWeightedNode
{
    fn compute_dp_table<'a>(
        &self,
        op: Ordering,
    ) -> (
        Vec<Vec<(TreeNodeWeight<'a, Self>, u32)>>,
        Vec<Vec<Vec<TreeNodeID<'a, Self>>>>,
        Vec<Vec<(TreeNodeWeight<'a, Self>, u32)>>,
        Vec<Vec<Vec<TreeNodeID<'a, Self>>>>,
    );

    fn compute_norm_min<'a>(
        &self,
    ) -> (
        Vec<Vec<(TreeNodeWeight<'a, Self>, u32)>>,
        Vec<Vec<Vec<TreeNodeID<'a, Self>>>>,
        Vec<Vec<(TreeNodeWeight<'a, Self>, u32)>>,
        Vec<Vec<Vec<TreeNodeID<'a, Self>>>>,
    ) {
        self.compute_dp_table(Ordering::Less)
    }

    fn compute_norm_max<'a>(
        &self,
    ) -> (
        Vec<Vec<(TreeNodeWeight<'a, Self>, u32)>>,
        Vec<Vec<Vec<TreeNodeID<'a, Self>>>>,
        Vec<Vec<(TreeNodeWeight<'a, Self>, u32)>>,
        Vec<Vec<Vec<TreeNodeID<'a, Self>>>>,
    ) {
        self.compute_dp_table(Ordering::Greater)
    }
}

pub trait UnrootedPhylogeneticDiversity: for<'a> RootedWeightedTree<'a>+ for <'a> Clusters<'a>
where 
    for <'a> <Self as RootedTree<'a>>::Node: RootedWeightedNode
{
    fn compute_dp_table<'a>(
        &self,
        op: Ordering,
    ) -> (
        Vec<Vec<(TreeNodeWeight<'a, Self>, u32)>>,
        Vec<Vec<Vec<TreeNodeID<'a, Self>>>>,
        Vec<Vec<(TreeNodeWeight<'a, Self>, u32)>>,
        Vec<Vec<Vec<TreeNodeID<'a, Self>>>>,
    );

    fn compute_norm_min<'a>(
        &self,
    ) -> (
        Vec<Vec<(TreeNodeWeight<'a, Self>, u32)>>,
        Vec<Vec<Vec<TreeNodeID<'a, Self>>>>,
        Vec<Vec<(TreeNodeWeight<'a, Self>, u32)>>,
        Vec<Vec<Vec<TreeNodeID<'a, Self>>>>,
    ) {
        self.compute_dp_table(Ordering::Less)
    }

    fn compute_norm_max<'a>(
        &self,
    ) -> (
        Vec<Vec<(TreeNodeWeight<'a, Self>, u32)>>,
        Vec<Vec<Vec<TreeNodeID<'a, Self>>>>,
        Vec<Vec<(TreeNodeWeight<'a, Self>, u32)>>,
        Vec<Vec<Vec<TreeNodeID<'a, Self>>>>,
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
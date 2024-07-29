use phylo::prelude::*;

use std::cmp::Ordering;

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

pub trait PhylogeneticDiversity: for<'a> RootedWeightedTree<'a>+ for <'a> Clusters<'a>
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

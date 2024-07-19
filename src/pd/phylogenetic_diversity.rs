use phylo::prelude::*;

use std::{collections::HashSet, fmt::Display, ops::IndexMut};
use num::{NumCast,Float,Signed,Zero};

#[derive(Debug, Clone)]
pub enum PDAttribute<T: Float + NumCast + Signed + Zero, U: Display> {
    PD(Vec<T>),
    EdgeCount(u32),
    Set(HashSet<U>),
}

pub enum PDNodeAttributeType {
    MinPD,
    MinNormPD,
    MinGenPD,
    MaxPD,
    MaxNormPD,
    MaxGenPD,
    EdgeCount,
    Set
}

pub trait PDNodeAttributes<T: Float + NumCast + Signed + Zero, U: Display>:
    IndexMut<PDNodeAttributeType, Output = PDAttribute<T, U>>
{
    fn reset(&mut self) {}
}

pub trait PhylogeneticDiversity<'a> 
where
    <Self::Tree as RootedTree<'a>>::Node: RootedTreeNode+RootedWeightedNode
{
    type Tree: RootedTree<'a>;

    fn precompute_minPDs(&mut self);

    fn get_tree(&self) -> &Self::Tree;

    fn compute_norm_min(
        &self,
    ) -> (
        Vec<Vec<(f32, u32)>>,
        Vec<Vec<Vec<usize>>>,
        Vec<Vec<(f32, u32)>>,
        Vec<Vec<Vec<usize>>>,
    );

    fn get_minPD(
        &self,
        num_taxa: usize,
    ) -> TreeNodeWeight<'a, Self::Tree>;

    fn get_min_genPD(
        &self,
    ) -> TreeNodeWeight<'a, Self::Tree>;

    fn get_min_genPD_set(
        &self,
    ) -> impl Iterator<Item = TreeNodeID<'a, Self::Tree>>;

    fn get_norm_minPD(
        &self,
        num_taxa: usize,
    ) -> TreeNodeWeight<'a, Self::Tree>;

    fn backtrack_min(
        &self,
        node_id: TreeNodeID<'a, Self::Tree>,
        num_taxa: usize,
        taxaset: &mut Vec<TreeNodeID<'a, Self::Tree>>,
    );

    fn get_minPD_taxa_set_node(
        &self,
        node_id: TreeNodeID<'a, Self::Tree>,
        num_taxa: usize,
    ) -> impl Iterator<Item = TreeNodeID<'a, Self::Tree>> {
        let mut taxa_set: Vec<TreeNodeID<'a, Self::Tree>> =
            vec![];
        self.backtrack_min(node_id, num_taxa, &mut taxa_set);
        taxa_set.into_iter()
    }

    fn get_minPD_taxa_set(
        &self,
        num_taxa: usize,
    ) -> impl Iterator<Item = TreeNodeID<'a, Self::Tree>> {
        let mut taxa_set: Vec<TreeNodeID<'a, Self::Tree>> =
            vec![];
        self.backtrack_min(self.get_tree().get_root_id(), num_taxa, &mut taxa_set);
        taxa_set.into_iter()
    }

    fn get_norm_minPD_taxa_set(
        &self,
        num_taxa: usize,
    ) -> impl Iterator<Item = TreeNodeID<'a, Self::Tree>>;

    fn get_minPD_node(
        &self,
        node_id: TreeNodeID<'a, Self::Tree>,
        num_taxa: usize,
    ) -> TreeNodeWeight<'a, Self::Tree>;

    fn precompute_maxPDs(&mut self);

    fn compute_norm_max(
        &self,
    ) -> (
        Vec<Vec<(f32, u32)>>,
        Vec<Vec<Vec<usize>>>,
        Vec<Vec<(f32, u32)>>,
        Vec<Vec<Vec<usize>>>,
    );

    fn get_maxPD(
        &self,
        num_taxa: usize,
    ) -> TreeNodeWeight<'a, Self::Tree>;

    fn get_max_genPD(
        &self,
    ) -> TreeNodeWeight<'a, Self::Tree>;

    fn get_max_genPD_set(
        &self,
    ) -> impl Iterator<Item = TreeNodeID<'a, Self::Tree>>;

    fn get_norm_maxPD(
        &self,
        num_taxa: usize,
    ) -> TreeNodeWeight<'a, Self::Tree>;

    fn backtrack_max(
        &self,
        node_id: TreeNodeID<'a, Self::Tree>,
        num_taxa: usize,
        taxaset: &mut Vec<TreeNodeID<'a, Self::Tree>>,
    );

    fn get_maxPD_taxa_set_node(
        &self,
        node_id: TreeNodeID<'a, Self::Tree>,
        num_taxa: usize,
    ) -> impl Iterator<Item = TreeNodeID<'a, Self::Tree>> {
        let mut taxa_set: Vec<TreeNodeID<'a, Self::Tree>> =
            vec![];
        self.backtrack_max(node_id, num_taxa, &mut taxa_set);
        taxa_set.into_iter()
    }

    fn get_maxPD_taxa_set(
        &self,
        num_taxa: usize,
    ) -> impl Iterator<Item = TreeNodeID<'a, Self::Tree>> {
        let mut taxa_set: Vec<TreeNodeID<'a, Self::Tree>> =
            vec![];
        self.backtrack_max(self.get_tree().get_root_id(), num_taxa, &mut taxa_set);
        taxa_set.into_iter()
    }

    fn get_norm_maxPD_taxa_set(
        &self,
        num_taxa: usize,
    ) -> impl Iterator<Item = TreeNodeID<'a, Self::Tree>>;

    fn get_maxPD_node(
        &self,
        node_id: TreeNodeID<'a, Self::Tree>,
        num_taxa: usize,
    ) -> TreeNodeWeight<'a, Self::Tree>;

}

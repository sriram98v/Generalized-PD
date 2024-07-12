use phylo::prelude::*;

pub trait PhylogeneticDiversity {
    type Tree: RootedTree<Node: RootedWeightedNode>;

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
    ) -> <<<Self as PhylogeneticDiversity>::Tree as RootedTree>::Node as RootedWeightedNode>::Weight;

    fn get_min_genPD(
        &self,
    ) -> <<<Self as PhylogeneticDiversity>::Tree as RootedTree>::Node as RootedWeightedNode>::Weight;

    fn get_min_genPD_set(
        &self,
    ) -> impl Iterator<Item = <<Self as PhylogeneticDiversity>::Tree as RootedTree>::NodeID>;

    fn get_norm_minPD(
        &self,
        num_taxa: usize,
    ) -> <<<Self as PhylogeneticDiversity>::Tree as RootedTree>::Node as RootedWeightedNode>::Weight;

    fn backtrack_min(
        &self,
        node_id: <<Self as PhylogeneticDiversity>::Tree as RootedTree>::NodeID,
        num_taxa: usize,
        taxaset: &mut Vec<<<Self as PhylogeneticDiversity>::Tree as RootedTree>::NodeID>,
    );

    fn get_minPD_taxa_set_node(
        &self,
        node_id: <<Self as PhylogeneticDiversity>::Tree as RootedTree>::NodeID,
        num_taxa: usize,
    ) -> impl Iterator<Item = <<Self as PhylogeneticDiversity>::Tree as RootedTree>::NodeID> {
        let mut taxa_set: Vec<<<Self as PhylogeneticDiversity>::Tree as RootedTree>::NodeID> =
            vec![];
        self.backtrack_min(node_id, num_taxa, &mut taxa_set);
        taxa_set.into_iter()
    }

    fn get_minPD_taxa_set(
        &self,
        num_taxa: usize,
    ) -> impl Iterator<Item = <<Self as PhylogeneticDiversity>::Tree as RootedTree>::NodeID> {
        let mut taxa_set: Vec<<<Self as PhylogeneticDiversity>::Tree as RootedTree>::NodeID> =
            vec![];
        self.backtrack_min(self.get_tree().get_root_id(), num_taxa, &mut taxa_set);
        taxa_set.into_iter()
    }

    fn get_norm_minPD_taxa_set(
        &self,
        num_taxa: usize,
    ) -> impl Iterator<Item = <<Self as PhylogeneticDiversity>::Tree as RootedTree>::NodeID>;

    fn get_minPD_node(
        &self,
        node_id: <<Self as PhylogeneticDiversity>::Tree as RootedTree>::NodeID,
        num_taxa: usize,
    ) -> <<<Self as PhylogeneticDiversity>::Tree as RootedTree>::Node as RootedWeightedNode>::Weight;

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
    ) -> <<<Self as PhylogeneticDiversity>::Tree as RootedTree>::Node as RootedWeightedNode>::Weight;

    fn get_max_genPD(
        &self,
    ) -> <<<Self as PhylogeneticDiversity>::Tree as RootedTree>::Node as RootedWeightedNode>::Weight;

    fn get_max_genPD_set(
        &self,
    ) -> impl Iterator<Item = <<Self as PhylogeneticDiversity>::Tree as RootedTree>::NodeID>;

    fn get_norm_maxPD(
        &self,
        num_taxa: usize,
    ) -> <<<Self as PhylogeneticDiversity>::Tree as RootedTree>::Node as RootedWeightedNode>::Weight;

    fn backtrack_max(
        &self,
        node_id: <<Self as PhylogeneticDiversity>::Tree as RootedTree>::NodeID,
        num_taxa: usize,
        taxaset: &mut Vec<<<Self as PhylogeneticDiversity>::Tree as RootedTree>::NodeID>,
    );

    fn get_maxPD_taxa_set_node(
        &self,
        node_id: <<Self as PhylogeneticDiversity>::Tree as RootedTree>::NodeID,
        num_taxa: usize,
    ) -> impl Iterator<Item = <<Self as PhylogeneticDiversity>::Tree as RootedTree>::NodeID> {
        let mut taxa_set: Vec<<<Self as PhylogeneticDiversity>::Tree as RootedTree>::NodeID> =
            vec![];
        self.backtrack_max(node_id, num_taxa, &mut taxa_set);
        taxa_set.into_iter()
    }

    fn get_maxPD_taxa_set(
        &self,
        num_taxa: usize,
    ) -> impl Iterator<Item = <<Self as PhylogeneticDiversity>::Tree as RootedTree>::NodeID> {
        let mut taxa_set: Vec<<<Self as PhylogeneticDiversity>::Tree as RootedTree>::NodeID> =
            vec![];
        self.backtrack_max(self.get_tree().get_root_id(), num_taxa, &mut taxa_set);
        taxa_set.into_iter()
    }

    fn get_norm_maxPD_taxa_set(
        &self,
        num_taxa: usize,
    ) -> impl Iterator<Item = <<Self as PhylogeneticDiversity>::Tree as RootedTree>::NodeID>;

    fn get_maxPD_node(
        &self,
        node_id: <<Self as PhylogeneticDiversity>::Tree as RootedTree>::NodeID,
        num_taxa: usize,
    ) -> <<<Self as PhylogeneticDiversity>::Tree as RootedTree>::Node as RootedWeightedNode>::Weight;
}

/// Binary-heap entry ordered only by priority.
/// The node is payload and does not affect ordering.
/// `Eq` and `Ord` compare the priority only.
pub struct PriorityEntry<P, N> {
    /// Value used for heap ordering.
    priority: P,
    /// Payload carried with the priority.
    /// It does not affect ordering.
    node: N,
}

impl<P, N> PriorityEntry<P, N> {
    /// Creates an entry ordered by `priority`.
    /// `node` is carried as payload.
    pub fn new(priority: P, node: N) -> Self {
        Self { priority, node }
    }

    /// Consumes the entry and returns the priority and node.
    pub fn into_parts(self) -> (P, N) {
        (self.priority, self.node)
    }

    // getters
    pub fn priority(&self) -> &P {
        &self.priority
    }

    pub fn node(&self) -> &N {
        &self.node
    }
}

impl<P: PartialEq, N> PartialEq for PriorityEntry<P, N> {
    fn eq(&self, other: &Self) -> bool {
        self.priority == other.priority
    }
}

impl<P: Eq, N> Eq for PriorityEntry<P, N> {}

impl<P: PartialOrd, N> PartialOrd for PriorityEntry<P, N> {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        self.priority.partial_cmp(&other.priority)
    }
}

impl<P: Ord, N> Ord for PriorityEntry<P, N> {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        self.priority.cmp(&other.priority)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::cmp::Ordering;
    use std::collections::BinaryHeap;

    /// Ordering and equality follow the priority and ignore the node.
    #[test]
    fn ordering_ignores_the_node() {
        let low = PriorityEntry::new(1, "later");
        let high = PriorityEntry::new(2, "earlier");
        let same_priority = PriorityEntry::new(2, "different");

        assert_eq!(*low.priority(), 1);
        assert_eq!(*low.node(), "later");
        assert_eq!(low.cmp(&high), Ordering::Less);
        assert!(high == same_priority);
        assert_eq!(high.cmp(&same_priority), Ordering::Equal);
    }

    /// A binary heap pops the highest priority first.
    #[test]
    fn binary_heap_pops_highest_priority_first() {
        let mut heap = BinaryHeap::new();
        heap.push(PriorityEntry::new(1, "low"));
        heap.push(PriorityEntry::new(3, "high"));
        heap.push(PriorityEntry::new(2, "mid"));

        assert_eq!(heap.pop().unwrap().into_parts(), (3, "high"));
        assert_eq!(heap.pop().unwrap().into_parts(), (2, "mid"));
        assert_eq!(heap.pop().unwrap().into_parts(), (1, "low"));
    }
}

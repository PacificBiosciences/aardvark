
use std::sync::OnceLock;

/// A lazily filled value that only affects equality when *both* values are set.
/// If either cache is empty, the two caches compare equal.
#[derive(Clone, Debug)]
pub struct Cache<T> {
    /// Value stored by the first successful `set`.
    value: OnceLock<T>,
}

impl<T> Cache<T> {
    /// Creates an empty cache.
    pub fn new() -> Self {
        Self { value: OnceLock::new() }
    }

    /// Returns the stored value when one has been set.
    pub fn get(&self) -> Option<&T> {
        self.value.get()
    }

    /// Stores `value` when the cache is still empty.
    /// A value that arrives after another caller stored one is returned unchanged.
    pub fn set(&self, value: T) -> Result<(), T> {
        self.value.set(value)
    }
}

impl<T> Default for Cache<T> {
    fn default() -> Self {
        Self::new()
    }
}

impl<T: PartialEq> PartialEq for Cache<T> {
    fn eq(&self, other: &Self) -> bool {
        match (self.get(), other.get()) {
            (Some(left), Some(right)) => left == right,
            _ => true,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// An empty cache compares equal to another empty cache and to a filled one.
    #[test]
    fn empty_cache_equals_any_cache() {
        let empty = Cache::new();
        let filled = Cache::new();
        filled.set(1).unwrap();
        let also_empty = Cache::<usize>::new();

        assert_eq!(empty, also_empty);
        assert_eq!(empty, filled);
    }

    /// Two filled caches compare their stored values.
    #[test]
    fn filled_caches_compare_stored_values() {
        let left = Cache::new();
        let same = Cache::new();
        let different = Cache::new();
        left.set(1).unwrap();
        same.set(1).unwrap();
        different.set(2).unwrap();

        assert_eq!(left, same);
        assert_ne!(left, different);
    }

    /// The first set is stored and a later set is returned.
    #[test]
    fn set_keeps_the_first_value() {
        let cache = Cache::<usize>::new();
        assert_eq!(cache.get().copied(), None);
        assert!(cache.set(7).is_ok());
        assert_eq!(cache.get().copied(), Some(7));
        assert_eq!(cache.set(9), Err(9));
        assert_eq!(cache.get().copied(), Some(7));
    }

    /// A clone copies a stored value, and cloning an empty cache stays empty.
    #[test]
    fn clone_copies_a_stored_value() {
        let filled = Cache::new();
        filled.set(4).unwrap();
        assert_eq!(filled.clone().get().copied(), Some(4));

        let empty: Cache<usize> = Cache::new();
        assert_eq!(empty.clone().get().copied(), None);
    }
}

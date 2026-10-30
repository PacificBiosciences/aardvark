
use itertools::Itertools;

/// Collapses consecutive identical bases (homopolymer compression).
/// For example, `AAACCCG` becomes `ACG`.
pub fn homopolymer_compress(sequence: &[u8]) -> Vec<u8> {
    sequence.iter().copied().dedup().collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_homopolymer_compress() {
        assert_eq!(homopolymer_compress(b""), b"");
        assert_eq!(homopolymer_compress(b"A"), b"A");
        assert_eq!(homopolymer_compress(b"ACGT"), b"ACGT");
        assert_eq!(homopolymer_compress(b"AAACCCG"), b"ACG");
        assert_eq!(homopolymer_compress(b"AAAA"), b"A");
    }
}

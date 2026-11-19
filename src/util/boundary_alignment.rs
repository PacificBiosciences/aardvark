
use anyhow::ensure;

use crate::util::sequence_alignment::{advance_block, horizontal_delta};

/// Match bitmasks for a growing sequence, one mask per distinct byte value.
/// Bit `i` of the mask for byte `b` is set when the sequence holds `b` at position `i`.
#[derive(Clone, Debug, Eq, Hash, PartialEq)]
struct SymbolMatches {
    /// Index into `masks` for each byte value, or `None` if it has not been seen
    slot_of: [Option<u8>; 256],
    /// One mask per byte value seen so far, each `len.div_ceil(64)` words long
    masks: Vec<Vec<u64>>,
    /// Number of symbols pushed so far
    len: usize,
}

impl Default for SymbolMatches {
    fn default() -> Self {
        Self {
            slot_of: [None; 256],
            masks: vec![],
            len: 0,
        }
    }
}

impl SymbolMatches {
    /// Appends one symbol to the sequence.
    /// # Arguments
    /// * `symbol` - the byte being appended
    fn push(&mut self, symbol: u8) {
        let block = self.len / 64;
        let bit = self.len % 64;

        // a new block starts, so every existing mask grows by one zero word
        if bit == 0 {
            for mask in self.masks.iter_mut() {
                mask.push(0);
            }
        }

        // a byte value seen for the first time gets a zero mask of the current size
        let slot = match self.slot_of[symbol as usize] {
            Some(slot) => slot as usize,
            None => {
                let slot = self.masks.len() as u8;
                self.slot_of[symbol as usize] = Some(slot);
                self.masks.push(vec![0; block + 1]);
                slot as usize
            }
        };
        self.masks[slot][block] |= 1u64 << bit;
        self.len += 1;
    }

    /// Returns the match mask for a byte value, or `None` if the sequence never holds it.
    /// # Arguments
    /// * `symbol` - the byte to look up
    fn get(&self, symbol: u8) -> Option<&[u64]> {
        self.slot_of[symbol as usize]
            .map(|slot| &self.masks[slot as usize][..])
    }
}

/// Exact edit-distance grid boundaries for two sequences that grow by appending.
/// Conceptually the grid `D[i][j]` is the edit distance between the first `i` symbols of the baseline and the first `j` symbols of the other sequence.
/// The baseline runs down the rows and the other sequence runs across the columns.
/// Only the last row and the last column are stored, as the +1/0/-1 steps between neighboring cells.
/// Appending to the other sequence extends the last column through the bit-parallel recurrence of Myers (1999) and Hyyrö (2001).
/// The references are in the `sequence_alignment` module comment.
/// Appending to the baseline does the same for the last row, with the roles of the sequences swapped.
/// Each appended symbol costs O(length / 64), so the whole grid is never rebuilt.
#[derive(Clone, Debug, Default, Eq, Hash, PartialEq)]
pub struct BoundaryTracker {
    /// Number of baseline symbols consumed, which is the last row of the grid
    rows: usize,
    /// Number of other-sequence symbols consumed, which is the last column of the grid
    cols: usize,

    /// Positive steps down the last column, bit `i` is `D[i+1][cols] - D[i][cols] == +1`
    column_positive: Vec<u64>,
    /// Negative steps down the last column, bit `i` is `D[i+1][cols] - D[i][cols] == -1`
    column_negative: Vec<u64>,

    /// Positive steps along the last row, bit `j` is `D[rows][j+1] - D[rows][j] == +1`
    row_positive: Vec<u64>,
    /// Negative steps along the last row, bit `j` is `D[rows][j+1] - D[rows][j] == -1`
    row_negative: Vec<u64>,

    /// Leading blocks of the last column that are known to be all -1 steps, which never change again
    column_settled: usize,
    /// Leading blocks of the last row that are known to be all -1 steps, which never change again
    row_settled: usize,
    /// First block of a run of all +1 steps that reaches the end of the last column; the block count when there is none
    column_tail: usize,
    /// First block of a run of all +1 steps that reaches the end of the last row; the block count when there is none
    row_tail: usize,

    /// Match masks over the baseline, used when a column is appended
    baseline_matches: SymbolMatches,
    /// Match masks over the other sequence, used when a row is appended
    other_matches: SymbolMatches,
}

impl BoundaryTracker {
    /// Brings the boundaries up to date with the current sequences.
    /// Both sequences must only have been appended to since the previous call.
    /// # Arguments
    /// * `baseline_seq` - the baseline sequence, down the rows
    /// * `other_seq` - the other sequence, across the columns
    /// # Errors
    /// * if either sequence is shorter than what was already consumed
    pub fn extend(&mut self, baseline_seq: &[u8], other_seq: &[u8]) -> anyhow::Result<()> {
        ensure!(baseline_seq.len() >= self.rows && other_seq.len() >= self.cols, "sequences can only be appended to");

        // new rows run across the columns that exist now
        for &symbol in baseline_seq[self.rows..].iter() {
            let (positive, negative) = advance_vector(
                self.other_matches.get(symbol), &mut self.row_positive, &mut self.row_negative, self.cols, &mut self.row_settled, &mut self.row_tail
            );
            push_step(&mut self.column_positive, &mut self.column_negative, self.rows, positive, negative, &mut self.column_tail);
            self.baseline_matches.push(symbol);
            self.rows += 1;
        }

        // new columns run down all the rows, including the new ones
        for &symbol in other_seq[self.cols..].iter() {
            let (positive, negative) = advance_vector(
                self.baseline_matches.get(symbol), &mut self.column_positive, &mut self.column_negative, self.rows, &mut self.column_settled, &mut self.column_tail
            );
            push_step(&mut self.row_positive, &mut self.row_negative, self.cols, positive, negative, &mut self.row_tail);
            self.other_matches.push(symbol);
            self.cols += 1;
        }
        Ok(())
    }

    /// Returns the smallest value on the last row or the last column of the grid.
    /// This is the edit distance at which at least one of the two sequences has been consumed completely.
    pub fn boundary_minimum(&self) -> usize {
        // the last column starts at D[0][cols] = cols, and the last row starts at D[rows][0] = rows
        let (column_minimum, column_end) = scan_steps(self.cols, &self.column_positive, &self.column_negative, self.rows);
        let (row_minimum, row_end) = scan_steps(self.rows, &self.row_positive, &self.row_negative, self.cols);
        debug_assert_eq!(column_end, row_end);
        column_minimum.min(row_minimum)
    }

    /// Returns the edit distance between the full baseline and the full other sequence.
    pub fn corner(&self) -> usize {
        scan_steps(self.cols, &self.column_positive, &self.column_negative, self.rows).1
    }
}

/// Appends one symbol to the sequence that `positive` and `negative` run along.
/// The vectors hold the steps along the last line of the grid, and are updated to the next line over.
/// Returns the step the new line takes at the far end of the vector as `(positive, negative)`.
/// # Arguments
/// * `matches` - match mask of the other sequence for the appended symbol, `None` if it never occurs
/// * `positive` - positive steps along the line, updated in place
/// * `negative` - negative steps along the line, updated in place
/// * `len` - number of steps along the line
/// * `settled` - leading blocks already known to be all -1 steps; raised to cover any newly settled ones
/// * `tail` - first block of a known run of all +1 steps that reaches the end; updated for the next call
fn advance_vector(matches: Option<&[u64]>, positive: &mut [u64], negative: &mut [u64], len: usize, settled: &mut usize, tail: &mut usize) -> (bool, bool) {
    // the first cell of every new line is one more than the cell before it, which is the +1 carry entering block 0
    if len == 0 {
        return (true, false);
    }

    let blocks = len.div_ceil(64);
    let last_block = blocks - 1;
    let last_bit_shift = (len - 1) % 64;
    let (mut carry_positive, mut carry_negative) = (1u64, 0u64);

    // a block whose steps are all -1 stays that way under a +1 carry, whatever the symbol, and passes the +1 on
    let mut known = (*settled).min(last_block);
    while known < last_block && positive[known] == 0 && negative[known] == u64::MAX {
        known += 1;
    }
    *settled = known;

    // more generally, a block stays the same under a +1 carry and passes it on when it has no +1 steps and every row matching the symbol steps by -1
    while known < last_block && positive[known] == 0 && matches.map_or(0, |mask| mask[known]) & !negative[known] == 0 {
        known += 1;
    }
    let settled = known;

    // likewise a block whose steps are all +1 stays that way under a -1 carry and passes the -1 on, so a run of them at the end can be skipped
    let real_bits = if len.is_multiple_of(64) { u64::MAX } else { (1u64 << (len % 64)) - 1 };
    let mut fixed_tail = (*tail).min(blocks);
    while fixed_tail > settled {
        let block_bits = if fixed_tail == blocks { real_bits } else { u64::MAX };
        let block = fixed_tail - 1;
        if positive[block] & block_bits != block_bits || negative[block] & block_bits != 0 {
            break;
        }
        fixed_tail -= 1;
    }

    // advance up to the tail, and stop there if the carry reaching it is -1, since the far end then steps by -1
    let middle_end = fixed_tail.min(last_block);
    (carry_positive, carry_negative) = advance_blocks(matches, positive, negative, settled, middle_end, (carry_positive, carry_negative));
    if fixed_tail < blocks && (carry_positive, carry_negative) == (0, 1) {
        *tail = fixed_tail;
        return (false, true);
    }

    // the tail, if any, is about to change, so no run is known after this call
    *tail = blocks;
    (carry_positive, carry_negative) = advance_blocks(matches, positive, negative, middle_end, last_block, (carry_positive, carry_negative));

    // the step at the far end is read from the last block before it is advanced
    let equality = matches.map_or(0, |mask| mask[last_block]);
    let (horizontal_positive, horizontal_negative) = horizontal_delta(equality | carry_negative, positive[last_block], negative[last_block]);
    let _ = advance_block(equality, &mut positive[last_block], &mut negative[last_block], carry_positive, carry_negative);
    ((horizontal_positive >> last_bit_shift) & 1 == 1, (horizontal_negative >> last_bit_shift) & 1 == 1)
}

/// Advances blocks `start..end` of a line by one symbol and returns the carry out of the last of them.
/// # Arguments
/// * `matches` - match mask of the other sequence for the appended symbol, `None` if it never occurs
/// * `positive` - positive steps along the line, updated in place
/// * `negative` - negative steps along the line, updated in place
/// * `start` - first block to advance
/// * `end` - one past the last block to advance
/// * `carry` - carry entering block `start` as `(positive, negative)`
fn advance_blocks(matches: Option<&[u64]>, positive: &mut [u64], negative: &mut [u64], start: usize, end: usize, carry: (u64, u64)) -> (u64, u64) {
    let (mut carry_positive, mut carry_negative) = carry;
    let (body_positive, body_negative) = (&mut positive[start..end], &mut negative[start..end]);
    match matches {
        Some(mask) => {
            for ((&equality, block_positive), block_negative) in mask[start..end].iter().zip(body_positive.iter_mut()).zip(body_negative.iter_mut()) {
                (carry_positive, carry_negative) = advance_block(equality, block_positive, block_negative, carry_positive, carry_negative);
            }
        },
        None => {
            for (block_positive, block_negative) in body_positive.iter_mut().zip(body_negative.iter_mut()) {
                (carry_positive, carry_negative) = advance_block(0, block_positive, block_negative, carry_positive, carry_negative);
            }
        }
    }
    (carry_positive, carry_negative)
}

/// Sets step `index` in a pair of step vectors, growing them by a word when needed.
/// Bits above the last real step are padding that the recurrence never reads downward, so they are overwritten here.
/// # Arguments
/// * `positive` - positive steps
/// * `negative` - negative steps
/// * `index` - the step to set
/// * `is_positive` - true if the step is +1
/// * `is_negative` - true if the step is -1
/// * `tail` - first block of the run of all +1 steps that reaches the end, updated for the new step
fn push_step(positive: &mut Vec<u64>, negative: &mut Vec<u64>, index: usize, is_positive: bool, is_negative: bool, tail: &mut usize) {
    let block = index / 64;
    let bit = index % 64;
    if block == positive.len() {
        positive.push(0);
        negative.push(0);
    }
    positive[block] = (positive[block] & !(1u64 << bit)) | ((is_positive as u64) << bit);
    negative[block] = (negative[block] & !(1u64 << bit)) | ((is_negative as u64) << bit);

    // a +1 step extends a run of +1 steps at the end, and anything else ends it
    if !is_positive || is_negative {
        *tail = block + 1;
    }
}

/// Walks a line of steps from its known starting value.
/// Returns the smallest value seen, including the start, and the value at the end.
/// # Arguments
/// * `start` - value of the first cell on the line
/// * `positive` - positive steps along the line
/// * `negative` - negative steps along the line
/// * `len` - number of steps along the line
fn scan_steps(start: usize, positive: &[u64], negative: &[u64], len: usize) -> (usize, usize) {
    let mut current = start as i64;
    let mut best = current;
    for block in 0..len.div_ceil(64) {
        // only the real steps of the last block count
        let bits = (len - block * 64).min(64);
        let mask = if bits == 64 { u64::MAX } else { (1u64 << bits) - 1 };
        let block_positive = positive[block] & mask;
        let block_negative = negative[block] & mask;

        // a block that cannot fall below the best so far only needs its net change
        if current - bits as i64 >= best {
            current += block_positive.count_ones() as i64 - block_negative.count_ones() as i64;
            continue;
        }

        // a block that never steps up is lowest at its end, and one that never steps down is lowest at its start
        if block_positive == 0 || block_negative == 0 {
            current += block_positive.count_ones() as i64 - block_negative.count_ones() as i64;
            best = best.min(current);
            continue;
        }
        for bit in 0..bits {
            current += ((block_positive >> bit) & 1) as i64 - ((block_negative >> bit) & 1) as i64;
            best = best.min(current);
        }
    }
    (best as usize, current as usize)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Builds the full grid `D[i][j]` the slow way.
    /// # Arguments
    /// * `baseline_seq` - the baseline sequence, down the rows
    /// * `other_seq` - the other sequence, across the columns
    fn full_grid(baseline_seq: &[u8], other_seq: &[u8]) -> Vec<Vec<usize>> {
        let mut grid = vec![vec![0usize; other_seq.len() + 1]; baseline_seq.len() + 1];
        for (i, row) in grid.iter_mut().enumerate() {
            row[0] = i;
        }
        for (j, cell) in grid[0].iter_mut().enumerate() {
            *cell = j;
        }
        for i in 1..=baseline_seq.len() {
            for j in 1..=other_seq.len() {
                let diagonal = grid[i - 1][j - 1] + usize::from(baseline_seq[i - 1] != other_seq[j - 1]);
                grid[i][j] = diagonal.min(grid[i - 1][j] + 1).min(grid[i][j - 1] + 1);
            }
        }
        grid
    }

    /// Returns the smallest value on the last row or last column of a full grid, and its bottom-right corner.
    /// # Arguments
    /// * `grid` - the full grid from `full_grid`
    fn grid_boundary(grid: &[Vec<usize>]) -> (usize, usize) {
        let last_row = grid.len() - 1;
        let last_col = grid[0].len() - 1;
        let row_minimum = grid[last_row].iter().copied().min().unwrap();
        let column_minimum = grid.iter().map(|row| row[last_col]).min().unwrap();
        (row_minimum.min(column_minimum), grid[last_row][last_col])
    }

    /// The tracker matches the full grid when both sequences grow in random chunks, including empty and block-border sizes.
    #[test]
    fn test_boundary_tracker_matches_full_grid() {
        let mut state: u64 = 0xD1B54A32D192ED03;
        let mut next = || { state ^= state << 13; state ^= state >> 7; state ^= state << 17; state };
        let alphabet = b"ACGTN";
        for case in 0..300 {
            let alphabet_size = if case % 3 == 0 { 2 } else { 4 + (case % 2) as u64 };
            let mut baseline: Vec<u8> = vec![];
            let mut other: Vec<u8> = vec![];
            let mut tracker = BoundaryTracker::default();

            // a fresh tracker already describes the empty grid
            assert_eq!(tracker.boundary_minimum(), 0);
            assert_eq!(tracker.corner(), 0);

            for _step in 0..8 {
                // chunk sizes straddle the 64-symbol block size, and either side may not grow
                let chunk = |value: u64| match value % 6 {
                    0 => 0,
                    1 => 1,
                    2 => 63,
                    3 => 64,
                    4 => 65,
                    _ => (value / 6 % 150) as usize,
                };
                for _ in 0..chunk(next()) {
                    baseline.push(alphabet[(next() % alphabet_size) as usize]);
                }
                for _ in 0..chunk(next()) {
                    other.push(alphabet[(next() % alphabet_size) as usize]);
                }
                tracker.extend(&baseline, &other).unwrap();
                let (expected_minimum, expected_corner) = grid_boundary(&full_grid(&baseline, &other));
                assert_eq!(tracker.boundary_minimum(), expected_minimum, "case {case} lengths {} {}", baseline.len(), other.len());
                assert_eq!(tracker.corner(), expected_corner, "case {case} lengths {} {}", baseline.len(), other.len());
            }
        }
    }

    /// The tracker matches the full grid when one sequence runs far ahead of the other, so whole blocks settle at both ends of each line.
    #[test]
    fn test_boundary_tracker_lopsided_growth() {
        let mut state: u64 = 0x05EE_D0FB_10C5;
        let mut next = || { state ^= state << 13; state ^= state >> 7; state ^= state << 17; state };
        for case in 0..12 {
            let alphabet: &[u8] = if case % 2 == 0 { b"AC" } else { b"ACGT" };
            let mut short: Vec<u8> = (0..(next() % 40) as usize).map(|_| alphabet[(next() % alphabet.len() as u64) as usize]).collect();
            let mut long: Vec<u8> = (0..1_200 + (next() % 200) as usize).map(|_| alphabet[(next() % alphabet.len() as u64) as usize]).collect();
            let mut tracker = BoundaryTracker::default();
            for _step in 0..6 {
                // either side may be the baseline, and the short side keeps growing a little at a time
                let (baseline, other) = if case % 4 < 2 { (&short, &long) } else { (&long, &short) };
                tracker.extend(baseline, other).unwrap();
                let (expected_minimum, expected_corner) = grid_boundary(&full_grid(baseline, other));
                assert_eq!(tracker.boundary_minimum(), expected_minimum, "case {case} lengths {} {}", baseline.len(), other.len());
                assert_eq!(tracker.corner(), expected_corner, "case {case} lengths {} {}", baseline.len(), other.len());
                short.extend((0..(next() % 130) as usize).map(|_| alphabet[(next() % alphabet.len() as u64) as usize]));
                long.extend((0..(next() % 70) as usize).map(|_| alphabet[(next() % alphabet.len() as u64) as usize]));
            }
        }
    }

    /// The tracker matches the full grid through many small, mostly one-sided appends, including from cloned trackers.
    #[test]
    fn test_boundary_tracker_small_lopsided_appends() {
        let mut state: u64 = 0xC0FF_EE15_600D;
        let mut next = || { state ^= state << 13; state ^= state >> 7; state ^= state << 17; state };
        for case in 0..16 {
            let alphabet: &[u8] = if case % 2 == 0 { b"AC" } else { b"ACGT" };
            let mut baseline: Vec<u8> = vec![];
            let mut other: Vec<u8> = vec![];
            let mut tracker = BoundaryTracker::default();
            for step in 0..60 {
                // one side races ahead for a while, then the other catches up
                let (fast, slow) = if (step / 20 + case) % 2 == 0 { (&mut baseline, &mut other) } else { (&mut other, &mut baseline) };
                fast.extend((0..(next() % 25) as usize).map(|_| alphabet[(next() % alphabet.len() as u64) as usize]));
                slow.extend((0..(next() % 3) as usize).map(|_| alphabet[(next() % alphabet.len() as u64) as usize]));

                // a clone carries the cached block runs along with the steps
                let mut grown = tracker.clone();
                grown.extend(&baseline, &other).unwrap();
                tracker = grown;
                if step % 5 == 4 {
                    let (expected_minimum, expected_corner) = grid_boundary(&full_grid(&baseline, &other));
                    assert_eq!(tracker.boundary_minimum(), expected_minimum, "case {case} step {step} lengths {} {}", baseline.len(), other.len());
                    assert_eq!(tracker.corner(), expected_corner, "case {case} step {step} lengths {} {}", baseline.len(), other.len());
                }
            }
        }
    }

    /// The tracker matches the full grid when one sequence is thousands of bases ahead, so the skipped head and tail cover many blocks.
    #[test]
    fn test_boundary_tracker_long_lopsided_growth() {
        let mut state: u64 = 0x57E5_5109_B10C;
        let mut next = || { state ^= state << 13; state ^= state >> 7; state ^= state << 17; state };
        let alphabet = b"ACGT";

        // exact block multiples and lengths just off them; the short side stays small so the full grid fits
        for (case, (long_len, short_len)) in [(4_096usize, 128usize), (4_159, 192), (8_192, 256), (8_257, 320)].into_iter().enumerate() {
            let long: Vec<u8> = (0..long_len).map(|_| alphabet[(next() % 4) as usize]).collect();
            let short: Vec<u8> = (0..short_len).map(|_| alphabet[(next() % 4) as usize]).collect();
            for swapped in [false, true] {
                let (baseline, other) = if swapped { (&short, &long) } else { (&long, &short) };
                let mut tracker = BoundaryTracker::default();

                // check the halfway append, then the rest, carrying the cached runs across a clone
                tracker.extend(&baseline[..baseline.len() / 2], &other[..other.len() / 2]).unwrap();
                let (expected_minimum, expected_corner) = grid_boundary(&full_grid(&baseline[..baseline.len() / 2], &other[..other.len() / 2]));
                assert_eq!(tracker.boundary_minimum(), expected_minimum, "case {case} swapped {swapped} halfway");
                assert_eq!(tracker.corner(), expected_corner, "case {case} swapped {swapped} halfway");
                let mut grown = tracker.clone();
                grown.extend(baseline, other).unwrap();
                tracker = grown;

                let (expected_minimum, expected_corner) = grid_boundary(&full_grid(baseline, other));
                assert_eq!(tracker.boundary_minimum(), expected_minimum, "case {case} swapped {swapped} lengths {} {}", baseline.len(), other.len());
                assert_eq!(tracker.corner(), expected_corner, "case {case} swapped {swapped} lengths {} {}", baseline.len(), other.len());

                // the long line ends in a run of +1 steps once the short sequence has been absorbed, and that run is what the tail skip leaves untouched
                let long_blocks = long_len.div_ceil(64);
                let tail = if swapped { tracker.row_tail } else { tracker.column_tail };
                assert!(tail + 32 <= long_blocks, "case {case} swapped {swapped} tail {tail} blocks {long_blocks}");
            }
        }

        // a deletion like region 189343, scaled down: the truth drops a long middle stretch of a 4 kbp query
        let reference: Vec<u8> = (0..4_096).map(|_| alphabet[(next() % 4) as usize]).collect();
        let query: Vec<u8> = reference.iter().enumerate().map(|(index, &base)| {
            if index % 256 == 0 { alphabet[(base as usize + 1) % 4] } else { base }
        }).collect();
        let mut tracker = BoundaryTracker::default();
        for position in [1_024usize, 2_048, 3_072, 4_096] {
            let truth_len = if position <= 256 { position } else { 256 + position.saturating_sub(3_072) };
            let truth = reference[..256].iter().chain(reference[3_072..].iter()).take(truth_len).copied().collect::<Vec<u8>>();
            tracker.extend(&truth, &query[..position]).unwrap();
            let (expected_minimum, expected_corner) = grid_boundary(&full_grid(&truth, &query[..position]));
            assert_eq!(tracker.boundary_minimum(), expected_minimum, "deletion at {position}");
            assert_eq!(tracker.corner(), expected_corner, "deletion at {position}");
        }
    }

    /// Times tracker appends shaped like region 189343. Ignored; run it by name in release.
    #[test]
    #[ignore]
    fn bench_tracker_extend() {
        let mut state: u64 = 0x1234_5678_9ABC_DEF0;
        let mut next = || { state ^= state << 13; state ^= state >> 7; state ^= state << 17; state };
        let alphabet = b"ACGT";
        let baseline: Vec<u8> = (0..50_000).map(|_| alphabet[(next() % 4) as usize]).collect();
        let other: Vec<u8> = (0..138_000).map(|_| alphabet[(next() % 4) as usize]).collect();

        // one long column append, then the same in small pieces, then rows appended across the long line
        let mut start = BoundaryTracker::default();
        start.extend(&baseline, &other[..50_000]).unwrap();
        for chunk in [88_000usize, 100] {
            let mut best = f64::MAX;
            for _ in 0..3 {
                let mut tracker = start.clone();
                let started = std::time::Instant::now();
                let mut end = 50_000;
                while end < other.len() {
                    end = (end + chunk).min(other.len());
                    tracker.extend(&baseline, &other[..end]).unwrap();
                }
                std::hint::black_box(tracker.boundary_minimum());
                best = best.min(started.elapsed().as_secs_f64() * 1_000.0);
            }
            println!("columns in chunks of {chunk:>6}: {best:8.1} ms");
        }
        let mut wide = BoundaryTracker::default();
        wide.extend(&baseline[..1_000], &other).unwrap();
        let mut best = f64::MAX;
        for _ in 0..3 {
            let mut tracker = wide.clone();
            let started = std::time::Instant::now();
            tracker.extend(&baseline[..20_000], &other).unwrap();
            std::hint::black_box(tracker.boundary_minimum());
            best = best.min(started.elapsed().as_secs_f64() * 1_000.0);
        }
        println!("19,000 rows across 138,000 columns: {best:8.1} ms");
        let mut best = f64::MAX;
        for _ in 0..200 {
            let started = std::time::Instant::now();
            std::hint::black_box(wide.boundary_minimum());
            std::hint::black_box(start.clone());
            best = best.min(started.elapsed().as_secs_f64() * 1_000_000.0);
        }
        println!("boundary_minimum + clone: {best:8.1} us");
    }

    /// Times one search branch shaped like region 189343: the truth deletes 88 kbp that the query keeps. Ignored; run it by name in release.
    #[test]
    #[ignore]
    fn bench_tracker_deletion() {
        let mut state: u64 = 0x0F1E_2D3C_4B5A_6978;
        let mut next = || { state ^= state << 13; state ^= state >> 7; state ^= state << 17; state };
        let alphabet = b"ACGT";

        // a 138 kbp reference; the truth deletes [1_100, 89_100) and the query carries a SNP every ~1 kbp
        let reference: Vec<u8> = (0..138_000).map(|_| alphabet[(next() % 4) as usize]).collect();
        let truth: Vec<u8> = reference[..1_100].iter().chain(reference[89_100..].iter()).copied().collect();
        let query: Vec<u8> = reference.iter().map(|&base| if next() % 1_000 == 0 { alphabet[((base as u64 + 1) % 4) as usize] } else { base }).collect();

        let mut best = f64::MAX;
        let mut checksum = 0;
        for _ in 0..3 {
            let mut tracker = BoundaryTracker::default();
            let started = std::time::Instant::now();
            for position in (100..=reference.len()).step_by(100) {
                let truth_len = if position <= 1_100 { position } else { 1_100 + position.saturating_sub(89_100) };
                tracker.extend(&truth[..truth_len], &query[..position]).unwrap();
                checksum += tracker.boundary_minimum();
            }
            best = best.min(started.elapsed().as_secs_f64() * 1_000.0);
        }
        println!("deletion branch: {best:8.1} ms (checksum {checksum})");
    }

    /// All 256 byte values can be tracked, not only the DNA alphabet.
    #[test]
    fn test_boundary_tracker_all_byte_values() {
        let baseline: Vec<u8> = (0..=255).collect();
        let mut other = baseline.clone();
        other.reverse();
        let mut tracker = BoundaryTracker::default();
        tracker.extend(&baseline, &other).unwrap();
        let (expected_minimum, expected_corner) = grid_boundary(&full_grid(&baseline, &other));
        assert_eq!(tracker.boundary_minimum(), expected_minimum);
        assert_eq!(tracker.corner(), expected_corner);
    }
}

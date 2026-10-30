//! Byte-level GFA scanning. Step lists are parsed straight off the bytes: a
//! chromosome GFA reaches 26 GB, and the serial read loop must not
//! UTF-8-validate it.

use flate2::read::MultiGzDecoder;
use std::fs::File;
use std::io::{self, BufRead, BufReader, Read, Seek, SeekFrom};

/// Sequential scan of a multi-GB file; BufReader's 8 KiB default would cost
/// millions of read(2) calls per pass.
const BUF: usize = 4 << 20;

/// Open plain or gzip/bgzip GFA, detected from the magic bytes.
/// MultiGzDecoder, not GzDecoder: bgzip writes many gzip members and
/// GzDecoder stops after the first.
pub fn open(path: &str) -> io::Result<Box<dyn BufRead>> {
    let mut f = File::open(path)?;
    let mut magic = [0u8; 2];
    let gz = f.read(&mut magic)? == 2 && magic == [0x1f, 0x8b];
    f.seek(SeekFrom::Start(0))?;
    Ok(if gz {
        Box::new(BufReader::with_capacity(BUF, MultiGzDecoder::new(f)))
    } else {
        Box::new(BufReader::with_capacity(BUF, f))
    })
}

pub fn trim_eol(b: &[u8]) -> &[u8] {
    let mut e = b.len();
    while e > 0 && (b[e - 1] == b'\n' || b[e - 1] == b'\r') {
        e -= 1;
    }
    &b[..e]
}

/// ASCII decimal; None on empty input, a non-digit, or overflow.
pub fn parse_uint(b: &[u8]) -> Option<usize> {
    if b.is_empty() {
        return None;
    }
    b.iter().try_fold(0usize, |v, &c| {
        c.is_ascii_digit().then_some(())?;
        v.checked_mul(10)?.checked_add((c - b'0') as usize)
    })
}

/// Node lengths indexed by `id - min_id`; ids absent from the GFA have length 0.
pub struct Segments {
    pub len: Vec<usize>,
    pub min_id: usize,
}

impl Segments {
    pub fn new(segs: &[(usize, usize)]) -> Self {
        let min_id = segs.iter().map(|s| s.0).min().unwrap_or(0);
        let span = segs.iter().map(|s| s.0 - min_id + 1).max().unwrap_or(0);
        let mut len = vec![0; span];
        for &(id, l) in segs {
            len[id - min_id] = l;
        }
        Segments { len, min_id }
    }
}

/// One P-line's steps: node ids, the bp prefix sum (`cum[i]` = bp before step
/// i, so `cum.len() == ids.len() + 1`), and `plus[i]` for a `+` step (filled
/// only when GAF output needs it).
pub struct Path {
    pub ids: Vec<usize>,
    pub cum: Vec<usize>,
    pub plus: Vec<bool>,
}

impl Path {
    /// Hot loop of the whole program: one integer parse and one random probe
    /// into `seg.len` per step (5.97e9 step-visits on HPRCv2.1 MC chr1), so
    /// it runs once per path, never once per bucket.
    pub fn parse(steps: &[u8], seg: &Segments, orient: bool) -> Self {
        let n = steps.iter().filter(|&&c| c == b',').count() + 1;
        let mut p = Path { ids: Vec::with_capacity(n), cum: Vec::with_capacity(n + 1), plus: Vec::new() };
        let mut acc = 0;
        p.cum.push(acc);
        for s in steps.split(|&c| c == b',') {
            let (id, o) = s.split_at(s.len() - 1);
            let id = parse_uint(id).expect("non-numeric GFA path step");
            acc += seg.len[id - seg.min_id];
            p.ids.push(id);
            p.cum.push(acc);
            if orient {
                p.plus.push(o == b"+");
            }
        }
        p
    }

    /// Turn this into the reverse traversal in place, with no re-parse: node
    /// order flips and `cum_rev[i] = total - cum_fwd[n - i]`. Exact because a
    /// node's length does not depend on traversal direction.
    pub fn reverse(&mut self) {
        let total = self.cum[self.cum.len() - 1];
        self.ids.reverse();
        self.plus.reverse();
        self.cum.reverse();
        for c in &mut self.cum {
            *c = total - *c;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn steps(ids: &[(usize, bool)]) -> Vec<u8> {
        let v: Vec<String> = ids.iter().map(|&(i, p)| format!("{i}{}", if p { '+' } else { '-' })).collect();
        v.join(",").into_bytes()
    }

    #[test]
    fn reverse_matches_reparse_of_reversed_steps() {
        let mut state: u64 = 0x2026_08_31;
        let mut next = move || {
            state = state.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
            (state >> 33) as usize
        };
        for _ in 0..500 {
            // Sparse, cactus-style ids far from zero.
            let min_id = if next() % 2 == 0 { 0 } else { 162_720_608 };
            let n_nodes = 1 + next() % 12;
            let segs: Vec<(usize, usize)> = (0..n_nodes).map(|i| (min_id + i, 1 + next() % 1024)).collect();
            let seg = Segments::new(&segs);
            let walk: Vec<(usize, bool)> = (0..1 + next() % 40).map(|_| (min_id + next() % n_nodes, next() % 2 == 0)).collect();
            let mut p = Path::parse(&steps(&walk), &seg, true);
            let fwd_total = p.cum[p.cum.len() - 1];
            p.reverse();
            let rev: Vec<_> = walk.iter().rev().copied().collect();
            let r = Path::parse(&steps(&rev), &seg, true);
            assert_eq!((&p.ids, &p.cum, &p.plus), (&r.ids, &r.cum, &r.plus));
            assert_eq!((p.cum[0], p.cum[p.cum.len() - 1]), (0, fwd_total));
        }
    }

    #[test]
    fn parse_uint_rejects_junk() {
        assert_eq!(parse_uint(b"162720608"), Some(162720608));
        assert_eq!(parse_uint(b""), None);
        assert_eq!(parse_uint(b"12a"), None);
        assert_eq!(parse_uint(b"99999999999999999999999"), None);
    }
}

//! One pass over the GFA: S-lines build the node-length table, then P-lines
//! are walked in parallel batches. Workers turn each (path, orientation)
//! bucket into coverage contributions (and GAF lines); a single-threaded
//! merge applies them in GFA order. That makes the output independent of
//! scheduling, and makes dedup winners match a sequential walk -- which
//! matters, because ~1.05% of dedup keys on PGGB graphs have more than one
//! distinct downstream walk.

use crate::gfa::{self, Path, Segments};
use crate::input::{invalid, Record};
use rayon::prelude::*;
use std::collections::{HashMap, HashSet};
use std::io::{self, Write};
use std::time::Instant;

pub struct Inputs<'a> {
    pub records: &'a [Record],
    pub starts: &'a [usize],
    pub names: &'a HashMap<String, usize>,
    pub dedup: bool,
    pub gaf: bool,
}

/// Per-record worker result; contributions live in the bucket's arena.
#[derive(Clone, Copy)]
struct Rec {
    key: u64,
    contrib_end: u32,
    /// The match ran off the end of the path: its coverage is still applied,
    /// but it gets no GAF row.
    overrun: bool,
}

struct Bucket {
    name: String,
    num_records: usize,
    recs: Vec<Rec>,
    /// (node_id - min_id, bp covered)
    contribs: Vec<(u32, u32)>,
    gaf: Vec<u8>,
    gaf_end: Vec<usize>,
    out_of_range: u64,
}

/// Move the step cursor forward until `cum[*i] <= bp < cum[*i + 1]`. Records
/// are sorted by path_bp within a bucket, so the cursor only ever advances:
/// O(steps + records) per bucket.
#[inline]
fn advance(bp: usize, cum: &[usize], i: &mut usize) {
    let last = cum.len() - 1;
    while *i + 1 < last && bp >= cum[*i + 1] {
        *i += 1;
    }
}

/// Dedup key, low to high bits: read_id 20 | read_st 9 | node_id - min_id 25
/// | offset 10. Strand is deliberately absent. Field widths are checked once
/// at startup (see `check_key_fields`).
#[inline]
fn key(r: &Record, rel_node: usize, offset: usize) -> u64 {
    r.read_id as u64 | (r.read_st as u64) << 20 | (rel_node as u64) << 29 | (offset as u64) << 54
}

fn check_key_fields(records: &[Record], seg: &Segments) {
    let max_id = records.iter().map(|r| r.read_id).max().unwrap_or(0);
    let max_st = records.iter().map(|r| r.read_st).max().unwrap_or(0);
    let max_len = seg.len.iter().copied().max().unwrap_or(0);
    assert!(max_id < 1 << 20, "read_id {max_id} exceeds the 20-bit dedup field");
    assert!(max_st < 1 << 9, "read_st {max_st} exceeds the 9-bit dedup field");
    assert!(seg.len.len() <= 1 << 25, "node-id range {} exceeds the 25-bit dedup field", seg.len.len());
    assert!(max_len <= 1024, "node length {max_len} exceeds the 10-bit dedup offset field");
}

/// Spread one match over the nodes from its start step. Returns (bp of the
/// nodes touched, bp left over when the path ran out).
#[inline]
fn traverse(ids: &[usize], seg: &Segments, offset: usize, match_len: usize, out: &mut Vec<(u32, u32)>) -> (usize, usize) {
    let (mut path_len, mut rem, mut skip) = (0, match_len, offset);
    if match_len == 0 {
        return (0, 0);
    }
    for &id in ids {
        let len = seg.len[id - seg.min_id];
        let c = rem.min(len - skip);
        path_len += len;
        rem -= c;
        skip = 0;
        out.push(((id - seg.min_id) as u32, c as u32));
        if rem == 0 {
            break;
        }
    }
    (path_len, rem)
}

fn walk_bucket(p: &Path, reverse: bool, name: String, records: &[Record], seg: &Segments, inp: &Inputs) -> Bucket {
    let mut b = Bucket {
        name,
        num_records: records.len(),
        recs: Vec::with_capacity(records.len()),
        contribs: Vec::new(),
        gaf: Vec::new(),
        gaf_end: Vec::new(),
        out_of_range: 0,
    };
    let total = p.cum[p.cum.len() - 1];
    let mut i = 0;
    for r in records {
        let bp = r.path_bp as usize;
        if bp >= total {
            b.out_of_range += 1;
            continue;
        }
        advance(bp, &p.cum, &mut i);
        let offset = bp - p.cum[i];
        let start = b.contribs.len();
        let (path_len, rem) = traverse(&p.ids[i..], seg, offset, r.match_len as usize, &mut b.contribs);
        b.recs.push(Rec {
            key: if inp.dedup { key(r, p.ids[i] - seg.min_id, offset) } else { 0 },
            contrib_end: b.contribs.len() as u32,
            overrun: rem > 0,
        });
        if inp.gaf {
            if rem == 0 {
                let g = &mut b.gaf;
                write!(g, "{}\t{}\t", r.read_id, r.read_st).unwrap();
                for k in i..i + (b.contribs.len() - start) {
                    g.push(if p.plus[k] != reverse { b'>' } else { b'<' });
                    write!(g, "{}", p.ids[k]).unwrap();
                }
                let m = r.match_len as usize;
                writeln!(g, "\t{path_len}\t{offset}\t{}\t{m}\t{}", offset + m, b.name).unwrap();
            }
            b.gaf_end.push(b.gaf.len());
        }
    }
    b
}

/// A P-line's forward and reverse buckets. The path is parsed once; the
/// reverse arrays are derived from it in place.
fn walk_pline(line: &[u8], seg: &Segments, inp: &Inputs) -> Vec<Bucket> {
    let mut f = gfa::trim_eol(line).split(|&c| c == b'\t').skip(1);
    let (Some(name), Some(steps)) = (f.next(), f.next()) else {
        eprintln!("Unable to parse GFA path line");
        return Vec::new();
    };
    let Ok(name) = std::str::from_utf8(name) else {
        eprintln!("Skipping path with a non-UTF-8 name");
        return Vec::new();
    };
    let Some(&idx) = inp.names.get(name) else {
        eprintln!("Skipping path: {name} not found");
        return Vec::new();
    };
    let s = &inp.starts[2 * idx..2 * idx + 3];
    if s[0] == s[2] {
        return Vec::new(); // no records in either orientation: skip the parse
    }
    let mut p = Path::parse(steps, seg, inp.gaf);
    let fwd = walk_bucket(&p, false, name.to_string(), &inp.records[s[0]..s[1]], seg, inp);
    p.reverse();
    let rev = walk_bucket(&p, true, format!("{name}_reverse"), &inp.records[s[1]..s[2]], seg, inp);
    vec![fwd, rev]
}

struct Merge<'w> {
    coverage: Vec<f64>,
    seen: Option<HashSet<u64>>,
    skipped: u64,
    entries: u64,
    gaf: Option<&'w mut dyn Write>,
}

impl Merge<'_> {
    /// Visiting buckets in GFA order and records in index order reproduces a
    /// sequential walk's dedup winners. The key is consumed before the
    /// overrun check, and overrun records keep their coverage.
    fn apply(&mut self, b: &Bucket) -> io::Result<()> {
        let (mut c0, mut overruns) = (0, 0);
        for (j, r) in b.recs.iter().enumerate() {
            let c = c0..r.contrib_end as usize;
            c0 = c.end;
            if self.seen.as_mut().is_some_and(|s| !s.insert(r.key)) {
                self.skipped += 1;
                continue;
            }
            for &(n, bp) in &b.contribs[c] {
                self.coverage[n as usize] += bp as f64;
            }
            if r.overrun {
                overruns += 1;
                continue;
            }
            self.entries += 1;
            if let Some(w) = self.gaf.as_mut() {
                w.write_all(&b.gaf[if j == 0 { 0 } else { b.gaf_end[j - 1] }..b.gaf_end[j]])?;
            }
        }
        if b.out_of_range > 0 {
            eprintln!("path {}: {} of {} records have path_bp past the path end; skipped", b.name, b.out_of_range, b.num_records);
        }
        if overruns > 0 {
            eprintln!("path {}: {overruns} matches run off the path end; coverage kept, no GAF row", b.name);
        }
        Ok(())
    }
}

/// Returns per-node coverage (bp, indexed by node_id - min_id) and the node table.
pub fn run(gfa_path: &str, inp: &Inputs, threads: usize, gaf: Option<&mut dyn Write>) -> io::Result<(Vec<f64>, Segments)> {
    let pool = rayon::ThreadPoolBuilder::new().num_threads(threads).build().map_err(|e| invalid(e.to_string()))?;
    // Small batches: peak memory is one batch of step arrays plus its
    // contributions, not the whole GFA.
    let batch = (threads * 2).max(2);
    let mut reader = gfa::open(gfa_path)?;
    let mut segs: Vec<(usize, usize)> = Vec::new();
    let mut seg: Option<Segments> = None;
    let mut m = Merge { coverage: Vec::new(), seen: None, skipped: 0, entries: 0, gaf };
    // P-lines are read straight into the batch arena; every other line is
    // read at the arena's tail and truncated away, so nothing is copied twice.
    let mut arena: Vec<u8> = Vec::new();
    let mut spans: Vec<(usize, usize)> = Vec::with_capacity(batch);
    let (t0, mut t_work, mut t_merge, mut lines) = (Instant::now(), 0.0, 0.0, 0u64);

    loop {
        let start = arena.len();
        let eof = reader.read_until(b'\n', &mut arena)? == 0;
        if !eof {
            lines += 1;
            match arena[start] {
                b'S' => {
                    // GFA lists S before P in every graph we build (vg convert);
                    // a late S-line would leave its node out of the length table.
                    assert!(seg.is_none(), "GFA S-line after the first P-line; S-lines must come first");
                    let mut f = gfa::trim_eol(&arena[start..]).split(|&c| c == b'\t').skip(1);
                    if let (Some(id), Some(seq)) = (f.next(), f.next()) {
                        segs.push((gfa::parse_uint(id).expect("non-numeric GFA segment id"), seq.len()));
                    }
                    arena.truncate(start);
                    continue;
                }
                b'P' => {
                    if seg.is_none() {
                        let s = Segments::new(&std::mem::take(&mut segs));
                        if inp.dedup {
                            check_key_fields(inp.records, &s);
                        }
                        m.coverage = vec![0.0; s.len.len()];
                        m.seen = inp.dedup.then(|| HashSet::with_capacity(inp.records.len() / 2 + 16));
                        eprintln!("PHASE_TIMING gfa_parse_s: {:.6}", t0.elapsed().as_secs_f64());
                        seg = Some(s);
                    }
                    spans.push((start, arena.len()));
                }
                _ => {
                    arena.truncate(start);
                    continue;
                }
            }
        }
        if spans.len() >= batch || (eof && !spans.is_empty()) {
            let (s, t) = (seg.as_ref().unwrap(), Instant::now());
            let a = &arena;
            let out: Vec<Bucket> = pool.install(|| spans.par_iter().flat_map_iter(|&(i, j)| walk_pline(&a[i..j], s, inp)).collect());
            t_work += t.elapsed().as_secs_f64();
            let t = Instant::now();
            for b in &out {
                m.apply(b)?;
            }
            t_merge += t.elapsed().as_secs_f64();
            spans.clear();
            arena.clear();
        }
        if eof {
            break;
        }
    }
    let seg = seg.unwrap_or_else(|| Segments::new(&segs));
    m.coverage.resize(seg.len.len(), 0.0); // GFA without P-lines
    let total = t0.elapsed().as_secs_f64();
    eprintln!("PHASE_TIMING path_walk_s: {total:.6}");
    eprintln!("PHASE_TIMING   walk_dispatch_s: {t_work:.6}   (parallel per-bucket work)");
    eprintln!("PHASE_TIMING   walk_merge_s: {t_merge:.6}   (ordered merge)");
    eprintln!("PHASE_TIMING   walk_read_s: {:.6}   (serial GFA scan, {lines} lines)", (total - t_work - t_merge).max(0.0));
    eprintln!("Total GAF entries: {}", m.entries);
    if let Some(seen) = &m.seen {
        let n = inp.records.len();
        eprintln!("Dedup: {} unique keys; {} duplicates skipped ({:.2}% of {n} records)", seen.len(), m.skipped, 100.0 * m.skipped as f64 / n as f64);
    }
    Ok((m.coverage, seg))
}

#[cfg(test)]
mod tests {
    use super::*;

    /// cum for node lengths [5, 7, 1, 3]; path length 16.
    const CUM: [usize; 5] = [0, 5, 12, 13, 16];

    #[test]
    fn cursor_lands_on_the_containing_step() {
        // (path_bp, step, offset): boundaries, a 1-bp node, the last bp.
        for (bp, step, off) in [(0, 0, 0), (4, 0, 4), (5, 1, 0), (6, 1, 1), (11, 1, 6), (12, 2, 0), (13, 3, 0), (15, 3, 2)] {
            let mut i = 0;
            advance(bp, &CUM, &mut i);
            assert_eq!((i, bp - CUM[i]), (step, off), "path_bp {bp}");
        }
    }

    #[test]
    fn cursor_is_monotonic_over_a_sorted_bucket() {
        let mut i = 0;
        for (bp, step) in [(0, 0), (1, 0), (1, 0), (5, 1), (5, 1), (6, 1), (12, 2), (13, 3), (15, 3)] {
            advance(bp, &CUM, &mut i);
            assert_eq!(i, step, "path_bp {bp}");
        }
    }

    #[test]
    fn traverse_splits_a_match_and_reports_overrun() {
        let seg = Segments::new(&[(10, 5), (11, 7), (12, 1), (13, 3)]);
        let ids = [10, 11, 12, 13];
        let mut out = Vec::new();
        assert_eq!(traverse(&ids, &seg, 3, 10, &mut out), (13, 0));
        assert_eq!(out, [(0, 2), (1, 7), (2, 1)]);
        out.clear();
        assert_eq!(traverse(&ids[2..], &seg, 0, 6, &mut out), (4, 2));
        assert_eq!(out, [(2, 1), (3, 3)]);
    }
}

//! find_mems outputs: 16-byte hit records, per-bucket start indices, path names.

use bytemuck::{Pod, Zeroable};
use std::collections::HashMap;
use std::fs::File;
use std::io::{self, BufRead, BufReader, Read};

/// One MEM hit from `_path_pos_v2.bin` (little-endian, 16 B), bucketed by
/// seq_id = 2 * path_idx + orientation and sorted by path_bp within a bucket.
/// node_id and offset are derived from path_bp by the walker. There is
/// deliberately no strand field: the BWT hit's is_rev is per-node and
/// independent of bucket parity; asserting the two equal once silently
/// dropped 1.4M valid records.
#[repr(C)]
#[derive(Clone, Copy, Pod, Zeroable)]
pub struct Record {
    pub path_bp: u32,
    pub match_len: u32,
    pub read_st: u32,
    pub read_id: u32,
}
const _: () = assert!(std::mem::size_of::<Record>() == 16);

pub fn invalid(msg: String) -> io::Error {
    io::Error::new(io::ErrorKind::InvalidData, msg)
}

/// Read straight into a calloc'd Vec<Record>: one allocation, Record-aligned.
pub fn read_records(path: &str) -> io::Result<Vec<Record>> {
    let mut f = File::open(path)?;
    let len = f.metadata()?.len() as usize;
    if len % 16 != 0 {
        return Err(invalid(format!("{path}: {len} bytes is not a multiple of 16; expected a _path_pos_v2.bin")));
    }
    let mut recs = bytemuck::zeroed_vec::<Record>(len / 16);
    f.read_exact(bytemuck::cast_slice_mut(&mut recs))?;
    Ok(recs)
}

/// Trimmed lines, skipping blanks and `#` comments.
fn data_lines(path: &str) -> io::Result<Vec<String>> {
    let mut out = Vec::new();
    for line in BufReader::new(File::open(path)?).lines() {
        let line = line?;
        let t = line.trim();
        if !t.is_empty() && !t.starts_with('#') {
            out.push(t.to_string());
        }
    }
    Ok(out)
}

/// `_seq_id_starts.out`: line k is the index of seq_id k's first record.
pub fn read_starts(path: &str) -> io::Result<Vec<usize>> {
    data_lines(path)?
        .iter()
        .map(|l| l.parse().map_err(|_| invalid(format!("{path}: bad bucket start {l:?}"))))
        .collect()
}

/// Path name -> path_idx (its line number in the `.paths` file).
pub fn read_path_names(path: &str) -> io::Result<HashMap<String, usize>> {
    Ok(data_lines(path)?.into_iter().enumerate().map(|(i, n)| (n, i)).collect())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn record_decodes_little_endian() {
        let bytes: [u8; 16] = [4, 3, 2, 1, 0x20, 0x4E, 0, 0, 7, 0, 0, 0, 0x2A, 0, 0, 0];
        let r: Record = bytemuck::pod_read_unaligned(&bytes);
        assert_eq!((r.path_bp, r.match_len, r.read_st, r.read_id), (0x01020304, 20_000, 7, 42));
    }
}

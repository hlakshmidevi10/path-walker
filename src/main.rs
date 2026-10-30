//! path-walker: project find_mems hits onto GFA paths as per-node coverage
//! (and optionally GAF), for MEM-Projection.

mod gfa;
mod input;
mod walk;

use clap::Parser;
use input::invalid;
use std::fs::File;
use std::io::{self, BufWriter, Write};
use std::time::Instant;

#[derive(Parser)]
#[command(version, about)]
struct Args {
    /// GFA graph, plain or gzip/bgzip; every S-line must precede the first P-line
    #[arg(long)]
    gfa: String,
    /// find_mems `<prefix>_path_pos_v2.bin`
    #[arg(long)]
    path_pos: String,
    /// find_mems `<prefix>_seq_id_starts.out`
    #[arg(long)]
    seq_id_starts: String,
    /// Path names, one per line (path_extract `.paths`); line k is path_idx k
    #[arg(long)]
    path_names: String,
    /// Write `<prefix>.gaf` and `<prefix>_coverage.csv`
    #[arg(long)]
    gaf_file_prefix: Option<String>,
    /// Write `<prefix>_coverage.csv` only
    #[arg(long)]
    coverage_prefix: Option<String>,
    /// Divide coverage by node length
    #[arg(short, long)]
    len_scale: bool,
    /// Count each (read_id, read_st, start node, offset) once; needed for find_mems --lightweight-tags output
    #[arg(long)]
    dedup_read_node: bool,
    /// Worker threads; output is byte-identical at any count
    #[arg(long, default_value_t = 8)]
    threads: usize,
}

fn run(a: Args) -> io::Result<()> {
    let t = Instant::now();
    let records = input::read_records(&a.path_pos)?;
    eprintln!("Loaded {} path_pos records", records.len());
    eprintln!("PHASE_TIMING records_load_s: {:.6}", t.elapsed().as_secs_f64());
    let starts = input::read_starts(&a.seq_id_starts)?;
    let names = input::read_path_names(&a.path_names)?;
    if starts.len() != 2 * names.len() + 1 || starts.last() != Some(&records.len()) {
        return Err(invalid(format!(
            "{} bucket starts ending at {:?} do not fit {} paths and {} records; were these files built from the same index?",
            starts.len(), starts.last(), names.len(), records.len())));
    }

    let prefix = a.gaf_file_prefix.as_deref().or(a.coverage_prefix.as_deref()).unwrap_or("output");
    let mut gaf = match &a.gaf_file_prefix {
        Some(p) => {
            let mut w = BufWriter::with_capacity(1 << 20, File::create(format!("{p}.gaf"))?);
            writeln!(w, "read_id\tread_st\tpath_str\tpath_len\tpath_st\tpath_end\tmatch_len\tpath_name")?;
            Some(w)
        }
        None => None,
    };
    let inp = walk::Inputs { records: &records, starts: &starts, names: &names, dedup: a.dedup_read_node, gaf: gaf.is_some() };
    let (coverage, seg) = walk::run(&a.gfa, &inp, a.threads.max(1), gaf.as_mut().map(|w| w as &mut dyn Write))?;
    if let Some(mut w) = gaf {
        w.flush()?;
    }

    let t = Instant::now();
    let path = format!("{prefix}_coverage.csv");
    let mut w = BufWriter::with_capacity(1 << 20, File::create(&path)?);
    writeln!(w, "node_id,node_coverage")?;
    for (i, &v) in coverage.iter().enumerate() {
        let v = if a.len_scale { v / seg.len[i] as f64 } else { v };
        writeln!(w, "{},{:.2}", seg.min_id + i, v)?;
    }
    w.flush()?;
    eprintln!("PHASE_TIMING coverage_write_s: {:.6}", t.elapsed().as_secs_f64());
    eprintln!("Coverage data written to: {path}");
    Ok(())
}

fn main() {
    if let Err(e) = run(Args::parse()) {
        eprintln!("path-walker: {e}");
        std::process::exit(1);
    }
}

# path-walker

Projects `find_mems` hits onto the paths of a pangenome GFA and writes per-node
coverage (optionally GAF). It is the step after `find_mems` in
[MEM-Projection](https://github.com/hlakshmidevi10/MEM-Projection). Forked from
[gafpack](https://github.com/pangenome/gafpack).

## Build

```bash
cargo build --release   # target/release/path-walker
```

## Usage

```bash
path-walker --gfa graph.gfa --path-pos mem_path_pos_v2.bin \
            --seq-id-starts mem_seq_id_starts.out --path-names graph.paths \
            --dedup-read-node --coverage-prefix mem      # writes mem_coverage.csv
```

| Option | |
|---|---|
| `--gfa` | Graph, plain or gzip/bgzip. Must share node ids with the GBZ the index was built from (MEM-Projection derives it with `vg convert -fW --no-translation`). All S-lines must precede the first P-line. |
| `--path-pos`, `--seq-id-starts` | `find_mems` outputs |
| `--path-names` | `path_extract` output; line k is path k |
| `--coverage-prefix P` | write `P_coverage.csv` |
| `--gaf-file-prefix P` | write `P.gaf` and `P_coverage.csv` |
| `--dedup-read-node` | count each (read, read start, node, offset) once; required for `find_mems --lightweight-tags` output |
| `-l, --len-scale` | divide coverage by node length |
| `--threads N` | default 8; output is byte-identical at any thread count |

## Formats

`_path_pos_v2.bin` is an array of 16-byte little-endian records
`u32 path_bp | u32 match_len | u32 read_st | u32 read_id`, grouped by
`seq_id = 2 × path_idx + orientation` and sorted by `path_bp` within each
group. `_seq_id_starts.out` has one line per `seq_id`, plus a final line: the
index of the group's first record. Node and offset are recovered by walking
`path_bp` along the path; the GAF strand is the orientation XOR the step's `+`/`-`.

`_coverage.csv`: `node_id,node_coverage`, with bp covered per node (mean depth with `-l`).

`.gaf`: tab-separated `read_id read_st path_str path_len path_st path_end match_len path_name`,
with a header line. Reverse-orientation rows carry the path name suffixed `_reverse`.

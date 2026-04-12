# BioBB DYNA Command Line Help

Generic pattern (each tool accepts optional `-c` / `--config` for YAML or JSON):

```text
<command> [-h] [-c CONFIG] <required arguments>
```

Configuration files can supply a `properties` dictionary; paths can also be passed on the command line as documented per tool.

---

## dccm

Computes covariance / correlation / LMI-style matrices from an MD trajectory (MDiGest).

### Get help

```bash
dccm -h
```

### Arguments (typical)

| Argument | Description |
|----------|-------------|
| `-i`, `--input_topology_path` | Topology (e.g. PDB) |
| `-f`, `--input_trajectory_path` | Trajectory (e.g. XTC, DCD) |
| `-o`, `--output_matrix_path` | Output `.npz` (contains array `dccm`) |
| `-c`, `--config` | Optional YAML/JSON config with `properties` |

### Example

```bash
dccm -c dccm_config.yml \
  -i topology.pdb \
  -f trajectory.xtc \
  -o matrix.npz
```

---

## create_graph

Builds a residue network (GraphML) from a correlation matrix NPZ and a topology.

### Get help

```bash
create_graph -h
```

### Arguments (typical)

| Argument | Description |
|----------|-------------|
| `-i`, `--input_matrix_path` | Input `.npz` (key `dccm`) |
| `-s`, `--input_topology_path` | Topology PDB |
| `-o`, `--output_matrix_path` | Output GraphML path (filename; content is graph) |
| `-ct`, `--input_contact_matrix_path` | Optional contact `.npz` |
| `-c`, `--config` | Optional YAML/JSON config |

### Example

```bash
create_graph -c graph_config.yml \
  -i matrix.npz \
  -s topology.pdb \
  -o network.graphml
```

---

## communities

Community detection on a GraphML graph (NetworkX).

### Get help

```bash
communities -h
```

### Arguments

| Argument | Description |
|----------|-------------|
| `-i`, `--input_graph_path` | Input GraphML |
| `-o`, `--output_communities_path` | Output JSON |
| `-c`, `--config` | Optional YAML/JSON (`properties`: `algorithm`, `resolution`, …) |

### Example

```bash
communities -c communities_config.yml \
  -i network.graphml \
  -o communities.json
```

---

## shortest_paths

Shortest paths between two residues (by `resid` stored on graph nodes).

### Get help

```bash
shortest_paths -h
```

### Arguments

| Argument | Description |
|----------|-------------|
| `-i`, `--input_graph_path` | Input GraphML |
| `-o`, `--output_shortest_paths_path` | Output JSON |
| `-c`, `--config` | Optional YAML/JSON (`properties`: `source_residue`, `sink_residue`, `num_paths`, …) |

### Example

```bash
shortest_paths -c paths_config.yml \
  -i network.graphml \
  -o paths.json
```

---

## JSON schemas

Machine-readable descriptions of I/O and `properties` for each block live under `biobb_dyna/json_schemas/` in the package (e.g. `dccm.json`, `create_graph.json`, `communities.json`, `shortest_paths.json`, `biobb_dyna.json`).

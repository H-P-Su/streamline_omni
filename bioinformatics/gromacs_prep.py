"""
GROMACS MD Preparation Tool
============================
Prepares a PDB file for molecular dynamics simulation with GROMACS.
Generates all required input files and a ready-to-run shell script.

Required packages:
    pip install streamlit biopython pdbfixer acpype

Note on pdbfixer / OpenMM:
    If `pip install pdbfixer` fails, use conda:
        conda install -c conda-forge openmm pdbfixer

Note on acpype / AmberTools:
    acpype requires antechamber from AmberTools for ligand parametrization.
    Install AmberTools first:
        conda install -c conda-forge ambertools
    Then:
        pip install acpype
"""

import io
import zipfile

import streamlit as st
from Bio.PDB import PDBParser, PDBIO, Select

# Optional: pdbfixer for structure repair
try:
    from pdbfixer import PDBFixer
    from openmm.app import PDBFile
    PDBFIXER_AVAILABLE = True
except ImportError:
    PDBFIXER_AVAILABLE = False

# ─── Constants ───────────────────────────────────────────────────────────────

STANDARD_RESIDUES = {
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
    "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
    "HID", "HIE", "HIP", "HSD", "HSE", "HSP",  # histidine variants
    "CYX", "CYM",                                # cysteine variants
    "ASH", "GLH",                                # protonated forms
    "ACE", "NME",                                # terminal caps
}

WATER_RESIDUES = {"HOH", "WAT", "TIP", "SOL", "H2O", "OHH"}

FORCE_FIELDS = {
    "AMBER99SB-ILDN": "amber99sb-ildn",
    "AMBER03":        "amber03",
    "CHARMM27":       "charmm27",
    "CHARMM36":       "charmm36-jul2022",
    "GROMOS96 54A7":  "gromos96",
    "OPLS-AA":        "oplsaa",
}

WATER_MODELS = {
    "TIP3P (recommended for AMBER)": "tip3p",
    "TIP4P":                         "tip4p",
    "SPC":                           "spc",
    "SPC/E":                         "spce",
}

BOX_TYPES = {
    "Dodecahedron (recommended)": "dodecahedron",
    "Cubic":                      "cubic",
    "Triclinic":                  "triclinic",
}

# ─── PDB Parsing ─────────────────────────────────────────────────────────────

def parse_pdb_info(pdb_string):
    """Parse PDB string and return summary info about chains, ligands, and waters."""
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("input", io.StringIO(pdb_string))

    protein_chains = []
    ligands = []
    water_count = 0

    for model in structure:
        for chain in model:
            chain_id = chain.get_id()
            aa_residues = []
            chain_waters = 0

            for residue in chain:
                resname = residue.get_resname().strip()
                hetflag = residue.get_id()[0]

                if resname in WATER_RESIDUES:
                    chain_waters += 1
                elif hetflag == " " or resname in STANDARD_RESIDUES:
                    aa_residues.append(resname)
                elif hetflag.startswith("H_"):
                    ligands.append({
                        "resname": resname,
                        "chain":   chain_id,
                        "resnum":  residue.get_id()[1],
                    })

            water_count += chain_waters
            if aa_residues:
                protein_chains.append({
                    "chain_id":   chain_id,
                    "n_residues": len(aa_residues),
                })

    return {
        "protein_chains": protein_chains,
        "ligands":        ligands,
        "water_count":    water_count,
        "structure":      structure,
        "n_atoms":        sum(1 for _ in structure.get_atoms()),
    }

# ─── PDB Cleaning ────────────────────────────────────────────────────────────

class ProteinSelect(Select):
    """Select only standard protein residues, excluding waters and alternate conformations."""
    def accept_residue(self, residue):
        resname = residue.get_resname().strip()
        hetflag = residue.get_id()[0]
        return (hetflag == " " or resname in STANDARD_RESIDUES) and resname not in WATER_RESIDUES

    def accept_atom(self, atom):
        return atom.get_altloc() in (" ", "A")


class LigandSelect(Select):
    """Select a single ligand residue by name, chain, and residue number."""
    def __init__(self, resname, chain, resnum):
        self.resname = resname
        self.chain   = chain
        self.resnum  = resnum

    def accept_residue(self, residue):
        return (
            residue.get_resname().strip() == self.resname
            and residue.get_parent().get_id() == self.chain
            and residue.get_id()[1] == self.resnum
        )

    def accept_atom(self, atom):
        return atom.get_altloc() in (" ", "A")


def structure_to_pdb_string(structure, select=None):
    """Serialise a BioPython structure to a PDB-format string."""
    output = io.StringIO()
    io_obj = PDBIO()
    io_obj.set_structure(structure)
    io_obj.save(output, select) if select else io_obj.save(output)
    return output.getvalue()


def fix_with_pdbfixer(pdb_string):
    """Use PDBFixer to repair structure (requires openmm + pdbfixer installed)."""
    notes = []
    fixer = PDBFixer(pdbfile=io.StringIO(pdb_string))
    fixer.findMissingResidues()
    fixer.findNonstandardResidues()

    if fixer.missingResidues:
        notes.append(f"Found {len(fixer.missingResidues)} missing residue segment(s)")

    fixer.replaceNonstandardResidues()
    fixer.removeHeterogens(keepWater=False)
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    notes.append("Added missing heavy atoms via PDBFixer")

    output = io.StringIO()
    PDBFile.writeFile(fixer.topology, fixer.positions, output)
    return output.getvalue(), notes

# ─── MDP File Generators ─────────────────────────────────────────────────────

def generate_em_mdp(nsteps=50000):
    return f"""; Energy Minimization
; Run parameters
integrator  = steep         ; Steepest descent
emtol       = 1000.0        ; Stop when max force < 1000 kJ/mol/nm
emstep      = 0.01          ; Initial step size (nm)
nsteps      = {nsteps}

; Neighbour searching
cutoff-scheme   = Verlet
ns_type         = grid
nstlist         = 1
rcoulomb        = 1.0
rvdw            = 1.0

; Electrostatics
coulombtype     = PME
pme_order       = 4
fourierspacing  = 0.16

; Output
nstxout     = 0
nstvout     = 0
nstenergy   = 500
nstlog      = 500
pbc         = xyz
"""


def generate_nvt_mdp(nsteps=50000, temp=300, dt=0.002):
    runtime_ps = nsteps * dt
    return f"""; NVT Equilibration ({runtime_ps:.0f} ps)
define          = -DPOSRES      ; Position restrain protein heavy atoms

; Run parameters
integrator      = md
nsteps          = {nsteps}
dt              = {dt}

; Output
nstxout         = 500
nstvout         = 500
nstenergy       = 500
nstlog          = 500

; Bond constraints
continuation            = no
constraint_algorithm    = lincs
constraints             = h-bonds
lincs_iter              = 1
lincs_order             = 4

; Neighbour searching
cutoff-scheme   = Verlet
ns_type         = grid
nstlist         = 10
rcoulomb        = 1.0
rvdw            = 1.0
DispCorr        = EnerPres

; Electrostatics
coulombtype     = PME
pme_order       = 4
fourierspacing  = 0.16

; Temperature coupling
tcoupl      = V-rescale
tc-grps     = Protein Non-Protein
tau_t       = 0.1   0.1
ref_t       = {temp}  {temp}

; Pressure coupling
pcoupl      = no

; Velocity generation
gen_vel     = yes
gen_temp    = {temp}
gen_seed    = -1
pbc         = xyz
"""


def generate_npt_mdp(nsteps=50000, temp=300, pressure=1.0, dt=0.002):
    runtime_ps = nsteps * dt
    return f"""; NPT Equilibration ({runtime_ps:.0f} ps)
define          = -DPOSRES      ; Position restrain protein heavy atoms

; Run parameters
integrator      = md
nsteps          = {nsteps}
dt              = {dt}

; Output
nstxout         = 500
nstvout         = 500
nstenergy       = 500
nstlog          = 500

; Bond constraints
continuation            = yes
constraint_algorithm    = lincs
constraints             = h-bonds
lincs_iter              = 1
lincs_order             = 4

; Neighbour searching
cutoff-scheme   = Verlet
ns_type         = grid
nstlist         = 10
rcoulomb        = 1.0
rvdw            = 1.0
DispCorr        = EnerPres

; Electrostatics
coulombtype     = PME
pme_order       = 4
fourierspacing  = 0.16

; Temperature coupling
tcoupl      = V-rescale
tc-grps     = Protein Non-Protein
tau_t       = 0.1   0.1
ref_t       = {temp}  {temp}

; Pressure coupling
pcoupl          = Parrinello-Rahman
pcoupltype      = isotropic
tau_p           = 2.0
ref_p           = {pressure}
compressibility = 4.5e-5
refcoord_scaling = com

; Velocity generation
gen_vel     = no
pbc         = xyz
"""


def generate_md_mdp(nsteps=500000, temp=300, pressure=1.0, dt=0.002):
    runtime_ns = nsteps * dt / 1000
    return f"""; Production MD ({runtime_ns:.2f} ns)

; Run parameters
integrator      = md
nsteps          = {nsteps}
dt              = {dt}

; Output (coordinates every 10 ps, energy every 10 ps)
nstxout             = 0
nstvout             = 0
nstenergy           = 5000
nstlog              = 5000
nstxout-compressed  = 5000
compressed-x-grps   = System

; Bond constraints
continuation            = yes
constraint_algorithm    = lincs
constraints             = h-bonds
lincs_iter              = 1
lincs_order             = 4

; Neighbour searching
cutoff-scheme   = Verlet
ns_type         = grid
nstlist         = 10
rcoulomb        = 1.0
rvdw            = 1.0
DispCorr        = EnerPres

; Electrostatics
coulombtype     = PME
pme_order       = 4
fourierspacing  = 0.16

; Temperature coupling
tcoupl      = V-rescale
tc-grps     = Protein Non-Protein
tau_t       = 0.1   0.1
ref_t       = {temp}  {temp}

; Pressure coupling
pcoupl          = Parrinello-Rahman
pcoupltype      = isotropic
tau_p           = 2.0
ref_p           = {pressure}
compressibility = 4.5e-5

; Velocity generation
gen_vel     = no
pbc         = xyz
"""

# ─── Shell Script Generator ──────────────────────────────────────────────────

def generate_run_script(info, settings):
    ff        = settings["force_field"]
    water     = settings["water_model"]
    box_type  = settings["box_type"]
    box_dist  = settings["box_distance"]
    salt_conc = settings["salt_concentration"]
    ligands   = info["ligands"]
    has_ligands = bool(ligands)

    unique_ligands = {}
    for lig in ligands:
        key = f"{lig['resname']}_{lig['chain']}_{lig['resnum']}"
        unique_ligands[key] = lig

    lines = [
        "#!/bin/bash",
        "# GROMACS MD Preparation Script",
        "# Generated by Streamline Omni — GROMACS Prep Tool",
        "#",
        "# Requirements:",
        "#   - GROMACS installed and in PATH  (https://www.gromacs.org)",
    ]
    if has_ligands:
        lines += [
            "#   - AmberTools installed (antechamber in PATH)",
            "#     conda install -c conda-forge ambertools",
            "#   - acpype installed",
            "#     pip install acpype",
        ]
    lines += [
        "#",
        "# Usage:",
        "#   bash run_gromacs.sh",
        "",
        "set -euo pipefail  # Exit on error, unset variable, or pipe failure",
        "",
        "echo '================================================'",
        "echo ' GROMACS MD Preparation'",
        "echo '================================================'",
        "",
    ]

    # ── Ligand parametrization ──────────────────────────────────────────────
    if has_ligands:
        lines += [
            "# ──────────────────────────────────────────────────────",
            "# STEP 0: Ligand parametrization with ACPYPE (GAFF2)",
            "# ──────────────────────────────────────────────────────",
            "echo '[0/8] Parametrizing ligands...'",
            "",
        ]
        for key, lig in unique_ligands.items():
            fname = f"ligand_{lig['resname']}.pdb"
            lines += [
                f"acpype -i {fname} -c bcc -a gaff2 -n 0",
                f"# Produces: {lig['resname']}.acpype/{lig['resname']}_GMX.itp",
                f"#           {lig['resname']}.acpype/{lig['resname']}_GMX.gro",
                "",
            ]

    # ── pdb2gmx ─────────────────────────────────────────────────────────────
    lines += [
        "# ──────────────────────────────────────────────────────",
        "# STEP 1: Generate protein topology (pdb2gmx)",
        "# ──────────────────────────────────────────────────────",
        "echo '[1/8] Running pdb2gmx...'",
        f"gmx pdb2gmx -f protein_clean.pdb -o protein_processed.gro \\",
        f"            -water {water} -ff {ff} -ignh",
        "",
    ]

    # ── Merge ligand topology ────────────────────────────────────────────────
    if has_ligands:
        lines += [
            "# ──────────────────────────────────────────────────────",
            "# STEP 1b: Merge ligand topology into topol.top",
            "# ──────────────────────────────────────────────────────",
            "echo '[1b/8] Merging ligand topologies...'",
            "python3 merge_topology.py",
            "",
            "# Combine protein and ligand GRO files into system.gro",
        ]
        gro_files = ["protein_processed.gro"] + [
            f"{lig['resname']}.acpype/{lig['resname']}_GMX.gro"
            for lig in unique_ligands.values()
        ]
        lines += [
            f"python3 combine_gro.py {' '.join(gro_files)} system.gro",
            "",
        ]
    else:
        lines += ["cp protein_processed.gro system.gro", ""]

    # ── editconf ─────────────────────────────────────────────────────────────
    lines += [
        "# ──────────────────────────────────────────────────────",
        "# STEP 2: Define simulation box (editconf)",
        "# ──────────────────────────────────────────────────────",
        "echo '[2/8] Defining simulation box...'",
        f"gmx editconf -f system.gro -o box.gro -c -d {box_dist} -bt {box_type}",
        "",
        "# ──────────────────────────────────────────────────────",
        "# STEP 3: Solvate system",
        "# ──────────────────────────────────────────────────────",
        "echo '[3/8] Solvating...'",
        "gmx solvate -cp box.gro -cs spc216.gro -o solv.gro -p topol.top",
        "",
        "# ──────────────────────────────────────────────────────",
        "# STEP 4: Add ions (neutralise + set salt concentration)",
        "# ──────────────────────────────────────────────────────",
        "echo '[4/8] Adding ions...'",
        "gmx grompp -f em.mdp -c solv.gro -p topol.top -o ions.tpr -maxwarn 2",
        f'echo "SOL" | gmx genion -s ions.tpr -o solv_ions.gro \\',
        f'    -p topol.top -pname NA -nname CL -neutral -conc {salt_conc}',
        "",
        "# ──────────────────────────────────────────────────────",
        "# STEP 5: Energy minimization",
        "# ──────────────────────────────────────────────────────",
        "echo '[5/8] Running energy minimization...'",
        "gmx grompp -f em.mdp -c solv_ions.gro -p topol.top -o em.tpr",
        "gmx mdrun -v -deffnm em",
        "",
        "# ──────────────────────────────────────────────────────",
        "# STEP 6: NVT equilibration",
        "# ──────────────────────────────────────────────────────",
        "echo '[6/8] Running NVT equilibration...'",
        "gmx grompp -f nvt.mdp -c em.gro -r em.gro -p topol.top -o nvt.tpr",
        "gmx mdrun -deffnm nvt",
        "",
        "# ──────────────────────────────────────────────────────",
        "# STEP 7: NPT equilibration",
        "# ──────────────────────────────────────────────────────",
        "echo '[7/8] Running NPT equilibration...'",
        "gmx grompp -f npt.mdp -c nvt.gro -r nvt.gro -t nvt.cpt -p topol.top -o npt.tpr",
        "gmx mdrun -deffnm npt",
        "",
        "# ──────────────────────────────────────────────────────",
        "# STEP 8: Prepare production MD TPR",
        "# ──────────────────────────────────────────────────────",
        "echo '[8/8] Preparing production MD...'",
        "gmx grompp -f md.mdp -c npt.gro -t npt.cpt -p topol.top -o md.tpr",
        "",
        "echo ''",
        "echo '================================================'",
        "echo ' System is ready for production MD!'",
        "echo ''",
        "echo ' Run on CPU:'",
        "echo '   gmx mdrun -v -deffnm md -ntmpi 1 -ntomp 8'",
        "echo ''",
        "echo ' Run with GPU:'",
        "echo '   gmx mdrun -v -deffnm md -ntmpi 1 -ntomp 8 -gpu_id 0'",
        "echo ''",
        "echo ' On an HPC cluster, submit md.tpr via your scheduler.'",
        "echo '================================================'",
    ]

    return "\n".join(lines)

# ─── Helper Script: merge_topology.py ────────────────────────────────────────

def generate_merge_topology_script(ligands):
    """Python helper script to insert ligand ITP entries into topol.top."""
    unique_resnames = list({lig["resname"] for lig in ligands})

    itp_includes = "\n".join(
        f'    f\'#include "{r}.acpype/{r}_GMX.itp"\\n\','
        for r in unique_resnames
    )
    mol_entries = "\n".join(
        f'    f"{r:<20} 1\\n",'
        for r in unique_resnames
    )

    return f'''#!/usr/bin/env python3
"""
Helper: merge ligand ITP files and molecule entries into topol.top.
Run this after pdb2gmx and acpype, before editconf.
"""

with open("topol.top", "r") as f:
    content = f.read()

# Insert ligand #include lines immediately before [ system ]
ligand_includes = (
{itp_includes}
)
content = content.replace("[ system ]", ligand_includes + "[ system ]")

# Append ligand entries to [ molecules ] section
ligand_molecules = (
{mol_entries}
)
content = content.rstrip() + "\\n" + "".join(ligand_molecules)

with open("topol.top", "w") as f:
    f.write(content)

print("topol.top updated with ligand ITP includes and molecule entries.")
'''

# ─── Helper Script: combine_gro.py ───────────────────────────────────────────

COMBINE_GRO_SCRIPT = '''#!/usr/bin/env python3
"""
Helper: combine multiple GROMACS GRO files into one.

Usage:
    python3 combine_gro.py file1.gro file2.gro [file3.gro ...] output.gro
"""

import sys

def read_gro(path):
    with open(path) as f:
        lines = f.readlines()
    title  = lines[0]
    natoms = int(lines[1].strip())
    atoms  = lines[2:2 + natoms]
    box    = lines[2 + natoms]
    return title, natoms, atoms, box

if len(sys.argv) < 4:
    print("Usage: combine_gro.py file1.gro file2.gro [...] output.gro")
    sys.exit(1)

inputs = sys.argv[1:-1]
output = sys.argv[-1]

all_atoms = []
total     = 0
box_line  = None
title     = None

for path in inputs:
    t, n, atoms, box = read_gro(path)
    if title is None:
        title    = t.strip()
        box_line = box
    all_atoms.extend(atoms)
    total += n

with open(output, "w") as f:
    f.write(f"{title}\\n")
    f.write(f"{total}\\n")
    for i, line in enumerate(all_atoms, 1):
        # Re-number atom index (columns 15-20)
        if len(line) >= 20:
            f.write(line[:15] + f"{i % 100000:5d}" + line[20:])
        else:
            f.write(line)
    f.write(box_line)

print(f"Combined {len(inputs)} GRO files → {total} atoms → {output}")
'''

# ─── README Generator ────────────────────────────────────────────────────────

def generate_readme(info, settings):
    ligands     = info["ligands"]
    has_ligands = bool(ligands)
    ns          = settings["md_steps"] * settings["dt"] / 1000

    lines = [
        "GROMACS MD Preparation Package",
        "=" * 50,
        "Generated by Streamline Omni — GROMACS Prep Tool",
        "",
        "SYSTEM SUMMARY",
        "-" * 30,
        f"Protein chains : {len(info['protein_chains'])}",
    ]
    for ch in info["protein_chains"]:
        lines.append(f"  Chain {ch['chain_id']}: {ch['n_residues']} residues")
    lines += [
        f"Ligands        : {len(ligands)}",
    ]
    for lig in ligands:
        lines.append(f"  {lig['resname']} (chain {lig['chain']}, residue {lig['resnum']})")
    lines += [
        f"Crystal waters : {info['water_count']} (removed)",
        f"Total atoms    : {info['n_atoms']}",
        "",
        "SIMULATION SETTINGS",
        "-" * 30,
        f"Force field    : {settings['force_field']}",
        f"Water model    : {settings['water_model']}",
        f"Box type       : {settings['box_type']}",
        f"Box distance   : {settings['box_distance']} nm",
        f"Salt conc      : {settings['salt_concentration']} M NaCl",
        f"Temperature    : {settings['temperature']} K",
        f"Pressure       : {settings['pressure']} bar",
        f"Timestep       : {settings['dt'] * 1000:.0f} fs",
        f"NVT            : {settings['nvt_steps']} steps  ({settings['nvt_steps'] * settings['dt']:.0f} ps)",
        f"NPT            : {settings['npt_steps']} steps  ({settings['npt_steps'] * settings['dt']:.0f} ps)",
        f"Production     : {settings['md_steps']} steps  ({ns:.2f} ns)",
        "",
        "FILES IN THIS PACKAGE",
        "-" * 30,
        "  protein_clean.pdb     Cleaned protein structure",
    ]
    if has_ligands:
        for resname in {lig["resname"] for lig in ligands}:
            lines.append(f"  ligand_{resname}.pdb       Extracted ligand structure")
        lines += [
            "  merge_topology.py     Merges ligand ITP files into topol.top",
            "  combine_gro.py        Combines GRO files",
        ]
    lines += [
        "  em.mdp                Energy minimization parameters",
        "  nvt.mdp               NVT equilibration parameters",
        "  npt.mdp               NPT equilibration parameters",
        "  md.mdp                Production MD parameters",
        "  run_gromacs.sh        Complete preparation script — RUN THIS",
        "",
        "HOW TO RUN",
        "-" * 30,
        "Requirements:",
        "  - GROMACS  https://www.gromacs.org",
    ]
    if has_ligands:
        lines += [
            "  - AmberTools (antechamber)",
            "      conda install -c conda-forge ambertools",
            "  - acpype",
            "      pip install acpype",
        ]
    lines += [
        "",
        "Steps:",
        "  1. Unzip this package to a clean working directory",
        "  2. Run:  bash run_gromacs.sh",
        "  3. When the script finishes, launch production MD:",
        "",
        "     CPU only:",
        "       gmx mdrun -v -deffnm md -ntmpi 1 -ntomp 8",
        "",
        "     With GPU:",
        "       gmx mdrun -v -deffnm md -ntmpi 1 -ntomp 8 -gpu_id 0",
        "",
        "     HPC cluster (example SLURM):",
        "       srun gmx_mpi mdrun -v -deffnm md",
        "",
        "Output files after production MD:",
        "  md.xtc   — compressed trajectory",
        "  md.edr   — energy file",
        "  md.log   — log file",
        "  md.gro   — final coordinates",
    ]
    return "\n".join(lines)

# ─── ZIP Packaging ───────────────────────────────────────────────────────────

def create_zip(files):
    """Pack a dict of {filename: str_content} into a zip and return bytes."""
    buf = io.BytesIO()
    with zipfile.ZipFile(buf, "w", zipfile.ZIP_DEFLATED) as zf:
        for fname, content in files.items():
            zf.writestr(fname, content)
    buf.seek(0)
    return buf.read()

# ─── Streamlit App ───────────────────────────────────────────────────────────

st.title("GROMACS MD Preparation")
st.caption(
    "Upload a PDB file to generate all input files for a GROMACS MD simulation. "
    "Download the zip, unzip it, and run `bash run_gromacs.sh` on any machine with GROMACS installed."
)

# ── Sidebar ──────────────────────────────────────────────────────────────────
with st.sidebar:
    st.header("Simulation Settings")

    ff_label    = st.selectbox("Force Field",  list(FORCE_FIELDS.keys()), index=0)
    water_label = st.selectbox("Water Model",  list(WATER_MODELS.keys()), index=0)
    box_label   = st.selectbox("Box Type",     list(BOX_TYPES.keys()),    index=0)
    box_distance = st.number_input(
        "Box edge distance (nm)", min_value=0.5, max_value=3.0, value=1.0, step=0.1
    )
    salt_conc = st.number_input(
        "Salt concentration (M NaCl)", min_value=0.0, max_value=1.0, value=0.15, step=0.01
    )

    st.subheader("MD Parameters")
    temperature = st.number_input("Temperature (K)", min_value=273, max_value=373, value=300, step=1)
    pressure    = st.number_input("Pressure (bar)",  min_value=0.5, max_value=2.0, value=1.0, step=0.1)
    dt_fs       = st.selectbox("Timestep (fs)", [1, 2, 4], index=1)

    st.subheader("Simulation Length")
    em_steps  = st.number_input("EM steps",         min_value=1000,  max_value=500000,   value=50000,  step=1000)
    nvt_steps = st.number_input("NVT steps",        min_value=1000,  max_value=5000000,  value=50000,  step=1000)
    npt_steps = st.number_input("NPT steps",        min_value=1000,  max_value=5000000,  value=50000,  step=1000)
    md_steps  = st.number_input("Production steps", min_value=10000, max_value=100000000, value=500000, step=10000)

    dt_ps = dt_fs / 1000
    st.caption(f"Production run = {md_steps * dt_ps / 1000:.2f} ns")

    if PDBFIXER_AVAILABLE:
        use_pdbfixer = st.checkbox("Fix missing atoms (PDBFixer)", value=True)
    else:
        use_pdbfixer = False
        st.info("pdbfixer not installed — structure repair disabled.\nInstall with: pip install pdbfixer")

settings = {
    "force_field":        FORCE_FIELDS[ff_label],
    "water_model":        WATER_MODELS[water_label],
    "box_type":           BOX_TYPES[box_label],
    "box_distance":       box_distance,
    "salt_concentration": salt_conc,
    "temperature":        int(temperature),
    "pressure":           pressure,
    "dt":                 dt_ps,
    "em_steps":           int(em_steps),
    "nvt_steps":          int(nvt_steps),
    "npt_steps":          int(npt_steps),
    "md_steps":           int(md_steps),
}

# ── File upload ───────────────────────────────────────────────────────────────
uploaded_file = st.file_uploader("Upload PDB file", type=["pdb"])

if uploaded_file is not None:
    pdb_string = uploaded_file.read().decode("utf-8")

    with st.spinner("Parsing PDB..."):
        info = parse_pdb_info(pdb_string)

    # System summary
    st.subheader("System Summary")
    col1, col2, col3 = st.columns(3)
    col1.metric("Protein Chains", len(info["protein_chains"]))
    col2.metric("Ligands",        len(info["ligands"]))
    col3.metric("Total Atoms",    info["n_atoms"])

    if info["protein_chains"]:
        with st.expander("Chain details"):
            for ch in info["protein_chains"]:
                st.write(f"Chain **{ch['chain_id']}**: {ch['n_residues']} residues")

    if info["ligands"]:
        with st.expander("Ligands detected"):
            for lig in info["ligands"]:
                st.write(f"**{lig['resname']}** — chain {lig['chain']}, residue {lig['resnum']}")

    if info["water_count"] > 0:
        st.info(f"{info['water_count']} crystallographic water(s) found — will be removed.")

    if not info["protein_chains"]:
        st.warning("No standard protein residues detected. Please check your PDB file.")

    st.divider()

    if st.button("Prepare MD Files", type="primary"):
        files = {}
        notes = []

        progress = st.progress(0, text="Starting...")

        # 1. Clean protein PDB
        progress.progress(10, text="Cleaning protein structure...")
        if use_pdbfixer:
            clean_pdb, fix_notes = fix_with_pdbfixer(pdb_string)
            notes.extend(fix_notes)
        else:
            clean_pdb = structure_to_pdb_string(info["structure"], ProteinSelect())
            notes.append("BioPython used for cleaning (pdbfixer not available)")
        files["protein_clean.pdb"] = clean_pdb

        # 2. Extract ligands
        progress.progress(25, text="Extracting ligands...")
        for lig in info["ligands"]:
            lig_pdb = structure_to_pdb_string(
                info["structure"],
                LigandSelect(lig["resname"], lig["chain"], lig["resnum"]),
            )
            files[f"ligand_{lig['resname']}.pdb"] = lig_pdb

        # 3. MDP files
        progress.progress(40, text="Generating MDP parameter files...")
        files["em.mdp"]  = generate_em_mdp(settings["em_steps"])
        files["nvt.mdp"] = generate_nvt_mdp(settings["nvt_steps"], settings["temperature"], settings["dt"])
        files["npt.mdp"] = generate_npt_mdp(settings["npt_steps"], settings["temperature"], settings["pressure"], settings["dt"])
        files["md.mdp"]  = generate_md_mdp(settings["md_steps"],  settings["temperature"], settings["pressure"], settings["dt"])

        # 4. Helper scripts for ligands
        progress.progress(55, text="Generating helper scripts...")
        if info["ligands"]:
            files["merge_topology.py"] = generate_merge_topology_script(info["ligands"])
            files["combine_gro.py"]    = COMBINE_GRO_SCRIPT

        # 5. Shell run script
        progress.progress(70, text="Generating run_gromacs.sh...")
        files["run_gromacs.sh"] = generate_run_script(info, settings)

        # 6. README
        progress.progress(85, text="Generating README...")
        files["README.txt"] = generate_readme(info, settings)

        # 7. Package
        progress.progress(95, text="Packaging zip...")
        zip_bytes = create_zip(files)

        progress.progress(100, text="Done!")
        st.success(f"Package ready — {len(files)} files.")

        if notes:
            with st.expander("Processing notes"):
                for note in notes:
                    st.write(f"- {note}")

        # Run command display
        st.subheader("Run Commands")
        st.caption("After running `bash run_gromacs.sh`, launch production MD:")
        st.code("gmx mdrun -v -deffnm md -ntmpi 1 -ntomp 8", language="bash")
        st.caption("With GPU:")
        st.code("gmx mdrun -v -deffnm md -ntmpi 1 -ntomp 8 -gpu_id 0", language="bash")

        with st.expander("Files in package"):
            for fname in files:
                st.write(f"- `{fname}`")

        pdb_stem = uploaded_file.name.removesuffix(".pdb")
        st.download_button(
            label="Download MD Package (.zip)",
            data=zip_bytes,
            file_name=f"{pdb_stem}_gromacs_md.zip",
            mime="application/zip",
            type="primary",
        )

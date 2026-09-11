"""Readers for CP-PAW input and output files.

Supported files:

* ``*.prot`` protocol files (plain or gzipped): total energies, k-point count,
  band count, chemical potential, linear-constraint values, the per-step
  minimisation trace and the final atom list with positions and forces.
* ``*.strc`` and ``*.strc_out`` structure files: lattice, atoms, species and the
  linear-constraint block, parsed with a small tokenizer for the CP-PAW
  ``!BLOCK KEY=value ... !END`` grammar.
* ``*.dos`` files written by ``paw_dos.x``: energy grid with occupied and
  unoccupied density of states.
* band files written by ``paw_bands.x``: k-path segments with band energies.

Everything is returned as plain dataclasses and NumPy arrays.
"""

from __future__ import annotations

import gzip
import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterator

import numpy as np

from pdhdiff.constants import BOHR_ANG

# --------------------------------------------------------------------------------------
# generic helpers
# --------------------------------------------------------------------------------------


def read_text(path: str | Path) -> str:
    """Read a text file, transparently decompressing ``.gz`` files.

    Args:
        path: File path. If the suffix is ``.gz`` the file is gunzipped on the fly.

    Returns:
        The file content as a string.
    """
    path = Path(path)
    if path.suffix == ".gz":
        with gzip.open(path, "rt", encoding="utf-8", errors="replace") as fh:
            return fh.read()
    return path.read_text(encoding="utf-8", errors="replace")


def _to_float(token: str) -> float:
    """Convert a Fortran-style number (``1.D-6`` allowed) to float."""
    return float(token.replace("D", "E").replace("d", "e"))


# --------------------------------------------------------------------------------------
# protocol files
# --------------------------------------------------------------------------------------

_ATOM_LINE = re.compile(
    r"^(?P<name>\S+)\s+\(\s*(?P<x>[-\d.]+),\s*(?P<y>[-\d.]+),\s*(?P<z>[-\d.]+)\)"
    r"\s+(?P<mass>[-\d.]+)\s+(?P<mpsi>[-\d.]+)\s+(?P<charge>[-\d.]+)"
    r"(?:\s+\(\s*(?P<fx>[-\d.]+),\s*(?P<fy>[-\d.]+),\s*(?P<fz>[-\d.]+)\))?"
)
_TOTAL_ENERGY = re.compile(r"^TOTAL ENERGY\s*:\s*([-\d.]+)\s*H", re.M)
_CHEM_POT = re.compile(r"^CHEMICAL POTENTIAL\.*:\s*([-\dE.+]+)\s*EV", re.M)
_NKPT = re.compile(r"^NUMBER OF K-POINTS\.*:\s*(\d+)", re.M)
_NBANDS = re.compile(r"^NUMBER OF BANDS\.*:\s*(\d+)", re.M)
_CUTOFF = re.compile(r"^WAVEFUNCTION PLANE WAVE CUTOFF\.*:\s*([-\d.]+)\s*RY", re.M)
_CONSTRAINT = re.compile(r"\(VAL/FORCE\):\(\s*([-\dE.+]+)/\s*([-\dE.+]+)\)")
_TRACE = re.compile(r"^!>\s+(\d+)\s+([-\d.]+)\s+(\d+)\s+([-\d.]+)\s+([-\d.]+)\s+([-\d.]+)", re.M)
_NITER = re.compile(r"^NUMBER OF ITERATIONS\.*:\s*(\d+)", re.M)
_WALL = re.compile(r"^ELAPSED WALLCLOCK TIME\.*:\s*(.+)$", re.M)
_FINISHED = re.compile(r"PROGRAM FINISHED")
_LATTICE = re.compile(r"^T([123])\[ANGSTROM\]=\s*([-\d.]+)\s+([-\d.]+)\s+([-\d.]+)", re.M)


@dataclass
class Atom:
    """One atom of the final CP-PAW atom list.

    Attributes:
        name: Atom label from the structure file, for example ``H_01`` or ``PD17``.
        position: Cartesian position in angstrom.
        mass: Mass in atomic mass units.
        charge: Charge used by the code in units of the elementary charge.
        force: Force in milli-Hartree per bohr, or ``None`` when the report has no forces.
    """

    name: str
    position: np.ndarray
    mass: float
    charge: float
    force: np.ndarray | None = None

    @property
    def element(self) -> str:
        """Element symbol derived from the label (``PD17`` gives ``Pd``)."""
        symbol = re.match(r"[A-Za-z]+", self.name).group(0)
        return symbol.capitalize()


@dataclass
class ProtocolSummary:
    """Quantities extracted from one CP-PAW protocol file.

    Attributes:
        path: Source file.
        n_kpoints: Number of k-points in the irreducible wedge used by the run.
        n_bands: Number of bands.
        cutoff_ry: Plane-wave cutoff for the wave functions in rydberg.
        total_energies: Every ``TOTAL ENERGY`` report in the file, in Hartree.
        chemical_potentials: Every chemical-potential report, in eV.
        constraint_values: ``(value, force)`` pairs of the linear constraint in bohr
            and Hartree per bohr, in the order they were printed.
        trace: Array with columns ``nfi, t_ps, T_K, ekin_psi, e_rho, econs`` for every
            ``!>`` line of the minimisation trace, all runs of the file concatenated.
        traces: The same trace split per run. A protocol file can hold several
            consecutive runs (CP-PAW appends to it), for example a wave-function
            optimisation followed by the restarted atomic relaxation.
        lattice: Lattice vectors as rows in angstrom.
        atoms: Final atom list.
        n_iterations: Number of time steps reported in the run-time report.
        wallclock: Elapsed wall-clock time string from the run-time report.
        finished: Whether the file ends with a ``PROGRAM FINISHED`` line.
    """

    path: Path
    n_kpoints: int | None
    n_bands: int | None
    cutoff_ry: float | None
    total_energies: list[float]
    chemical_potentials: list[float]
    constraint_values: list[tuple[float, float]]
    trace: np.ndarray
    lattice: np.ndarray | None
    traces: list[np.ndarray] = field(default_factory=list)
    atoms: list[Atom] = field(default_factory=list)
    n_iterations: int | None = None
    wallclock: str | None = None
    finished: bool = False

    @property
    def final_energy(self) -> float:
        """Last reported total energy in Hartree."""
        if not self.total_energies:
            raise ValueError(f"no TOTAL ENERGY found in {self.path}")
        return self.total_energies[-1]

    @property
    def n_runs(self) -> int:
        """Number of consecutive CP-PAW runs recorded in the file."""
        return max(1, len(self.traces))

    @property
    def last_trace(self) -> np.ndarray:
        """Minimisation trace of the last run in the file."""
        return self.traces[-1] if self.traces else self.trace

    @property
    def final_constraint_value(self) -> float | None:
        """Last reported value of the linear constraint in bohr, if any."""
        return self.constraint_values[-1][0] if self.constraint_values else None

    @property
    def fermi_energy_ev(self) -> float | None:
        """Last reported chemical potential in eV."""
        return self.chemical_potentials[-1] if self.chemical_potentials else None

    def atom(self, name: str) -> Atom:
        """Return the atom with a given label."""
        for atom in self.atoms:
            if atom.name == name:
                return atom
        raise KeyError(name)

    def positions(self) -> np.ndarray:
        """Cartesian positions of all atoms in angstrom, shape ``(n_atoms, 3)``."""
        return np.array([a.position for a in self.atoms])

    def forces(self) -> np.ndarray:
        """Forces of all atoms in milli-Hartree per bohr, shape ``(n_atoms, 3)``."""
        return np.array([a.force if a.force is not None else np.zeros(3) for a in self.atoms])


def _parse_atomlist_block(block: str) -> tuple[np.ndarray | None, list[Atom]]:
    """Parse one ATOMLIST REPORT block (text starting at the title line)."""
    lattice = np.zeros((3, 3))
    for m in _LATTICE.finditer(block[:2000]):
        lattice[int(m.group(1)) - 1] = [float(m.group(i)) for i in (2, 3, 4)]
    atoms: list[Atom] = []
    started = False
    for line in block.splitlines():
        if line.startswith("NAME"):
            started = True
            continue
        if not started:
            continue
        if not line.strip():
            break
        m = _ATOM_LINE.match(line.strip())
        if m is None:
            break
        pos = np.array([float(m.group(k)) for k in ("x", "y", "z")])
        force = None
        if m.group("fx") is not None:
            force = np.array([float(m.group(k)) for k in ("fx", "fy", "fz")])
        atoms.append(
            Atom(
                name=m.group("name"),
                position=pos,
                mass=float(m.group("mass")),
                charge=float(m.group("charge")),
                force=force,
            )
        )
    if not np.any(lattice):
        lattice = None
    return lattice, atoms


def _parse_atomlist(text: str) -> tuple[np.ndarray | None, list[Atom]]:
    """Parse the last ATOMLIST REPORT block of a protocol."""
    idx = text.rfind("ATOMLIST REPORT")
    if idx < 0:
        return None, []
    return _parse_atomlist_block(text[idx:])


@dataclass
class AtomListReport:
    """One of the periodic ``ATOMLIST REPORT`` blocks of a run, with forces.

    CP-PAW prints an energy report and an atom list every ``IPRINT`` steps and at
    the end of the run. The atom list carries the positions ``R(0)`` and the forces
    used to propagate the atoms in that step (the ``FORCE`` array of the atoms
    object, printed before it is reset by the switch to the next step).

    Attributes:
        nfi: Time-step counter of the last ``!>`` trace line before the report.
        total_energy_h: The ``TOTAL ENERGY`` of the accompanying energy report.
        atoms: Atoms with positions in angstrom and forces in milli-Hartree per bohr.
    """

    nfi: int
    total_energy_h: float
    atoms: list[Atom]

    def atom(self, name: str) -> Atom:
        """Return the atom with a given label."""
        for a in self.atoms:
            if a.name == name:
                return a
        raise KeyError(name)


def parse_reports(path: str | Path, last_run_only: bool = True) -> list[AtomListReport]:
    """Return every atom-list report that carries forces, in order.

    Args:
        path: Protocol file, plain or gzipped.
        last_run_only: Restrict to the last run recorded in the file (the atomic
            relaxation in the two-stage protocols of this project).

    Returns:
        List of :class:`AtomListReport`. Reports without forces (atoms frozen or
        velocities reset, as at the start of a run) are skipped.
    """
    text = read_text(path)
    if last_run_only:
        text = text.split("PROGRAM STARTED")[-1]
    reports = []
    for m in re.finditer(r"ATOMLIST REPORT", text):
        before = text[: m.start()]
        energies = _TOTAL_ENERGY.findall(before)
        steps = re.findall(r"^!>\s+(\d+)", before, re.M)
        if not energies:
            continue
        _, atoms = _parse_atomlist_block(text[m.start() : m.start() + 20000])
        if not atoms or atoms[0].force is None:
            continue
        reports.append(
            AtomListReport(
                nfi=int(steps[-1]) if steps else 0,
                total_energy_h=float(energies[-1]),
                atoms=atoms,
            )
        )
    return reports


def parse_prot(path: str | Path) -> ProtocolSummary:
    """Parse a CP-PAW protocol file.

    Args:
        path: Path to a ``.prot`` or ``.prot.gz`` file.

    Returns:
        A :class:`ProtocolSummary` with energies, constraint values, the minimisation
        trace and the final atom list.
    """
    path = Path(path)
    text = read_text(path)
    energies = [float(m.group(1)) for m in _TOTAL_ENERGY.finditer(text)]
    mus = [_to_float(m.group(1)) for m in _CHEM_POT.finditer(text)]
    nk = _NKPT.search(text)
    nb = _NBANDS.search(text)
    cut = _CUTOFF.search(text)
    constraints = [
        (_to_float(m.group(1)), _to_float(m.group(2))) for m in _CONSTRAINT.finditer(text)
    ]
    traces = []
    for chunk in re.split(r"PROGRAM STARTED", text)[1:] or [text]:
        rows = [[float(v) for v in m.groups()] for m in _TRACE.finditer(chunk)]
        if rows:
            traces.append(np.array(rows))
    trace = np.vstack(traces) if traces else np.zeros((0, 6))
    lattice, atoms = _parse_atomlist(text)
    niter_all = _NITER.findall(text)
    niter = niter_all[-1] if niter_all else None
    wall = _WALL.search(text)
    return ProtocolSummary(
        path=path,
        n_kpoints=int(nk.group(1)) if nk else None,
        n_bands=int(nb.group(1)) if nb else None,
        cutoff_ry=float(cut.group(1)) if cut else None,
        total_energies=energies,
        chemical_potentials=mus,
        constraint_values=constraints,
        trace=trace,
        lattice=lattice,
        traces=traces,
        atoms=atoms,
        n_iterations=int(niter) if niter else None,
        wallclock=wall.group(1).strip() if wall else None,
        finished=bool(_FINISHED.search(text)),
    )


# --------------------------------------------------------------------------------------
# structure files
# --------------------------------------------------------------------------------------

_TOKEN = re.compile(
    r"(?P<block>![A-Za-z_][A-Za-z0-9_]*)"
    r"|(?P<key>[A-Za-z_][A-Za-z0-9_\[\]/().+-]*=)"
    r"|(?P<str>'[^']*')"
    r"|(?P<word>[^\s']+)"
)


@dataclass
class Block:
    """A ``!NAME ... !END`` block of a CP-PAW input file.

    Attributes:
        name: Block name without the leading exclamation mark, upper case.
        entries: Mapping of ``KEY`` (without ``=``) to the list of values that followed it.
        children: Nested blocks in file order.
    """

    name: str
    entries: dict[str, list[str | float]] = field(default_factory=dict)
    children: list["Block"] = field(default_factory=list)

    def find(self, name: str) -> list["Block"]:
        """All direct children with a given name (case insensitive)."""
        return [c for c in self.children if c.name == name.upper()]

    def first(self, name: str) -> "Block":
        """First direct child with a given name."""
        found = self.find(name)
        if not found:
            raise KeyError(f"block !{name.upper()} not found in !{self.name}")
        return found[0]

    def get(self, key: str, default=None):
        """Values of a key of this block, or ``default`` when the key is absent."""
        return self.entries.get(key.upper(), default)

    def scalar(self, key: str, default=None):
        """First value of a key, or ``default``."""
        vals = self.get(key)
        return vals[0] if vals else default


def _tokens(text: str) -> Iterator[tuple[str, str]]:
    for m in _TOKEN.finditer(text):
        kind = m.lastgroup
        yield kind, m.group(kind)


def parse_blocks(text: str) -> Block:
    """Parse CP-PAW block syntax into a tree.

    Args:
        text: Content of a ``.strc``, ``.cntl``, ``.strc_out`` or similar file.

    Returns:
        A root :class:`Block` named ``ROOT`` whose children are the top-level blocks.
    """
    root = Block("ROOT")
    stack = [root]
    current_key: str | None = None
    for kind, tok in _tokens(text):
        if kind == "block":
            name = tok[1:].upper()
            if name == "END":
                stack.pop()
                current_key = None
            elif name == "EOB":
                break
            else:
                blk = Block(name)
                stack[-1].children.append(blk)
                stack.append(blk)
                current_key = None
        elif kind == "key":
            current_key = tok[:-1].upper()
            stack[-1].entries.setdefault(current_key, [])
        elif kind == "str":
            if current_key is not None:
                stack[-1].entries[current_key].append(tok.strip("'"))
        else:
            if current_key is None:
                continue
            try:
                stack[-1].entries[current_key].append(_to_float(tok))
            except ValueError:
                stack[-1].entries[current_key].append(tok)
    return root


@dataclass
class StructureInput:
    """Content of a CP-PAW structure file.

    Attributes:
        lunit_ang: Length unit of the file in angstrom (the lattice constant here).
        lattice: Lattice vectors as rows, in units of ``lunit_ang``.
        names: Atom labels in file order.
        positions: Atom positions in units of ``lunit_ang``, shape ``(n_atoms, 3)``.
        species: Species name of every atom, derived from the label prefix.
        kpoint_density: Value of ``!KPOINTS R``.
        empty_bands: Value of ``!OCCUPATIONS EMPTY``.
        linear_constraints: List of mappings ``atom label -> coefficient vector`` for
            every ``!CONSTRAINTS!LINEAR`` block.
        constraint_move: The ``MOVE`` flag of the first linear constraint, if any.
    """

    lunit_ang: float
    lattice: np.ndarray
    names: list[str]
    positions: np.ndarray
    species: list[str]
    kpoint_density: float | None
    empty_bands: int | None
    linear_constraints: list[dict[str, np.ndarray]]
    constraint_move: bool | None = None

    @property
    def cell_ang(self) -> np.ndarray:
        """Lattice vectors as rows in angstrom."""
        return self.lattice * self.lunit_ang

    def cartesian(self) -> np.ndarray:
        """Cartesian positions in angstrom."""
        return self.positions * self.lunit_ang

    def index(self, name: str) -> int:
        """Index of the atom with a given label."""
        return self.names.index(name)


def parse_strc(path: str | Path) -> StructureInput:
    """Parse a CP-PAW structure file (``.strc`` or ``.strc_out``).

    Args:
        path: Path to the structure file.

    Returns:
        A :class:`StructureInput`.
    """
    root = parse_blocks(read_text(path))
    strc = root.first("STRUCTURE")
    generic = strc.first("GENERIC")
    if generic.get("LUNIT[AA]"):
        lunit = float(generic.scalar("LUNIT[AA]"))
    elif generic.get("LUNIT"):
        lunit = float(generic.scalar("LUNIT")) * BOHR_ANG
    else:
        lunit = BOHR_ANG
    lattice = np.array(strc.first("LATTICE").get("T"), dtype=float).reshape(3, 3)
    names, positions, species = [], [], []
    for atom in strc.find("ATOM"):
        name = str(atom.scalar("NAME"))
        names.append(name)
        positions.append([float(v) for v in atom.get("R")])
        sp = atom.scalar("SP")
        if sp is None:
            sp = re.match(r"[A-Za-z_]+", name).group(0).rstrip("_")
        species.append(str(sp).rstrip("_"))
    kp = strc.find("KPOINTS")
    occ = strc.find("OCCUPATIONS")
    constraints: list[dict[str, np.ndarray]] = []
    move = None
    for cblock in strc.find("CONSTRAINTS"):
        for lin in cblock.find("LINEAR"):
            if move is None and lin.get("MOVE"):
                move = str(lin.scalar("MOVE")).upper().startswith("T")
            entry = {}
            for atom in lin.find("ATOM"):
                entry[str(atom.scalar("NAME"))] = np.array(atom.get("R"), dtype=float)
            constraints.append(entry)
    return StructureInput(
        lunit_ang=lunit,
        lattice=lattice,
        names=names,
        positions=np.array(positions, dtype=float),
        species=species,
        kpoint_density=float(kp[0].scalar("R")) if kp and kp[0].get("R") else None,
        empty_bands=int(occ[0].scalar("EMPTY")) if occ and occ[0].get("EMPTY") else None,
        linear_constraints=constraints,
        constraint_move=move,
    )


# --------------------------------------------------------------------------------------
# density of states and band structure
# --------------------------------------------------------------------------------------


@dataclass
class DosCurve:
    """A density-of-states curve from ``paw_dos.x``.

    Attributes:
        energy_ev: Energy grid in eV (absolute, as written by the code).
        occupied: Occupied part of the DOS in states per eV.
        unoccupied: Unoccupied part of the DOS in states per eV.
    """

    energy_ev: np.ndarray
    occupied: np.ndarray
    unoccupied: np.ndarray

    @property
    def total(self) -> np.ndarray:
        """Occupied plus unoccupied DOS."""
        return self.occupied + self.unoccupied


def parse_dos(path: str | Path) -> DosCurve:
    """Read a ``.dos`` file written by ``paw_dos.x`` for plotting with Grace.

    The file stores two closed polygons on the same energy grid: column 2 is the
    full density of states and column 3 the occupied part, both traversed once
    with zeros (the baseline) and once with the negative curve values. The
    function folds the two passes onto the unique energy grid and returns
    positive curves.

    Args:
        path: Path to a ``.dos`` or ``.dos.gz`` file.

    Returns:
        A :class:`DosCurve` on a strictly increasing energy grid.
    """
    data = np.loadtxt(Path(path))
    energy, inverse = np.unique(data[:, 0], return_inverse=True)
    total = np.zeros_like(energy)
    occupied = np.zeros_like(energy)
    np.maximum.at(total, inverse, np.abs(data[:, 1]))
    np.maximum.at(occupied, inverse, np.abs(data[:, 2]))
    occupied = np.minimum(occupied, total)
    return DosCurve(energy_ev=energy, occupied=occupied, unoccupied=total - occupied)


@dataclass
class BandSegment:
    """One straight line of a band-structure path.

    Attributes:
        k_start: First k-point in relative coordinates.
        k_end: Last k-point in relative coordinates.
        x: Cumulative path coordinate of every point of the segment.
        energies_ev: Band energies, shape ``(n_k, n_bands)``, absolute eV.
    """

    k_start: np.ndarray
    k_end: np.ndarray
    x: np.ndarray
    energies_ev: np.ndarray


def parse_bands(path: str | Path) -> list[BandSegment]:
    """Read a band file written by ``paw_bands.x`` with several appended line blocks.

    Args:
        path: Path to the ``.dat`` band file.

    Returns:
        List of :class:`BandSegment` in file order. The path coordinate ``x`` is
        cumulative over segments as written by the code.
    """
    segments: list[BandSegment] = []
    k_start = k_end = None
    rows: list[list[float]] = []

    def flush() -> None:
        if rows:
            arr = np.array(rows)
            segments.append(
                BandSegment(
                    k_start=np.array(k_start if k_start is not None else [0, 0, 0]),
                    k_end=np.array(k_end if k_end is not None else [0, 0, 0]),
                    x=arr[:, 0],
                    energies_ev=arr[:, 1:],
                )
            )

    for line in read_text(path).splitlines():
        s = line.strip()
        if not s:
            continue
        if s.startswith("#"):
            if "DATA FOR LINE BLOCK" in s:
                flush()
                rows = []
                k_start = k_end = None
            elif "FIRST K-POINT IN RELATIVE" in s:
                k_start = [float(v) for v in s.split()[-3:]]
            elif "LAST  K-POINT IN RELATIVE" in s or "LAST K-POINT IN RELATIVE" in s:
                k_end = [float(v) for v in s.split()[-3:]]
            continue
        rows.append([float(v) for v in s.split()])
    flush()
    return segments

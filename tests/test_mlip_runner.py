"""Tests for the MLIP stage runner (subprocess target).

Only the argument parser and the constraint-application path are tested here;
they do not require an MLIP package. End-to-end optimization is the user's
live smoke test (needs an MLIP env).
"""

import json
import os
import shutil
import tempfile

from ase import Atoms
from ase.constraints import Hookean

from asetools.workflow import mlip_runner


class TestParser:
    def test_parses_required_and_optional_args(self):
        argv = [
            "--structure", "CONTCAR",
            "--mlip", "uma-s-1p2",
            "--optimizer", "bfgs",
            "--fmax", "0.10",
            "--max-steps", "300",
            "--device", "auto",
            "--name", "OPT_MLIP",
        ]
        args = mlip_runner.build_parser().parse_args(argv)
        assert args.structure == "CONTCAR"
        assert args.mlip == "uma-s-1p2"
        assert args.optimizer == "bfgs"
        assert args.fmax == 0.10
        assert args.max_steps == 300
        assert args.device == "auto"
        assert args.name == "OPT_MLIP"
        assert args.relax_cell is False
        assert args.constraints_json is None

    def test_relax_cell_is_a_flag(self):
        args = mlip_runner.build_parser().parse_args(
            ["--structure", "P", "--mlip", "m", "--name", "N", "--relax-cell"]
        )
        assert args.relax_cell is True

    def test_constraint_args_parsed(self):
        args = mlip_runner.build_parser().parse_args(
            [
                "--structure", "P", "--mlip", "m", "--name", "N",
                "--constraints-json", "pairs.json",
                "--spring-constant", "30.0",
                "--distance-factor", "1.2",
            ]
        )
        assert args.constraints_json == "pairs.json"
        assert args.spring_constant == 30.0
        assert args.distance_factor == 1.2


class TestConstraintApplication:
    def setup_method(self):
        self.dir = tempfile.mkdtemp()
        self.cwd = os.getcwd()
        os.chdir(self.dir)

    def teardown_method(self):
        os.chdir(self.cwd)
        shutil.rmtree(self.dir)

    def _atoms(self):
        return Atoms(
            "H4",
            positions=[[0, 0, 0], [1, 0, 0], [2, 0, 0], [3, 0, 0]],
            cell=[10, 10, 10],
        )

    def test_applies_hookean_from_json(self):
        with open("pairs.json", "w") as f:
            json.dump({"pairs": [[0, 1], [2, 3]]}, f)
        args = mlip_runner.build_parser().parse_args(
            [
                "--structure", "P", "--mlip", "m", "--name", "N",
                "--constraints-json", "pairs.json",
                "--spring-constant", "20.0",
            ]
        )
        atoms = self._atoms()
        mlip_runner.apply_constraints_from_args(atoms, args)
        hookeans = [c for c in atoms.constraints if isinstance(c, Hookean)]
        assert len(hookeans) == 2

    def test_no_constraints_json_leaves_atoms_unconstrained(self):
        args = mlip_runner.build_parser().parse_args(
            ["--structure", "P", "--mlip", "m", "--name", "N"]
        )
        atoms = self._atoms()
        mlip_runner.apply_constraints_from_args(atoms, args)
        assert not atoms.constraints

"""MLIP stage runner — the subprocess target for ``engine: mlip`` stages.

The workflow manager runs in the VASP/asetools job env, which cannot host an
MLIP torch stack (see mlip_platform ADR 0001: per-MLIP envs). An MLIP stage is
therefore executed as a subprocess in a separate MLIP env, invoking this module:

    <mlip-env-python> -m asetools.workflow.mlip_runner \
        --structure CONTCAR --mlip uma-s-1p2 --optimizer bfgs \
        --fmax 0.10 --max-steps 300 --device auto --name OPT_MLIP \
        [--relax-cell] [--uma-task omat] [--mace-head omat_pbe] \
        [--constraints-json hookean_c_pairs.json --spring-constant 20.0 \
         --distance-factor 1.134]

It reads the current structure, re-applies the stage's constraints from the same
JSON the VASP stages use (``asetools`` is a dependency of ``mlip_platform``, so
``ConstraintManager`` is available here), attaches the MLIP calculator, runs the
optimization, and writes a ``CONTCAR`` (plus a ``CONTCAR_{name}`` backup). The
process exits non-zero if the optimization did not converge so the manager
refuses to mark the stage done.

Heavy MLIP imports (``mlip_platform``, torch) are deferred into :func:`main` so
the parser and the constraint helper can be imported and tested without an MLIP
env.
"""

import argparse
import logging
import shutil
import sys

from ase.io import read, write

from asetools.workflow.constraints import ConstraintManager

logger = logging.getLogger(__name__)


def build_parser() -> argparse.ArgumentParser:
    """Build the command-line parser for the MLIP stage runner."""
    p = argparse.ArgumentParser(
        prog="asetools.workflow.mlip_runner",
        description="Run an MLIP geometry optimization for a workflow stage.",
    )
    p.add_argument("--structure", required=True, help="Input structure file to relax.")
    p.add_argument("--name", required=True, help="Stage name (used for the backup suffix).")
    p.add_argument("--mlip", required=True, help="MLIP tag passed to setup_calculator.")
    p.add_argument("--optimizer", default="bfgs", help="ASE optimizer name.")
    p.add_argument("--fmax", type=float, default=0.05, help="Force convergence (eV/Ang).")
    p.add_argument("--max-steps", type=int, default=200, help="Maximum optimizer steps.")
    p.add_argument("--device", default="auto", help="Compute device (auto/cuda/cpu).")
    p.add_argument("--relax-cell", action="store_true", help="Also relax the cell.")
    p.add_argument("--uma-task", default="omat", help="UMA task head (UMA tags only).")
    p.add_argument("--mace-head", default="omat_pbe", help="MACE-MH head (MACE-MH tags only).")
    # Constraint re-application (optional)
    p.add_argument("--constraints-json", default=None, help="Hookean pairs JSON file.")
    p.add_argument(
        "--spring-constant", type=float, default=20.0, help="Hookean spring constant (eV/Ang^2)."
    )
    p.add_argument(
        "--distance-factor", type=float, default=None, help="Override the covalent-radii factor."
    )
    return p


def apply_constraints_from_args(atoms, args) -> None:
    """Re-apply the stage's constraints (if any) from the JSON config.

    Uses the same ``ConstraintManager`` path the VASP stages use, so Hookean
    springs are rebuilt identically inside the MLIP env. No-op if no
    ``--constraints-json`` was given.
    """
    if not args.constraints_json:
        return
    constraint_config = {
        "type": "hookean",
        "config_file": args.constraints_json,
        "spring_constant": args.spring_constant,
    }
    if args.distance_factor is not None:
        constraint_config["distance_factor"] = args.distance_factor
    ConstraintManager().apply_stage_constraints(atoms, constraint_config)


def main(argv=None) -> int:
    """Entry point. Returns 0 on convergence, 1 otherwise."""
    logging.basicConfig(level=logging.INFO, format="%(message)s")
    args = build_parser().parse_args(argv)

    # Deferred heavy imports: only available inside an MLIP env.
    from mlip_platform.cli.utils import setup_calculator
    from mlip_platform.core.optimize import run_optimization

    logger.info(f"MLIP stage '{args.name}': reading {args.structure}")
    atoms = read(args.structure)

    apply_constraints_from_args(atoms, args)

    setup_calculator(
        atoms,
        mlip=args.mlip,
        uma_task=args.uma_task,
        device=args.device,
        mace_head=args.mace_head,
    )

    converged = run_optimization(
        atoms,
        optimizer=args.optimizer,
        fmax=args.fmax,
        max_steps=args.max_steps,
        relax_cell=args.relax_cell,
        output_dir=".",
        model_name=args.mlip,
    )

    # Write the handoff structure as CONTCAR (so the next VASP stage's
    # load_structure picks it up) plus a stage-suffixed backup.
    write("CONTCAR", atoms, format="vasp")
    shutil.copy("CONTCAR", f"CONTCAR_{args.name}")
    logger.info(f"MLIP stage '{args.name}': wrote CONTCAR (converged={converged})")

    return 0 if converged else 1


if __name__ == "__main__":
    sys.exit(main())

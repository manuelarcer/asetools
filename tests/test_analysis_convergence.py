"""Tests for analysis/vasp.py convergence checking using synthetic OUTCAR content."""

import os
import tempfile

from asetools.analysis.vasp import check_outcar_convergence


def _write_outcar(content: str) -> str:
    """Write content to a temp file and return the path."""
    fd, path = tempfile.mkstemp(suffix="_OUTCAR")
    with os.fdopen(fd, "w") as f:
        f.write(content)
    return path


class TestCheckOutcarConvergence:
    def test_missing_file(self):
        converged, _vasp = check_outcar_convergence("/nonexistent/path/OUTCAR")
        assert converged is False
        assert _vasp == ""

    def test_single_point_converged(self):
        content = """\
 vasp.6.3.0 19Jan22 (build Jul 10 2022 14:26:11)
   IBRION =       -1
   NSW    =        0
General timing and accounting informations for this job:
"""
        path = _write_outcar(content)
        try:
            converged, _vasp = check_outcar_convergence(path)
            assert converged is True
            assert _vasp == "vasp6"
        finally:
            os.unlink(path)

    def test_optimization_converged(self):
        content = """\
 vasp.6.4.1
   IBRION =        2
   NSW    =      100
 reached required accuracy - Loss
General timing and accounting informations for this job:
"""
        path = _write_outcar(content)
        try:
            converged, _vasp = check_outcar_convergence(path)
            assert converged is True
        finally:
            os.unlink(path)

    def test_optimization_not_converged(self):
        content = """\
 vasp.5.4.4
   IBRION =        2
   NSW    =      100
General timing and accounting informations for this job:
"""
        path = _write_outcar(content)
        try:
            converged, _vasp = check_outcar_convergence(path)
            assert converged is False
            assert _vasp == "vasp5"
        finally:
            os.unlink(path)

    def test_ase_optimizer_ibrion_neg1_nsw_positive(self):
        """IBRION=-1 with NSW>0 is ASE-controlled optimization."""
        content = """\
 vasp.6.3.0
   IBRION =       -1
   NSW    =       50
General timing and accounting informations for this job:
"""
        path = _write_outcar(content)
        try:
            converged, _vasp = check_outcar_convergence(path)
            assert converged is True
        finally:
            os.unlink(path)

    def test_ibrion1_optimization(self):
        content = """\
 vasp.6.3.0
   IBRION =        1
   NSW    =       50
 reached required accuracy
General timing and accounting informations for this job:
"""
        path = _write_outcar(content)
        try:
            converged, _ = check_outcar_convergence(path)
            assert converged is True
        finally:
            os.unlink(path)

    def test_ibrion3_optimization_not_converged(self):
        content = """\
 vasp.6.3.0
   IBRION =        3
   NSW    =       50
"""
        path = _write_outcar(content)
        try:
            converged, _ = check_outcar_convergence(path)
            assert converged is False
        finally:
            os.unlink(path)

    def test_verbose_output(self, capsys):
        content = """\
 vasp.6.3.0
   IBRION =        2
   NSW    =       50
 reached required accuracy
General timing and accounting informations for this job:
"""
        path = _write_outcar(content)
        try:
            converged, _ = check_outcar_convergence(path, verbose=True)
            assert converged is True
            captured = capsys.readouterr()
            assert "CONVERGED" in captured.out
        finally:
            os.unlink(path)

    def test_empty_file(self):
        path = _write_outcar("")
        try:
            converged, _vasp = check_outcar_convergence(path)
            assert converged is False
            assert _vasp == ""
        finally:
            os.unlink(path)

"""Tests for top-level package imports and lazy loading."""

import pytest


class TestLazyImports:
    """Test that lazy subpackage loading works correctly."""

    def test_version(self):
        import asetools

        assert hasattr(asetools, "__version__")
        assert isinstance(asetools.__version__, str)

    @pytest.mark.parametrize(
        "subpackage",
        [
            "analysis",
            "electronic",
            "structure",
            "electrochemistry",
            "thermodynamics",
            "database",
            "plotting",
            "cli",
            "parsers",
        ],
    )
    def test_subpackage_import(self, subpackage):
        """Each subpackage should be accessible via asetools.<name>."""
        import asetools

        mod = getattr(asetools, subpackage)
        assert mod is not None

    def test_invalid_attr_raises(self):
        import asetools

        with pytest.raises(AttributeError):
            _ = asetools.nonexistent_module

    def test_thermodynamics_reexports(self):
        """Backward-compatible thermodynamics re-exports."""
        import asetools

        # These should be directly importable if they exist
        thermo = asetools.thermodynamics
        assert thermo is not None

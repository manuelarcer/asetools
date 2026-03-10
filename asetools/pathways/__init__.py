"""Reaction pathway tools (NEB, dimer)."""


def __getattr__(name):
    _neb_names = {
        'extract_neb_data', 'plot_nebs', 'redistribute_images_evenly',
        'interpolate_neb_images', 'check_atomic_distances',
        'check_neb_images_sanity', 'setup_neb_calculation',
    }
    _dimer_names = {
        'read_modecar', 'write_modecar', 'generate_displacement_vector',
        'setup_dimer_atoms', 'check_dimer_convergence',
        'extract_saddle_point_info', 'save_dimer_trajectory',
        'validate_dimer_kwargs',
    }
    if name in _neb_names:
        from . import neb
        return getattr(neb, name)
    if name in _dimer_names:
        from . import dimer
        return getattr(dimer, name)
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")


__all__ = [
    'extract_neb_data',
    'plot_nebs',
    'redistribute_images_evenly',
    'interpolate_neb_images',
    'check_atomic_distances',
    'check_neb_images_sanity',
    'setup_neb_calculation',
    'read_modecar',
    'write_modecar',
    'generate_displacement_vector',
    'setup_dimer_atoms',
    'check_dimer_convergence',
    'extract_saddle_point_info',
    'save_dimer_trajectory',
    'validate_dimer_kwargs',
]

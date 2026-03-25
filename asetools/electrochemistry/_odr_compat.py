"""Compatibility layer for ODR fitting.

Prefers the standalone ``odrpack`` package (recommended replacement for
``scipy.odr`` which is deprecated as of SciPy 1.17).  Falls back to
``scipy.odr`` with the deprecation warning suppressed when ``odrpack``
is not installed.
"""

from __future__ import annotations

import warnings
from typing import Any, Callable

import numpy as np

try:
    from odrpack import OdrResult, odr_fit

    _HAS_ODRPACK = True
except ImportError:
    _HAS_ODRPACK = False
    with warnings.catch_warnings():
        warnings.filterwarnings(
            "ignore",
            message=r".*scipy\.odr.*deprecated.*",
            category=DeprecationWarning,
        )
        from scipy import odr as _scipy_odr


class OdrFitResult:
    """Unified result object wrapping either odrpack.OdrResult or scipy.odr.Output.

    Exposes ``.beta`` (fitted parameters) which is the primary attribute
    used by the electrochemistry module.
    """

    def __init__(self, wrapped: Any) -> None:
        self._wrapped = wrapped

    @property
    def beta(self) -> np.ndarray:
        return np.asarray(self._wrapped.beta)

    def __repr__(self) -> str:
        return f"OdrFitResult(beta={self.beta})"


def fit_odr(
    func: Callable[..., np.ndarray],
    x: np.ndarray,
    y: np.ndarray,
    beta0: list[float],
) -> OdrFitResult:
    """Run ODR fit using the best available backend.

    Parameters
    ----------
    func : callable
        Model function with signature ``f(beta, x)`` (scipy convention).
        Internally adapted to ``f(x, beta)`` for odrpack.
    x, y : array-like
        Observed data.
    beta0 : list of float
        Initial parameter guesses.

    Returns
    -------
    OdrFitResult
        Result with ``.beta`` attribute.
    """
    if _HAS_ODRPACK:
        # odrpack uses f(x, beta)
        result = odr_fit(
            lambda xd, beta: func(beta, xd),
            np.asarray(x, dtype=np.float64),
            np.asarray(y, dtype=np.float64),
            np.asarray(beta0, dtype=np.float64),
        )
        return OdrFitResult(result)
    else:
        model = _scipy_odr.Model(func)
        data = _scipy_odr.Data(x, y)
        odr_obj = _scipy_odr.ODR(data, model, beta0=beta0)
        output = odr_obj.run()
        return OdrFitResult(output)


def fit_polynomial_odr(
    x: np.ndarray,
    y: np.ndarray,
    order: int,
) -> OdrFitResult:
    """Fit a standard polynomial via ODR.

    Equivalent to ``scipy.odr.polynomial(order)`` — no fixed constant.

    Parameters
    ----------
    x, y : array-like
        Data to fit.
    order : int
        Polynomial order.

    Returns
    -------
    OdrFitResult
        Result with ``.beta`` containing polynomial coefficients.
    """
    if _HAS_ODRPACK:
        # Build polynomial model: beta[0] + beta[1]*x + ... + beta[order]*x^order
        def poly_model(xd: np.ndarray, beta: np.ndarray) -> np.ndarray:
            result = np.zeros_like(xd)
            for i, b in enumerate(beta):
                result += b * np.power(xd, i)
            return result

        result = odr_fit(
            poly_model,
            np.asarray(x, dtype=np.float64),
            np.asarray(y, dtype=np.float64),
            np.ones(order + 1, dtype=np.float64),
        )
        return OdrFitResult(result)
    else:
        model = _scipy_odr.polynomial(order)
        data = _scipy_odr.Data(x, y)
        odr_obj = _scipy_odr.ODR(data, model)
        output = odr_obj.run()
        return OdrFitResult(output)


def is_odr_result(obj: Any) -> bool:
    """Check if an object is an ODR fit result (any backend)."""
    if isinstance(obj, OdrFitResult):
        return True
    if _HAS_ODRPACK and isinstance(obj, OdrResult):
        return True
    return not _HAS_ODRPACK and isinstance(obj, _scipy_odr.Output)

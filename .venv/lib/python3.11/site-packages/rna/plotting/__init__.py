"""
This module is built to be able to use multiple plotting backends in the same
way.
No matter, which backend you use, plotting should work with the steps shown in
the following
Examples:
    >>> import rna
    >>> import tfields

    Build a cube
    >>> cube = tfields.Mesh3D.grid(
    ...         (-1, 1, 2),
    ...         (-1, 1, 2),
    ...         (-5, -3, 2))
    >>> cube.maps[3].fields = [tfields.Tensors(list(range(len(cube.maps[3]))))]
    >>> map_index = 0
    >>> scalars = cube.maps[3].fields[map_index]

    Tell rna to use the pyqtgraph backend. Default is matplotlib
    >>> rna.plotting.use('pyqtgraph')

    Switching back to default for this test. The code below stays the same if you keep pyqtgraph.
    >>> rna.plotting.use('matplotlib')

    Get an axes element
    >>> axes = rna.plotting.add_subplot(dim=3)

    Plot the cube
    >>> artist = rna.plotting.plot_mesh(axes, cube, cube.maps[3], color=scalars,
    ...     cmap='plasma')

    # Plot a colorbar
    >>> cbar = rna.plotting.set_colorbar(axes, artist)  # TODO: Broken in pyqtgraph atm

    Show the app and start the loop. Commented for doctests
    >>> rna.plotting.show()

"""
import typing
import numpy as np
from . import base
from . import colors
from .core import (
    Api,
    ApiBackend,
    plot_signature,
    Figure,
    PlotKwargsFunctions,
)


# TODO: remove EXPOSED_BACKEND_METHOD in favor of check for ApiBackend member
EXPOSED_BACKEND_METHODS = [
    "axes_dim",
    "add_subplot",
    "gca",
    "is_axes",
    "show",
    "clear",
    "plot_mesh",
    "plot_tensor_map",
    "plot_array",
    "plot_tensor",
    "save",
    "set_colorbar",
    "set_style",
    "set_labels",
    "set_legend",
    "set_aspect_equal",
    "set_limits",
]


def use(backend: typing.Union[str, ApiBackend]):
    """
    Activate a backend of you choice. Non specific methods will be exposed to
    the namespace of rna.plotting . Further methods are exposed to
    rna.plotting.backend
    Args:
        newbackend (str):
            The name of the backend to use.
    Examples:
        >>> import rna
        >>> rna.plotting.use('matplotlib')
        >>> mpl_show = rna.plotting.show
        >>> rna.plotting.use('pyqtgraph')
        >>> pqg_show = rna.plotting.show
        >>> assert callable(mpl_show)
        >>> assert callable(pqg_show)
        >>> assert mpl_show is not pqg_show
        >>> rna.plotting.use('matplotlib')
    """
    Api.STRATEGY_DEFAULT = backend


def __getattr__(name: str) -> typing.Any:
    """
    Forward to generic Api methods

    Note:
        Possibly only on 3.7 +
        Maybe look closer to
        https://stackoverflow.com/questions/2447353/getattr-on-a-module
    """
    if name in EXPOSED_BACKEND_METHODS:
        backend = Api().backend
        try:
            return getattr(backend, name)
        except AttributeError as err:
            raise NotImplementedError(
                "Backend {backend} requires implementation of {name}".format(**locals())
            ) from err
    elif not name.startswith("_") and name in dir(PlotKwargsFunctions):
        return getattr(PlotKwargsFunctions, name)
    else:
        raise AttributeError("%r has no attribute %r" % ("rna.plotting", name))


def figsize(width_pt, scale=1.0, ratio=None):
    """
    Args:
        width_pt (float): figure width width in pt. Can be scaled by
            Get textwidth from latex using '\the\textwidth'
        scale (float)
        ratio (float | None):
            float: ratio of figure height / figure width
            None: golden ratio
    """
    inches_per_pt = 1.0 / 72.27  # Convert pt to inch
    if ratio is None:
        # golden ratio
        ratio = (np.sqrt(5.0) - 1.0) / 2.0
    fig_width = width_pt * inches_per_pt * scale  # fig_width in inches
    fig_height = fig_width * ratio  # height in inches
    fig_size = (fig_width, fig_height)
    return fig_size


use("matplotlib")

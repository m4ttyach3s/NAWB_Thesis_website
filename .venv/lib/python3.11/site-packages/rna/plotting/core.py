"""
Abstraction of plotting with arbitrary plotting apis (backends)
"""
import typing
import warnings
import abc
import pickle
import rna.pattern.backend
import rna
import rna.pattern.decorator


def _get_axes_from_input(backend: "ApiBackend", *data, **kwargs):
    """
    Args:
        backend: backend providing add_subplot, axes_dim and is_axes callable
    """
    # 'ax' is the standard in matplotlib, seaborn etc. Sometimes passed, e.g. seaborn.Pairgrid.map...
    ax = kwargs.pop("ax", None)
    # 'axes' is self explanatory as a variable while ax is too short
    axes = kwargs.pop("axes", None)
    if ax and axes:
        raise ValueError("Axes defined twice (ax and axes)")
    else:
        axes = axes or ax
    dim = kwargs.pop("dim", None)

    if data:
        if backend.is_axes(data[0]):
            if axes is not None and axes is not data[0]:
                raise ValueError("Conflicting axes passed as arg and kwarg")
            axes = data[0]
            data = data[1:]

    if axes and dim:
        if dim != backend.axes_dim(axes):
            raise ValueError(
                "Axes dimension and dim argument in conflict:" f"{backend.dim} vs {dim}"
            )
    if axes is None:
        axes = backend.add_subplot(dim=dim)
    return axes, data, kwargs


# pylint:disable=invalid-name,too-few-public-methods
class plot_signature(rna.pattern.decorator.Decorator):
    """
    This decorator enforces a signature for plot methods and functions.
    The wrapped callable must then satisfy the following signature:

    >>> import rna.plotting
    >>> @rna.plotting.plot_signature
    ... def plot_something(axes, *args, **kwargs):
    ...     artist = axes.plot(*args, **kwargs)
    ...     # usually it is a good idea to return the artist(s) but this is not enforced
    ...     return artist, axes

    And given exemplary inputs

    >>> data = ([1, 2, 3], [4, 5, 6])
    >>> kwargs = {'color': 'r'}

    it can be called lazily (without explicit axes object)

    >>> artist, axes = plot_something(*data, **kwargs)

    where the dimension of the generated axes can be controlled by the `dim` kwarg

    >>> artist, axes_2 = plot_something(*data, dim=2, **kwargs)

    The returned axes are new ones

    >>> assert axes is not axes_2

    Or more pythonic (expicit is better than implicit) in the following ways:

    >>> artist, axes_3 = plot_something(axes, *data, **kwargs)
    >>> artist, axes_4 = plot_something(*data, axes=axes, **kwargs)
    >>> artist, axes_5 = plot_something(*data, ax=axes, **kwargs)

    The returned axes are the ones given

    >>> assert axes is axes_3 is axes_4 is axes_5

    """

    def _wrap(self, this, func, *args, **kwargs):
        if isinstance(this, ApiBackend):
            backend = this
        else:
            backend = Api().backend
        axes, data, kwargs = _get_axes_from_input(backend, *args, **kwargs)
        return func(axes, *data, **kwargs)


class PlotKwargsFunctions:
    """
    Collection of staticmethod functions acting on kwargs passed to a plotting method.
    This provides a standardized signature of plotting functions.
    """

    @staticmethod
    def retrieve(
        kwargs: dict,
        *attrs: str,
        default: typing.Any = None,
        keep: typing.Union[bool, typing.List[str]] = True,
        pop: typing.Optional[typing.List[str]] = None,
        deprecated: typing.Optional[typing.List[str]] = None,
    ):
        """
        Try to find a value in kwargs. Loop through attrs keys and return the first valid
        occurance.
        This can be used e.g. to allow old style arguments where another attribute name should be
        established.

        Args:
            kwargs:
            *attrs: attributes to use as keys for lookup
            default: Default value to return if no value was found. This is not set.
            keep: If False pop the keys from kwargs. If True, keep all.
                If list of attrs -> keep these attributes, remove all others.
            pop: List of attrs to pop. Disjunct with keep so keep is set to False automatically.
            depcreated:

        Examples:
            >>> import rna

            In the simplest case it behaves like get with None as default
            >>> kwargs = dict(a=3)
            >>> rna.plotting.retrieve(kwargs, 'b') is None
            True
            >>> kwargs
            {'a': 3}

            In the second to simplest case it behaves like pop with None as default
            >>> rna.plotting.retrieve(kwargs, 'a', keep=False)
            3
            >>> kwargs
            {}

            Multiple lookups are allowed in a chain
            >>> kwargs = {}
            >>> rna.plotting.retrieve(kwargs, 'edgecolor', 'c', 'color', default='k')
            'k'

            >>> kwargs = dict(color='b', edgecolor='r')
            >>> rna.plotting.retrieve(kwargs, 'color', 'c')
            'b'
            >>> rna.plotting.retrieve(kwargs, 'c', 'color', deprecated=['c'])
            'b'
            >>> rna.plotting.retrieve(kwargs, 'edgecolor', 'c', 'color')
            'r'
            >>> assert 'color' in kwargs
            >>> assert 'edgecolor' in kwargs

            Pop values if no longer required in kwargs:
            >>> rna.plotting.retrieve(kwargs, 'edgecolor', 'c', 'color', keep=False)
            'r'
            >>> assert len(kwargs) == 0

            Pop values except those in the list keep if no longer required in kwargs:
            >>> kwargs = dict(color='b', edgecolor='r', facecolor='g')
            >>> rna.plotting.retrieve(kwargs, 'edgecolor', 'c', 'color', pop=['edgecolor'])
            'r'
            >>> kwargs
            {'color': 'b', 'facecolor': 'g'}
            >>> assert 'color' in kwargs
            >>> assert len(kwargs) == 2

            >>> rna.plotting.retrieve(kwargs, 'facecolor', 'color', keep=False)
            'g'
            >>> assert len(kwargs) == 0

        """
        if pop and (isinstance(keep, list) or not keep):
            raise ValueError("pop and keep arguments are exclusory.")

        found = False
        for attr in attrs:
            if attr not in kwargs:
                continue

            # attr present. pop or get?
            if pop and attr in pop:
                tmp = kwargs.pop(attr)
            elif not keep or (isinstance(keep, list) and attr not in keep):
                tmp = kwargs.pop(attr)
            else:
                tmp = kwargs.get(attr)

            if not found:
                # first occurance
                value = tmp

            found = True

            if not (deprecated or pop or not keep or isinstance(keep, list)):
                # no need to continue
                break

            if deprecated and attr in deprecated:
                warnings.warn(
                    "Attribute {attrs[0]} is deprecated.".format(**locals()),
                    DeprecationWarning,
                )

        if found:
            return value

        return default

    @staticmethod
    def get_norm_args(kwargs, vmin_default=0, vmax_default=1, cmap_default=None):
        """
        Examples:
            >>> import rna
            >>> from rna.plotting.backends.matplotlib import ApiBackendMatplotlib

            >>> rna.plotting.set_style()
            >>> ApiBackendMatplotlib.get_norm_args(dict(vmin=2, vmax=42))
            ('viridis', 2, 42)
        """
        if cmap_default is None:
            import matplotlib.pyplot as plt  # pylint: disable = import-outside-toplevel

            cmap_default = plt.rcParams["image.cmap"]
        cmap = kwargs.get("cmap", cmap_default)
        vmin = kwargs.get("vmin", vmin_default)
        vmax = kwargs.get("vmax", vmax_default)
        if vmin is None:
            vmin = vmin_default
        if vmax is None:
            vmax = vmax_default
        return cmap, vmin, vmax

    @staticmethod
    def pop_norm_args(kwargs, **defaults):
        """
        Pop vmin, vmax and cmap from plot_kwargs
        Args:
            **defaults:
                see get_norm_args method
        """
        cmap, vmin, vmax = ApiBackend.get_norm_args(kwargs, **defaults)
        kwargs.pop("cmap", None)
        kwargs.pop("vmin", None)
        kwargs.pop("vmax", None)
        return cmap, vmin, vmax

    @staticmethod
    def pop_xyz_index(kwargs: dict):
        """
        Args:
            kwargs: contains optional integer keys x_index, y_index and z_index
                Thes indices refer to the index of some data object which holds data of x, y and
                z axis respectively
        """
        x_index = kwargs.pop("x_index", 0)
        y_index = kwargs.pop("y_index", 1)
        z_index = kwargs.pop("z_index", None)
        if z_index is None:
            z_index = {1, 2, 3}
            z_index.difference_update()
            z_index = z_index.pop()
        return x_index, y_index, z_index

    @staticmethod
    def format_colors(kwargs, colors, fmt="rgba", length=None, dtype=None):
        """
        format colors according to fmt argument
        Args:
            kwargs: dict, potentially containing the color argument
            colors (list/one value of rgba tuples/int/float/str): This argument
                will be interpreted as color
            fmt (str): rgba | hex | norm
            length (int/None): if not None: correct colors lenght
            dtype (np.dtype): output data type

        Returns:
            colors in fmt
        """
        cmap, vmin, vmax = ApiBackend.get_norm_args(
            kwargs, cmap_default=None, vmin_default=None, vmax_default=None
        )
        from rna.plotting.colors import to_colors

        return to_colors(
            colors, fmt, length=length, vmin=vmin, vmax=vmax, cmap=cmap, dtype=dtype
        )


class ApiBackend(
    rna.pattern.backend.Backend,
    rna.polymorphism.Plottable,
    rna.polymorphism.Storable,
    PlotKwargsFunctions,
):
    """
    Base class for a backend implementation for a Api subclass

    The plotting api backend is thought of as unified interface for plot arguments taking care of
        * interdependent arguments (axes/dimension, cmap/vmin/vmax, x_index,...)
        * color formating (format_colors)

    processing kwargs for plotting functions and providing easy
    access to axes, dimension and plotting method as well as indices
    for array choice (x..., y..., z_index)

    Examples:
        >>> import rna
        >>> from rna.plotting.backends.matplotlib import ApiBackendMatplotlib

        Dimension arguemnt (dim)
        >>> api = rna.plotting.Api()
        >>> api.backend = 'matplotlib'  # which is default anyway

        As a functionality of the Backend/Strategy pattern the backend is instantiated at
        access time
        >>> backend = api.backend
        >>> assert not isinstance(backend, str)
        >>> assert isinstance(backend, ApiBackendMatplotlib)

        Default dimension of the axis is 2
        >>> backend.axes_dim(backend.add_subplot())
        2

        You can select the correct axes object with the dimension key
        >>> axes = ApiBackendMatplotlib.add_subplot(dim=3)
        >>> backend.axes_dim(axes)
        3

        Switching the backend keeps the top level behaviour the same
        >>> api = rna.plotting.Api()
        >>> api.backend = 'pyqtgraph'
        >>> api.backend.axes_dim(api.backend.add_subplot(dim=2))
        2

    """

    STRATEGY_MODULE_BASE = "rna.plotting.backends"

    @plot_signature
    def plot_array(self, *args, **kwargs):
        warnings.warn("Use plot_tensor instead.", DeprecationWarning)
        return self.plot_tensor(*args, **kwargs)  # pylint: disable=no-member

    @staticmethod
    @abc.abstractmethod
    def axes_dim(axes):
        """
        Retrieve the current axes object which resembles the canvas you plot on to

        Args:
            dim: dimension of the axes object
        """

    @staticmethod
    @abc.abstractmethod
    def gca(dim: typing.Optional[int] = None, **kwargs):
        """
        Retrieve the current axes object which resembles the canvas you plot on to

        Args:
            dim: dimension of the axes object
        """

    @staticmethod
    @abc.abstractmethod
    def is_axes(obj: typing.Any) -> bool:
        """
        Check if the given object is an axes instance of this backend
        """

    @abc.abstractmethod
    def show(self):
        """
        Render the current api and promt it to display.
        """


class Api(
    rna.pattern.backend.Frontend, rna.polymorphism.Plottable, rna.polymorphism.Storable
):
    """
    Base class of a generic api.

    Examples:
        >>> import numpy as np
        >>> from rna.plotting import Api
        >>> api = Api()
        >>> api.backend  # doctest +ELLIPSIS
        <rna.plotting.backends.matplotlib.ApiBackendMatplotlib object at ...>

        To use plot, implement _plot
        >>> class Demo(Api):
        ...     def _plot(self, **kwargs):
        ...         self.set_style()
        ...         labels = kwargs.pop("labels", ["x (m)", "y (m)", "z (m)"])
        ...         ax = self.add_subplot()
        ...         artist = self.plot_tensor(ax, self.data, **kwargs)
        ...         self.set_labels(ax, *labels)
        ...         return artist
        >>> demo = Demo()

        >>> artists = demo.plot(np.array([[1,1], [2,2]]), color='r')
        >>> _ = demo.set_legend(artists)

        # saving with a data format stores the data and allows recovery
        >>> path = "./tmp.pkl"
        >>> demo.save(path)
        >>> demo_restored = demo.load(path)

        >>> artists = demo_restored.plot()

        # rendering the api
        >>> demo_restored.show()

        # saving with an api format stores the Api
        >>> path_demo = "./tmp.png"
        >>> demo.save(path_demo)

        The backend can be easily exchanged:
        >>> demo.backend="plotly"
    """

    STRATEGY_TYPE = ApiBackend
    STRATEGY_DEFAULT = "matplotlib"

    def __init__(self, data=None):
        # Make sure you pass nothing to backend
        super().__init__()
        self.data = data

    def plot(self, data=None, **kwargs):
        """
        Args:
            data: data object. If you loaded this Api from data file, this is not required
                data is set to self.data which is then available in _plot
            **kwargs: passed to _plot
        """
        if data is None and self.data is None:
            # "No data given (first argument)")
            pass
        elif data is not None:
            self.data = data
        return self._plot(**kwargs)

    def _plot(self, **kwargs):
        raise NotImplementedError(
            "{self} requires implementation of _plot(**kwargs)".format(**locals())
        )

    def save(self, *args, **kwargs):
        """
        Data is saved from Api class directly.
        Plots are specific to the backend so _save_<ext> must be implemented there. Also they have
        no load method.
        Data has priority if same format is present in Api and backend. Use Api.backend.save
        if you want to save backend also/only in this format.
        """
        try:
            # data saving
            if self.data is None:
                raise ValueError(
                    "{self}.data is None! You forgot to set the data attribute in plotting."
                )
            super().save(*args, **kwargs)  # pylint: disable=no-member
        except NotImplementedError:
            # figure saving
            self.backend.save(*args, **kwargs)

    def _save_pkl(self, path, **kwargs):
        kwargs.setdefault("protocol", pickle.HIGHEST_PROTOCOL)
        with open(path, "wb") as handle:
            pickle.dump(self.data, handle, **kwargs)

    @classmethod
    def _load_pkl(cls, path, **kwargs):
        with open(path, "rb") as handle:
            data = pickle.load(handle, **kwargs)
        obj = cls()
        obj.data = data
        return obj

    def __getattr__(self, name: str) -> typing.Any:
        """
        Forward to generic Api methods either defined by EXPOSED_BACKEND_METHODS or part of
        ApiBackend.

        Note:
            Possibly only on 3.7 +
            Maybe look closer to
            https://stackoverflow.com/questions/2447353/getattr-on-a-module
        """
        if name in rna.plotting.EXPOSED_BACKEND_METHODS:
            backend = (
                self.backend
            )  # CAREFUL: could break the clever logic of the strategy pattern
            try:
                return getattr(backend, name)
            except AttributeError as err:
                raise NotImplementedError(
                    "Backend {backend} requires implementation of {name}".format(
                        **locals()
                    )
                ) from err
        if hasattr(ApiBackend, name):
            return getattr(self.backend, name)
        return self.__getattribute__(name)  # raise AttributeError


class Figure(Api):
    """
    The Figure backend should be derived from
    """

    pass

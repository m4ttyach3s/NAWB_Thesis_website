"""
Color conversion
"""
from itertools import cycle
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt


FMT_DTYPES = {
    # 'hex': np.dtype([('hex', '<U7')]),
    "rgba": np.dtype([("r", "<f8"), ("g", "<f8"), ("b", "<f8"), ("a", "<f8")]),
    "rgb": np.dtype([("r", "<f8"), ("g", "<f8"), ("b", "<f8")]),
}


def color_fmt(colors):
    """
    Determine the color representation type
    Definition confirms to https://matplotlib.org/api/colors_api.html

    Args:
        colors (various)

    Returns:
        str | None: color representation type
            'str' -> predefined colors e.g. 'r', 'g', 'k', ...
            'gray' -> grayscale value e.g. '0.5'
            'norm' -> arbitrary numbers e.g. [1,2,3.42]
            'rgba' -> tuples of (red, green, blue, transparency)
            'rgb' -> tuples or (red, green, blue)
            None -> fmt unknown

    """
    if isinstance(colors, np.ndarray):
        for fmt, dtype in FMT_DTYPES.items():
            if dtype == colors.dtype:
                return fmt

    if not hasattr(colors, "__iter__") or isinstance(colors, (str, tuple)):  # NOQA
        # colors is just one element
        colors = [colors]
    if len(colors) == 0:
        # empty colors
        return None
    color = colors[0]
    norm_types = (int, float, np.float32, np.float64, np.integer)

    if isinstance(color, str):
        try:
            col = eval(color)
        except Exception:
            pass
        else:
            if isinstance(col, norm_types):
                return "gray"
        if color.startswith("#"):
            return "hex"
        return "str"
    if isinstance(color, tuple):
        if len(color) == len(FMT_DTYPES["rgba"]):
            return "rgba"
        if len(color) == len(FMT_DTYPES["rgb"]):
            return "rgb"
    if isinstance(color, norm_types):
        return "norm"
    return None


def truncate_colormap(cmap, vmin=0.0, vmax=1.0, n=100):
    """
    Truncate colormap cmap at vmin and vmax.
    """
    cmap = plt.get_cmap(cmap)
    new_cmap = mpl.colors.LinearSegmentedColormap.from_list(
        "trunc({n},{a:.2f},{b:.2f})".format(n=cmap.name, a=vmin, b=vmax),
        cmap(np.linspace(vmin, vmax, n)),
    )
    return new_cmap


def to_colors(
    colors, fmt="rgba", length=None, vmin=None, vmax=None, cmap=None, dtype=None
):
    """
    format colors according to fmt argument
    Args:
        colors (List(tuple|int|float|str)): rgb|rgba|string representation|hex
            interpreted as color
        fmt (str): one of rgba, rgb, hex or norm
        length (int/None): if not None: correct colors lenght
        vmin (float)
        vmax (float)
        cmap (str|mpl.cm.ColorMap)
        dtype (np.dtype): output data type

    Returns:
        colors in fmt

    Examples:
        >>> import rna
        >>> import numpy as np

        >>> rna.plotting.colors.to_colors(1, cmap='gray')
        array((1., 1., 1., 1.),
              dtype=[('r', '<f8'), ('g', '<f8'), ('b', '<f8'), ('a', '<f8')])
        >>> rna.plotting.colors.to_colors(1, fmt='rgb', length=2, cmap='gray')
        array([(1., 1., 1.), (1., 1., 1.)],
              dtype=[('r', '<f8'), ('g', '<f8'), ('b', '<f8')])
        >>> rna.plotting.colors.to_colors([1], cmap='gray', vmin=0)
        array([(1., 1., 1., 1.)],
              dtype=[('r', '<f8'), ('g', '<f8'), ('b', '<f8'), ('a', '<f8')])
        >>> rna.plotting.colors.to_colors('k', cmap='gray', vmin=0, vmax=0)
        array((0., 0., 0., 1.),
              dtype=[('r', '<f8'), ('g', '<f8'), ('b', '<f8'), ('a', '<f8')])
        >>> rna.plotting.colors.to_colors(
        ...     ['xkcd:red'], fmt='rgb', length=3, cmap='gray', vmin=0, vmax=0
        ... )  # doctest: +ELLIPSIS
        array([(0.89803922, 0., 0.), (0.898..., 0., 0.), (0.898..., 0., 0.)],
              dtype=[('r', '<f8'), ('g', '<f8'), ('b', '<f8')])
        >>> rna.plotting.colors.to_colors(1, cmap='gray', length=2, fmt='hex')
        array(['#ffffff', '#ffffff'], dtype='<U7')

        Keeps the shape
        >>> rna.plotting.colors.to_colors(np.array([[[1,1],[1,1]]]),
        ...     cmap='gray', vmin=0)
        array([[[(1., 1., 1., 1.), (1., 1., 1., 1.)],
                [(1., 1., 1., 1.), (1., 1., 1., 1.)]]],
              dtype=[('r', '<f8'), ('g', '<f8'), ('b', '<f8'), ('a', '<f8')])

    """
    # mapping also scalars to lists
    has_iter = True
    if not hasattr(colors, "__iter__") or isinstance(colors, (str, tuple)):
        # colors is just one element
        has_iter = False
        colors = [colors]

    # work on a flat color array and reshape afterwards
    colors = np.array(colors)
    shape = list(colors.shape)
    colors = np.array(colors.flat)

    # Assume all colors to be of the same fmt
    orig_fmt = color_fmt(colors)

    # check cmap, vmin, vmax for norm conversions
    if orig_fmt == "norm" or fmt == "norm":
        if orig_fmt == "norm":
            if cmap is None:
                raise TypeError("converting from norm always requires cmap")
            if vmin is None:
                vmin = min(colors) if len(set(colors)) > 1 else min(0, min(colors))
            if vmax is None:
                vmax = max(colors) if len(set(colors)) > 1 else max(1, max(colors))

            # already cut initial norm values
            colors[colors < vmin] = vmin
            colors[colors > vmax] = vmax
        else:
            if vmin is None or vmax is None or cmap is None:
                raise TypeError(
                    "converting to norm always requires cmap, " "vmin and vmax."
                )

    # conversion
    if fmt == orig_fmt:
        # already correct
        pass
    elif fmt == "norm":
        # route everything through rgba_to_norm
        colors_rgba = to_colors(colors, "rgba")
        colors = rgba_to_norm(colors_rgba, cmap, vmin, vmax)
    elif fmt == "rgb":
        colors_rgba = to_colors(colors, "rgba", vmin=vmin, vmax=vmax, cmap=cmap)
        colors = colors_rgba[["r", "g", "b"]]
    elif fmt == "rgba":
        if orig_fmt in ("str", "hex"):
            colors = [mpl.colors.to_rgba(color) for color in colors]
        elif orig_fmt == "rgb":
            colors = [tuple(c) + (1,) for c in colors]
        elif orig_fmt == "norm":
            colors = norm_to_rgba(colors, vmin=vmin, vmax=vmax, cmap=cmap)
        else:
            raise NotImplementedError(
                f"Conversion from '{orig_fmt}' to rgba not implemented."
            )
    elif fmt == "hex":
        colors_rgba = to_colors(colors, "rgba", vmin=vmin, vmax=vmax, cmap=cmap)
        colors = [mpl.colors.to_hex(tuple(c)) for c in colors_rgba]
    else:
        raise NotImplementedError(
            f"Color fmt '{fmt}' not implemented."
        )

    # set correct length
    if length is not None:
        # just one colors value given
        if len(colors) != length:
            if len(colors) != 1:
                raise ValueError("Can not correct color length")
            colors = [colors[0] for i in range(length)]
            # colors = np.repeat(colors, length, axis=0)  # this does not cast
            # to structured array properly
            shape[0] = length
    elif not has_iter:
        shape = shape[1:]
        colors = colors[0]

    # add fmt info to array dtype
    if fmt in FMT_DTYPES:
        if dtype is None:
            dtype = FMT_DTYPES[fmt]
        else:
            shape.append(len(FMT_DTYPES[fmt]))
    elif hasattr(colors, "dtype"):
        if dtype is None:
            dtype = colors.dtype
        else:
            colors = np.array([c.tolist() for c in colors])

    colors = np.array(colors, dtype=dtype)

    shape = tuple(shape)
    colors = colors.reshape(shape)

    return colors


def norm_to_rgba(scalars, cmap=None, vmin=None, vmax=None):
    """
    retrieve the rgba colors for a list of scalars

    Returns:
        List(tuple)

    Examples:
        >>> import rna
        >>> colors = rna.plotting.colors.norm_to_rgba([0, 1, 1.5, 2, 3, 4],
        ...                                           vmin=1, vmax=2)
        >>> colors = np.array(colors)
        >>> colors = colors[:, :3]  # strip a for test
        >>> assert all(colors[0] == colors[1])
        >>> assert all(colors[1] != colors[2])
        >>> assert all(colors[2] != colors[3])
        >>> assert all(colors[3] == colors[4])
        >>> assert all(colors[4] == colors[5])

    """
    color_map = plt.get_cmap(cmap)
    norm = mpl.colors.Normalize(vmin=vmin, vmax=vmax)
    colors = color_map([norm(s) for s in scalars])
    colors = [tuple(c) for c in colors]
    return colors


def rgba_to_norm(
    colors, cmap, vmin, vmax, max_distance=1e-6, normalizer=mpl.colors.Normalize
):
    """
    Inverse 'norm_to_rgba'
    Reconstruct the numeric values (0 - 1) of given
    Args:
        colors (list or rgba tuple)
        cmap (matplotlib colormap)
        vmin (float)
        vmax (float)
        max_distance (float): raise warning if the distance in color space is
            larger than max_distance
    """
    # any to np.ndarray
    colors = np.array(colors)
    # above does not convert struct dtypes. -> struct to float array
    colors = colors.view(np.float64).reshape(colors.shape + (-1,))

    rnge = np.linspace(vmin, vmax, 256)
    norm = normalizer(vmin, vmax)
    if isinstance(cmap, str):
        cmap = mpl.cm.get_cmap(cmap)
    mapvals = cmap(norm(rnge))[:, :4]  # there are 4 channels: r,g,b,a
    scalars = []
    for color in colors:
        distance = np.sum((mapvals - color) ** 2, axis=1)
        min_distance = np.min(distance)
        if min_distance > max_distance:
            import warnings

            warnings.warn(
                f"Distance in color space is large ({min_distance:.6f})"
            )
        scalars.append(rnge[np.argmin(distance)])
    scalars = np.array(scalars)
    return scalars


def colormap(seq):
    """
    Args:
        seq (iterable): a sequence of floats and RGB-tuples.
            The floats should be increasing and in the interval (0,1).
    Returns:
        LinearSegmentedColormap
    """
    seq = [(None,) * 3, 0.0] + list(seq) + [1.0, (None,) * 3]
    cdict = {"red": [], "green": [], "blue": []}
    for i, item in enumerate(seq):
        if isinstance(item, float):
            r1, g1, b1 = seq[i - 1]
            r2, g2, b2 = seq[i + 1]
            cdict["red"].append([item, r1, r2])
            cdict["green"].append([item, g1, g2])
            cdict["blue"].append([item, b1, b2])
    return mpl.colors.LinearSegmentedColormap("CustomMap", cdict)


def color_cycle(cmap=None, n=None):
    """
    Args:
        cmap (matplotlib colormap): e.g. plt.cm.coolwarm
        n (int): number of colors. Only used for cmap argument

    Examples:
        >>> color_cycle('jet', 3)  # doctest: +ELLIPSIS
        <itertools.cycle object at 0x...>
    """
    if cmap:
        if n is None:
            raise AttributeError("color_cycle with cmap requires n")
        if n == 0:
            # At least give one color (easier to debug)
            n = 1
        colors = to_colors(np.linspace(0, 1, n), cmap=cmap, vmin=0, vmax=1, fmt="hex")
    else:
        colors = list([color["color"] for color in mpl.rcParams["axes.prop_cycle"]])
    return cycle(colors)


if __name__ == "__main__":
    import doctest

    doctest.testmod()

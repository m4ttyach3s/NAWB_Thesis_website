"""
Matplotlib implementation of an interchangeable backend for plotting
"""
import typing
import logging
import os
import contextlib
import re
import functools
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.animation as manimation
import matplotlib.collections as mcollections
import matplotlib.dates as dates
import matplotlib.patches as mplpatch
import more_itertools
import mpl_toolkits.mplot3d as plt3d
import mpl_toolkits.axes_grid1 as pltg1
import rna
import numpy as np
import warnings


# ANIMATION

try:
    FFMPEGWRITER = manimation.writers["ffmpeg"]
    FFMPEG = True
except Exception as err:
    logging.getLogger().warning(str(err))
    FFMPEGWRITER = object
    FFMPEG = False


class Writer(FFMPEGWRITER):
    """
    FFMPEGWRITER child

    Args:
        **kwargs:
            title (str)
            comment (str)
            fps (int)

    Examples:
        >>> import numpy as np
        >>> import rna
        >>> from rna.plotting.backends.matplotlib import Writer
        >>> axes = rna.plotting.add_subplot()
        >>> writer = Writer()
        >>> pictures = np.array([[[1,1], [0,0]],
        ...                      [[0.5,1], [0,0.5]]])
        >>> with writer.saving():
        ...     for s in pictures:
        ...         artist = axes.imshow(s.T)
        ...         writer.grab_frame()
        ...         axes.clear()
    """

    def __init__(self, **kwargs):
        if FFMPEG:
            # set Defaults
            title = kwargs.pop("title", "")
            comment = kwargs.pop("comment", "")
            artist = kwargs.pop("artist", "rna plotting module")
            fps = kwargs.pop("fps", 15)
            if len(kwargs) > 0:
                raise ValueError("Some Attributes are not used.")

            # set up writer
            metadata = dict(title=title, artist=artist, comment=comment)
            super(Writer, self).__init__(fps=fps, metadata=metadata)

    @contextlib.contextmanager
    def saving(self, fig=None, path=None, dpi=300, **kwargs):
        """
        Overwritinig to FFMPEGWRITER.saving
        setting defaults though
        """
        log = logging.getLogger()
        if FFMPEG:
            if fig is None:
                fig = plt.gcf()
            if path is None:
                path = "/tmp/lastMovie.mp4"
            else:
                path = rna.path.resolve(path)
            log.info("Write movie to %s", path)
            with super(Writer, self).saving(fig, path, dpi, **kwargs):
                yield
        else:
            log.warning("FFMPEG Writer not available. Nothing is written!")


class set_zoomable(object):
    """
    Left click a colorbar and release in order to zoom.
    Upper and lower 5% of the colorbar will zoom out.
    """

    def __init__(self, cbar):
        self._cbar = cbar
        self._v_press = None
        self._v_release = None

        artist = self._cbar.mappable
        self._press_connection_id = artist.axes.figure.canvas.mpl_connect(
            "button_press_event", functools.partial(self.on_press)
        )
        self._release_connection_id = artist.axes.figure.canvas.mpl_connect(
            "button_release_event", functools.partial(self.on_release)
        )

    def on_press(self, event):
        if event.inaxes is self._cbar.ax:
            self._v_press = event.ydata

    def on_release(self, event):
        if event.inaxes is self._cbar.ax:
            self._v_release = event.ydata
            self.update_v_min_max()

    def update_v_min_max(self):
        # sort press and release event
        if self._v_press > self._v_release:
            x_up = self._v_press
            x_low = self._v_release
        if self._v_press < self._v_release:
            x_up = self._v_release
            x_low = self._v_press

        # zoom out if in 5% margin
        if x_up > 0.95:
            x_up = 1.5
        if x_low < 0.05:
            x_low = -0.5

        artist = self._cbar.mappable
        vmin, vmax = artist.get_clim()
        v_range = vmax - vmin
        vmax = vmin + x_up * v_range
        vmin = vmin + x_low * v_range
        artist.set_clim(vmin, vmax)
        artist.axes.figure.canvas.draw()


# Backend


class ApiBackendMatplotlib(rna.plotting.ApiBackend):
    """
    # PlotManager "method" not required any more
    """

    EXPLICIT_SAVE_METHODS = {
        ext: "_plt_savefig"
        for ext in plt.gcf().canvas.get_supported_filetypes()
        if ext not in ("pgf")
    }

    def show(self):
        plt.show()

    @staticmethod
    def axes_dim(axes):
        """
        Returns int: axes dimension
        """
        if axes.name == "3d":
            return 3
        return 2

    @staticmethod
    def add_subplot(dim=None, **kwargs):
        if dim is not None and "projection" not in kwargs:
            if dim == 2:
                kwargs["projection"] = "rectilinear"
            elif dim == 3:
                kwargs["projection"] = "3d"
            else:
                raise NotImplementedError(
                    "Dimension {dim} not implemented".format(**locals())
                )
        return plt.gcf().add_subplot(**kwargs)

    @staticmethod
    def gca(dim=None, **kwargs):
        """
        Forwarding to plt.add_subplot but translating the dimension to projection
        correct dimension

        """
        # this is deprecated, use axes instead
        warnings.warn("gca is deprecated, use add_suplot instead")
        return ApiBackendMatplotlib.add_subplot(dim=dim, **kwargs)
        # breakpoint()
        # if dim is not None and "projection" not in kwargs:
        #     if dim == 2:

        #         kwargs["projection"] = "rectilinear"
        #         # the line above would already be good ...
        #         # but there is a bug in matplotlib which ignores setting the default
        #         axes = plt.gca()
        #         if ApiBackendMatplotlib.axes_dim(axes) != dim:
        #             return plt.gcf().add_subplot(**kwargs)
        #         else:
        #             return axes
        #     elif dim == 3:
        #         kwargs["projection"] = "3d"
        #     else:
        #         raise NotImplementedError(
        #             "Dimension {dim} not implemented".format(**locals())
        #         )
        # axes = plt.gca()
        # return axes

    @staticmethod
    def is_axes(obj):
        return isinstance(obj, plt.Axes)

    def set_labels(self, axes, *labels):
        """
        Set x, y and possibly z labels
        """
        axes.set_xlabel(labels[0])
        axes.set_ylabel(labels[1])
        if self.axes_dim(axes) == 3:
            axes.set_zlabel(labels[2])
        elif self.axes_dim(axes) != 2:
            raise ValueError(
                "Axis dimensions other than 2 or 3 is not supported for set_labels."
            )

    @rna.plotting.plot_signature
    def set_colorbar(
        self,
        axes,
        artist,
        label=None,
        divide=True,
        position="right",
        size="2%",
        pad=0.05,
        labelpad=None,
        zoom=False,
        **kwargs
    ):
        """
        Note:
            Bug found in matplotlib:
                when calling axes.clear(), the colorbar has to be removed by hand,
                since it will not be removed by clear.
            >> import rna
            >> axes = rna.plotting.add_subplot()
            >> artist = ...
            >> cbar = rna.plotting.set_colorbar(axes, artist)
            >> rna.plotting.save(...)
            >> cbar.remove()  # THIS IS IMPORTANT. Otherwise you will have problems
            # at the next call of savefig.
            >> axes.clear()
        """
        ticks_position = "default"
        label_position = "bottom"
        labelpad = 30 if labelpad is None else labelpad
        if position == "right":
            rotation = 270
        elif position == "left":
            rotation = 90
        elif position == "top":
            rotation = 0
            ticks_position = "top"
            label_position = "top"
            labelpad = 5
        elif position == "bottom":
            rotation = 0
        # colorbar
        if divide and self.axes_dim(axes) == 2:
            divider = pltg1.make_axes_locatable(axes)
            axes = divider.append_axes(position, size=size, pad=pad)
            cbar = plt.colorbar(artist, cax=axes, **kwargs)
        else:
            cbar = plt.colorbar(artist, ax=axes, **kwargs)
        cbar.ax.xaxis.set_ticks_position(ticks_position)
        cbar.ax.xaxis.set_label_position(label_position)
        cbar.ax.tick_params(axis="x", which="major", pad=0)

        # label
        if label is None:
            art_label = artist.get_label()
            if art_label:
                label = art_label
        if label is not None:
            cbar.set_label(label, rotation=rotation, labelpad=labelpad)

        if zoom:
            set_zoomable(cbar)

        return cbar

    @staticmethod
    def clear(axes):
        axes.clear()

    @staticmethod
    def upgrade_style(style, source, dest=None):
        """
        Copy a style file at <origionalFilePath> to the <dest> which is the
        foreseen local matplotlib rc dir by default
        The style will be name <style>.mplstyle
        Args:
            style (str): name of style
            source (str): full path to mplstyle file to use
            dest (str): local directory to copy the file to. Matpotlib has to
                search this directory for mplstyle files!
        Examples:
            >>> import rna
            >>> import os
            >>> rna.plotting.use('matplotlib')
            >>> from rna.plotting.backends.matplotlib import ApiBackendMatplotlib
            >>> ApiBackendMatplotlib.upgrade_style(
            ...     'rna',
            ...     os.path.join(os.path.dirname(rna.plotting.__file__),
            ...                  'mplstyles/rna.mplstyle'))

        """
        if dest is None:
            dest = mpl.get_configdir()
        style_extension = "mplstyle"
        path = rna.path.resolve(os.path.join(dest, style + "." + style_extension))
        source = rna.path.resolve(source)
        rna.path.cp(source, path)

    def set_style(
        self,
        style: str = "rna",
        dest: typing.Optional[str] = None,
        upgrade: bool = False,
    ):
        """
        Set the matplotlib style

        Args:
            style: name of style
            dest: mplstydirectory that stores the mplstyle. if None, use default
                maplotlib destination
            upgrade: if True copy the style to the source destination (dest)
        """
        if dest is None:
            dest = mpl.get_configdir()

        style_extension = "mplstyle"
        path = rna.path.resolve(os.path.join(dest, style + "." + style_extension))
        if style in mpl.style.available and not upgrade:
            plt.style.use(style)
        elif os.path.exists(path) and not upgrade:
            plt.style.use(path)
        else:
            log = logging.getLogger()
            source = rna.path.resolve(
                os.path.join(
                    os.path.join(os.path.dirname(__file__), "mplstyles"),
                    style + "." + style_extension,
                )
            )
            if os.path.exists(source):
                log.warning("I will copy the style %s to %s.", style, dest)
                try:
                    self.upgrade_style(style, source, dest)
                    self.set_style(style)
                except Exception:
                    log.error("Could not set style")
            else:
                log.error(
                    "Could not set style %s. Probably you would want to"
                    "call rna.plotting.upgrade_style(<style>, "
                    "<path to mplstyle file that should be copied>)"
                    "once",
                    path,
                )

    @contextlib.contextmanager
    def _prepared_figure(self, kwargs):
        """
        Args:
            path
            *args: joined behind path
            **kwargs:
                axes (Axes): specify which axes to save. Higher priority than fig
                fig (Api): specifying which figure should be saved. If none gcf
                width (float): figure width in pts.
                scale (float): Fraction of the width which you wish
                    the figure to occupy
                ratio (float): height to width ratio


        """
        # catch figure from axes or fig
        axes = kwargs.pop("axes", None)
        if axes is None:
            fig_default = plt.gcf()
            axes = self.add_subplot()
            if fig_default is None:
                raise ValueError("fig_default may not be None")
        else:
            fig_default = axes.figure
        fig = kwargs.pop("fig", fig_default)

        # resize the figure to match a given width (in pt).
        width = kwargs.pop("width", None)
        scale = kwargs.pop("scale", 1.0)
        orig_size = fig.get_size_inches()
        orig_ratio = orig_size[1] / orig_size[0]
        ratio = kwargs.pop("ratio", orig_ratio)
        if width is not None:
            fig.set_size_inches(*rna.plotting.figsize(width, scale=scale, ratio=ratio))

        # crop the plot down based on the extents of the artists in the plot
        kwargs["bbox_inches"] = kwargs.pop("bbox_inches", "tight")
        if kwargs["bbox_inches"] == "tight":
            extra_artists = None
            for ax in fig.get_axes():
                first_label = ax.get_legend_handles_labels()[0] or None
                if first_label:
                    if not extra_artists:
                        extra_artists = []
                    if isinstance(first_label, list):
                        extra_artists.extend(first_label)
                    else:
                        extra_artists.append(first_label)
            kwargs["bbox_extra_artists"] = kwargs.pop(
                "bbox_extra_artists", extra_artists
            )

        # set current figure
        plt.figure(fig.number)

        yield

        # reset the figure size
        if width is not None:
            fig.set_size_inches(*orig_size)

    def _plt_savefig(self, path, **kwargs):
        with self._prepared_figure(kwargs):
            plt.savefig(path, **kwargs)

    def _save_tikz(self, path, **kwargs):
        width = kwargs.get("width", None)
        scale = kwargs.get("scale", 1.0)

        with self._prepared_figure(kwargs):
            import matplotlib2tikz as tikz

            if width is None:
                figurewidth = r"\textwidth"
            else:
                figurewidth = "{0}pt".format(scale * width)
            tikz.save(path, figurewidth=figurewidth)

    def _save_pgf(self, path, **kwargs):
        """
        Args:
            **kwargs:
                post_process (bool): switch post-processing of pgfs on and of
                substitutes (list of tuple of str):
                scaling (str): Default is \\igraph@width
                - further are passed to savefig and _pre_save
        """
        post_process = kwargs.pop("post_process", False)
        substitutes = kwargs.pop("substitutes", [(r"\sffamily", "")])
        scaling = kwargs.pop("scaling", r"\igraph@width")

        with self._prepared_figure(kwargs):
            plt.savefig(path, **kwargs)

            # pgf post-processing works but matpltolib pgf backend writes pgftext left
            # adjusted and creates right adjusted text by moving the text in x
            # this makes it impossible to rescale text correctly. If one knew, the
            # original was right adjusted a correction could be applied. This way
            # however, we cann not do anything.
            if post_process and path.endswith("pgf"):
                # use it like \igraph[width=0.5\columnwidth]{path_to_the.pgf}
                with open(path, "r") as f:
                    content = []
                    exp = r"(\\pgftext.+)(\\fontsize{.+)\\selectfont(.+)}%"
                    exp_width = (
                        r"\\pgfpathrectangle{\\pgfpointorigin}"
                        r"{\\pgfqpoint{(.+)}{(.+)}%"
                    )
                    for line in f:
                        res = re.match(exp, line)
                        res_width = re.match(exp_width, line)
                        if res_width:
                            pgf_path_size = res_width.groups()[0]
                            # factor 72.27 (latex definition differs from pdf
                            # definition) for conversion from in to pt
                            pgf_path_size = (
                                r"(" + pgf_path_size.replace("in", "*72.27") + r")"
                            )
                        if res:
                            groups = res.groups()
                            line = (
                                groups[0]
                                + r"\selectfont \scalebox{\fpeval{(%s/(%s))}}{"
                                % (pgf_path_size, scaling)
                                + groups[2]
                                + r"}}%"
                                + "\n"
                            )
                        for x, y in substitutes:
                            line = line.replace(x, y)
                        content.append(line)
                with open(path, "w") as f:
                    f.writelines(content)

    @rna.plotting.plot_signature
    def plot_function(self, axes, fun, **kwargs):
        """
        Args:
            axes (matplotlib.Axis) object
            fun: callable
            **kwargs:
                num: degree of detail (third argument of linspace)
                futher are forwarded to plot

        Returns:
            Artist or list of Artists (imitating the axes.scatter/plot behaviour).
            Better Artist not list of Artists
        """
        xmin, xmax = kwargs.pop("xmin", 0), kwargs.pop("xmax", 1)
        num = kwargs.pop("num", 100)
        vals = np.linspace(xmin, xmax, num)
        args = (vals, map(fun, vals))
        artist = axes.plot(*args, **kwargs)
        return artist

    @rna.plotting.plot_signature
    def plot_errorbar(self, axes, points, errors_up, errors_down=None, **kwargs):
        """
        Args:
            axes (matplotlib.Axis) object

        Returns:
            Artist or list of Artists (imitating the axes.scatter/plot behaviour).
            Better Artist not list of Artists
        """
        po = rna.plotting.PlotManager(kwargs)
        po.setdefault("marker", "_")

        if errors_down is None:
            errors_down = errors_up

        artists = []

        # plot errorbars
        for i in range(len(points)):
            artists.append(
                axes.plot(
                    [points[i, 0] + errors_up[i, 0], points[i, 0] - errors_down[i, 0]],
                    [points[i, 1], points[i, 1]],
                    [points[i, 2], points[i, 2]],
                    **kwargs
                )
            )
            artists.append(
                axes.plot(
                    [points[i, 0], points[i, 0]],
                    [points[i, 1] + errors_up[i, 1], points[i, 1] - errors_down[i, 1]],
                    [points[i, 2], points[i, 2]],
                    **kwargs
                )
            )
            artists.append(
                axes.plot(
                    [points[i, 0], points[i, 0]],
                    [points[i, 1], points[i, 1]],
                    [points[i, 2] + errors_up[i, 2], points[i, 2] - errors_down[i, 2]],
                    **kwargs
                )
            )

        return artists

    @rna.plotting.plot_signature
    def plot_tensor(
        self,
        axes,
        points: np.ndarray,
        field: typing.Optional[np.ndarray] = None,
        **kwargs
    ):
        """
        Args:
            points: base vectors
            field: direction field
            **kwargs: passed to PlotManager
        """
        x_index, y_index, z_index = self.pop_xyz_index(kwargs)
        points = np.array(points)
        if field is not None:
            field = np.array(field)

        if field is not None and len(field.shape) == 2 and field.shape[1] == 1:
            # scalar
            field = field.reshape(len(points))

        if field is None or len(field.shape) == 1:
            # scalar
            if self.axes_dim(axes) == 2:
                points = [points[:, x_index], points[:, y_index]]
            else:
                points = [points[:, x_index], points[:, y_index], points[:, z_index]]
            color = kwargs.pop("color", field)
            artists = axes.scatter(*points, c=color, **kwargs)
        elif len(field.shape) == 2:
            color = kwargs.get("color", None)
            if rna.plotting.colors.color_fmt(color) == "norm":
                kwargs.pop("color")
                kwargs.setdefault("norm", mpl.colors.Normalize)
                norm = kwargs.pop("norm")
                if issubclass(norm, mpl.colors.Normalize):
                    norm = norm(color.min(), color.max())
                kwargs["norm"] = norm
                args = (color,)
            else:
                kwargs.pop("norm", None)
                args = tuple()

            # vector
            if points is None:
                points = np.full(field.shape, 0.0)
            if self.axes_dim(axes) == 3:
                artists = axes.quiver(
                    points[:, x_index],
                    points[:, y_index],
                    points[:, z_index],
                    field[:, x_index],
                    field[:, y_index],
                    field[:, z_index],
                    *args,
                    **kwargs
                )
            elif self.axes_dim(axes) == 2:
                artists = axes.quiver(
                    points[:, x_index],
                    points[:, y_index],
                    field[:, x_index],
                    field[:, y_index],
                    *args,
                    **kwargs
                )
            else:
                raise NotImplementedError("Dimension != 2|3")
        else:
            raise NotImplementedError("Only Scalars and Vectors implemented")
        return artists

    @rna.plotting.plot_signature
    def plot_tensor_map(self, axes, tensors, map_, **kwargs):
        """
        Args:
            tensors: np.ndarray
            map_: indices of the tensors to connect to a graph

        Examples:
            >>> import rna
            >>> import tfields
            >>> tm = tfields.TensorMaps(
            ...     [[0, 0, 0], [0, 1, 0], [1, 0, 0], [1, 1, 0], [0.5, 0.5, 0]]
            ... )
            >>> tm.maps[2] = tfields.Tensors([[0, 1], [0, 2], [1, 3], [2, 3]])
            >>> artist = tm.plot(map=2, dim=2, color="r", linewidths=4)
            >>> rna.plotting.show()
        """
        x_index, y_index, _ = self.pop_xyz_index(kwargs)
        tensors = np.array(tensors)
        map_ = np.array(map_)
        if map_.shape[0] == 0:
            logging.warning("No map to plot")
            return None
        if max(map_.flat) > tensors.shape[0]:
            raise ValueError("Some indices point to non existing tensors.")

        dim = self.axes_dim(axes)
        if dim == 2:
            color = self.retrieve(
                kwargs, "edgecolor", "color", default=None, keep=False
            )
            if color is not None:
                kwargs["edgecolor"] = color
            tensors = np.transpose([tensors[:, x_index], tensors[:, y_index]])
            lines = np.array([tensors[edge] for edge in map_])
            artist = mcollections.LineCollection(lines, **kwargs)
            axes.add_collection(artist)
        else:
            raise NotImplementedError("")
        return artist

    @rna.plotting.plot_signature
    def plot_mesh(self, axes, vertices, faces, **kwargs):
        """
        Args:
            axes (matplotlib axes)
            vertices (np.ndarray): vertices of the mesh (shape: (n, 3))
            faces (np.ndarray): faces of the mesh (indices of vertices) (shape: (m, 3))
            x_index (int)
            y_index (int)
            z_index (int)
            edgecolor (color)
            color (color): if given, use this color for faces in 2D
            cmap
            vmin
            vmax
        """
        x_index, y_index, z_index = self.pop_xyz_index(kwargs)
        vertices = np.array(vertices)
        faces = np.array(faces)
        if faces.shape[0] == 0:
            logging.warning("No faces to plot")
            return None
        if max(faces.flat) > vertices.shape[0]:
            raise ValueError("Some faces point to non existing vertices.")

        dim = self.axes_dim(axes)
        if dim == 2:
            full = True
            import tfields

            mesh = tfields.Mesh3D(vertices, faces=faces)
            transformation = np.eye(4)
            transformation[:3, :3] = kwargs.pop("rotation", np.eye(3))
            transformation[:3, 3] = kwargs.pop("translation", np.array([0, 0, 0]))

            # rotation
            mesh[:] = np.einsum("...i,t...i->t...", transformation[:3, :3], mesh)
            # translation
            mesh += transformation[:3, 3]

            facecolors = self.retrieve(
                kwargs, "facecolors", "color", default=0, keep=False
            )
            if full:
                # implementation that will sort the triangles by z_index
                centroids = mesh.centroids()
                axis_indices = [0, 1, 2]
                axis_indices.pop(axis_indices.index(x_index))
                axis_indices.pop(axis_indices.index(y_index))
                z_index = axis_indices[0]
                zs = centroids[:, z_index]
                try:
                    # no numbers
                    iter(facecolors)
                    # no strings
                    if isinstance(facecolors, str):
                        raise TypeError()
                except TypeError:
                    zs, faces = tfields.lib.util.multi_sort(zs, faces)
                else:
                    zs, faces, facecolors = tfields.lib.util.multi_sort(
                        zs, faces, facecolors
                    )

                n_faces_initial = len(faces)
            else:
                # cut away "back sides" implementation
                direction_vector = np.array([1.0, 1.0, 1.0])
                direction_vector[x_index] = 0.0
                direction_vector[y_index] = 0.0
                norm_vectors = mesh.triangles().norms()
                dot_product = np.dot(norm_vectors, direction_vector)
                n_faces_initial = len(faces)
                faces = faces[dot_product > 0]

            vertices = [mesh[:, x_index], mesh[:, y_index]]

            # sort out color arguments
            facecolors = self.format_colors(
                kwargs, facecolors, fmt="norm", length=n_faces_initial
            )
            if not full:
                facecolors = facecolors[dot_product > 0]

            artist = axes.tripcolor(
                *vertices, triangles=faces, facecolors=facecolors, **kwargs
            )
        elif dim == 3:
            label = kwargs.pop("label", None)
            color = self.retrieve(
                kwargs, "color", "c", "facecolors", default="grey", keep=False
            )
            from rna.plotting.colors import color_fmt

            orig_fmt = color_fmt(color)
            if orig_fmt == "norm":
                orig_color = self.format_colors(
                    kwargs, color, fmt="norm", length=len(faces)
                )

            color = self.format_colors(
                kwargs, color, fmt="rgba", length=len(faces), dtype=float
            )

            nan_mask = np.isnan(color)
            if nan_mask.any():
                logging.warning("nan found in colors. Removing corresponding faces!")
                color = color[~nan_mask]
                faces = faces[~nan_mask]

            edgecolor = kwargs.pop("edgecolor", None)
            alpha = kwargs.pop("alpha", None)
            self.pop_norm_args(kwargs)

            triangles = np.array([vertices[face] for face in faces])
            artist = plt3d.art3d.Poly3DCollection(triangles, **kwargs)
            axes.add_collection3d(artist)

            if edgecolor is not None:
                artist.set_edgecolor(edgecolor)
                artist.set_facecolors(color)
            else:
                artist.set_color(color)
            if orig_fmt == "norm":
                # sets vmin, vmax to min, max of val
                artist.set_array(np.array(orig_color))
            if alpha is not None:
                artist.set_alpha(alpha)

            # for some reason auto-scale does not work
            self.autoscale_3d(axes, array=vertices)

            # legend lables do not work at all as an argument
            if label:
                artist.set_label(label)

            # when plotting the legend edgecolors/facecolors2d are needed
            artist._edgecolors2d = None
            artist._facecolors2d = None

        else:
            raise NotImplementedError("Dimension != 2|3")

        return artist

    @rna.plotting.plot_signature
    def plot_sphere(self, axes, point, radius, **kwargs):
        # Make data
        u = np.linspace(0, 2 * np.pi, 100)
        v = np.linspace(0, np.pi, 100)
        x = point[0] + radius * np.outer(np.cos(u), np.sin(v))
        y = point[1] + radius * np.outer(np.sin(u), np.sin(v))
        z = point[2] + radius * np.outer(np.ones(np.size(u)), np.cos(v))

        # Plot the surface
        return axes.plot_surface(x, y, z, **kwargs)

    @rna.plotting.plot_signature
    def plot_plane(self, axes, point, normal, **kwargs):
        def rotation_matrix(d):
            sin_angle = np.linalg.norm(d)
            if sin_angle == 0:
                return np.identity(3)
            d /= sin_angle
            eye = np.eye(3)
            ddt = np.outer(d, d)
            skew = np.array(
                [[0, d[2], -d[1]], [-d[2], 0, d[0]], [d[1], -d[0], 0]], dtype=np.float64
            )

            M = ddt + np.sqrt(1 - sin_angle**2) * (eye - ddt) + sin_angle * skew
            return M

        def pathpatch_2d_to_3d(pathpatch, z, normal):
            if type(normal) is str:  # Translate strings to normal vectors
                index = "xyz".index(normal)
                normal = np.roll((1.0, 0, 0), index)

            normal /= np.linalg.norm(normal)  # Make sure the vector is normalised
            path = pathpatch.get_path()  # Get the path and associated transform
            trans = pathpatch.get_patch_transform()

            path = trans.transform_path(path)  # Apply the transform

            pathpatch.__class__ = plt3d.art3d.PathPatch3D  # Change the class
            pathpatch._code3d = path.codes  # Copy the codes
            pathpatch._facecolor3d = pathpatch.get_facecolor  # Get the face color

            verts = path.vertices  # Get the vertices in 2D

            d = np.cross(normal, (0, 0, 1))  # Obtain the rotation vector
            M = rotation_matrix(d)  # Get the rotation matrix

            pathpatch._segment3d = np.array(
                [np.dot(M, (x, y, 0)) + (0, 0, z) for x, y in verts]
            )

        def pathpatch_translate(pathpatch, delta):
            pathpatch._segment3d += delta

        kwargs["alpha"] = kwargs.pop("alpha", 0.5)
        patch = mplpatch.Circle((0, 0), **kwargs)
        axes.add_patch(patch)
        pathpatch_2d_to_3d(patch, z=0, normal=normal)
        pathpatch_translate(patch, (point[0], point[1], point[2]))

    @rna.plotting.plot_signature
    def autoscale_3d(self, axes, array=None, xLim=None, yLim=None, zLim=None):
        if array is not None:
            xMin, yMin, zMin = array.min(axis=0)
            xMax, yMax, zMax = array.max(axis=0)
            xLim = (xMin, xMax)
            yLim = (yMin, yMax)
            zLim = (zMin, zMax)
        xLimAxis = axes.get_xlim()
        yLimAxis = axes.get_ylim()
        zLimAxis = axes.get_zlim()

        if not False:
            # not empty axes
            xMin = min(xLimAxis[0], xLim[0])
            yMin = min(yLimAxis[0], yLim[0])
            zMin = min(zLimAxis[0], zLim[0])
            xMax = max(xLimAxis[1], xLim[1])
            yMax = max(yLimAxis[1], yLim[1])
            zMax = max(zLimAxis[1], zLim[1])
        axes.set_xlim([xMin, xMax])
        axes.set_ylim([yMin, yMax])
        axes.set_zlim([zMin, zMax])

    @rna.plotting.plot_signature
    def set_legend(self, axes, *artists, **kwargs):
        """
        Convenience method to set a legend from multiple artists to an axes.
        Args:
            **kwargs
                table (bool): if True, labels containing ',' are mapped to table
                table_title (str): value of the table entry top left - only active
                    if table
        Examples:
            >>> import rna
            >>> import matplotlib.pyplot as plt
            >>> import pylab

            >>> fig = plt.figure()
            >>> ax = fig.add_subplot(111)

            >>> im1 = ax.plot(range(10), pylab.randn(10), "r--",
            ...              label=(r"$i = 1$,$j = 1$"))
            >>> im2 = ax.plot(range(10), pylab.randn(10), "g--",
            ...              label=(r"$i = 1$,$j = 2$"))
            >>> im3 = ax.plot(range(10), pylab.randn(10), "b--",
            ...              label=(r"$i = 1$,$j = 3$"))
            >>> im4 = ax.plot(range(10), pylab.randn(10), "r.",
            ...              label=(r"$i = 2$,$j = 1$"))
            >>> im5 = ax.plot(range(10), pylab.randn(10), "g.",
            ...              label=(r"$i = 2$,$j = 2$"))
            >>> im6 = ax.plot(range(10), pylab.randn(10), "b.",
            ...              label=(r"$i = 2$,$j = 3$"))
            >>> im7 = ax.plot(range(10), pylab.randn(10), "r^",
            ...              label=(r"$i = 3$,$j = 1$"))
            >>> im8 = ax.plot(range(10), pylab.randn(10), "g^",
            ...              label=(r"$i = 3$,$j = 2$"))
            >>> im9 = ax.plot(range(10), pylab.randn(10), "b^",
            ...              label=(r"$i = 3$,$j = 3$"))
            >>> handles = [im1, im2, im3, im4, im5, im6, im7, im8, im9]

            >>> legend = rna.plotting.set_legend(*handles, table=True, axes=ax)

            >>> rna.plotting.show()
        """
        table = kwargs.pop("table", False)
        labels = kwargs.pop("labels", None)
        ncol = kwargs.pop("ncol", None)

        handles = []
        for artist in artists:
            if isinstance(artist, list):
                handles.extend(more_itertools.collapse(artist, base_type=tuple))
            elif isinstance(artist, tuple):
                # Dont know any more when this occurs
                tuple_handle = tuple()
                for sub_artist in artist:
                    if isinstance(sub_artist, list):
                        tuple_handle += (sub_artist[0],)
                    else:
                        tuple_handle += (sub_artist,)
                handles.append(tuple_handle)
            else:
                handles.append(artist)

        if labels is None and any([isinstance(h, tuple) for h in handles]):
            labels = []
            for h in handles:
                if isinstance(h, tuple):
                    sub_labels = [sub_h.get_label() for sub_h in h]
                    label = sub_labels[0]
                    for sub_label in sub_labels[1:]:
                        if not sub_label.startswith("_"):
                            label += sub_label
                else:
                    label = h.get_label()
                labels.append(label)

        if table and labels is None and ncol is None:
            table_title = kwargs.pop("table_title", "")
            labels = np.array([h.get_label() for h in handles])
            labels = [label.split(",") for label in labels]
            captions_i = []
            captions_j = []
            for label in labels:
                if label[0] not in captions_i:
                    captions_i.append(label[0])
                if label[1] not in captions_j:
                    captions_j.append(label[1])

            shape = (len(captions_i), len(captions_j))

            # initialize
            shape = np.array(shape)
            handles = np.array(handles)

            # create blank rectangle
            extra = mplpatch.Rectangle(
                (0, 0), 1, 1, fc="w", fill=False, edgecolor="none", linewidth=0
            )

            # Create organized list containing all handles for table.
            # Extra represent empty space
            handles_table = np.full(shape + 1, extra, dtype=object)
            for handle, label in zip(handles, labels):
                i = captions_i.index(label[0])
                j = captions_j.index(label[1])
                if handles_table[i + 1, j + 1] != extra:
                    raise ValueError("Duplicate label {label}".format(**locals()))
                handles_table[i + 1, j + 1] = handle

            # Define the label captions
            labels_table = np.full(shape + 1, "", dtype="S80")
            labels_table[0, 0] = table_title
            labels_table[0, 1:] = captions_j
            labels_table[1:, 0] = captions_i
            labels_table = labels_table.astype(str)

            handles = list(handles_table.flat)
            labels = list(labels_table.flat)
            kwargs["ncol"] = shape[0] + 1
            # negative numbers in handletextpad move to the right
            kwargs["handletextpad"] = kwargs.pop("handletextpad", -1.5)
            kwargs["columnspacing"] = kwargs.pop("columnspacing", 1.5)

        return axes.legend(handles=handles, labels=labels, **kwargs)

    def set_axis_off(self, axes):
        """
        Turn off axis display
        """
        if self.axes_dim(axes) == 2:
            axes.set_axis_off()
        else:
            axes._axis3don = False

    def set_aspect_equal(self, axes):
        """Fix equal aspect bug for 3D plots."""

        if self.axes_dim(axes) == 2:
            axes.set_aspect("equal")
            return

        xlim = axes.get_xlim3d()
        ylim = axes.get_ylim3d()
        zlim = axes.get_zlim3d()

        from numpy import mean

        xmean = mean(xlim)
        ymean = mean(ylim)
        zmean = mean(zlim)

        plot_radius = max(
            [
                abs(lim - mean_)
                for lims, mean_ in ((xlim, xmean), (ylim, ymean), (zlim, zmean))
                for lim in lims
            ]
        )

        axes.set_xlim3d([xmean - plot_radius, xmean + plot_radius])
        axes.set_ylim3d([ymean - plot_radius, ymean + plot_radius])
        axes.set_zlim3d([zmean - plot_radius, zmean + plot_radius])

    def set_limits(self, axes, *limits):
        """
        Args:
            *limits: list of 2-tuples containing min and max limit of axes
        """
        dim = self.axes_dim(axes)
        for coord, limit in zip(["x", "y", "z"][:dim], limits):
            if dim == 2:
                set_method = getattr(axes, "set_" + coord + "lim")
            elif dim == 3:
                set_method = getattr(axes, "set_" + coord + "limd3")

            set_method(*limit)  # limit = xmin, xmax

    @rna.plotting.plot_signature
    def set_formatter(
        self, axes, sub_axes=None, formatter=dates.DateFormatter("%d-%m-%y")
    ):
        """
        Examples:
            >>> from rna.plotting.backends.matplotlib import ApiBackendMatplotlib
            >>> axes = ApiBackendMatplotlib().add_subplot()
            >>> ApiBackendMatplotlib().set_formatter(axes)
        """
        if sub_axes is None:
            sub_axes = axes.xaxis
        sub_axes.set_major_formatter(formatter)


class ScientificFormatter(mpl.ticker.ScalarFormatter):
    """
    Examples:
        >> cbar = rna.plotting.set_colorbar(
        ...     axes, artist,
        ...     label=r"$q_c\;(MW/m^2)$",  # NOQA
        ...     format=rna.plotting.backend.ScientificFormatter(
        ...         None, useMathText=False))

    Args:
        oom (int): order of magnitued on the axes
        **kwargs: forwarded to ScalarFormatter

    """

    def __init__(self, oom=None, **kwargs):
        self._oom = oom
        super().__init__(**kwargs)

    def _set_orderOfMagnitude(self, oom):  # pylint: disable=invalid-name
        exp = int(np.log10(oom))
        if self._oom is not None:
            oom = self._oom
        else:  # Default: -3, 0, 3, 6, ...
            oom = exp - (exp % 3)
        self.orderOfMagnitude = oom

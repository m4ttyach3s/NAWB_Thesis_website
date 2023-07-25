"""
Notes on mapping to matplotlib language:
References:
    * https://pyqtgraph.readthedocs.io/en/latest/plotting.html

Examples:
    $ python -m pyqtgraph.examples

Philosophy:
?       -> QMainWindow

Figure  ->
            LayoutWidget (3d) -> QWidget
            pg.GraphicsLayoutWidget (2D) -> GraphicsWidget -> (GraphicsItem, QGraphicsWidget)
                    -> WINDOW
Axes    -> PlotItem (2D) / GLViewWidget (3D)
Artist  -> GraphicsItem
"""
import typing
import logging
import sys
import warnings
import numpy as np
import pyqtgraph as pg
import pyqtgraph.opengl
import pyqtgraph.exporters
import rna


pg.setConfigOption("background", "w")
pg.setConfigOption("foreground", "w")
pg.setConfigOptions(antialias=True)


global WINDOW, VIEW
WINDOW, VIEW = None, None


def axes_dim(axes):
    if isinstance(axes, pg.PlotItem):
        return 2
    if isinstance(axes, pg.opengl.GLViewWidget):
        return 3
    raise NotImplementedError("Dimension of {axes} not mapped".format(**locals()))


def add_subplot(dim=None, **kwargs):
    global WINDOW, VIEW
    if dim is None:
        dim = 2

    if WINDOW is None:
        pg.mkQApp()

    if dim == 2:
        VIEW = pg.PlotItem()  # -> GraphicsWidget -> (GraphicsItem, QgraphicsWidget)
        WINDOW = pg.GraphicsLayoutWidget(title="rna 2D plot")  # show=False
        WINDOW.addItem(VIEW)  # , row, col, rowspan, colspan)
    elif dim == 3:
        VIEW = pg.opengl.GLViewWidget()
        WINDOW = pg.LayoutWidget()
        WINDOW.addWidget(VIEW)
    else:
        raise NotImplementedError()

    WINDOW.resize(800, 800)

    return VIEW


def gca(dim=None):
    """
    Returns a 2D or 3D plot view object depending on the value of the `dim` parameter.
    The returned object is of type `pg.PlotItem` for 2D views, or `pg.opengl.GLViewWidget` for 3D views.
    """
    global WINDOW, VIEW

    initiate = False

    if VIEW is not None:
        if dim is None:
            return VIEW
        if axes_dim(VIEW) == dim:
            return VIEW
        else:
            initiate = True

    if dim is None:
        dim = 2

    if WINDOW is None:
        pg.mkQApp()
        initiate = True

    if dim == 2:
        VIEW = pg.PlotItem()  # -> GraphicsWidget -> (GraphicsItem, QgraphicsWidget)
        if initiate:
            WINDOW = pg.GraphicsLayoutWidget(title="rna 2D plot")  # show=False
        WINDOW.addItem(VIEW)  # , row, col, rowspan, colspan)

    elif dim == 3:
        VIEW = pg.opengl.GLViewWidget()
        if initiate:
            WINDOW = pg.LayoutWidget()
        WINDOW.addWidget(VIEW)
    else:
        raise NotImplementedError()

    if initiate:
        WINDOW.resize(800, 800)

    return VIEW


def show():
    global WINDOW, VIEW
    WINDOW.show()

    if (sys.flags.interactive != 1) or not hasattr(pg.Qt.QtCore, "PYQT_VERSION"):
        # Protect exectuion: https://stackoverflow.com/a/59691688/6233888
        pg.Qt.QtGui.QApplication.instance().exec_()


class GraphicsViewOverlay2D(pg.GraphicsView):
    """
    GraphicsView subclass that uses GLViewWidget as its canvas.
    This allows 2D graphics to be overlaid on a 3D background.
    """

    def __init__(self, axes):
        """
        Args:
            axes (gl.GLViewWidget): view upon which to copy all events
        """
        self.glView = axes
        pg.GraphicsView.__init__(self, background=None)
        self.setStyleSheet("background: transparent")
        self.setViewport(self.glView)

    def mousePos(self, event):
        """
        Distribute mouse events to both widgets
        """
        pg.GraphicsView.mousePos(self, event)
        self.glView.mousePos(event)

    def paintEvent(self, event):
        """
        Distribute paint events to both widgets
        """
        self.glView.paintEvent(event)
        return pg.GraphicsView.paintEvent(self, event)

    def mousePressEvent(self, event):
        """
        Distribute mouse events to both widgets
        """
        pg.GraphicsView.mousePressEvent(self, event)
        self.glView.mousePressEvent(event)

    def mouseMoveEvent(self, event):
        """
        Distribute mouse events to both widgets
        """
        pg.GraphicsView.mouseMoveEvent(self, event)
        self.glView.mouseMoveEvent(event)

    def mouseReleaseEvent(self, event):
        """
        Distribute mouse events to both widgets
        """
        pg.GraphicsView.mouseReleaseEvent(self, event)
        self.glView.mouseReleaseEvent(event)


class ApiBackendPyqtgraph(rna.plotting.ApiBackend):
    EXPLICIT_SAVE_METHODS = {ext: "_frame_save" for ext in ("png")}

    def show(self):
        show()

    @staticmethod
    def clear(axes):
        global WINDOW, VIEW
        VIEW.setParent(None)
        # WINDOW.removeWidget(VIEW)
        VIEW = None

    @staticmethod
    def add_subplot(dim=None, **kwargs):
        return add_subplot(dim=dim, **kwargs)

    @staticmethod
    def gca(dim=None, **kwargs):
        warnings.warn("gca is deprecated, use add_suplot instead")
        return gca(dim=dim, **kwargs)

    @staticmethod
    def is_axes(obj):
        return isinstance(obj, (pg.PlotItem, pg.opengl.GLViewWidget))

    @staticmethod
    def axes_dim(axes):
        return axes_dim(axes)

    @staticmethod
    def set_colorbar(axes, artist, label=None, position="right", **kwargs):
        # Make a 2D color bar using the same ColorMap
        overlay_view = GraphicsViewOverlay2D(axes)
        kwargs["size"] = kwargs.pop("size", (20, 200))
        kwargs["offset"] = kwargs.pop("offset", (15, -25))
        colorbar_legend = Colorbar(artist, label=label, **kwargs)
        overlay_view.addItem(colorbar_legend)
        return overlay_view

    @staticmethod
    def get_cmap(cm_name):
        """
        Use matplotlib conventions to retrieve a pg.ColorMap

        Examples:
            >>> import rna
            >>> rna.plotting.use('pyqtgraph')
            >>> from rna.plotting.backends.pyqtgraph import ApiBackendPyqtgraph
            >>> cm = ApiBackendPyqtgraph.get_cmap('viridis')
            >>> rna.plotting.use('matplotlib')  # reset for test suite to start at default
        """
        import matplotlib.pyplot as plt

        plt_map = plt.get_cmap(cm_name)
        colors = plt_map.colors
        colors = [c + [1.0] for c in colors]
        positions = np.linspace(0, 1, len(colors))
        cmap = pg.ColorMap(positions, colors)
        return cmap

    @staticmethod
    def _frame_save(path, **kwargs):
        axes = kwargs.get("axes", None)
        if axes is None:
            axes = gca()
        axes.grabFrameBuffer().save(path, **kwargs)
        # global WINDOW
        # exporter = pg.exporters.ImageExporter(WINDOW.plotItem)
        # exporter.parameters().update(kwargs)
        # exporter.export(path)

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
        raise NotImplementedError(
            "Put a lot of stuff from plot_tensor to ApiBackend and implement"
        )

    @rna.plotting.plot_signature
    def plot_mesh(self, axes, vertices, faces, **kwargs):
        # TODO-1 drawEdges -> edgecolor
        kwargs.setdefault("drawEdges", False)
        kwargs.setdefault("smooth", False)
        kwargs.setdefault("shader", None)
        scalars = kwargs.pop("color", None)
        if scalars is None:
            scalars = [0] * len(faces)
        cmap, vmin, vmax = self.get_norm_args(
            kwargs,
            vmin_default=min(scalars),
            vmax_default=max(scalars),
            cmap_default="viridis",
        )

        colors = self.format_colors(kwargs, scalars, fmt="rgba", length=len(faces))
        colors = np.array([list(c) for c in colors])

        mesh_data = pg.opengl.MeshData(
            vertexes=np.array(vertices), faces=np.array(faces), faceColors=colors
        )
        artist = pg.opengl.GLMeshItem(meshdata=mesh_data, **kwargs)
        artist.colors = colors
        artist.vmin = vmin
        artist.vmax = vmax
        artist.cmap = cmap
        artist.setGLOptions("opaque")

        axes.addItem(artist)
        return artist


class Colorbar(pg.GradientLegend):
    def __init__(self, artist, *args, **kwargs):
        self.label = kwargs.pop("label")
        self.format = kwargs.pop("format", "{0:.2f}")
        self.vmin, self.vmax = artist.vmin, artist.vmax
        ticks = kwargs.pop("ticks", np.linspace(self.vmin, self.vmax, 4))
        super(Colorbar, self).__init__(*args, **kwargs)
        logging.error("Problem with colorbar not yet solved")
        # color_map = get_cmap(artist.cmap)
        # self.setGradient(color_map.getGradient())
        self.set_tick_positions(ticks)

    def set_tick_positions(self, tick_positions):
        if self.vmin > min(tick_positions) or self.vmax < max(tick_positions):
            raise ValueError("tick positions must lie in (vmin, vmax)")

        self.setLabels(
            dict(
                [
                    (
                        self.format.format(v),
                        1 / (self.vmax - self.vmin) * v
                        - self.vmin / (self.vmax - self.vmin),
                    )  # NOQA
                    for v in tick_positions
                ]
            )
        )

    def paint(self, p, opt, widget):
        pg.UIGraphicsItem.paint(self, p, opt, widget)
        rect = self.boundingRect()  # Bounds of visible area in scene coords.
        unit = self.pixelSize()  # Size of one view pixel in scene coords.
        if unit[0] is None:
            return

        # determine max width of all labels
        labelWidth = 0
        labelHeight = 0
        for k in self.labels:
            b = p.boundingRect(
                pg.Qt.QtCore.QRectF(0, 0, 0, 0),
                pg.Qt.QtCore.Qt.AlignLeft | pg.Qt.QtCore.Qt.AlignVCenter,
                str(k),
            )
            labelWidth = max(labelWidth, b.width())
            labelHeight = max(labelHeight, b.height())

        b_tick = p.boundingRect(
            pg.Qt.QtCore.QRectF(0, 0, 0, 0),
            pg.Qt.QtCore.Qt.AlignCenter | pg.Qt.QtCore.Qt.AlignVCenter,
            "-",
        )
        b_label = p.boundingRect(
            pg.Qt.QtCore.QRectF(0, 0, 0, 0),
            pg.Qt.QtCore.Qt.AlignCenter | pg.Qt.QtCore.Qt.AlignVCenter,
            str(self.label),
        )

        labelWidth *= unit[0]
        labelHeight *= unit[1]

        textPadding = 2  # in px

        if self.offset[0] < 0:
            x3 = rect.right() + unit[0] * self.offset[0]
            x2 = x3 - labelWidth - unit[0] * textPadding * 2
            x1 = x2 - unit[0] * self.size[0]
        else:
            x1 = rect.left() + unit[0] * self.offset[0]
            x2 = x1 + unit[0] * self.size[0]
            x3 = x2 + labelWidth + unit[0] * textPadding * 2
        if self.offset[1] < 0:
            y2 = rect.top() - unit[1] * self.offset[1]
            y1 = y2 + unit[1] * self.size[1]
        else:
            y1 = rect.bottom() - unit[1] * self.offset[1]
            y2 = y1 - unit[1] * self.size[1]
        self.b = [x1, x2, x3, y1, y2, labelWidth]

        # Have to scale painter so that text and gradients are correct size.
        p.scale(unit[0], unit[1])

        # Draw color bar
        self.gradient.setStart(0, y1 / unit[1])
        self.gradient.setFinalStop(0, y2 / unit[1])
        p.setBrush(self.gradient)
        rect = pg.Qt.QtCore.QRectF(
            pg.Qt.QtCore.QPointF(x1 / unit[0], y1 / unit[1]),
            pg.Qt.QtCore.QPointF(x2 / unit[0], y2 / unit[1]),
        )
        p.drawRect(rect)

        # draw labels
        p.setPen(pg.Qt.QtGui.QPen(pg.Qt.QtGui.QColor(0, 0, 0)))
        tx = x2 + unit[0] * textPadding
        lh = labelHeight / unit[1]
        for k in self.labels:
            y = y1 + self.labels[k] * (y2 - y1)
            p.drawText(
                pg.Qt.QtCore.QRectF(
                    tx / unit[0] + b_tick.width(), y / unit[1] - lh / 2.0, 1000, lh
                ),
                pg.Qt.QtCore.Qt.AlignLeft | pg.Qt.QtCore.Qt.AlignVCenter,
                str(k),
            )
            p.drawText(
                pg.Qt.QtCore.QRectF(
                    x2,
                    y / unit[1]
                    - b_tick.height() / 2.0
                    - 2,  # the -2 is a hack to move the ticks up  # NOQA
                    b_tick.width(),
                    b_tick.height(),
                ),
                pg.Qt.QtCore.Qt.AlignLeft | pg.Qt.QtCore.Qt.AlignVCenter,
                "-",
            )
        p.save()
        p.translate(x3 + b_label.height(), y1 + 0.5 * (y2 - y1) - b_label.width() / 2)
        p.rotate(90)
        p.drawText(0.0, 0.0, str(self.label))
        p.restore()

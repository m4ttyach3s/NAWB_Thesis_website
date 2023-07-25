"""
Various polymorphisms which can be used as mixins
"""
import pathlib
import typing
import logging
import io

import rna.path


class Storable:
    """
    Polymorphism with abstract standardized save method and load factory.

    Examples:
        >>> import rna
        >>> class MyStorable(rna.polymorphism.Storable):
        ...     def _save_dummy(self, path):
        ...         print('Saving to', path)
        ...     @classmethod
        ...     def _load_dummy(cls, path):
        ...         print('Loading from', path)
        >>> obj = MyStorable()
        >>> obj.save("/base/path", "to", "obj.dummy", create_parents=False)
        Saving to /base/path/to/obj.dummy
        >>> obj2 = MyStorable.load("/base/path", "to", "obj.dummy")
        Loading from /base/path/to/obj.dummy

    """

    # mapping from extension[str] to method name [str]
    EXPLICIT_SAVE_METHODS: typing.Optional[dict] = None

    def save(
        self,
        path: typing.Union[str, io.IOBase],
        *args: [str],
        create_parents: bool = True,
        extension: typing.Optional[typing.Union[str, typing.List[str]]] = None,
        **kwargs
    ) -> None:
        """
        Saving by redirecting to the correct save method depending on path.
        Implement _save_{extension}(self, path, **kwargs) for saving to your extension of choice

        Args:
            path (str | buffer)
            *args: joined with path
            create_parents: create the parent dir(s) if required
            **kwargs:
                extension (str | List[str]): explicit extension required if path is buffer,
                    if list of extensions save as all extension formats.
                ... remaining:forwarded to extension specific method
        """
        # get the extension
        if isinstance(extension, list):
            for ext in extension:
                print(ext)
                self.save(
                    path, *args, create_parents=create_parents, extension=ext, **kwargs
                )
            return

        if isinstance(path, (str, pathlib.Path)):
            path = rna.path.resolve(path, *args)
            if extension is None:
                extension = rna.path.extension(path)
            if not path.endswith(extension):
                path += "." + extension
        if not extension:
            raise ValueError(
                "Path {path} requires extension for auto rooting.".format(**locals())
            )

        # get the save method
        if (
            self.EXPLICIT_SAVE_METHODS is not None
            and extension
            in self.EXPLICIT_SAVE_METHODS  # pylint:disable=unsupported-membership-test
        ):
            # pylint:disable=unsubscriptable-object
            method_name = self.EXPLICIT_SAVE_METHODS[extension]
        else:
            method_name = "_save_" + extension
        try:
            save_method = getattr(self, method_name)
        except AttributeError as err:
            raise NotImplementedError(
                "Save method {method_name} for extension: {extension} required.".format(
                    **locals()
                )
            ) from err

        if create_parents:
            rna.path.mkdir(path)

        log = logging.getLogger(__name__)
        log.info("Saving to path '%s' with method %s.", path, save_method)
        save_method(path, **kwargs)

    @classmethod
    def load(cls, path: typing.Union[str, io.IOBase], *args: [str], **kwargs):
        """
        Instantiate an object of this class from file.
        Implement _load_{extension}(cls, path, *args) -> obj for loading from extension of choice.

        Args:
            path (str or buffer)
            *args: joined with path
            **kwargs:
                extension (str): only needed if path is buffer
                ... remaining:forwarded to extension specific method
        """
        if isinstance(path, (str, pathlib.Path)):
            # pylint: disable=possibly-unused-variable
            path = rna.path.resolve(path, *args)
            extension = kwargs.pop("extension", pathlib.Path(path).suffix.lstrip("."))
        else:
            extension = kwargs.get("extension")

        if not extension:
            raise ValueError("Path requires extension for auto rooting.")

        try:
            load_method = getattr(cls, "_load_{extension}".format(**locals()))
        except AttributeError as err:
            raise NotImplementedError(
                "Load method '_load_{extension}' for extension: {extension} required.".format(
                    **locals()
                )
            ) from err

        return load_method(path, **kwargs)

    def reload(self, path):
        """
        Update the object from path. This requires the inheritor of this polymorphism
        to provide the 'def __update__(self, other):' method
        """
        if not hasattr(self, "__update__"):
            raise NotImplementedError("__update__ required")
        other = type(self).load(path)
        self.___update__(other)

    def _save_pickle(self, path):
        import pickle  # pylint:disable = import-outside-toplevel

        with open(path, "wb") as out_file:
            pickle.dump(self, out_file)

    @classmethod
    def _load_pickle(cls, path):
        import pickle  # pylint:disable = import-outside-toplevel

        with open(path, "rb") as in_file:
            obj = pickle.load(in_file)
        return obj

    def _save_dill(self, path):
        import dill  # pylint:disable = import-outside-toplevel

        with open(path, "wb") as out_file:
            dill.dump(self, out_file)

    @classmethod
    def _load_dill(cls, path):
        import dill  # pylint:disable = import-outside-toplevel

        with open(path, "rb") as in_file:
            obj = dill.load(in_file)
        return obj


class Plottable:  # pylint: disable=too-few-public-methods
    """
    Polymorphism for object which can be plotted
    The polymorphisms lies in the method name *plot' and the suggested kwargs for this method.
    """

    def plot(self, *args, **kwargs):
        """
        Args:
            **kwargs:
                axes: Most fundamental object of the plotting backend e.g. matplotlib.pyplot.Axes
                    Ideally the plot method implements multiple backends which are automatically
                    selected on the basis of the axes object.
                color: color to plot this object. This is to unify the various different signatures
                    of e.g. matplotlib
                cmap (str | colormap object): colorbar map from scalar to color
                vmin: minimum color value for color bar range
                vmax: maximum color value for color bar range
                label(str): label of plotting object
        """
        raise NotImplementedError()

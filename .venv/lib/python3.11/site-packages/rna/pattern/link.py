import typing
from rna.pattern.switch import Switch


NOTSET = object()  #


# TODO-0(@dboe): Allthough this is nice, the default state of exposed should be True
# TODO-0(@dboe): the global switch variable must be on the same namespace as the Linker class,
# otherwise you can not use the Linker patter twice in the same interpreter
class exposed(Switch):
    """
    This switch is related to the Linker pattern and - if enabled - exposes the Link behind
    a Linker attribute.
    """

    ENABLED = False


class Reference:
    # TODO-1(@dboe): Rename _references _linkers
    # TODO-1(@dboe): Think about garbage collection. If e.g. states are saved in references
    # I guess they are not removed via garbage collection...
    def __init__(self, references: typing.List[typing.Tuple[typing.Any, str]] = None):
        if references is None:
            references = []
        self._references = references


class Link:
    """
    Base class of an attribute in the Linker class that acts as a symlink to a reference.

    Args:
        value: value of the attribute.
        ref: reference containing the value
        fget: callable, following fget(ref) -> value or fget() -> value if ref is INDEPENDENT
        cache: if true, store the value after fget and return the stored value after that.
    """

    def __init__(
        self,
        ref: typing.Union[Reference, typing.Any] = NOTSET,
        fget: typing.Callable = NOTSET,
        value: typing.Any = NOTSET,
        cache=False,
    ):
        self.value = value
        self.ref = ref
        self.fget = fget
        self.cache = cache
        if value is NOTSET:
            assert fget is not None


class LinkNotFoundError(TypeError):
    """
    To be raised when a Link type is not expected found
    """

    pass


class Linker:
    """
    Linkers are handling the attribute acces of Links by __getattribute__.
    Linker attributes reveal the ref with the exposed switch set to True

    Examples:

        >>> import dataclasses
        >>> from rna.pattern.link import Linker, Link, Reference, exposed

        >>> @dataclasses.dataclass
        ... class Wout(Reference):
        ...     path: str = None
        ...     data_1: any = None
        ...     data_2: any = None
        ...
        ...     def __post_init__(self):
        ...         super().__init__()

        >>> wout = Wout(path="wout_asdf.nc", data_1=21, data_2=42)

        >>> @dataclasses.dataclass
        ... class Equilibrium(Linker):
        ...     data_1: typing.Union[Link, float] = Link(
        ...         ref=wout, fget=lambda wo: wo.data_1
        ...     )
        ...     data_2: typing.Union[Link, float] = Link(
        ...         ref=wout, fget=lambda wo: wo.data_2
        ...     )
        ...     flux_surfaces: typing.Union[Link, float] = Link(value=None)
        ...     dummy: any = "any value"

        >>> equi_explicit = Equilibrium(flux_surfaces=42)
        >>> assert equi_explicit.flux_surfaces == 42

        >>> equi = Equilibrium(
        ...     flux_surfaces=Link(
        ...         "./flux_surfaces_123.txt -> 123",
        ...         lambda ref: int(ref[-3:])
        ...     )
        ... )
        >>> assert equi.flux_surfaces == 123

        Linker attributes reveal the ref with the exposed switch set to True

        >>> with exposed(True):
        ...     assert "./flux_surfaces" in equi.flux_surfaces.ref
        >>> assert equi.flux_surfaces == 123

        You can also use the 'get_link' method to do the same

        >>> assert equi.data_2 == 42
        >>> link = equi.get_link("data_2")
        >>> assert link.ref is wout

        Or in short use the 'get_ref' method to directly access the reference

        >>> ref = equi.get_ref("data_2")
        >>> assert ref is wout

        It raises LinkNotFoundError if no Link instance is found

        >>> equi.get_link("dummy")  # doctest: +ELLIPSIS
        Traceback (most recent call last):
        ...
        rna.pattern.link.LinkNotFoundError: dummy
        >>> equi.dummy
        'any value'

        You can find all linkers refering to the ref by the _references attribute.
        This is only possible, if the referenced object is an instance of Reference.

        >>> assert isinstance(wout, Reference)

        The first two entries in _references stam from equi_explicit, the last two from equi

        >>> assert wout._references[2][0] is equi
        >>> assert wout._references[2][1] == "data_1"

        The cache attribute allows to - if true - store the call value of fget (at the expense
        of not being up to date with the ref). -> Only use the cache=True value if you are
        sure that the reference object will never change the referred value.

        >>> equi_cached = Equilibrium(
        ...     data_2=Link(
        ...         ref=wout, fget=lambda wo: wo.data_2, cache=True
        ...     )
        ... )
        >>> assert equi_cached.data_2 == 42
        >>> wout.data_2 = 12
        >>> assert equi_cached.data_2 == 42
        >>> assert equi.data_2 == 12
    """

    def __getattribute__(self, name):
        val = super().__getattribute__(name)
        if isinstance(val, Link):
            if exposed.enabled():
                return val
            if not val.cache or val.value is NOTSET:
                if val.ref is NOTSET:
                    value = val.fget()
                else:
                    value = val.fget(val.ref)
                if val.cache:
                    val.value = value
            else:
                value = val.value
            return value
        return val

    def __setattr__(self, name, val):
        if isinstance(val, Link):
            if isinstance(val.ref, Reference):
                val.ref._references.append((self, name))
        super().__setattr__(name, val)

    def __repr__(self) -> str:
        with exposed(True):
            repr = super().__repr__()
        return repr

    def get_link(self, attr_name: str) -> Link:
        """
        Returns:
            exposed link under the given attr_name
        """
        with exposed(True):
            val = getattr(self, attr_name)
        if not isinstance(val, Link):
            raise LinkNotFoundError(attr_name)
        return val

    def get_ref(self, attr_name: str) -> Reference:
        """
        Returns:
            reference found under attribute_names
        """
        link = self.get_link(attr_name)
        return link.ref

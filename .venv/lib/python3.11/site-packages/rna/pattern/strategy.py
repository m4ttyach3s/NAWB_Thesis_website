"""
Strategy design pattern, see https://refactoring.guru/design-patterns/strategy

Examples:
    >>> from abc import abstractmethod
    >>> from rna.pattern.strategy import Strategy, Context

    Further user abstraction on top
    >>> class Backend(Strategy):
    ...     @abstractmethod
    ...     def all_backends_should_have_this(self):
    ...         pass

    >>> class Node(Context):
    ...     @property
    ...     def backend(self):
    ...         return self.strategy
    ...
    ...     @backend.setter
    ...     def backend(self, backend):
    ...         self.strategy = backend

    High level description of the still abstact Strategy but the final Context class
    Setting STRATEGY_TYPE defines the Abstract Strategy class every strategy should derive of.
    >>> class MyBackend(Backend):
    ...     @abstractmethod
    ...     def required(self):
    ...         pass
    >>> class MyNode(Node):
    ...     STRATEGY_TYPE = MyBackend
    ...     def compute(self):
    ...         return self.backend.required()

    Finally an implementation of two strategy workers
    >>> class MyBackendVersion1(MyBackend):
    ...     def all_backends_should_have_this(self):
    ...         pass
    ...     def required(self):
    ...         return 42
    >>> class MyBackendVersion2(MyBackend):
    ...     def all_backends_should_have_this(self):
    ...         pass
    ...     def required(self):
    ...         return 21

    >>> a = MyNode(strategy=MyBackendVersion1)
    >>> a.compute()
    42
    >>> a.backend = MyBackendVersion2
    >>> a.compute()
    21

    The backend needs to be a subclass of MyNode.STRATEGY_TYPE
    >>> a.strategy = Backend  # doctest: +ELLIPSIS
    Traceback (most recent call last):
    ...
    TypeError: Strategy must be of type <class '...MyBackend'>, not <class '....Backend'>

    Reset the strategy by setting it to None
    >>> a.backend = None

    Context.strategy automatically generates a strategy from STRATEGY_DEFAULT if set.
    >>> MyNode.STRATEGY_DEFAULT = MyBackendVersion1
    >>> a.compute()
    42
"""
from abc import ABC
import typing
import importlib
import inspect
import sys


class Strategy(ABC):
    """
    The Strategy interface declares operations common to all supported versions of some
    algorithm.

    The Context uses this interface to call the algorithm defined by Concrete Strategies.

    Attributes:
        STRATEGY_MODULE_BASE: If not None given, this is used as a base namespace for the
            from_module method.
        STRATEGY_MODULE_NAME: last part of namespace for the from_module_method.
    """

    STRATEGY_MODULE_BASE: typing.Optional[str] = None  # e.g. "w7x.simulation.backends"
    STRATEGY_MODULE_NAME: typing.Optional[typing.Union[str, int]] = None  # e.g. "flt"

    @classmethod
    def inheritors(cls) -> set:
        """
        Returns:
            flat subclasses and subsub...classes
        """
        subclasses = set()
        work = [cls]
        while work:
            parent = work.pop()
            for child in parent.__subclasses__():
                if child not in subclasses:
                    subclasses.add(child)
                    work.append(child)
        return subclasses

    @classmethod
    def from_module(
        cls, strategy: str, *args, auto_import=True, **kwargs
    ) -> "Strategy":
        """
        Examples:
            >>> from rna.pattern.strategy import Strategy
            >>> backend = Strategy.from_module('rna.pattern.backend')
            >>> from rna.pattern.backend import Backend
            >>> type(backend) is Backend
            True
        """
        # importing the backend should add a subclass of self
        module = ""
        if cls.STRATEGY_MODULE_BASE is not None:
            module = cls.STRATEGY_MODULE_BASE + "."

        module += strategy

        if cls.STRATEGY_MODULE_NAME is not None:
            if isinstance(cls.STRATEGY_MODULE_NAME, str):
                module += "." + cls.STRATEGY_MODULE_NAME
            else:
                if cls.__module__ == "__main__":
                    parts = (
                        sys.modules[cls.__module__].__file__.rstrip(".py").split("/")
                    )
                else:
                    parts = cls.__module__.split(".")

                if isinstance(cls.STRATEGY_MODULE_NAME, int):
                    module += "." + parts[cls.STRATEGY_MODULE_NAME]
                elif isinstance(cls.STRATEGY_MODULE_NAME, slice):
                    module += "." + ".".join(parts[cls.STRATEGY_MODULE_NAME])
        for tpe in cls.inheritors():
            if module in tpe.__module__:
                return tpe(*args, **kwargs)

        if auto_import:
            importlib.import_module(module)
            return cls.from_module(strategy, auto_import=False, *args, **kwargs)
        raise ValueError("No subclass of {cls} in module {module}".format(**locals()))


# pylint: disable=too-few-public-methods
class Context:
    """
    The Context defines the interface of interest to clients.

    Args:
        *args: forwarded to strategy (on strategy initialization time - see :meth:`strategy`)
        **kwargs: forwarded to strategy (on strategy initialization time - see :meth:`strategy`)
    """

    STRATEGY_TYPE: Strategy = None  # base class of all allowed strategies
    STRATEGY_DEFAULT: typing.Union[typing.Type[Strategy], str] = None

    def __init__(
        self,
        *args,
        strategy: typing.Optional[
            typing.Union[Strategy, typing.Type[Strategy], str]
        ] = None,
        **kwargs,
    ) -> None:
        """
        Usually, the Context accepts a strategy through the constructor, but
        also provides a setter to change it at runtime.
        """
        self._strategy_init_values = (args, kwargs)
        self.strategy = strategy

    @property
    def strategy(self) -> Strategy:
        """
        The Context maintains a reference to one of the Strategy objects. The
        Context does not know the concrete class of a strategy. It should work
        with all strategies via the Strategy interface.

        Actually instantiate self._strategy with self._strategy_init_values if self._strategy is
        a type instance or string.

        Returns:
            self._strategy: the newly instantiated strategy attribute.
        """
        (
            args,
            kwargs,
        ) = self._strategy_init_values  # TODO use setters all over here. Why not?
        if self._strategy is None:
            if self.STRATEGY_DEFAULT is None:
                return None
            # pylint:disable=not-callable
            if isinstance(self.STRATEGY_DEFAULT, str):
                self._strategy = self.STRATEGY_TYPE.from_module(
                    self.STRATEGY_DEFAULT, *args, **kwargs
                )
            else:
                self._strategy = self.STRATEGY_DEFAULT(*args, **kwargs)
        elif inspect.isclass(self._strategy) and issubclass(self._strategy, Strategy):
            self._strategy = self._strategy(*args, **kwargs)
        elif isinstance(self._strategy, str):
            self._strategy = self.STRATEGY_TYPE.from_module(
                self._strategy, *args, **kwargs
            )
        return self._strategy

    @strategy.setter
    def strategy(
        self,
        strategy: typing.Union[Strategy, typing.Type[Strategy], str, typing.Type[None]],
    ) -> None:
        """
        Specialty of this setter is that it does also accept uninitialized classes.
        The getter will take care of instantiation. The idea is to not produce any overhead
        until the strategy is actually required.
        """
        if inspect.isclass(strategy) and not issubclass(strategy, self.STRATEGY_TYPE):
            raise TypeError(
                "Strategy must be of type {self.STRATEGY_TYPE}, not {strategy}".format(
                    **locals()
                )
            )
        self._strategy = strategy

    # The below complicates error search too much.
    # def __getattr__(self, attr: str):
    #     """
    #     Forwarding to strategy namespace defined by abstractmethod.

    #     Examples:
    #         The Context class exposes the abstract strategy methods.
    #         >>> a.required()
    #         42
    #     """
    #     if hasattr(self.STRATEGY_TYPE, attr):
    #         value = getattr(self.STRATEGY_TYPE, attr)
    #         if getattr(value, "__isabstractmethod__", False):
    #             # Parent strategy defines this attribute as exposed!
    #             if self.strategy is None:
    #                 raise AttributeError(
    #                     "Request to forward '{self.__class__.__name__}' object attribute {attr}"
    #                     " to un-set strategy.".format(**locals())
    #                 )
    #             return getattr(self.strategy, attr)
    #     raise AttributeError(
    #         "'{self.__class__.__name__}' object has no attribute '{attr}'".format(
    #             **locals()
    #         )
    #     )

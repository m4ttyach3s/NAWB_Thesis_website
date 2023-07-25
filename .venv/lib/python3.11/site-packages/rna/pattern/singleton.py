"""
Strategy design pattern, see https://refactoring.guru/design-patterns/singleton
"""


class SingletonMeta(type):
    """
    Meta class for creation of only one instance
    """

    _instances = {}

    def __call__(cls, *args, **kwargs):
        """
        Possible changes to the value of the `__init__` argument do not affect
        the returned instance.
        """
        if cls not in cls._instances:
            instance = super().__call__(*args, **kwargs)
            cls._instances[cls] = instance
        return cls._instances[cls]


class Singleton(metaclass=SingletonMeta):  # pylint: disable=too-few-public-methods
    """
    The Singleton class can be implemented in different ways in Python. Some
    possible methods include: base class, decorator, metaclass. We will use the
    metaclass because it is best suited for this purpose.

    Examples:
        >>> s1 = Singleton()
        >>> s2 = Singleton()
        >>> assert id(s1) == id(s2)
    """

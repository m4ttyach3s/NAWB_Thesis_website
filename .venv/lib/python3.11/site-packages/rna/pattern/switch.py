class SwitchMeta(type):
    """
    This switch behaves very different from typical classes. You can not instantiate it (calling
    it returns the cls). Everything is stored under ENABLED attribute.
    """

    def __call__(cls, enabled=True):
        """Actually only used for context manager"""
        if not hasattr(cls, "CONTEXT_VALUE"):
            cls.CONTEXT_VALUE = []
        cls.CONTEXT_VALUE.append(enabled)
        return cls

    def __enter__(cls):
        if not hasattr(cls, "MEMORY"):
            cls.MEMORY = []
        cls.MEMORY.append(cls.ENABLED)
        cls.ENABLED = cls.CONTEXT_VALUE.pop(-1)
        return cls

    def __exit__(cls, exc_type, exc_value, exc_traceback):
        cls.ENABLED = cls.MEMORY.pop(-1)


class Switch(metaclass=SwitchMeta):
    """
    Global switch for the use of distributed computing
    """

    ENABLED = True

    @classmethod
    def enabled(cls) -> bool:
        """
        Returns value of this global like variable
        """
        return cls.ENABLED

    @classmethod
    def enable(cls):
        """
        Set value to True
        """
        cls.ENABLED = True

    @classmethod
    def disable(cls):
        """
        Set value to False
        """
        cls.ENABLED = False

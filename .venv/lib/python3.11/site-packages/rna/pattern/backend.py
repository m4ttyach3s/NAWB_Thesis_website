"""
Alias for strategy pattern. Aliasing makes the use case more understandable in some cases.
"""
from rna.pattern import strategy


class Backend(strategy.Strategy):
    """
    Actual Worker
    """


class Frontend(strategy.Context):
    """
    Client side, deligating work to Backends
    """

    @property
    def backend(self):
        """
        Alias for strategy
        """
        return self.strategy

    @backend.setter
    def backend(self, backend):
        self.strategy = backend

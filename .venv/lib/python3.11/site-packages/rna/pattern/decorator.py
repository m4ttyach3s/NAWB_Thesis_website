"""
Decorator base class providing the same behaviour for both

    @Decorator
    def ..

and

    @Decorator(my_option=True)
    def ..

Additionally it can be applied to methods, function, staticmethods, ... without thinking
"""
import abc
import functools
import inspect

# TODO(@dboe): look at https://github.com/mrocklin/multipledispatch/blob/master/multipledispatch/core.py
# Maybe this can replace the currently complicated scheme


class DecoratorMeta(abc.ABCMeta):
    """
    Meta class for creation of only one instance
    """

    def __call__(cls, *args, **kwargs):
        """
        Allow two signatures:
            1)
                @Decorator
                def fun... / class Class...

                -> args = (fun,)
                -> kwargs = {}
            2)
                @Decorator(3, option=24)
                def fun... / class Class...

                -> args = (3,)
                -> kwargs = {'option': 24}

        We subdivide the signatures further depending on the wrapped item:
            x.1) function
                @Decorator[(...)]
                def fun(*args, **kwargs)

            x.2) staticmethod
                class ...
                    @staticmethod
                    @Decorator[(...)]
                    def fun(*args, **kwargs)

            x.3) method
                class ...
                    @Decorator[(...)]
                    def fun(self, *args, **kwargs)

            x.4) classmethod
                class ...
                    @classmethod
                    @Decorator[(...)]
                    def fun(cls, *args, **kwargs)


        Note:
            To pass a single callable to init instead of wrapping it, add an Ellipsis (...)
        """
        func = None
        if args:
            if hasattr(cls, "FLEXIBLE") and not cls.FLEXIBLE:
                # Enforce Signature 2!
                args = args
            elif args[-1] is Ellipsis:
                # Signature 2
                args = args[:-1]
            elif callable(args[0]) and not args[1:] and not kwargs:
                # Signature 1
                func = args[0]
                args = args[1:]
            # Else Signature 2

        instance = cls.__new__(cls, *args, **kwargs)
        instance._func = func  # pylint:disable=attribute-defined-outside-init
        instance._obj = None  # pylint:disable=attribute-defined-outside-init
        instance._wrapped = None  # pylint:disable=attribute-defined-outside-init

        if instance._func is not None:
            # Signature 1
            functools.update_wrapper(instance, func)

        cls.__init__(instance, *args, **kwargs)

        return instance


# pylint:disable=access-member-before-definition,attribute-defined-outside-init
class Decorator(metaclass=DecoratorMeta):
    """
    Abstract decorator base class providing the same behaviour for two signatures:
        1)
            @DecoratorBase
            def fun...

            -> args = (fun,)
            -> kwargs = {}
        2)
            @DecoratorBase(3, option=24)
            def fun...

            -> args = 3
            -> kwargs = {'option': 24}

    Additionally it can be applied to methods, function, staticmethods and classmethods without
    thinking. Additionally functools.wraps is automatically applied to wrapped method. We
    recommend to not use it in addition.

    The Decorator class provides easy wrapping with only :py:meth:`_wrap` to be implemented.
    """

    def __get__(self, obj, obj_type=None):
        # Signature 2: self._func is a method!
        self._obj = obj
        return self

    def __call__(self, *args, **kwargs):
        """
        Note:
            - `__call__` is called at different times for signatures 1 and 2.
            - this is intentioanlly as flat as possible (as opposed to nested calles) in
                order to reduce the depth of the stack trace. Dont be clever.

        """
        #################
        # Determine base signature and complete self._func
        #################

        signature = (1,)
        if self._func is None:
            # Signature 2
            signature = (2,)
            assert len(args) == 1  # possibly remove
            assert not kwargs
            self._func = args[0]
            functools.update_wrapper(self, self._func)

        #################
        # Return Cache
        #################

        if self._wrapped is None:
            #################
            # Determine signature
            #################

            # first time request of wrapped property
            cache = True
            method = False
            method_requires_partial = False
            if self._obj:
                # Signature 2.3: _func is a method
                method = True
                method_requires_partial = True
                cache = False
                # NOTE: The return directly here is crucial because self._obj might change.
                #       If you do not return here, you will always return the method of the
                #       first object requested. That causes horror.
            else:
                inspect_signature = inspect.Signature.from_callable(self._func)
                first_arg_name = next(iter(inspect_signature.parameters), None)
                if first_arg_name in ("self", "cls"):
                    # self: Signature 1.3: _func is a method      but _obj is None
                    # cls:  Signature x.3: _func is a classmethod
                    method = True

            #################
            # Build wrapper
            #################

            if method:

                def wrap(this, *args, **kwargs):
                    func = functools.partial(self._func, this)
                    return self._wrap(this, func, *args, **kwargs)

                if method_requires_partial:
                    wrap = functools.partial(wrap, self._obj)
            else:

                def wrap(*args, **kwargs):
                    return self._wrap(None, self._func, *args, **kwargs)

            wrap = functools.wraps(self._func)(wrap)
            if cache:
                self._wrapped = wrap
        else:
            wrap = self._wrapped

        #################
        # Call
        #################

        if signature[0] == 1:
            return wrap(*args, **kwargs)
        return wrap

    @abc.abstractmethod
    def _wrap(self, this, func, *args, **kwargs):
        """
        This is a generic wrapper method for both methods and classes. Special about this is that
        it needs to have.

        Args:
            this: obj if func is method or classmethod then 'this' corresponds to 'self' or 'cls'
                in the wrapped method correspondingly.
                If this is None, func is a function or staticmethod.
            func: func is the function to be wrapped and will be passed as an explicit kwarg
                NOTE: When calling the function, do not use self as the first attribute even if
                func is a method. This is handled internally.
            *args: Arguments passed to decorated function. Can be passed to func.
            **kwargs: Kwargs passed to decorated function. Can be passed to func.

        Examples:
            >>> from rna.pattern.decorator import Decorator
            >>> # pylint:disable=too-few-public-methods,invalid-name
            >>> class add_args_multiply_kwargs(Decorator):
            ...     def __init__(self, att=None):
            ...         self.att = att
            ...
            ...     def _wrap(
            ...         self, this, func, *args, **kwargs
            ...     ):
            ...         if this is None:  # we are wrapping a function!
            ...             args = args + tuple([sum(args)])
            ...         else:  # we are wrapping a method
            ...             args = args + tuple([this.fast_sum(args)])
            ...         kwargs["res"] = func(*args, **kwargs)
            ...         assert self.att is None or self.att == 42  # we have access to decorator
            ...         return args, kwargs

            >>> class C:
            ...     def fast_sum(self, args):
            ...         return sum(args)
            ...
            ...     @add_args_multiply_kwargs
            ...     def multiply_kwargs(self, *args, **kwargs):
            ...         val = 1
            ...         for v in kwargs.values():
            ...             val *= v
            ...         return val

            >>> @add_args_multiply_kwargs
            ... def multiply_kwargs(*args, **kwargs):
            ...     val = 1
            ...     for v in kwargs.values():
            ...         val *= v
            ...     return val

            >>> C().multiply_kwargs(3, 3, 4, a=2, b=3, c=10)
            ((3, 3, 4, 10), {'a': 2, 'b': 3, 'c': 10, 'res': 60})
            >>> multiply_kwargs(3, 3, 4, a=2, b=3, c=10)
            ((3, 3, 4, 10), {'a': 2, 'b': 3, 'c': 10, 'res': 60})
        """

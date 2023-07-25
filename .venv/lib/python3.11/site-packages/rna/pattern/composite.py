"""
See [here](https://refactoring.guru/design-patterns/composite).

Examples:
    >>> import rna
    >>> from rna.pattern.composite import Component, Leaf, Composite
    >>> class MyLeaf(Leaf):
    ...     def operation(self) -> str:
    ...         return "Leaf"

    >>> class MyComposite(Composite):
    ...     def operation(self) -> str:
    ...         results = []
    ...         for child in self._children:
    ...             results.append(child.operation())
    ...         return f"Branch({'+'.join(results)})"

    The client code works with all of the components via the base interface.
    >>> def client_code(component: Component) -> None:
    ...     print(f"RESULT: {component.operation()}", end="")

    Thanks to the fact that the child-management operations are declared in the
    base Component class, the client code can work with any component, simple or
    complex, without depending on their concrete classes.
    >>> def client_code2(component1: Component, component2: Component) -> None:
    ...     if component1.is_composite():
    ...         component1.add(component2)
    ...     print(f"RESULT: {component1.operation()}", end="")


    This way the client code can support the simple leaf components...
    >>> simple = MyLeaf()

    Client: I've got a simple component:
    >>> client_code(simple)
    RESULT: Leaf

    ...as well as the complex composites.
    >>> tree = MyComposite()

    >>> branch1 = MyComposite()
    >>> branch1.add(MyLeaf())
    >>> branch1.add(MyLeaf())

    >>> branch2 = MyComposite()
    >>> branch2.add(MyLeaf())

    >>> tree.add(branch1)
    >>> tree.add(branch2)

    Client: Now I've got a composite tree:
    >>> client_code(tree)
    RESULT: Branch(Branch(Leaf+Leaf)+Branch(Leaf))

    Client: I don't need to check the components classes even when managing the tree:
    >>> client_code2(tree, simple)
    RESULT: Branch(Branch(Leaf+Leaf)+Branch(Leaf)+Leaf)

    This pattern is generically implemented such that you can subclass Component and then
    auto-generate Leave and Composite classes from this.
    >>> from abc import abstractmethod
    >>> class Graphic(Component):
    ...     @abstractmethod
    ...     def render(self) -> None:
    ...         pass
    >>> class CompositeGraphic(Graphic.Composite):
    ...     def render(self):
    ...         for child in self._children:
    ...             child.render()
    >>> class Ellipse(Graphic.Leaf):
    ...     def __init__(self, name):
    ...         self.name = name
    ...         super().__init__()
    ...
    ...     def render(self):
    ...         print("Ellipse", self.name)
    >>> ellipse1 = Ellipse("1")
    >>> ellipse2 = Ellipse("2")
    >>> ellipse3 = Ellipse("3")
    >>> ellipse4 = Ellipse("4")
    >>> graphic1 = CompositeGraphic()
    >>> graphic2 = CompositeGraphic()
    >>> graphic1.add(ellipse1)
    >>> graphic1.add(ellipse2)
    >>> graphic1.add(ellipse3)
    >>> graphic2.add(ellipse4)
    >>> graphic = CompositeGraphic()
    >>> graphic.add(graphic1)
    >>> graphic.add(graphic2)
    >>> graphic.render()
    Ellipse 1
    Ellipse 2
    Ellipse 3
    Ellipse 4

    The parent attribute of a leaf is a reference to the composite instance one level above in the
    tree.
    >>> graphic_parent = CompositeGraphic()
    >>> child_1 = Ellipse("1")
    >>> child_2 = Ellipse("1")

    Initially Composite and Leaf instances' parents are None
    >>> assert graphic_parent.parent is None
    >>> assert child_1.parent is None

    When adding the leaf to the composite their relationship is set and the leaf has a parent
    >>> graphic_parent.add(child_1)
    >>> assert child_1.parent is graphic_parent
    >>> assert graphic_parent.parent is None

    This allows vertical naviagtion in the tree:
    >>> graphic_parent.add(child_2)
    >>> assert child_1.parent.get_leaves()[1] is child_2
"""
import typing
from abc import ABCMeta, abstractmethod


def get_leaves(
    children: typing.List[typing.Union["Leaf", "Composite"]],
    flat: bool = False,
    apply: typing.Optional[typing.Callable] = None,
    apply_composite: typing.Optional[typing.Callable] = None,
    _position: typing.Tuple[int] = (),
    _index: int = 0,
    **kwargs,
):
    """
    Recursively run through children and join component leafs

    Args:
        children: Composite tree
        flat: flatten nested sub components
        apply: function to apply to a leaf before appending to the results
        apply_composite: function to apply to a list of leaves before appending / extending to
            results
        _position: recursion argument. position of leaf/composite from the top of the tree
        _index: recursion argument. Index of last leaf in the flat tree.
        **kwargs: passed down to both apply and apply_composite

    Examples:
        >>> import rna.pattern.composite

        Build a nested tree
        >>> class L(rna.pattern.composite.Leaf):
        ...     pass
        >>> structure = rna.pattern.composite.Composite()
        >>> structure.add(L())
        >>> structure.add(L())
        >>> structure.add(L())
        >>> mid_structure = rna.pattern.composite.Composite()
        >>> mid_structure.add(L())
        >>> mid_structure.add(L())
        >>> inner_structure = rna.pattern.composite.Composite()
        >>> inner_structure.add(L())
        >>> mid_structure.add(inner_structure)
        >>> structure.add(mid_structure)

        Get the flat list of leafs
        >>> inner_structure.get_leaves()  # Doctest:+ELLIPSIS
        [<rna.pattern.composite.L object at 0x...>]
        >>> mid_structure.get_leaves(apply=lambda leaf, **kwargs: "L")
        ['L', 'L', ['L']]
        >>> mid_structure.get_leaves(flat=True, apply=lambda leaf, **kwargs: "L")
        ['L', 'L', 'L']
        >>> structure.get_leaves(apply=lambda leaf, **kwargs: "L")
        ['L', 'L', 'L', ['L', 'L', ['L']]]
        >>> structure.get_leaves(flat=True, apply=lambda leaf, **kwargs: "L")
        ['L', 'L', 'L', 'L', 'L', 'L']

        Since the structure of this function is occuring often slightly modified, you can
        use the apply function (as done before for the test):
        >>> def apply(leaf, **kwargs):
        ...     index = kwargs.get("index")
        ...     return index
        >>> structure.get_leaves(apply=apply)
        [0, 1, 2, [3, 4, [5]]]

        You have acces to the _index kwarg (see above) and the _position kwarg:
        >>> def apply(leaf, **kwargs):
        ...     pos = kwargs.get("position")
        ...     return pos
        >>> structure.get_leaves(apply=apply)
        [(0,), (1,), (2,), [(3, 0), (3, 1), [(3, 2, 0)]]]

        The apply_composite function can be used to manipulate composites
        >>> def apply_composite(leaves, **kwargs):
        ...     comp = rna.pattern.composite.Composite()
        ...     for l in leaves:
        ...         comp.add(l)
        ...     return comp
        >>> replica = rna.pattern.composite.Composite()
        >>> for thing in structure.get_leaves(apply_composite=apply_composite):
        ...     replica.add(thing)
        >>> replica.get_leaves(apply=apply)
        [(0,), (1,), (2,), [(3, 0), (3, 1), [(3, 2, 0)]]]
    """
    leaves = []
    for i, component in enumerate(children):
        _curr_position = _position + (i,)
        if component.is_composite():
            lfs = component.get_leaves(
                flat=flat,
                apply=apply,
                apply_composite=apply_composite,
                _position=_curr_position,
                _index=_index,
                **kwargs,
            )
            if apply_composite is not None:
                lfs = apply_composite(
                    lfs,
                    component=component,
                    flat=flat,
                    position=_curr_position,
                    index=_index,
                    **kwargs,
                )
            if flat:
                leaves.extend(lfs)
            else:
                leaves.append(lfs)
        else:
            if apply is not None:
                component = apply(
                    component, position=_curr_position, index=_index, **kwargs
                )
            leaves.append(component)
            _index += 1
    return leaves


class ComponentMeta(ABCMeta):
    """
    Meta class exposing class attributes as mixin generators
    """

    # pylint: disable=invalid-name
    def Subclass(cls, name: str, mixins: typing.Type, namespace: dict = None):
        """
        Auto-generate a sub class from this class together with mixins
        Args:
            name: name of new class
        """
        if namespace is None:
            namespace = {}
        return type(name, (cls,) + mixins, namespace)

    @property
    def Leaf(cls) -> typing.Type:  # pylint: disable=invalid-name
        """
        Auto-generate a Leaf mixin from this class.
        The Leaf class represents the end objects of a composition. A leaf can't
        have any children.  Usually, it's the Leaf objects that do the actual work, whereas
        Composite objects only delegate to their sub-components.
        """

        def __init__(self, *args, **kwargs) -> None:
            cls.__init__(self, *args, **kwargs)
            # super(cls, self).__init__(*args, **kwargs) this does not work

        def add(self, component: cls):  # pylint: disable=no-self-use
            """
            Leafes can not add
            """
            raise TypeError(f"Trying to add component {component} to Leaf")

        def remove(self, component: cls):  # pylint: disable=no-self-use
            """
            Leafes can not remove
            """
            raise TypeError(f"Trying to remove component {component} from Leaf")

        return cls.Subclass("Leaf", (), dict(__init__=__init__, add=add, remove=remove))

    @property
    def Composite(cls) -> typing.Type:  # pylint: disable=invalid-name
        """
        Auto-generate a Composite subclass from this class.
        The Composite class represents the complex components that may have children. Usually, the
        Composite objects delegate the actual work to their children and then "sums-up" the result.
        """

        def __init__(self, *args, **kwargs) -> None:
            self._children: typing.List[cls] = []
            cls.__init__(self, *args, **kwargs)
            # super(cls, self).__init__(*args, **kwargs) this does not work

        def add(self, component: cls) -> None:
            """
            Add component by appending it to children
            """
            self._children.append(component)  # pylint:disable=protected-access
            component.parent = self

        def remove(self, component: cls) -> None:
            """
            Remove component by appending it to children
            """
            self._children.remove(component)  # pylint:disable=protected-access
            component.parent = None

        def is_composite(self) -> bool:  # pylint: disable=unused-argument
            """
            Whether a component can bear children.
            """
            return True

        def get_leaves_(self, *args, **kwargs) -> bool:
            """
            Forwarding to :meth:get_leaves
            """
            return get_leaves(
                self._children, *args, **kwargs
            )  # pylint:disable=protected-access

        return cls.Subclass(
            "Composite",
            (),
            dict(
                __init__=__init__,
                add=add,
                remove=remove,
                is_composite=is_composite,
                get_leaves=get_leaves_,
            ),
        )


class Component(metaclass=ComponentMeta):
    """
    The base Component class declares common operations for both simple and complex objects of a
    composition.
    """

    def __init__(self):
        self._parent: typing.Optional["Component"] = None
        super().__init__()

    @property
    def parent(self) -> typing.Optional["Component"]:
        """
        Returns:
            parent component. None, if self is the root
        """
        return self._parent

    @parent.setter
    def parent(self, parent: "Component"):
        """
        Optionally, the base Component can declare an interface for setting and
        accessing a parent of the component in a tree structure. It can also
        provide some default implementation for these methods.
        """
        self._parent = parent

    def get_root(self) -> "Composite":
        """
        Return the root (i.e. the top node) of the composition
        """
        root = self
        while root.parent is not None:
            root = root.parent
        return root

    # In some cases, it is beneficial to define the child-management
    # operations right in the base Component class. This way, you won't need to
    # expose any concrete component classes to the client code, even during the
    # object tree assembly.

    @abstractmethod
    def add(self, component: "Component") -> None:
        """
        Polymorphism for adding components to a subclass.
        """

    @abstractmethod
    def remove(self, component: "Component") -> None:
        """
        Polymorphism for removing components to a subclass.
        """

    def is_composite(self) -> bool:  # pylint: disable=no-self-use
        """
        Whether a component can bear children.
        """
        return False


class Leaf(Component.Leaf):
    """
    Implementation for pattern namespace. If you want to subclass Component, use auto-generated
    Component.Leaf type
    """


class Composite(Component.Composite):
    """
    Implementation for pattern namespace. If you want to subclass Component, use auto-generated
    Component.Composite type
    """

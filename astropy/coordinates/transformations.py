# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""
This module contains a general framework for defining graphs of transformations
between coordinates, suitable for either spatial coordinates or more generalized
coordinate systems.

The fundamental idea is that each class is a node in the transformation graph,
and transitions from one node to another are defined as functions (or methods)
wrapped in transformation objects.

This module also includes more specific transformation classes for
celestial/spatial coordinate frames, generally focused around matrix-style
transformations that are typically how the algorithms are defined.
"""


import heapq
import inspect
import subprocess
from warnings import warn

from abc import ABCMeta, abstractmethod
from collections import defaultdict
from contextlib import suppress, contextmanager
from inspect import signature

import numpy as np

from astropy import units as u
from astropy.utils.exceptions import AstropyWarning

from .matrix_utilities import matrix_product

__all__ = ['TransformGraph', 'CoordinateTransform', 'FunctionTransform',
           'BaseAffineTransform', 'AffineTransform',
           'StaticMatrixTransform', 'DynamicMatrixTransform',
           'FunctionTransformWithFiniteDifference', 'CompositeTransform']


def frame_attrs_from_set(frame_set):
    """
    A `dict` of all the attributes of all frame classes in this
    `TransformGraph`.

    Broken out of the class so this can be called on a temporary frame set to
    validate new additions to the transform graph before actually adding them.
    """
    result = {}

    for frame_cls in frame_set:
        result.update(frame_cls.frame_attributes)

    return result


def frame_comps_from_set(frame_set):
    """
    A `set` of all component names every defined within any frame class in
    this `TransformGraph`.

    Broken out of the class so this can be called on a temporary frame set to
    validate new additions to the transform graph before actually adding them.
    """
    result = set()

    for frame_cls in frame_set:
        rep_info = frame_cls._frame_specific_representation_info
        for mappings in rep_info.values():
            for rep_map in mappings:
                result.update([rep_map.framename])

    return result


class TransformGraph:
    """
    A graph representing the paths between coordinate frames.
    """

    def __init__(self):
        self._graph = defaultdict(dict)
        self.invalidate_cache()  # generates cache entries

    @property
    def _cached_names(self):
        if self._cached_names_dct is None:
            self._cached_names_dct = dct = {}
            for c in self.frame_set:
                nm = getattr(c, 'name', None)
                if nm is not None:
                    if not isinstance(nm, list):
                        nm = [nm]
                    for name in nm:
                        dct[name] = c

        return self._cached_names_dct

    @property
    def frame_set(self):
        """
        A `set` of all the frame classes present in this `TransformGraph`.
        """
        if self._cached_frame_set is None:
            self._cached_frame_set = set()
            for a in self._graph:
                self._cached_frame_set.add(a)
                for b in self._graph[a]:
                    self._cached_frame_set.add(b)

        return self._cached_frame_set.copy()

    @property
    def frame_attributes(self):
        """
        A `dict` of all the attributes of all frame classes in this
        `TransformGraph`.
        """
        if self._cached_frame_attributes is None:
            self._cached_frame_attributes = frame_attrs_from_set(self.frame_set)

        return self._cached_frame_attributes

    @property
    def frame_component_names(self):
        """
        A `set` of all component names every defined within any frame class in
        this `TransformGraph`.
        """
        if self._cached_component_names is None:
            self._cached_component_names = frame_comps_from_set(self.frame_set)

        return self._cached_component_names

    def invalidate_cache(self):
        """
        Invalidates the cache that stores optimizations for traversing the
        transform graph.  This is called automatically when transforms
        are added or removed, but will need to be called manually if
        weights on transforms are modified inplace.
        """
        self._cached_names_dct = None
        self._cached_frame_set = None
        self._cached_frame_attributes = None
        self._cached_component_names = None
        self._shortestpaths = {}
        self._composite_cache = {}

    def add_transform(self, fromsys, tosys, transform):
        """
        Add a new coordinate transformation to the graph.

        Parameters
        ----------
        fromsys : class
            The coordinate frame class to start from.
        tosys : class
            The coordinate frame class to transform into.
        transform : `CoordinateTransform`
            The transformation object. Typically a `CoordinateTransform` object,
            although it may be some other callable that is called with the same
            signature.

        Raises
        ------
        TypeError
            If ``fromsys`` or ``tosys`` are not classes or ``transform`` is
            not callable.
        """

        if not inspect.isclass(fromsys):
            raise TypeError('fromsys must be a class')
        if not inspect.isclass(tosys):
            raise TypeError('tosys must be a class')
        if not callable(transform):
            raise TypeError('transform must be callable')

        frame_set = self.frame_set.copy()
        frame_set.add(fromsys)
        frame_set.add(tosys)

        # Now we check to see if any attributes on the proposed frames override
        # *any* component names, which we can't allow for some of the logic in
        # the SkyCoord initializer to work
        attrs = set(frame_attrs_from_set(frame_set).keys())
        comps = frame_comps_from_set(frame_set)

        invalid_attrs = attrs.intersection(comps)
        if invalid_attrs:
            invalid_frames = set()
            for attr in invalid_attrs:
                if attr in fromsys.frame_attributes:
                    invalid_frames.update([fromsys])

                if attr in tosys.frame_attributes:
                    invalid_frames.update([tosys])

            raise ValueError("Frame(s) {} contain invalid attribute names: {}"
                             "\nFrame attributes can not conflict with *any* of"
                             " the frame data component names (see"
                             " `frame_transform_graph.frame_component_names`)."
                             .format(list(invalid_frames), invalid_attrs))

        self._graph[fromsys][tosys] = transform
        self.invalidate_cache()

    def remove_transform(self, fromsys, tosys, transform):
        """
        Removes a coordinate transform from the graph.

        Parameters
        ----------
        fromsys : class or None
            The coordinate frame *class* to start from. If `None`,
            ``transform`` will be searched for and removed (``tosys`` must
            also be `None`).
        tosys : class or None
            The coordinate frame *class* to transform into. If `None`,
            ``transform`` will be searched for and removed (``fromsys`` must
            also be `None`).
        transform : callable or None
            The transformation object to be removed or `None`.  If `None`
            and ``tosys`` and ``fromsys`` are supplied, there will be no
            check to ensure the correct object is removed.
        """
        if fromsys is None or tosys is None:
            if not (tosys is None and fromsys is None):
                raise ValueError('fromsys and tosys must both be None if either are')
            if transform is None:
                raise ValueError('cannot give all Nones to remove_transform')

            # search for the requested transform by brute force and remove it
            for a in self._graph:
                agraph = self._graph[a]
                for b in agraph:
                    if agraph[b] is transform:
                        del agraph[b]
                        fromsys = a
                        break

                # If the transform was found, need to break out of the outer for loop too
                if fromsys:
                    break
            else:
                raise ValueError(f'Could not find transform {transform} in the graph')

        else:
            if transform is None:
                self._graph[fromsys].pop(tosys, None)
            else:
                curr = self._graph[fromsys].get(tosys, None)
                if curr is transform:
                    self._graph[fromsys].pop(tosys)
                else:
                    raise ValueError('Current transform from {} to {} is not '
                                     '{}'.format(fromsys, tosys, transform))

        # Remove the subgraph if it is now empty
        if self._graph[fromsys] == {}:
            self._graph.pop(fromsys)

        self.invalidate_cache()

    def find_shortest_path(self, fromsys, tosys):
        """
        Computes the shortest distance along the transform graph from
        one system to another.

        Parameters
        ----------
        fromsys : class
            The coordinate frame class to start from.
        tosys : class
            The coordinate frame class to transform into.

        Returns
        -------
        path : list of class or None
            The path from ``fromsys`` to ``tosys`` as an in-order sequence
            of classes.  This list includes *both* ``fromsys`` and
            ``tosys``. Is `None` if there is no possible path.
        distance : float or int
            The total distance/priority from ``fromsys`` to ``tosys``.  If
            priorities are not set this is the number of transforms
            needed. Is ``inf`` if there is no possible path.
        """

        inf = float('inf')

        # special-case the 0 or 1-path
        if tosys is fromsys:
            if tosys not in self._graph[fromsys]:
                # Means there's no transform necessary to go from it to itself.
                return [tosys], 0
        if tosys in self._graph[fromsys]:
            # this will also catch the case where tosys is fromsys, but has
            # a defined transform.
            t = self._graph[fromsys][tosys]
            return [fromsys, tosys], float(t.priority if hasattr(t, 'priority') else 1)

        # otherwise, need to construct the path:

        if fromsys in self._shortestpaths:
            # already have a cached result
            fpaths = self._shortestpaths[fromsys]
            if tosys in fpaths:
                return fpaths[tosys]
            else:
                return None, inf

        # use Dijkstra's algorithm to find shortest path in all other cases

        nodes = []
        # first make the list of nodes
        for a in self._graph:
            if a not in nodes:
                nodes.append(a)
            for b in self._graph[a]:
                if b not in nodes:
                    nodes.append(b)

        if fromsys not in nodes or tosys not in nodes:
            # fromsys or tosys are isolated or not registered, so there's
            # certainly no way to get from one to the other
            return None, inf

        edgeweights = {}
        # construct another graph that is a dict of dicts of priorities
        # (used as edge weights in Dijkstra's algorithm)
        for a in self._graph:
            edgeweights[a] = aew = {}
            agraph = self._graph[a]
            for b in agraph:
                aew[b] = float(agraph[b].priority if hasattr(agraph[b], 'priority') else 1)

        # entries in q are [distance, count, nodeobj, pathlist]
        # count is needed because in py 3.x, tie-breaking fails on the nodes.
        # this way, insertion order is preserved if the weights are the same
        q = [[inf, i, n, []] for i, n in enumerate(nodes) if n is not fromsys]
        q.insert(0, [0, -1, fromsys, []])

        # this dict will store the distance to node from ``fromsys`` and the path
        result = {}

        # definitely starts as a valid heap because of the insert line; from the
        # node to itself is always the shortest distance
        while len(q) > 0:
            d, orderi, n, path = heapq.heappop(q)

            if d == inf:
                # everything left is unreachable from fromsys, just copy them to
                # the results and jump out of the loop
                result[n] = (None, d)
                for d, orderi, n, path in q:
                    result[n] = (None, d)
                break
            else:
                result[n] = (path, d)
                path.append(n)
                if n not in edgeweights:
                    # this is a system that can be transformed to, but not from.
                    continue
                for n2 in edgeweights[n]:
                    if n2 not in result:  # already visited
                        # find where n2 is in the heap
                        for i in range(len(q)):
                            if q[i][2] == n2:
                                break
                        else:
                            raise ValueError('n2 not in heap - this should be impossible!')

                        newd = d + edgeweights[n][n2]
                        if newd < q[i][0]:
                            q[i][0] = newd
                            q[i][3] = list(path)
                            heapq.heapify(q)

        # cache for later use
        self._shortestpaths[fromsys] = result
        return result[tosys]

    def get_transform(self, fromsys, tosys):
        """
        Generates and returns the `CompositeTransform` for a transformation
        between two coordinate systems.

        Parameters
        ----------
        fromsys : class
            The coordinate frame class to start from.
        tosys : class
            The coordinate frame class to transform into.

        Returns
        -------
        trans : `CompositeTransform` or None
            If there is a path from ``fromsys`` to ``tosys``, this is a
            transform object for that path.   If no path could be found, this is
            `None`.

        Notes
        -----
        This function always returns a `CompositeTransform`, because
        `CompositeTransform` is slightly more adaptable in the way it can be
        called than other transform classes. Specifically, it takes care of
        intermediate steps of transformations in a way that is consistent with
        1-hop transformations.

        """
        if not inspect.isclass(fromsys):
            raise TypeError('fromsys is not a class')
        if not inspect.isclass(tosys):
            raise TypeError('tosys is not a class')

        path, distance = self.find_shortest_path(fromsys, tosys)

        if path is None:
            return None

        transforms = []
        currsys = fromsys
        for p in path[1:]:  # first element is fromsys so we skip it
            transforms.append(self._graph[currsys][p])
            currsys = p

        fttuple = (fromsys, tosys)
        if fttuple not in self._composite_cache:
            comptrans = CompositeTransform(transforms, fromsys, tosys,
                                           register_graph=False)
            self._composite_cache[fttuple] = comptrans
        return self._composite_cache[fttuple]

    def lookup_name(self, name):
        """
        Tries to locate the coordinate class with the provided alias.

        Parameters
        ----------
        name : str
            The alias to look up.

        Returns
        -------
        `BaseCoordinateFrame` subclass
            The coordinate class corresponding to the ``name`` or `None` if
            no such class exists.
        """

        return self._cached_names.get(name, None)

    def get_names(self):
        """
        Returns all available transform names. They will all be
        valid arguments to `lookup_name`.

        Returns
        -------
        nms : list
            The aliases for coordinate systems.
        """
        return list(self._cached_names.keys())

    def to_dot_graph(self, priorities=True, addnodes=[], savefn=None,
                     savelayout='plain', saveformat=None, color_edges=True):
        """
        Converts this transform graph to the graphviz_ DOT format.

        Optionally saves it (requires `graphviz`_ be installed and on your path).

        .. _graphviz: http://www.graphviz.org/

        Parameters
        ----------
        priorities : bool
            If `True`, show the priority values for each transform.  Otherwise,
            the will not be included in the graph.
        addnodes : sequence of str
            Additional coordinate systems to add (this can include systems
            already in the transform graph, but they will only appear once).
        savefn : None or str
            The file name to save this graph to or `None` to not save
            to a file.
        savelayout : str
            The graphviz program to use to layout the graph (see
            graphviz_ for details) or 'plain' to just save the DOT graph
            content. Ignored if ``savefn`` is `None`.
        saveformat : str
            The graphviz output format. (e.g. the ``-Txxx`` option for
            the command line program - see graphviz docs for details).
            Ignored if ``savefn`` is `None`.
        color_edges : bool
            Color the edges between two nodes (frames) based on the type of
            transform. ``FunctionTransform``: red, ``StaticMatrixTransform``:
            blue, ``DynamicMatrixTransform``: green.

        Returns
        -------
        dotgraph : str
            A string with the DOT format graph.
        """

        nodes = []
        # find the node names
        for a in self._graph:
            if a not in nodes:
                nodes.append(a)
            for b in self._graph[a]:
                if b not in nodes:
                    nodes.append(b)
        for node in addnodes:
            if node not in nodes:
                nodes.append(node)
        nodenames = []
        invclsaliases = dict([(f, [k for k, v in self._cached_names.items() if v == f])
                              for f in self.frame_set])
        for n in nodes:
            if n in invclsaliases:
                aliases = '`\\n`'.join(invclsaliases[n])
                nodenames.append('{0} [shape=oval label="{0}\\n`{1}`"]'.format(n.__name__, aliases))
            else:
                nodenames.append(n.__name__ + '[ shape=oval ]')

        edgenames = []
        # Now the edges
        for a in self._graph:
            agraph = self._graph[a]
            for b in agraph:
                transform = agraph[b]
                pri = transform.priority if hasattr(transform, 'priority') else 1
                color = trans_to_color[transform.__class__] if color_edges else 'black'
                edgenames.append((a.__name__, b.__name__, pri, color))

        # generate simple dot format graph
        lines = ['digraph AstropyCoordinateTransformGraph {']
        lines.append('graph [rankdir=LR]')
        lines.append('; '.join(nodenames) + ';')
        for enm1, enm2, weights, color in edgenames:
            labelstr_fmt = '[ {0} {1} ]'

            if priorities:
                priority_part = f'label = "{weights}"'
            else:
                priority_part = ''

            color_part = f'color = "{color}"'

            labelstr = labelstr_fmt.format(priority_part, color_part)
            lines.append(f'{enm1} -> {enm2}{labelstr};')

        lines.append('')
        lines.append('overlap=false')
        lines.append('}')
        dotgraph = '\n'.join(lines)

        if savefn is not None:
            if savelayout == 'plain':
                with open(savefn, 'w') as f:
                    f.write(dotgraph)
            else:
                args = [savelayout]
                if saveformat is not None:
                    args.append('-T' + saveformat)
                proc = subprocess.Popen(args, stdin=subprocess.PIPE,
                                        stdout=subprocess.PIPE,
                                        stderr=subprocess.PIPE)
                stdout, stderr = proc.communicate(dotgraph)
                if proc.returncode != 0:
                    raise OSError('problem running graphviz: \n' + stderr)

                with open(savefn, 'w') as f:
                    f.write(stdout)

        return dotgraph

    def to_networkx_graph(self):
        """
        Converts this transform graph into a networkx graph.

        .. note::
            You must have the `networkx <https://networkx.github.io/>`_
            package installed for this to work.

        Returns
        -------
        nxgraph : ``networkx.Graph``
            This `TransformGraph` as a `networkx.Graph <https://networkx.github.io/documentation/stable/reference/classes/graph.html>`_.
        """
        import networkx as nx

        nxgraph = nx.Graph()

        # first make the nodes
        for a in self._graph:
            if a not in nxgraph:
                nxgraph.add_node(a)
            for b in self._graph[a]:
                if b not in nxgraph:
                    nxgraph.add_node(b)

        # Now the edges
        for a in self._graph:
            agraph = self._graph[a]
            for b in agraph:
                transform = agraph[b]
                pri = transform.priority if hasattr(transform, 'priority') else 1
                color = trans_to_color[transform.__class__]
                nxgraph.add_edge(a, b, weight=pri, color=color)

        return nxgraph

    def transform(self, transcls, fromsys, tosys, priority=1, **kwargs):
        """
        A function decorator for defining transformations.

        .. note::
            If decorating a static method of a class, ``@staticmethod``
            should be  added *above* this decorator.

        Parameters
        ----------
        transcls : class
            The class of the transformation object to create.
        fromsys : class
            The coordinate frame class to start from.
        tosys : class
            The coordinate frame class to transform into.
        priority : float or int
            The priority if this transform when finding the shortest
            coordinate transform path - large numbers are lower priorities.

        Additional keyword arguments are passed into the ``transcls``
        constructor.

        Returns
        -------
        deco : function
            A function that can be called on another function as a decorator
            (see example).

        Notes
        -----
        This decorator assumes the first argument of the ``transcls``
        initializer accepts a callable, and that the second and third
        are ``fromsys`` and ``tosys``. If this is not true, you should just
        initialize the class manually and use `add_transform` instead of
        using this decorator.

        Examples
        --------
        ::

            graph = TransformGraph()

            class Frame1(BaseCoordinateFrame):
               ...

            class Frame2(BaseCoordinateFrame):
                ...

            @graph.transform(FunctionTransform, Frame1, Frame2)
            def f1_to_f2(f1_obj):
                ... do something with f1_obj ...
                return f2_obj

        """
        def deco(func):
            # this doesn't do anything directly with the transform because
            # ``register_graph=self`` stores it in the transform graph
            # automatically
            transcls(func, fromsys, tosys, priority=priority,
                     register_graph=self, **kwargs)
            return func
        return deco

    def _add_merged_transform(self, fromsys, tosys, *furthersys, priority=1):
        """
        Add a single-step transform that encapsulates a multi-step transformation path,
        using the transforms that already exist in the graph.

        The created transform internally calls the existing transforms.  If all of the
        transforms are affine, the merged transform is
        `~astropy.coordinates.transformations.DynamicMatrixTransform` (if there are no
        origin shifts) or `~astropy.coordinates.transformations.AffineTransform`
        (otherwise).  If at least one of the transforms is not affine, the merged
        transform is
        `~astropy.coordinates.transformations.FunctionTransformWithFiniteDifference`.

        This method is primarily useful for defining loopback transformations
        (i.e., where ``fromsys`` and the final ``tosys`` are the same).

        Parameters
        ----------
        fromsys : class
            The coordinate frame class to start from.
        tosys : class
            The coordinate frame class to transform to.
        *furthersys : class
            Additional coordinate frame classes to transform to in order.
        priority : number
            The priority of this transform when finding the shortest
            coordinate transform path - large numbers are lower priorities.

        Notes
        -----
        Even though the created transform is a single step in the graph, it
        will still internally call the constituent transforms.  Thus, there is
        no performance benefit for using this created transform.

        For Astropy's built-in frames, loopback transformations typically use
        `~astropy.coordinates.ICRS` to be safe.  Tranforming through an inertial
        frame ensures that changes in observation time and observer
        location/velocity are properly accounted for.

        An error will be raised if a direct transform between ``fromsys`` and
        ``tosys`` already exist.
        """
        frames = [fromsys, tosys, *furthersys]
        lastsys = frames[-1]
        full_path = self.get_transform(fromsys, lastsys)
        transforms = [self.get_transform(frame_a, frame_b)
                      for frame_a, frame_b in zip(frames[:-1], frames[1:])]
        if None in transforms:
            raise ValueError(f"This transformation path is not possible")
        if len(full_path.transforms) == 1:
            raise ValueError(f"A direct transform for {fromsys.__name__}->{lastsys.__name__} already exists")

        self.add_transform(fromsys, lastsys,
                           CompositeTransform(transforms, fromsys, lastsys,
                                              priority=priority)._as_single_transform())

    @contextmanager
    def impose_finite_difference_dt(self, dt):
        """
        Context manager to impose a finite-difference time step on all applicable transformations

        For each transformation in this transformation graph that has the attribute
        ``finite_difference_dt``, that attribute is set to the provided value.  The only standard
        transformation with this attribute is
        `~astropy.coordinates.transformations.FunctionTransformWithFiniteDifference`.

        Parameters
        ----------
        dt : `~astropy.units.Quantity` ['time'] or callable
            If a quantity, this is the size of the differential used to do the finite difference.
            If a callable, should accept ``(fromcoord, toframe)`` and return the ``dt`` value.
        """
        key = 'finite_difference_dt'
        saved_settings = []

        try:
            for to_frames in self._graph.values():
                for transform in to_frames.values():
                    if hasattr(transform, key):
                        old_setting = (transform, key, getattr(transform, key))
                        saved_settings.append(old_setting)
                        setattr(transform, key, dt)
            yield
        finally:
            for setting in saved_settings:
                setattr(*setting)


# <-------------------Define the builtin transform classes-------------------->

class CoordinateTransform(metaclass=ABCMeta):
    """
    An object that transforms a coordinate from one system to another.
    Subclasses must implement `__call__` with the provided signature.
    They should also call this superclass's ``__init__`` in their
    ``__init__``.

    Parameters
    ----------
    fromsys : `~astropy.coordinates.BaseCoordinateFrame` subclass
        The coordinate frame class to start from.
    tosys : `~astropy.coordinates.BaseCoordinateFrame` subclass
        The coordinate frame class to transform into.
    priority : float or int
        The priority if this transform when finding the shortest
        coordinate transform path - large numbers are lower priorities.
    register_graph : `TransformGraph` or None
        A graph to register this transformation with on creation, or
        `None` to leave it unregistered.
    """

    def __init__(self, fromsys, tosys, priority=1, register_graph=None):
        if not inspect.isclass(fromsys):
            raise TypeError('fromsys must be a class')
        if not inspect.isclass(tosys):
            raise TypeError('tosys must be a class')

        self.fromsys = fromsys
        self.tosys = tosys
        self.priority = float(priority)

        if register_graph:
            # this will do the type-checking when it adds to the graph
            self.register(register_graph)
        else:
            if not inspect.isclass(fromsys) or not inspect.isclass(tosys):
                raise TypeError('fromsys and tosys must be classes')

        self.overlapping_frame_attr_names = overlap = []
        if (hasattr(fromsys, 'get_frame_attr_names') and
                hasattr(tosys, 'get_frame_attr_names')):
            # the if statement is there so that non-frame things might be usable
            # if it makes sense
            for from_nm in fromsys.frame_attributes.keys():
                if from_nm in tosys.frame_attributes.keys():
                    overlap.append(from_nm)

    def register(self, graph):
        """
        Add this transformation to the requested Transformation graph,
        replacing anything already connecting these two coordinates.

        Parameters
        ----------
        graph : `TransformGraph` object
            The graph to register this transformation with.
        """
        graph.add_transform(self.fromsys, self.tosys, self)

    def unregister(self, graph):
        """
        Remove this transformation from the requested transformation
        graph.

        Parameters
        ----------
        graph : a TransformGraph object
            The graph to unregister this transformation from.

        Raises
        ------
        ValueError
            If this is not currently in the transform graph.
        """
        graph.remove_transform(self.fromsys, self.tosys, self)

    @abstractmethod
    def __call__(self, fromcoord, toframe):
        """
        Does the actual coordinate transformation from the ``fromsys`` class to
        the ``tosys`` class.

        Parameters
        ----------
        fromcoord : `~astropy.coordinates.BaseCoordinateFrame` subclass instance
            An object of class matching ``fromsys`` that is to be transformed.
        toframe : object
            An object that has the attributes necessary to fully specify the
            frame.  That is, it must have attributes with names that match the
            keys of the dictionary that ``tosys.get_frame_attr_names()``
            returns. Typically this is of class ``tosys``, but it *might* be
            some other class as long as it has the appropriate attributes.

        Returns
        -------
        tocoord : `BaseCoordinateFrame` subclass instance
            The new coordinate after the transform has been applied.
        """


class FunctionTransform(CoordinateTransform):
    """
    A coordinate transformation defined by a function that accepts a
    coordinate object and returns the transformed coordinate object.

    Parameters
    ----------
    func : callable
        The transformation function. Should have a call signature
        ``func(formcoord, toframe)``. Note that, unlike
        `CoordinateTransform.__call__`, ``toframe`` is assumed to be of type
        ``tosys`` for this function.
    fromsys : class
        The coordinate frame class to start from.
    tosys : class
        The coordinate frame class to transform into.
    priority : float or int
        The priority if this transform when finding the shortest
        coordinate transform path - large numbers are lower priorities.
    register_graph : `TransformGraph` or None
        A graph to register this transformation with on creation, or
        `None` to leave it unregistered.

    Raises
    ------
    TypeError
        If ``func`` is not callable.
    ValueError
        If ``func`` cannot accept two arguments.


    """

    def __init__(self, func, fromsys, tosys, priority=1, register_graph=None):
        if not callable(func):
            raise TypeError('func must be callable')

        with suppress(TypeError):
            sig = signature(func)
            kinds = [x.kind for x in sig.parameters.values()]
            if (len(x for x in kinds if x == sig.POSITIONAL_ONLY) != 2 and
                    sig.VAR_POSITIONAL not in kinds):
                raise ValueError('provided function does not accept two arguments')

        self.func = func

        super().__init__(fromsys, tosys, priority=priority,
                         register_graph=register_graph)

    def __call__(self, fromcoord, toframe):
        res = self.func(fromcoord, toframe)
        if not isinstance(res, self.tosys):
            raise TypeError(f'the transformation function yielded {res} but '
                            f'should have been of type {self.tosys}')
        if fromcoord.data.differentials and not res.data.differentials:
            warn("Applied a FunctionTransform to a coordinate frame with "
                 "differentials, but the FunctionTransform does not handle "
                 "differentials, so they have been dropped.", AstropyWarning)
        return res


class FunctionTransformWithFiniteDifference(FunctionTransform):
    r"""
    A coordinate transformation that works like a `FunctionTransform`, but
    computes velocity shifts based on the finite-difference relative to one of
    the frame attributes.  Note that the transform function should *not* change
    the differential at all in this case, as any differentials will be
    overridden.

    When a differential is in the from coordinate, the finite difference
    calculation has two components. The first part is simple the existing
    differential, but re-orientation (using finite-difference techniques) to
    point in the direction the velocity vector has in the *new* frame. The
    second component is the "induced" velocity.  That is, the velocity
    intrinsic to the frame itself, estimated by shifting the frame using the
    ``finite_difference_frameattr_name`` frame attribute a small amount
    (``finite_difference_dt``) in time and re-calculating the position.

    Parameters
    ----------
    finite_difference_frameattr_name : str or None
        The name of the frame attribute on the frames to use for the finite
        difference.  Both the to and the from frame will be checked for this
        attribute, but only one needs to have it. If None, no velocity
        component induced from the frame itself will be included - only the
        re-orientation of any existing differential.
    finite_difference_dt : `~astropy.units.Quantity` ['time'] or callable
        If a quantity, this is the size of the differential used to do the
        finite difference.  If a callable, should accept
        ``(fromcoord, toframe)`` and return the ``dt`` value.
    symmetric_finite_difference : bool
        If True, the finite difference is computed as
        :math:`\frac{x(t + \Delta t / 2) - x(t + \Delta t / 2)}{\Delta t}`, or
        if False, :math:`\frac{x(t + \Delta t) - x(t)}{\Delta t}`.  The latter
        case has slightly better performance (and more stable finite difference
        behavior).

    All other parameters are identical to the initializer for
    `FunctionTransform`.

    """

    def __init__(self, func, fromsys, tosys, priority=1, register_graph=None,
                 finite_difference_frameattr_name='obstime',
                 finite_difference_dt=1*u.second,
                 symmetric_finite_difference=True):
        super().__init__(func, fromsys, tosys, priority, register_graph)
        self.finite_difference_frameattr_name = finite_difference_frameattr_name
        self.finite_difference_dt = finite_difference_dt
        self.symmetric_finite_difference = symmetric_finite_difference

    @property
    def finite_difference_frameattr_name(self):
        return self._finite_difference_frameattr_name

    @finite_difference_frameattr_name.setter
    def finite_difference_frameattr_name(self, value):
        if value is None:
            self._diff_attr_in_fromsys = self._diff_attr_in_tosys = False
        else:
            diff_attr_in_fromsys = value in self.fromsys.frame_attributes
            diff_attr_in_tosys = value in self.tosys.frame_attributes
            if diff_attr_in_fromsys or diff_attr_in_tosys:
                self._diff_attr_in_fromsys = diff_attr_in_fromsys
                self._diff_attr_in_tosys = diff_attr_in_tosys
            else:
                raise ValueError('Frame attribute name {} is not a frame '
                                 'attribute of {} or {}'.format(value,
                                                                self.fromsys,
                                                                self.tosys))
        self._finite_difference_frameattr_name = value

    def __call__(self, fromcoord, toframe):
        from .representation import (CartesianRepresentation,
                                     CartesianDifferential)

        supcall = self.func
        if fromcoord.data.differentials:
            # this is the finite difference case

            if callable(self.finite_difference_dt):
                dt = self.finite_difference_dt(fromcoord, toframe)
            else:
                dt = self.finite_difference_dt
            halfdt = dt/2

            from_diffless = fromcoord.realize_frame(fromcoord.data.without_differentials())
            reprwithoutdiff = supcall(from_diffless, toframe)

            # first we use the existing differential to compute an offset due to
            # the already-existing velocity, but in the new frame
            fromcoord_cart = fromcoord.cartesian
            if self.symmetric_finite_difference:
                fwdxyz = (fromcoord_cart.xyz +
                          fromcoord_cart.differentials['s'].d_xyz*halfdt)
                fwd = supcall(fromcoord.realize_frame(CartesianRepresentation(fwdxyz)), toframe)
                backxyz = (fromcoord_cart.xyz -
                           fromcoord_cart.differentials['s'].d_xyz*halfdt)
                back = supcall(fromcoord.realize_frame(CartesianRepresentation(backxyz)), toframe)
            else:
                fwdxyz = (fromcoord_cart.xyz +
                          fromcoord_cart.differentials['s'].d_xyz*dt)
                fwd = supcall(fromcoord.realize_frame(CartesianRepresentation(fwdxyz)), toframe)
                back = reprwithoutdiff
            diffxyz = (fwd.cartesian - back.cartesian).xyz / dt

            # now we compute the "induced" velocities due to any movement in
            # the frame itself over time
            attrname = self.finite_difference_frameattr_name
            if attrname is not None:
                if self.symmetric_finite_difference:
                    if self._diff_attr_in_fromsys:
                        kws = {attrname: getattr(from_diffless, attrname) + halfdt}
                        from_diffless_fwd = from_diffless.replicate(**kws)
                    else:
                        from_diffless_fwd = from_diffless
                    if self._diff_attr_in_tosys:
                        kws = {attrname: getattr(toframe, attrname) + halfdt}
                        fwd_frame = toframe.replicate_without_data(**kws)
                    else:
                        fwd_frame = toframe
                    fwd = supcall(from_diffless_fwd, fwd_frame)

                    if self._diff_attr_in_fromsys:
                        kws = {attrname: getattr(from_diffless, attrname) - halfdt}
                        from_diffless_back = from_diffless.replicate(**kws)
                    else:
                        from_diffless_back = from_diffless
                    if self._diff_attr_in_tosys:
                        kws = {attrname: getattr(toframe, attrname) - halfdt}
                        back_frame = toframe.replicate_without_data(**kws)
                    else:
                        back_frame = toframe
                    back = supcall(from_diffless_back, back_frame)
                else:
                    if self._diff_attr_in_fromsys:
                        kws = {attrname: getattr(from_diffless, attrname) + dt}
                        from_diffless_fwd = from_diffless.replicate(**kws)
                    else:
                        from_diffless_fwd = from_diffless
                    if self._diff_attr_in_tosys:
                        kws = {attrname: getattr(toframe, attrname) + dt}
                        fwd_frame = toframe.replicate_without_data(**kws)
                    else:
                        fwd_frame = toframe
                    fwd = supcall(from_diffless_fwd, fwd_frame)
                    back = reprwithoutdiff

                diffxyz += (fwd.cartesian - back.cartesian).xyz / dt

            newdiff = CartesianDifferential(diffxyz)
            reprwithdiff = reprwithoutdiff.data.to_cartesian().with_differentials(newdiff)
            return reprwithoutdiff.realize_frame(reprwithdiff)
        else:
            return supcall(fromcoord, toframe)


class BaseAffineTransform(CoordinateTransform):
    """Base class for common functionality between the ``AffineTransform``-type
    subclasses.

    This base class is needed because ``AffineTransform`` and the matrix
    transform classes share the ``__call__()`` method, but differ in how they
    generate the affine parameters.  ``StaticMatrixTransform`` passes in a
    matrix stored as a class attribute, and both of the matrix transforms pass
    in ``None`` for the offset. Hence, user subclasses would likely want to
    subclass this (rather than ``AffineTransform``) if they want to provide
    alternative transformations using this machinery.
    """

    def _apply_transform(self, fromcoord, matrix, offset):
        from .representation import (UnitSphericalRepresentation,
                                     CartesianDifferential,
                                     SphericalDifferential,
                                     SphericalCosLatDifferential,
                                     RadialDifferential)

        data = fromcoord.data
        has_velocity = 's' in data.differentials

        # Bail out if no transform is actually requested
        if matrix is None and offset is None:
            return data

        # list of unit differentials
        _unit_diffs = (SphericalDifferential._unit_differential,
                       SphericalCosLatDifferential._unit_differential)
        unit_vel_diff = (has_velocity and
                         isinstance(data.differentials['s'], _unit_diffs))
        rad_vel_diff = (has_velocity and
                        isinstance(data.differentials['s'], RadialDifferential))

        # Some initial checking to short-circuit doing any re-representation if
        # we're going to fail anyways:
        if isinstance(data, UnitSphericalRepresentation) and offset is not None:
            raise TypeError("Position information stored on coordinate frame "
                            "is insufficient to do a full-space position "
                            "transformation (representation class: {})"
                            .format(data.__class__))

        elif (has_velocity and (unit_vel_diff or rad_vel_diff) and
              offset is not None and 's' in offset.differentials):
            # Coordinate has a velocity, but it is not a full-space velocity
            # that we need to do a velocity offset
            raise TypeError("Velocity information stored on coordinate frame "
                            "is insufficient to do a full-space velocity "
                            "transformation (differential class: {})"
                            .format(data.differentials['s'].__class__))

        elif len(data.differentials) > 1:
            # We should never get here because the frame initializer shouldn't
            # allow more differentials, but this just adds protection for
            # subclasses that somehow skip the checks
            raise ValueError("Representation passed to AffineTransform contains"
                             " multiple associated differentials. Only a single"
                             " differential with velocity units is presently"
                             " supported (differentials: {})."
                             .format(str(data.differentials)))

        # If the representation is a UnitSphericalRepresentation, and this is
        # just a MatrixTransform, we have to try to turn the differential into a
        # Unit version of the differential (if no radial velocity) or a
        # sphericaldifferential with zero proper motion (if only a radial
        # velocity) so that the matrix operation works
        if (has_velocity and isinstance(data, UnitSphericalRepresentation) and
                not unit_vel_diff and not rad_vel_diff):
            # retrieve just velocity differential
            unit_diff = data.differentials['s'].represent_as(
                data.differentials['s']._unit_differential, data)
            data = data.with_differentials({'s': unit_diff})  # updates key

        # If it's a RadialDifferential, we flat-out ignore the differentials
        # This is because, by this point (past the validation above), we can
        # only possibly be doing a rotation-only transformation, and that
        # won't change the radial differential. We later add it back in
        elif rad_vel_diff:
            data = data.without_differentials()

        # Convert the representation and differentials to cartesian without
        # having them attached to a frame
        rep = data.to_cartesian()
        diffs = dict([(k, diff.represent_as(CartesianDifferential, data))
                      for k, diff in data.differentials.items()])
        rep = rep.with_differentials(diffs)

        # Only do transform if matrix is specified. This is for speed in
        # transformations that only specify an offset (e.g., LSR)
        if matrix is not None:
            # Note: this applies to both representation and differentials
            rep = rep.transform(matrix)

        # TODO: if we decide to allow arithmetic between representations that
        # contain differentials, this can be tidied up
        if offset is not None:
            newrep = (rep.without_differentials() +
                      offset.without_differentials())
        else:
            newrep = rep.without_differentials()

        # We need a velocity (time derivative) and, for now, are strict: the
        # representation can only contain a velocity differential and no others.
        if has_velocity and not rad_vel_diff:
            veldiff = rep.differentials['s']  # already in Cartesian form

            if offset is not None and 's' in offset.differentials:
                veldiff = veldiff + offset.differentials['s']

            newrep = newrep.with_differentials({'s': veldiff})

        if isinstance(fromcoord.data, UnitSphericalRepresentation):
            # Special-case this because otherwise the return object will think
            # it has a valid distance with the default return (a
            # CartesianRepresentation instance)

            if has_velocity and not unit_vel_diff and not rad_vel_diff:
                # We have to first represent as the Unit types we converted to,
                # then put the d_distance information back in to the
                # differentials and re-represent as their original forms
                newdiff = newrep.differentials['s']
                _unit_cls = fromcoord.data.differentials['s']._unit_differential
                newdiff = newdiff.represent_as(_unit_cls, newrep)

                kwargs = dict([(comp, getattr(newdiff, comp))
                               for comp in newdiff.components])
                kwargs['d_distance'] = fromcoord.data.differentials['s'].d_distance
                diffs = {'s': fromcoord.data.differentials['s'].__class__(
                    copy=False, **kwargs)}

            elif has_velocity and unit_vel_diff:
                newdiff = newrep.differentials['s'].represent_as(
                    fromcoord.data.differentials['s'].__class__, newrep)
                diffs = {'s': newdiff}

            else:
                diffs = newrep.differentials

            newrep = newrep.represent_as(fromcoord.data.__class__)  # drops diffs
            newrep = newrep.with_differentials(diffs)

        elif has_velocity and unit_vel_diff:
            # Here, we're in the case where the representation is not
            # UnitSpherical, but the differential *is* one of the UnitSpherical
            # types. We have to convert back to that differential class or the
            # resulting frame will think it has a valid radial_velocity. This
            # can probably be cleaned up: we currently have to go through the
            # dimensional version of the differential before representing as the
            # unit differential so that the units work out (the distance length
            # unit shouldn't appear in the resulting proper motions)

            diff_cls = fromcoord.data.differentials['s'].__class__
            newrep = newrep.represent_as(fromcoord.data.__class__,
                                         diff_cls._dimensional_differential)
            newrep = newrep.represent_as(fromcoord.data.__class__, diff_cls)

        # We pulled the radial differential off of the representation
        # earlier, so now we need to put it back. But, in order to do that, we
        # have to turn the representation into a repr that is compatible with
        # having a RadialDifferential
        if has_velocity and rad_vel_diff:
            newrep = newrep.represent_as(fromcoord.data.__class__)
            newrep = newrep.with_differentials(
                {'s': fromcoord.data.differentials['s']})

        return newrep

    def __call__(self, fromcoord, toframe):
        params = self._affine_params(fromcoord, toframe)
        newrep = self._apply_transform(fromcoord, *params)
        return toframe.realize_frame(newrep)

    @abstractmethod
    def _affine_params(self, fromcoord, toframe):
        pass


class AffineTransform(BaseAffineTransform):
    """
    A coordinate transformation specified as a function that yields a 3 x 3
    cartesian transformation matrix and a tuple of displacement vectors.

    See `~astropy.coordinates.builtin_frames.galactocentric.Galactocentric` for
    an example.

    Parameters
    ----------
    transform_func : callable
        A callable that has the signature ``transform_func(fromcoord, toframe)``
        and returns: a (3, 3) matrix that operates on ``fromcoord`` in a
        Cartesian representation, and a ``CartesianRepresentation`` with
        (optionally) an attached velocity ``CartesianDifferential`` to represent
        a translation and offset in velocity to apply after the matrix
        operation.
    fromsys : class
        The coordinate frame class to start from.
    tosys : class
        The coordinate frame class to transform into.
    priority : float or int
        The priority if this transform when finding the shortest
        coordinate transform path - large numbers are lower priorities.
    register_graph : `TransformGraph` or None
        A graph to register this transformation with on creation, or
        `None` to leave it unregistered.

    Raises
    ------
    TypeError
        If ``transform_func`` is not callable

    """

    def __init__(self, transform_func, fromsys, tosys, priority=1,
                 register_graph=None):

        if not callable(transform_func):
            raise TypeError('transform_func is not callable')
        self.transform_func = transform_func

        super().__init__(fromsys, tosys, priority=priority,
                         register_graph=register_graph)

    def _affine_params(self, fromcoord, toframe):
        return self.transform_func(fromcoord, toframe)


class StaticMatrixTransform(BaseAffineTransform):
    """
    A coordinate transformation defined as a 3 x 3 cartesian
    transformation matrix.

    This is distinct from DynamicMatrixTransform in that this kind of matrix is
    independent of frame attributes.  That is, it depends *only* on the class of
    the frame.

    Parameters
    ----------
    matrix : array-like or callable
        A 3 x 3 matrix for transforming 3-vectors. In most cases will
        be unitary (although this is not strictly required). If a callable,
        will be called *with no arguments* to get the matrix.
    fromsys : class
        The coordinate frame class to start from.
    tosys : class
        The coordinate frame class to transform into.
    priority : float or int
        The priority if this transform when finding the shortest
        coordinate transform path - large numbers are lower priorities.
    register_graph : `TransformGraph` or None
        A graph to register this transformation with on creation, or
        `None` to leave it unregistered.

    Raises
    ------
    ValueError
        If the matrix is not 3 x 3

    """

    def __init__(self, matrix, fromsys, tosys, priority=1, register_graph=None):
        if callable(matrix):
            matrix = matrix()
        self.matrix = np.array(matrix)

        if self.matrix.shape != (3, 3):
            raise ValueError('Provided matrix is not 3 x 3')

        super().__init__(fromsys, tosys, priority=priority,
                         register_graph=register_graph)

    def _affine_params(self, fromcoord, toframe):
        return self.matrix, None


class DynamicMatrixTransform(BaseAffineTransform):
    """
    A coordinate transformation specified as a function that yields a
    3 x 3 cartesian transformation matrix.

    This is similar to, but distinct from StaticMatrixTransform, in that the
    matrix for this class might depend on frame attributes.

    Parameters
    ----------
    matrix_func : callable
        A callable that has the signature ``matrix_func(fromcoord, toframe)`` and
        returns a 3 x 3 matrix that converts ``fromcoord`` in a cartesian
        representation to the new coordinate system.
    fromsys : class
        The coordinate frame class to start from.
    tosys : class
        The coordinate frame class to transform into.
    priority : float or int
        The priority if this transform when finding the shortest
        coordinate transform path - large numbers are lower priorities.
    register_graph : `TransformGraph` or None
        A graph to register this transformation with on creation, or
        `None` to leave it unregistered.

    Raises
    ------
    TypeError
        If ``matrix_func`` is not callable

    """

    def __init__(self, matrix_func, fromsys, tosys, priority=1,
                 register_graph=None):
        if not callable(matrix_func):
            raise TypeError('matrix_func is not callable')
        self.matrix_func = matrix_func

        super().__init__(fromsys, tosys, priority=priority,
                         register_graph=register_graph)

    def _affine_params(self, fromcoord, toframe):
        return self.matrix_func(fromcoord, toframe), None


class CompositeTransform(CoordinateTransform):
    """
    A transformation constructed by combining together a series of single-step
    transformations.

    Note that the intermediate frame objects are constructed using any frame
    attributes in ``toframe`` or ``fromframe`` that overlap with the intermediate
    frame (``toframe`` favored over ``fromframe`` if there's a conflict).  Any frame
    attributes that are not present use the defaults.

    Parameters
    ----------
    transforms : sequence of `CoordinateTransform` object
        The sequence of transformations to apply.
    fromsys : class
        The coordinate frame class to start from.
    tosys : class
        The coordinate frame class to transform into.
    priority : float or int
        The priority if this transform when finding the shortest
        coordinate transform path - large numbers are lower priorities.
    register_graph : `TransformGraph` or None
        A graph to register this transformation with on creation, or
        `None` to leave it unregistered.
    collapse_static_mats : bool
        If `True`, consecutive `StaticMatrixTransform` will be collapsed into a
        single transformation to speed up the calculation.

    """

    def __init__(self, transforms, fromsys, tosys, priority=1,
                 register_graph=None, collapse_static_mats=True):
        super().__init__(fromsys, tosys, priority=priority,
                         register_graph=register_graph)

        if collapse_static_mats:
            transforms = self._combine_statics(transforms)

        self.transforms = tuple(transforms)

    def _combine_statics(self, transforms):
        """
        Combines together sequences of `StaticMatrixTransform`s into a single
        transform and returns it.
        """
        newtrans = []
        for currtrans in transforms:
            lasttrans = newtrans[-1] if len(newtrans) > 0 else None

            if (isinstance(lasttrans, StaticMatrixTransform) and
                    isinstance(currtrans, StaticMatrixTransform)):
                combinedmat = matrix_product(currtrans.matrix, lasttrans.matrix)
                newtrans[-1] = StaticMatrixTransform(combinedmat,
                                                     lasttrans.fromsys,
                                                     currtrans.tosys)
            else:
                newtrans.append(currtrans)
        return newtrans

    def __call__(self, fromcoord, toframe):
        curr_coord = fromcoord
        for t in self.transforms:
            # build an intermediate frame with attributes taken from either
            # `toframe`, or if not there, `fromcoord`, or if not there, use
            # the defaults
            # TODO: caching this information when creating the transform may
            # speed things up a lot
            frattrs = {}
            for inter_frame_attr_nm in t.tosys.get_frame_attr_names():
                if hasattr(toframe, inter_frame_attr_nm):
                    attr = getattr(toframe, inter_frame_attr_nm)
                    frattrs[inter_frame_attr_nm] = attr
                elif hasattr(fromcoord, inter_frame_attr_nm):
                    attr = getattr(fromcoord, inter_frame_attr_nm)
                    frattrs[inter_frame_attr_nm] = attr

            curr_toframe = t.tosys(**frattrs)
            curr_coord = t(curr_coord, curr_toframe)

        # this is safe even in the case where self.transforms is empty, because
        # coordinate objects are immutable, so copying is not needed
        return curr_coord

    def _as_single_transform(self):
        """
        Return an encapsulated version of the composite transform so that it appears to
        be a single transform.

        The returned transform internally calls the constituent transforms.  If all of
        the transforms are affine, the merged transform is
        `~astropy.coordinates.transformations.DynamicMatrixTransform` (if there are no
        origin shifts) or `~astropy.coordinates.transformations.AffineTransform`
        (otherwise).  If at least one of the transforms is not affine, the merged
        transform is
        `~astropy.coordinates.transformations.FunctionTransformWithFiniteDifference`.
        """
        # Create a list of the transforms including flattening any constituent CompositeTransform
        transforms = [t if not isinstance(t, CompositeTransform) else t._as_single_transform()
                      for t in self.transforms]

        if all([isinstance(t, BaseAffineTransform) for t in transforms]):
            # Check if there may be an origin shift
            fixed_origin = all([isinstance(t, (StaticMatrixTransform, DynamicMatrixTransform))
                                for t in transforms])

            # Dynamically define the transformation function
            def single_transform(from_coo, to_frame):
                if from_coo.is_equivalent_frame(to_frame):  # loopback to the same frame
                    return None if fixed_origin else (None, None)

                # Create a merged attribute dictionary for any intermediate frames
                # For any attributes shared by the "from"/"to" frames, the "to" frame takes
                #   precedence because this is the same choice implemented in __call__()
                merged_attr = {name: getattr(from_coo, name)
                               for name in from_coo.frame_attributes}
                merged_attr.update({name: getattr(to_frame, name)
                                    for name in to_frame.frame_attributes})

                affine_params = (None, None)
                # Step through each transform step (frame A -> frame B)
                for i, t in enumerate(transforms):
                    # Extract the relevant attributes for frame A
                    if i == 0:
                        # If frame A is actually the initial frame, preserve its attributes
                        a_attr = {name: getattr(from_coo, name)
                                  for name in from_coo.frame_attributes}
                    else:
                        a_attr = {k: v for k, v in merged_attr.items()
                                  if k in t.fromsys.frame_attributes}

                    # Extract the relevant attributes for frame B
                    b_attr = {k: v for k, v in merged_attr.items()
                              if k in t.tosys.frame_attributes}

                    # Obtain the affine parameters for the transform
                    # Note that we insert some dummy data into frame A because the transformation
                    #   machinery requires there to be data present.  Removing that limitation
                    #   is a possible TODO, but some care would need to be taken because some affine
                    #   transforms have branching code depending on the presence of differentials.
                    next_affine_params = t._affine_params(t.fromsys(from_coo.data, **a_attr),
                                                          t.tosys(**b_attr))

                    # Combine the affine parameters with the running set
                    affine_params = _combine_affine_params(affine_params, next_affine_params)

                # If there is no origin shift, return only the matrix
                return affine_params[0] if fixed_origin else affine_params

            # The return type depends on whether there is any origin shift
            transform_type = DynamicMatrixTransform if fixed_origin else AffineTransform
        else:
            # Dynamically define the transformation function
            def single_transform(from_coo, to_frame):
                if from_coo.is_equivalent_frame(to_frame):  # loopback to the same frame
                    return to_frame.realize_frame(from_coo.data)
                return self(from_coo, to_frame)

            transform_type = FunctionTransformWithFiniteDifference

        return transform_type(single_transform, self.fromsys, self.tosys, priority=self.priority)


def _combine_affine_params(params, next_params):
    """
    Combine two sets of affine parameters.

    The parameters for an affine transformation are a 3 x 3 Cartesian
    transformation matrix and a displacement vector, which can include an
    attached velocity.  Either type of parameter can be ``None``.
    """
    M, vec = params
    next_M, next_vec = next_params

    # Multiply the transformation matrices if they both exist
    if M is not None and next_M is not None:
        new_M = next_M @ M
    else:
        new_M = M if M is not None else next_M

    if vec is not None:
        # Transform the first displacement vector by the second transformation matrix
        if next_M is not None:
            vec = vec.transform(next_M)

        # Calculate the new displacement vector
        if next_vec is not None:
            if 's' in vec.differentials and 's' in next_vec.differentials:
                # Adding vectors with velocities takes more steps
                # TODO: Add support in representation.py
                new_vec_velocity = vec.differentials['s'] + next_vec.differentials['s']
                new_vec = vec.without_differentials() + next_vec.without_differentials()
                new_vec = new_vec.with_differentials({'s': new_vec_velocity})
            else:
                new_vec = vec + next_vec
        else:
            new_vec = vec
    else:
        new_vec = next_vec

    return new_M, new_vec


# map class names to colorblind-safe colors
trans_to_color = {}
trans_to_color[AffineTransform] = '#555555'  # gray
trans_to_color[FunctionTransform] = '#783001'  # dark red-ish/brown
trans_to_color[FunctionTransformWithFiniteDifference] = '#d95f02'  # red-ish
trans_to_color[StaticMatrixTransform] = '#7570b3'  # blue-ish
trans_to_color[DynamicMatrixTransform] = '#1b9e77'  # green-ish

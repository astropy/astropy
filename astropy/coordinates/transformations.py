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

from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import heapq
import inspect
import subprocess

from abc import ABCMeta, abstractmethod
from collections import defaultdict

import numpy as np

from ..utils.compat import ignored
from ..extern import six


__all__ = ['TransformGraph', 'CoordinateTransform', 'FunctionTransform',
           'StaticMatrixTransform', 'DynamicMatrixTransform', 'CompositeTransform']


class TransformGraph(object):
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
                    dct[nm] = c

        return self._cached_names_dct

    @property
    def frame_set(self):
        """
        A `set` of all the frame classes present in this `TransformGraph`.
        """
        if self._cached_frame_set is None:
            self._cached_frame_set = frm_set = set()
            for a in self._graph:
                frm_set.add(a)
                for b in self._graph[a]:
                    frm_set.add(b)

        return self._cached_frame_set.copy()

    def invalidate_cache(self):
        """
        Invalidates the cache that stores optimizations for traversing the
        transform graph.  This is called automatically when transforms
        are added or removed, but will need to be called manually if
        weights on transforms are modified inplace.
        """
        self._cached_names_dct = None
        self._cached_frame_set = None
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
        transform : CoordinateTransform or similar callable
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
        if not six.callable(transform):
            raise TypeError('transform must be callable')

        self._graph[fromsys][tosys] = transform
        self.invalidate_cache()

    def remove_transform(self, fromsys, tosys, transform):
        """
        Removes a coordinate transform from the graph.

        Parameters
        ----------
        fromsys : class or `None`
            The coordinate frame *class* to start from. If `None`,
            ``transform`` will be searched for and removed (``tosys`` must
            also be `None`).
        tosys : class or `None`
            The coordinate frame *class* to transform into. If `None`,
            ``transform`` will be searched for and removed (``fromsys`` must
            also be `None`).
        transform : callable or `None`
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
                    if b is transform:
                        del agraph[b]
                        break
            else:
                raise ValueError('Could not find transform {0} in the '
                                 'graph'.format(transform))

        else:
            if transform is None:
                self._graph[fromsys].pop(tosys, None)
            else:
                curr = self._graph[fromsys].get(tosys, None)
                if curr is transform:
                    self._graph[fromsys].pop(tosys)
                else:
                    raise ValueError('Current transform from {0} to {1} is not '
                                     '{2}'.format(fromsys, tosys, transform))
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
        path : list of classes or `None`
            The path from ``fromsys`` to ``tosys`` as an in-order sequence
            of classes.  This list includes *both* ``fromsys`` and
            ``tosys``. Is `None` if there is no possible path.
        distance : number
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

        #otherwise, need to construct the path:

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
        trans : `CompositeTransform` or `None`
            If there is a path from ``fromsys`` to ``tosys``, this is a
            transform object for that path.   If no path could be found, this is
            `None`.

        Notes
        -----
        This function always returns a `CompositeTransform`, because
        `CompositeTransform` is slightly more adaptable in the way it can be
        called than other transform classes. Specifically, it takes care of
        inetermediate steps of transformations in a way that is consistent with
        1-hop transformations.

        """
        if not inspect.isclass(fromsys):
            raise TypeError('fromsys is not a class')
        if not inspect.isclass(fromsys):
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
        coordcls
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
        return list(six.iterkeys(self._cached_names))

    def to_dot_graph(self, priorities=True, addnodes=[], savefn=None,
                     savelayout='plain', saveformat=None):
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
        savefn : `None` or str
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
        invclsaliases = dict([(v, k) for k, v in six.iteritems(self._cached_names)])
        for n in nodes:
            if n in invclsaliases:
                nodenames.append('{0} [shape=oval label="{0}\\n`{1}`"]'.format(n.__name__, invclsaliases[n]))
            else:
                nodenames.append(n.__name__ + '[ shape=oval ]')

        edgenames = []
        # Now the edges
        for a in self._graph:
            agraph = self._graph[a]
            for b in agraph:
                pri = agraph[b].priority if hasattr(agraph[b], 'priority') else 1
                edgenames.append((a.__name__, b.__name__, pri))

        # generate simple dot format graph
        lines = ['digraph AstropyCoordinateTransformGraph {']
        lines.append('; '.join(nodenames) + ';')
        for enm1, enm2, weights in edgenames:
            labelstr = '[ label = "{0}" ]'.format(weights) if priorities else ''
            lines.append('{0} -> {1}{2};'.format(enm1, enm2, labelstr))
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
                    raise IOError('problem running graphviz: \n' + stderr)

                with open(savefn, 'w') as f:
                    f.write(stdout)

        return dotgraph

    def to_networkx_graph(self):
        """
        Converts this transform graph into a networkx graph.

        .. note::
            You must have the `networkx <http://networkx.lanl.gov/>`_
            package installed for this to work.

        Returns
        -------
        nxgraph : `networkx.Graph <http://networkx.lanl.gov/reference/classes.graph.html>`_
            This `TransformGraph` as a `networkx.Graph`_.
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
                pri = agraph[b].priority if hasattr(agraph[b], 'priority') else 1
                nxgraph.add_edge(a, b, weight=pri)

        return nxgraph

    def transform(self, transcls, fromsys, tosys, priority=1):
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
        priority : number
            The priority if this transform when finding the shortest
            coordinate tranform path - large numbers are lower priorities.

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
            # this doesn't do anything directly with the trasnform because
            # ``register_graph=self`` stores it in the transform graph
            # automatically
            transcls(func, fromsys, tosys, priority=priority,
                     register_graph=self)
            return func
        return deco


#<--------------------Define the builtin transform classes--------------------->

@six.add_metaclass(ABCMeta)
class CoordinateTransform(object):
    """
    An object that transforms a coordinate from one system to another.
    Subclasses must implement `__call__` with the provided signature.
    They should also call this superclass's ``__init__`` in their
    ``__init__``.

    Parameters
    ----------
    fromsys : class
        The coordinate frame class to start from.
    tosys : class
        The coordinate frame class to transform into.
    priority : number
        The priority if this transform when finding the shortest
        coordinate tranform path - large numbers are lower priorities.
    register_graph : `TransformGraph` or `None`
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
            #the if statement is there so that non-frame things might be usable
            #if it makes sense
            for from_nm in fromsys.get_frame_attr_names():
                if from_nm in tosys.get_frame_attr_names():
                    overlap.append(from_nm)

    def register(self, graph):
        """
        Add this transformation to the requested Transformation graph,
        replacing anything already connecting these two coordinates.

        Parameters
        ----------
        graph : a TransformGraph object
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
        fromcoord : fromsys object
            An object of class matching ``fromsys`` that is to be transformed.
        toframe : object
            An object that has the attributes necessary to fully specify the
            frame.  That is, it must have attributes with names that match the
            keys of the dictionary that ``tosys.get_frame_attr_names()``
            returns. Typically this is of class ``tosys``, but it *might* be
            some other class as long as it has the appropriate attributes.

        Returns
        -------
        tocoord : tosys object
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
    priority : number
        The priority if this transform when finding the shortest
        coordinate tranform path - large numbers are lower priorities.
    register_graph : `TransformGraph` or `None`
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
        from inspect import getargspec

        if not six.callable(func):
            raise TypeError('func must be callable')

        with ignored(TypeError):
            # TypeError raised for things getargspec can't process.  We'll trust
            # the transform designer knows what they're doing, though, because
            # sometimes this is fine..
            argspec = getargspec(func)
            if (len(argspec[0]) - len(argspec[3]) != 2) and not argspec[1]:
                raise ValueError('provided function does not accept two arguments')

        self.func = func

        super(FunctionTransform, self).__init__(fromsys, tosys,
            priority=priority, register_graph=register_graph)

    def __call__(self, fromcoord, toframe):
        res = self.func(fromcoord, toframe)
        if not isinstance(res, self.tosys):
            raise TypeError('the transformation function yielded {0} but '
                'should have been of type {1}'.format(res, self.tosys))
        return res


class StaticMatrixTransform(CoordinateTransform):
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
    priority : number
        The priority if this transform when finding the shortest
        coordinate tranform path - large numbers are lower priorities.
    register_graph : `TransformGraph` or `None`
        A graph to register this transformation with on creation, or
        `None` to leave it unregistered.

    Raises
    ------
    ValueError
        If the matrix is not 3 x 3

    """
    def __init__(self, matrix, fromsys, tosys, priority=1, register_graph=None):
        if six.callable(matrix):
            matrix = matrix()
        self.matrix = np.array(matrix)

        if self.matrix.shape != (3, 3):
            raise ValueError('Provided matrix is not 3 x 3')

        super(StaticMatrixTransform, self).__init__(fromsys, tosys,
            priority=priority, register_graph=register_graph)

    def __call__(self, fromcoord, toframe):
        from .representation import CartesianRepresentation, \
                                    UnitSphericalRepresentation

        xyz = fromcoord.represent_as(CartesianRepresentation).xyz
        v = xyz.reshape((3, xyz.size // 3))
        v2 = np.dot(np.asarray(self.matrix), v)
        subshape = xyz.shape[1:]
        x = v2[0].reshape(subshape)
        y = v2[1].reshape(subshape)
        z = v2[2].reshape(subshape)

        newrep = CartesianRepresentation(x, y, z)
        if fromcoord.data.__class__ == UnitSphericalRepresentation:
            #need to special-case this because otherwise the new class will
            #think it has a valid distance
            newrep = newrep.represent_as(UnitSphericalRepresentation)

        frameattrs = dict([(attrnm, getattr(fromcoord, attrnm))
                           for attrnm in self.overlapping_frame_attr_names])

        return toframe.realize_frame(newrep, **frameattrs)


class DynamicMatrixTransform(CoordinateTransform):
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
    priority : number
        The priority if this transform when finding the shortest
        coordinate tranform path - large numbers are lower priorities.
    register_graph : `TransformGraph` or `None`
        A graph to register this transformation with on creation, or
        `None` to leave it unregistered.

    Raises
    ------
    TypeError
        If ``matrix_func`` is not callable

    """
    def __init__(self, matrix_func, fromsys, tosys, priority=1,
                 register_graph=None):
        if not six.callable(matrix_func):
            raise TypeError('matrix_func is not callable')
        self.matrix_func = matrix_func

        super(DynamicMatrixTransform, self).__init__(fromsys, tosys,
            priority=priority, register_graph=register_graph)

    def __call__(self, fromcoord, toframe):
        from .representation import CartesianRepresentation, \
                                    UnitSphericalRepresentation

        xyz = fromcoord.represent_as(CartesianRepresentation).xyz
        v = xyz.reshape((3, xyz.size // 3))
        v2 = np.dot(np.asarray(self.matrix_func(fromcoord, toframe)), v)
        subshape = xyz.shape[1:]
        x = v2[0].reshape(subshape)
        y = v2[1].reshape(subshape)
        z = v2[2].reshape(subshape)

        newrep = CartesianRepresentation(x, y, z)
        if fromcoord.data.__class__ == UnitSphericalRepresentation:
            #need to special-case this because otherwise the new class will
            #think it has a valid distance
            newrep = newrep.represent_as(UnitSphericalRepresentation)

        return toframe.realize_frame(newrep)


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
    transforms : sequence of `CoordinateTransform` objects
        The sequence of transformations to apply.
    fromsys : class
        The coordinate frame class to start from.
    tosys : class
        The coordinate frame class to transform into.
    priority : number
        The priority if this transform when finding the shortest
        coordinate tranform path - large numbers are lower priorities.
    register_graph : `TransformGraph` or `None`
        A graph to register this transformation with on creation, or
        `None` to leave it unregistered.
    collapse_static_mats : bool
        If `True`, consecutive `StaticMatrixTransform` will be collapsed into a
        single transformation to speed up the calculation.

    """
    def __init__(self, transforms, fromsys, tosys, priority=1,
                 register_graph=None, collapse_static_mats=True):
        super(CompositeTransform, self).__init__(fromsys, tosys,
                                                 priority=priority,
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
                combinedmat = np.dot(lasttrans.matrix, currtrans.matrix)
                newtrans[-1] = StaticMatrixTransform(combinedmat,
                                                     lasttrans.fromsys,
                                                     currtrans.tosys)
            else:
                newtrans.append(currtrans)
        return newtrans

    def __call__(self, fromcoord, toframe):
        curr_coord = fromcoord
        for t in self.transforms:
            #build an intermediate frame with attributes taken from either
            #`fromframe`, or if not there, `toframe`, or if not there, use
            #the defaults
            #TODO: caching this information when creating the transform may
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

        # this is safe even in the case enere self.transforms is empty, because
        # coordinate objects are immutible, so copying is not needed
        return curr_coord

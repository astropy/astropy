# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""
This module contains the framework for transforming points from
one coordinate system to another (e.g. equatorial to galactic). The
implementation is actually in individual coordinates in the
`builtin_systems` module, while this module provides the framework and
related utilities.
"""
from abc import ABCMeta, abstractmethod

import numpy as np

__all__ = ['StaticMatrixTransform', 'FunctionTransform',
           'DynamicMatrixTransform', 'CompositeStaticMatrixTransform',
           'static_transform_matrix', 'transform_function',
           'dynamic_transform_matrix', 'coordinate_alias'
          ]


class TransformGraph(object):
    """
    A graph representing the paths between coordinate systems.
    """

    def __init__(self):
        from collections import defaultdict

        self._graph = defaultdict(dict)
        self._clsaliases = {}

        self.invalidate_cache()  # generates cache entries

    def add_transform(self, fromsys, tosys, transform):
        """
        Add a new coordinate transformation to the graph.

        Parameters
        ----------
        fromsys : class
            The coordinate system *class* to start from
        tosys : class
            The coordinate system *class* to transform to
        transform : callable
            The transformation object. Should have call parameters compatible
            with `CoordinateTransform`.

        Raises
        ------
        TypeError
            If `fromsys` or `tosys` are not classes or `transform` is
            not callable.
        """
        from inspect import isclass

        if not isclass(fromsys):
            raise TypeError('fromsys must be a class')
        if not isclass(tosys):
            raise TypeError('tosys must be a class')
        if not callable(transform):
            raise TypeError('transform must be callable')

        self._graph[fromsys][tosys] = transform
        self.invalidate_cache()

    def remove_transform(self, fromsys, tosys, transform):
        """
        Removes a coordinate transform from the graph.

        Parameters
        ----------
        fromsys : class or None
            The coordinate system *class* to start from. If None,
            `transform` will be searched for and removed (`tosys` must
            also be None).
        tosys : class or None
            The coordinate system *class* to transform into. If None,
            `transform` will be searched for and removed (`fromsys` must
            also be None).
        transform : callable or None
            The transformation object to be removed or None.  If None
            and `tosys` and `fromsys` are supplied, there will be no
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
            The starting coordinate system.
        tosys : class
            The starting coordinate system.

        Returns
        -------
        path : list of classes or None
            The path from `fromsys` to `tosys` as an in-order sequence
            of classes.  This list includes *both* `fromsys` and
            `tosys`. Is None if there is no possible path.
        distance : number
            The total distance/priority from `fromsys` to `tosys`.  If
            priorities are not set this is the number of trasnforms
            needed. Is `inf` if there is no possible path.
        """
        import heapq

        inf = float('inf')

        # special-case the 0-path and 1-path
        if tosys is fromsys:
            return [tosys], 0
        elif tosys in self._graph[fromsys]:
            t = self._graph[fromsys][tosys]
            return [fromsys, tosys], float(t.priority if hasattr(t, 'priority') else 1)

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

        # this dict will store the distance to node from `fromsys` and the path
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

    def invalidate_cache(self):
        """
        Invalidates the cache that stores optimizations for traversing the
        transform cache.  This is called automatically when transforms
        are added or removed, but will need to be called manually if
        weights on transforms are modified inplace.
        """
        self._shortestpaths = {}

    # TODO: cache composites so they don't need to be generated every time?
    def get_transform(self, fromsys, tosys):
        """
        Determines or generates a transformation between two coordinate
        systems.

        Parameters
        ----------
        fromsys : class
            The coordinate system *class* to start from
        tosys : class
            The coordinate system *class* to transform into.

        Returns
        -------
        trans : `CoordinateTransform` or None
            If there is a path from `fromsys` to `tosys`, this is a transform
            object for that path.  If None, no path could be found.
        """
        if tosys in self._graph[fromsys]:
            return self._graph[fromsys][tosys]
        else:
            path, distance = self.find_shortest_path(fromsys, tosys)

            if path is None:
                return None

            transforms = []
            currsys = fromsys
            for p in path[1:]:  # first element is fromsys so we skip it
                transforms.append(self._graph[currsys][p])
                currsys = p

            # TODO: collapse "runs" of statics?
            if all([isinstance(p, StaticMatrixTransform) for p in path]):
                return CompositeStaticMatrixTransform(fromsys, tosys, transforms, register=False)
            else:
                return CompositeTransform(fromsys, tosys, transforms, register=False)

    def add_coord_name(self, name, coordcls):
        """
        Adds an alias for a coordinate, primarily for allowing
        attribute-style access of coordinate transformations (e.g.,
        ``coordasgal = coord.galactic``).

        Parameters
        ----------
        name : str
            The alias for the coordinate class. Should be a valid
            python identifier.
        coordcls : class
            The class object to be referenced by this name.

        Raises
        ------
        ValueError
            If `coordcls` already has a name assigned.
        """
        if coordcls in self._clsaliases.values():
            idx = self._clsaliases.values().index(coordcls)
            oldnm = self._clsaliases.keys()[idx]
            msg = 'Coordinate class {0} already has a name: {1}'
            raise ValueError(msg.format(coordcls, oldnm))
        self._clsaliases[name] = coordcls

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
            The coordinate class corresponding to the `name` or None if
            no such class exists.
        """
        return self._clsaliases.get(name, None)

    def get_aliases(self):
        """
        Returns all available transform aliases. They will all be
        valid arguments to `lookup_name`.

        Returns
        -------
        nms : list
            The aliases for coordinate systems.
        """
        return self._clsaliases.keys()

    def to_dot_graph(self, priorities=True, addnodes=[], savefn=None,
                     savelayout='plain', saveformat=None):
        """
        Converts this transform graph to the graphviz_ DOT format, and
        optionally saves it (requires graphviz_ be installed and on your
        path).

        Parameters
        ----------
        priorities : bool
            If True, show the priority values for each transform.  Otherwise,
            the will not be included in the graph.
        addnodes : sequence of str
            Additional coordinate systems to add (this can include systems
            already in the transform graph, but they will only appear once).
        savefn : None or str
            The file name to save this graph to or None to not save
            to a file.
        savelayout : str
            The graphviz program to use to layout the graph (see
            graphviz_ for details) or 'plain' to just save the DOT graph
            content. Ignored if `savefn` is None.
        saveformat : str
            The graphviz output format. (e.g. the ``-Txxx`` option for
            the command line program - see graphviz docs for details).
            Ignored if `savefn` is None.

        Returns
        -------
        dotgraph : str
            A string with the DOT format graph.



        .. _graphviz: http://www.graphviz.org/
        """
        from subprocess import Popen, PIPE

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
        invclsaliases = dict([(v, k) for k, v in self._clsaliases.iteritems()])
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
                proc = Popen(args, stdin=PIPE, stdout=PIPE, stderr=PIPE)
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
        nxgraph : `networkx.Graph`
            This `TransformGraph` as a `networkx.Graph`.
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


# The primary transform graph for astropy coordinates
master_transform_graph = TransformGraph()


class CoordinateTransform(object):
    """
    An object that transforms a coordinate from one system to another.
    Subclasses must implement `__call__` with the provided signature.
    They should also call this superclass's `__init__` in their
    `__init__`.
    """

    __metaclass__ = ABCMeta

    def __init__(self, fromsys, tosys, register=True):
        from inspect import isclass

        self.fromsys = fromsys
        self.tosys = tosys

        if register:
            # this will do the type-checking
            self.register()
        else:
            if not isclass(fromsys) or not isclass(tosys):
                raise TypeError('fromsys and tosys must be classes')

    def register(self):
        """
        Add this transformation to the master transformation graph, replacing
        anything already connecting these two coordinates.
        """
        master_transform_graph.add_transform(self.fromsys, self.tosys, self)

    def unregister(self):
        """
        Remove this transformation to the master transformation graph.

        Raises
        ------
        ValueError
            If this is not currently in the transform graph.
        """
        master_transform_graph.remove_transform(self.fromsys, self.tosys, self)

    @abstractmethod
    def __call__(self, fromcoord):
        """
        Accepts the provided coordinate object and yields a new coordinate
        object with the transform applied.
        """


# TODO: array: specify in the docs how arrays should be dealt with
class FunctionTransform(CoordinateTransform):
    """
    A coordinate transformation defined by a function that simply
    accepts a coordinate object and returns the transformed coordinate
    object.

    Parameters
    ----------
    fromsys : class
        The coordinate system *class* to start from.
    tosys : class
        The coordinate system *class* to transform into.
    func : callable
        The transformation function.
    copyobstime : bool
        If True (default) the value of the  `_obstime` attribute will be copied
        to the newly-produced coordinate.

    priority : number
        The priority if this transform when finding the shortest
        coordinate tranform path - large numbers are lower priorities.
    register : bool
        Determines if this transformation will be registered in the
        astropy master transform graph.

    Raises
    ------
    TypeError
        If `func` is not callable.
    ValueError
        If `func` cannot accept one argument.


    """
    def __init__(self, fromsys, tosys, func, copyobstime=True, priority=1,
                 register=True):
        from inspect import getargspec

        if not callable(func):
            raise TypeError('func must be callable')

        try:
            argspec = getargspec(func)
            if (len(argspec[0]) - len(argspec[3]) != 1) and not argspec[1]:
                raise ValueError('provided function does not accept a single argument')
        except TypeError:
            pass  # hopefully this person knows what they're doing...

        self.func = func
        self.priority = priority
        self.copyobstime = copyobstime

        super(FunctionTransform, self).__init__(fromsys, tosys)

    def __call__(self, fromcoord):
        res = self.func(fromcoord)
        if not isinstance(res, self.tosys):
            raise TypeError('the transformation function yielded {0} but '
                'should have been of type {1}'.format(res, self.tosys))

        if self.copyobstime:
            # copy over the obstime
            if hasattr(fromcoord, '_obstime') and hasattr(res, '_obstime'):
                res._obstime = fromcoord._obstime

        return res


class StaticMatrixTransform(CoordinateTransform):
    """
    A coordinate transformation defined as a 3 x 3 cartesian
    transformation matrix.

    Parameters
    ----------
    fromsys : class
        The coordinate system *class* to start from.
    tosys : class
        The coordinate system *class* to transform into.
    matrix: array-like
        A 3 x 3 matrix for transforming 3-vectors. In most cases will
        be unitary (although this is not strictly required).
    priority : number
        The priority if this transform when finding the shortest
        coordinate tranform path - large numbers are lower priorities.

    Raises
    ------
    ValueError
        If the matrix is not 3 x 3

    """
    def __init__(self, fromsys, tosys, matrix, priority=1, register=True):
        self.matrix = np.array(matrix)
        if self.matrix.shape != (3, 3):
            raise ValueError('Provided matrix is not 3 x 3')
        self.priority = priority
        super(StaticMatrixTransform, self).__init__(fromsys, tosys)

    # TODO: array: this needs some extra bits to do the broadcasting right
    def __call__(self, fromcoord):
        v = [fromcoord.x, fromcoord.y, fromcoord.z]
        x, y, z = np.dot(np.asarray(self.matrix), v)
        unit = None if fromcoord.distance is None else fromcoord.distance._unit
        result = self.tosys(x=x, y=y, z=z, unit=unit)

        # copy over the observation time
        if hasattr(fromcoord, '_obstime') and hasattr(result, '_obstime'):
            result._obstime = fromcoord._obstime

        return result


class CompositeStaticMatrixTransform(StaticMatrixTransform):
    """
    A `MatrixTransform` constructed by combining a sequence of matricies
    together.  See `MatrixTransform` for syntax details.

    Parameters
    ----------
    fromsys : class
        The coordinate system *class* to start from.
    tosys : class
        The coordinate system *class* to transform into.
    matrices: sequence of array-like
        A sequence of 3 x 3 cartesian transformation matricies.
    priority : number
        The priority if this transform when finding the shortest
        coordinate tranform path - large numbers are lower priorities.

    """
    def __init__(self, fromsys, tosys, matricies, priority=1, register=True):
        self.matricies = [np.array(m) for m in matricies]
        for m in matricies:
            if m.shape != (3, 3):
                raise ValueError('One of the provided matrices is not 3 x 3')

        matrix = np.array(self.matricies[0])
        if len(self.matricies) > 1:
            for m in self.matricies[1:]:
                matrix = np.dot(np.asarray(self.matrix), m)

        super(CompositeStaticMatrixTransform, self).__init__(self, fromsys,
            tosys, matrix, priority)


class DynamicMatrixTransform(CoordinateTransform):
    """
    A coordinate transformation specified as a function that yields a
    3 x 3 cartesian transformation matrix.

    Parameters
    ----------
    fromsys : class
        The coordinate system *class* to start from.
    tosys : class
        The coordinate system *class* to transform into.
    matrix_func: callable
        A callable that accepts a coordinate object and yields the 3 x 3
        matrix that converts it to the new coordinate system.
    priority : number
        The priority if this transform when finding the shortest
        coordinate tranform path - large numbers are lower priorities.

    Raises
    ------
    TypeError
        If `matrix_func` is not callable

    """
    def __init__(self, fromsys, tosys, matrix_func, priority=1, register=True):
        if not callable(matrix_func):
            raise TypeError('matrix_func is not callable')
        self.matrix_func = matrix_func
        self.priority = priority
        super(DynamicMatrixTransform, self).__init__(fromsys, tosys, register)

    # TODO: array: this needs some extra bits to do the broadcasting right
    def __call__(self, fromcoord):
        v = [fromcoord.x, fromcoord.y, fromcoord.z]
        x, y, z = np.dot(np.asarray(self.matrix_func(fromcoord)), v)
        unit = None if fromcoord.distance is None else fromcoord.distance._unit
        result = self.tosys(x=x, y=y, z=z, unit=unit)

        # copy over the observation time
        if hasattr(fromcoord, '_obstime') and hasattr(result, '_obstime'):
            result._obstime = fromcoord._obstime

        return result


class CompositeTransform(CoordinateTransform):
    """
    A `MatrixTransform` constructed by combining a sequence of matricies
    together.  See `MatrixTransform` for syntax details.

    Parameters
    ----------
    fromsys : class
        The coordinate system *class* to start from.
    tosys : class
        The coordinate system *class* to transform into.
    transforms: sequence of `CoordinateTransform`s
        A sequence of transformations to apply in sequence.
    priority : number
        The priority if this transform when finding the shortest
        coordinate tranform path - large numbers are lower priorities.

    """
    def __init__(self, fromsys, tosys, transforms, priority=1, register=True):
        self.transforms = transforms
        super(CompositeTransform, self).__init__(fromsys, tosys, register)

    def __call__(self, fromcoord):
        coord = fromcoord
        for t in self.transforms:
            coord = t(coord)
        return coord


#<------------function decorators for actual practical use--------------------->
def transform_function(fromsys, tosys, copyobstime=True, priority=1):
    """
    A function decorator for defining transformations between coordinate
    systems.

    .. note::
        If decorating a static method of a class, ``@staticmethod``
        should be  added *above* this decorator.

    Parameters
    ----------
    fromsys : class
        The coordinate system this function starts from.
    tosys : class
        The coordinate system this function results in.
    copyobstime : bool
        If True (default) the value of the  `_obstime` attribute will be
        copied to the newly-produced coordinate.
    priority : number
        The priority if this transform when finding the shortest
        coordinate tranform path - large numbers are lower priorities.

    """
    def deco(func):
        # this doesn't do anything directly with the trasnform because
        #``register=True`` stores it in the transform graph automatically
        FunctionTransform(fromsys, tosys, func, copyobstime=copyobstime,
                          priority=priority, register=True)
        return func
    return deco


def static_transform_matrix(fromsys, tosys, priority=1):
    """
    A function decorator for defining transformations between coordinate
    systems using a matrix.

    The decorated function should accept *no* arguments and yield a
    3 x 3 matrix.

    .. note::
        If decorating a static method of a class, ``@staticmethod``
        should be  added *above* this decorator.

    Parameters
    ----------
    fromsys : class
        The coordinate system this function starts from.
    tosys : class
        The coordinate system this function results in.
    priority : number
        The priority if this transform when finding the shortest
        coordinate tranform path - large numbers are lower priorities.
    """
    def deco(matfunc):
        StaticMatrixTransform(fromsys, tosys, matfunc(), priority, register=True)
        return matfunc
    return deco


def dynamic_transform_matrix(fromsys, tosys, priority=1):
    """
    A function decorator for defining transformations between coordinate
    systems using a function that yields a matrix.

    The decorated function should accept a single argument, the
    coordinate object to be transformed, and should return a 3 x 3
    matrix.

    .. note::
        If decorating a static method of a class, ``@staticmethod``
        should be  added *above* this decorator.

    Parameters
    ----------
    fromsys : class
        The coordinate system this function starts from.
    tosys : class
        The coordinate system this function results in.
    priority : number
        The priority if this transform when finding the shortest
        coordinate tranform path - large numbers are lower priorities.
    """
    def deco(matfunc):
        DynamicMatrixTransform(fromsys, tosys, matfunc, priority, register=True)
        return matfunc
    return deco


def coordinate_alias(name, coordcls=None):
    """
    Gives a short name to this coordinate system, allowing other coordinate
    objects to convert to this one using attribute-style access.

    Parameters
    ----------
    name : str
        The short alias to use for this coordinate class. Should be a
        valid python identifier.
    coordcls : class or None
        Either the coordinate class to register or None to use this as a
        decorator.

    Examples
    --------
    For use with a class already defined, do::

        coordinate_alias('fancycoords', MyFancyCoordinateClass)

    To use as a decorator, do::

        @coordiante_alias('fancycoords')
        class MyFancyCoordinateClass(SphericalCoordinatesBase):
            ...

    """
    if coordcls is None:
        def deco(cls):
            master_transform_graph.add_coord_name(name, cls)
            return cls
        return deco
    else:
        master_transform_graph.add_coord_name(name, coordcls)

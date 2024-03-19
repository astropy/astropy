# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""The class that implements the framework for defining graphs of transformations.

The fundamental idea is that each class is a node in the transformation graph,
and transitions from one node to another are defined as functions (or methods)
wrapped in transformation objects.
"""

import heapq
import subprocess
from collections import defaultdict
from contextlib import contextmanager
from itertools import pairwise

from astropy.coordinates.transformations.affine import (
    AffineTransform,
    DynamicMatrixTransform,
    StaticMatrixTransform,
)
from astropy.coordinates.transformations.composite import CompositeTransform
from astropy.coordinates.transformations.function import (
    FunctionTransform,
    FunctionTransformWithFiniteDifference,
)
from astropy.utils import lazyproperty

__all__ = ["TransformGraph"]


# map class names to colorblind-safe colors
trans_to_color = {}
trans_to_color[AffineTransform] = "#555555"  # gray
trans_to_color[FunctionTransform] = "#783001"  # dark red-ish/brown
trans_to_color[FunctionTransformWithFiniteDifference] = "#d95f02"  # red-ish
trans_to_color[StaticMatrixTransform] = "#7570b3"  # blue-ish
trans_to_color[DynamicMatrixTransform] = "#1b9e77"  # green-ish


def frame_attrs_from_set(frame_set):
    """
    A `dict` of all the attributes of all frame classes in this
    `~astropy.coordinates.TransformGraph`.

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
    this `~astropy.coordinates.TransformGraph`.

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

    @lazyproperty
    def _cached_names(self):
        dct = {}
        for c in self.frame_set:
            if (nm := getattr(c, "name", None)) is not None:
                dct |= dict.fromkeys(nm if isinstance(nm, list) else [nm], c)
        return dct

    @property
    def frame_set(self):
        """
        A `set` of all the frame classes present in this TransformGraph.
        """
        if self._cached_frame_set is None:
            self._cached_frame_set = set()
            for a in self._graph:
                self._cached_frame_set.add(a)
                for b in self._graph[a]:
                    self._cached_frame_set.add(b)

        return self._cached_frame_set.copy()

    @lazyproperty
    def frame_attributes(self):
        """
        A `dict` of all the attributes of all frame classes in this TransformGraph.
        """
        return frame_attrs_from_set(self.frame_set)

    @property
    def frame_component_names(self):
        """
        A `set` of all component names every defined within any frame class in
        this TransformGraph.
        """
        return frame_comps_from_set(self.frame_set)

    def invalidate_cache(self):
        """
        Invalidates the cache that stores optimizations for traversing the
        transform graph.  This is called automatically when transforms
        are added or removed, but will need to be called manually if
        weights on transforms are modified inplace.
        """
        del self._cached_names
        self._cached_frame_set = None
        del self.frame_attributes
        self._shortestpaths = {}
        self._composite_cache = {}

    def add_transform(self, fromsys, tosys, transform):
        """Add a new coordinate transformation to the graph.

        Parameters
        ----------
        fromsys : class
            The coordinate frame class to start from.
        tosys : class
            The coordinate frame class to transform into.
        transform : `~astropy.coordinates.CoordinateTransform`
            The transformation object. Typically a
            `~astropy.coordinates.CoordinateTransform` object, although it may
            be some other callable that is called with the same signature.

        Raises
        ------
        TypeError
            If ``fromsys`` or ``tosys`` are not classes or ``transform`` is
            not callable.

        """
        if not isinstance(fromsys, type):
            raise TypeError("fromsys must be a class")
        if not isinstance(tosys, type):
            raise TypeError("tosys must be a class")
        if not callable(transform):
            raise TypeError("transform must be callable")

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

            raise ValueError(
                f"Frame(s) {list(invalid_frames)} contain invalid attribute names:"
                f" {invalid_attrs}\nFrame attributes can not conflict with *any* of"
                " the frame data component names (see"
                " `frame_transform_graph.frame_component_names`)."
            )

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
            if tosys is not fromsys:
                raise ValueError("fromsys and tosys must both be None if either are")
            if transform is None:
                raise ValueError("cannot give all Nones to remove_transform")

            # search for the requested transform by brute force and remove it
            for a, agraph in self._graph.items():
                for b, bgraph in agraph.items():
                    if bgraph is transform:
                        del agraph[b]
                        fromsys = a
                        break

                # If the transform was found, need to break out of the outer for loop too
                if fromsys:
                    break
            else:
                raise ValueError(f"Could not find transform {transform} in the graph")

        elif transform is None or self._graph[fromsys].get(tosys, None) is transform:
            self._graph[fromsys].pop(tosys, None)
        else:
            raise ValueError(
                f"Current transform from {fromsys} to {tosys} is not {transform}"
            )

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
        inf = float("inf")

        # special-case the 0 or 1-path
        if tosys is fromsys and tosys not in self._graph[fromsys]:
            # Means there's no transform necessary to go from it to itself.
            return [tosys], 0
        if tosys in self._graph[fromsys]:
            # this will also catch the case where tosys is fromsys, but has
            # a defined transform.
            t = self._graph[fromsys][tosys]
            return [fromsys, tosys], float(getattr(t, "priority", 1))

        # otherwise, need to construct the path:

        if fromsys in self._shortestpaths:
            # already have a cached result
            return self._shortestpaths[fromsys].get(tosys, (None, inf))

        # use Dijkstra's algorithm to find shortest path in all other cases

        # We store nodes as `dict` keys because differently from `list` uniqueness is
        # guaranteed and differently from `set` insertion order is preserved.
        nodes = {}
        for node, node_graph in self._graph.items():
            nodes[node] = None
            nodes |= {node: None for node in node_graph}

        if fromsys not in nodes or tosys not in nodes:
            # fromsys or tosys are isolated or not registered, so there's
            # certainly no way to get from one to the other
            return None, inf

        edgeweights = {}
        # construct another graph that is a dict of dicts of priorities
        # (used as edge weights in Dijkstra's algorithm)
        for a, graph in self._graph.items():
            edgeweights[a] = {b: float(getattr(graph[b], "priority", 1)) for b in graph}

        # entries in q are [distance, count, nodeobj, pathlist]
        # count is needed because in py 3.x, tie-breaking fails on the nodes.
        # this way, insertion order is preserved if the weights are the same
        q = [[0, -1, fromsys, []]]
        q.extend([inf, i, n, []] for i, n in enumerate(nodes) if n is not fromsys)

        # this dict will store the distance to node from ``fromsys`` and the path
        result = {}

        # definitely starts as a valid heap because of the insert line; from the
        # node to itself is always the shortest distance
        while q:
            d, orderi, n, path = heapq.heappop(q)

            if d == inf:
                # everything left is unreachable from fromsys, just copy them to
                # the results and jump out of the loop
                result[n] = (None, d)
                for d, orderi, n, path in q:
                    result[n] = (None, d)
                break
            result[n] = (path, d)
            path.append(n)
            if n not in edgeweights:
                # this is a system that can be transformed to, but not from.
                continue
            for n2 in edgeweights[n]:
                if n2 not in result:  # already visited
                    # find where n2 is in the heap
                    for q_elem in q:
                        if q_elem[2] == n2:
                            if (newd := d + edgeweights[n][n2]) < q_elem[0]:
                                q_elem[0] = newd
                                q_elem[3] = list(path)
                                heapq.heapify(q)
                            break
                    else:
                        raise ValueError("n2 not in heap - this should be impossible!")

        # cache for later use
        self._shortestpaths[fromsys] = result
        return result[tosys]

    def get_transform(self, fromsys, tosys):
        """Generates and returns the CompositeTransform for a transformation
        between two coordinate systems.

        Parameters
        ----------
        fromsys : class
            The coordinate frame class to start from.
        tosys : class
            The coordinate frame class to transform into.

        Returns
        -------
        trans : `~astropy.coordinates.CompositeTransform` or None
            If there is a path from ``fromsys`` to ``tosys``, this is a
            transform object for that path.   If no path could be found, this is
            `None`.

        Notes
        -----
        A `~astropy.coordinates.CompositeTransform` is always returned, because
        `~astropy.coordinates.CompositeTransform` is slightly more adaptable in
        the way it can be called than other transform classes. Specifically, it
        takes care of intermediate steps of transformations in a way that is
        consistent with 1-hop transformations.

        """
        if not isinstance(fromsys, type):
            raise TypeError("fromsys is not a class")
        if not isinstance(tosys, type):
            raise TypeError("tosys is not a class")

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
            comptrans = CompositeTransform(
                transforms, fromsys, tosys, register_graph=False
            )
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

    def to_dot_graph(
        self,
        priorities=True,
        addnodes=[],
        savefn=None,
        savelayout="plain",
        saveformat=None,
        color_edges=True,
    ):
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
        savelayout : {"plain", "dot", "neato", "fdp", "sfdp", "circo", "twopi", "nop", "nop2", "osage", "patchwork"}
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
        # We store nodes as `dict` keys because differently from `list` uniqueness is
        # guaranteed and differently from `set` insertion order is preserved.
        nodes = {}
        for node, node_graph in self._graph.items():
            nodes[node] = None
            nodes |= {node: None for node in node_graph}
        nodes |= {node: None for node in addnodes}
        nodenames = []
        invclsaliases = {
            f: [k for k, v in self._cached_names.items() if v == f]
            for f in self.frame_set
        }
        for n in nodes:
            if n in invclsaliases:
                aliases = "`\\n`".join(invclsaliases[n])
                nodenames.append(
                    f'{n.__name__} [shape=oval label="{n.__name__}\\n`{aliases}`"]'
                )
            else:
                nodenames.append(n.__name__ + "[ shape=oval ]")

        edgenames = []
        # Now the edges
        for a, agraph in self._graph.items():
            for b, transform in agraph.items():
                pri = getattr(transform, "priority", 1)
                color = trans_to_color[transform.__class__] if color_edges else "black"
                edgenames.append((a.__name__, b.__name__, pri, color))

        # generate simple dot format graph
        lines = ["digraph AstropyCoordinateTransformGraph {"]
        lines.append("graph [rankdir=LR]")
        lines.append("; ".join(nodenames) + ";")
        for enm1, enm2, weights, color in edgenames:
            labelstr_fmt = "[ {0} {1} ]"
            priority_part = f'label = "{weights}"' if priorities else ""
            color_part = f'color = "{color}"'
            labelstr = labelstr_fmt.format(priority_part, color_part)
            lines.append(f"{enm1} -> {enm2}{labelstr};")
        lines.append("")
        lines.append("overlap=false")
        lines.append("}")
        dotgraph = "\n".join(lines)

        if savefn is not None:
            if savelayout == "plain":
                with open(savefn, "w") as f:
                    f.write(dotgraph)
            # Options from https://graphviz.org/docs/layouts/
            elif savelayout in (
                "dot",
                "neato",
                "fdp",
                "sfdp",
                "circo",
                "twopi",
                "nop",
                "nop2",
                "osage",
                "patchwork",
            ):
                args = [savelayout]
                if saveformat is not None:
                    args.append("-T" + saveformat)
                proc = subprocess.Popen(
                    args,
                    stdin=subprocess.PIPE,
                    stdout=subprocess.PIPE,
                    stderr=subprocess.PIPE,
                )
                stdout, stderr = proc.communicate(dotgraph)
                if proc.returncode != 0:
                    raise OSError("problem running graphviz: \n" + stderr)

                with open(savefn, "w") as f:
                    f.write(stdout)
            else:
                raise NotImplementedError(f'savelayout="{savelayout}" is not supported')

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
            This `~astropy.coordinates.TransformGraph` as a
            `networkx.Graph <https://networkx.github.io/documentation/stable/reference/classes/graph.html>`_.
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
                pri = transform.priority if hasattr(transform, "priority") else 1
                color = trans_to_color[transform.__class__]
                nxgraph.add_edge(a, b, weight=pri, color=color)

        return nxgraph

    def transform(self, transcls, fromsys, tosys, priority=1, **kwargs):
        """A function decorator for defining transformations.

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
        initializer accepts a callable, and that the second and third are
        ``fromsys`` and ``tosys``. If this is not true, you should just
        initialize the class manually and use
        `~astropy.coordinates.TransformGraph.add_transform` instead of this
        decorator.

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
            transcls(
                func, fromsys, tosys, priority=priority, register_graph=self, **kwargs
            )
            return func

        return deco

    def _add_merged_transform(self, fromsys, tosys, *furthersys, priority=1):
        """
        Add a single-step transform that encapsulates a multi-step transformation path,
        using the transforms that already exist in the graph.

        The created transform internally calls the existing transforms.  If all of the
        transforms are affine, the merged transform is
        `~astropy.coordinates.DynamicMatrixTransform` (if there are no
        origin shifts) or `~astropy.coordinates.AffineTransform`
        (otherwise).  If at least one of the transforms is not affine, the merged
        transform is
        `~astropy.coordinates.FunctionTransformWithFiniteDifference`.

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
        `~astropy.coordinates.ICRS` to be safe.  Transforming through an inertial
        frame ensures that changes in observation time and observer
        location/velocity are properly accounted for.

        An error will be raised if a direct transform between ``fromsys`` and
        ``tosys`` already exist.
        """
        frames = [fromsys, tosys, *furthersys]
        lastsys = frames[-1]
        full_path = self.get_transform(fromsys, lastsys)
        transforms = [
            self.get_transform(frame_a, frame_b)
            for frame_a, frame_b in pairwise(frames)
        ]
        if None in transforms:
            raise ValueError("This transformation path is not possible")
        if len(full_path.transforms) == 1:
            raise ValueError(
                f"A direct transform for {fromsys.__name__}->{lastsys.__name__} already"
                " exists"
            )

        self.add_transform(
            fromsys,
            lastsys,
            CompositeTransform(
                transforms, fromsys, lastsys, priority=priority
            )._as_single_transform(),
        )

    @contextmanager
    def impose_finite_difference_dt(self, dt):
        """
        Context manager to impose a finite-difference time step on all applicable transformations.

        For each transformation in this transformation graph that has the attribute
        ``finite_difference_dt``, that attribute is set to the provided value.  The only standard
        transformation with this attribute is
        `~astropy.coordinates.FunctionTransformWithFiniteDifference`.

        Parameters
        ----------
        dt : `~astropy.units.Quantity` ['time'] or callable
            If a quantity, this is the size of the differential used to do the finite difference.
            If a callable, should accept ``(fromcoord, toframe)`` and return the ``dt`` value.
        """
        key = "finite_difference_dt"
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

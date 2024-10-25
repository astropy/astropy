from itertools import permutations

import numpy as np

__all__ = ["auto_assign_coord_positions"]


def auto_assign_coord_positions(ax):
    """
    Given a ``WCSAxes`` instance, automatically update any dynamic tick, tick
    label and axis label positions.

    This function operates in-place on the axes and assumes that
    ``_update_ticks`` has already been called on all the ``CoordinateHelper``
    instances.
    """
    # Since ticks, tick labels and axis labels can all be auto or fixed, we need
    # a few rules to decide in what order to process things:
    #
    # * We should deal with axis labels last - if they are set to auto they
    #   should be made to match the tick labels.
    # * If ticks and tick labels are both set to automatic, then we first
    #   process tick labels and at the end we can make ticks match tick labels.
    # * We should use existing fixed tick label assignments to exclude certain
    #   spines.
    # * Fixed tick positions should not exclude other labels from being placed
    #   on that same spine since multiple different ticks might appear on the
    #   same axis.
    # * If ticks are not auto and tick labels are auto then tick labels should be
    #   placed at a spine containing ticks.

    # Start off by simplifying some cases that don't require algorithmic
    # positioning. First, if tick labels are shown at fixed positions, we
    # should just adjust any auto ticks and axis labels to match.
    for coords in ax._all_coords:
        for coord in coords:
            if "#" not in (pos := coord.get_ticklabel_position()):
                if "#" in coord.get_ticks_position():
                    coord.set_ticks_position(pos + ["#"])
                if "#" in coord.get_axislabel_position():
                    coord.set_axislabel_position(pos + ["#"])

    # At this point, all coordinates requiring automatic placement have
    # tick labels requiring automatic placement (any coordinates with fixed
    # tick label placement are dealt with). We construct a list of these
    # coordinates and also keep track of spines on which tick labels are
    # already fixed.

    auto_coords = []
    already_used = []
    for coords in ax._all_coords:
        for coord in coords:
            pos = coord.get_ticklabel_position()
            if "#" in pos:
                auto_coords.append(coord)
            else:
                already_used += list(pos)

    # If there are no more coordinates with auto settings, we are all done
    if len(auto_coords) == 0:
        return

    # Extract the spines for the frame
    spines = coords.frame._spine_auto_position_order

    # Construct a new list of spines taking into account excluded ones
    spines = "".join(s for s in spines if s not in already_used)

    # We create an iterable of different assignments of spines to
    # coords, where empty string means the coordinate will not be shown
    # on any axis.
    if len(auto_coords) > len(spines):
        spines = spines + "".join(" " * (len(auto_coords) - len(spines)))
    options = permutations(spines, r=len(auto_coords))

    # We now loop through and check different assignments based on the options
    # and try and find the one that maximises the number of tick labels in the
    # overall plot.
    n_tick_max = -1
    best_option = None
    for option in options:
        # Check if option is consistent with any fixed tick positions - that
        # is, if ticks were explicitly requested on axes say b and t then we
        # shouldn't automatically put tick labels on say l or r.
        consistent = True
        for c, s in zip(auto_coords, option):
            pos = c.get_ticks_position()
            if "#" not in pos and s not in pos:
                consistent = False
                break
        if not consistent:
            continue

        # Determine the number of tick labels on each axis
        n_on_each = {s: len(c._ticks.world[s]) for c, s in zip(auto_coords, option)}

        # Determine the total number of tick labels
        n_tick = sum(n_on_each.values())

        # Determine a sorted version of this list by spine axis order
        n_on_each_sorted = [n_on_each.get(s, 0) for s in spines]

        # We should ideally not have empty spines, so if there are any non-zero
        # values in n_on_each then we should add a penalty for every 0 value
        # that occurs before the last non-zero value
        idx = np.nonzero(n_on_each_sorted)[0]
        if len(idx) > 0:
            n_tick -= n_on_each_sorted[: idx[-1]].count(0)

        # Keep track of the best option so far
        if n_tick > n_tick_max:
            n_tick_max = n_tick
            best_option = option

    # Finalize assignments
    for coord, spine in zip(auto_coords, best_option):
        if "#" in coord.get_ticks_position():
            coord.set_ticks_position([spine, "#"])
        if "#" in coord.get_ticklabel_position():
            coord.set_ticklabel_position([spine, "#"])
        if "#" in coord.get_axislabel_position():
            coord.set_axislabel_position([spine, "#"])

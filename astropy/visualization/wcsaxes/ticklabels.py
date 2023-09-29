# Licensed under a 3-clause BSD style license - see LICENSE.rst
import warnings
from collections import defaultdict

import numpy as np
from matplotlib import rcParams
from matplotlib.artist import allow_rasterization
from matplotlib.text import Text

from astropy.utils.decorators import deprecated_renamed_argument
from astropy.utils.exceptions import AstropyDeprecationWarning

from .frame import RectangularFrame


def sort_using(X, Y):
    return [x for (y, x) in sorted(zip(Y, X))]


class TickLabels(Text):
    def __init__(self, frame, *args, **kwargs):
        self.clear()
        self._frame = frame
        super().__init__(*args, **kwargs)
        self.set_clip_on(True)
        self.set_visible_axes("all")
        self.set_pad(rcParams["xtick.major.pad"])
        self._exclude_overlapping = False

        # Mapping from axis > list[bounding boxes]
        self._axis_bboxes = defaultdict(list)

        # Stale if either xy positions haven't been calculated, or if
        # something changes that requires recomputing the positions
        self._stale = True

        # Check rcParams
        if "color" not in kwargs:
            self.set_color(rcParams["xtick.color"])

        if "size" not in kwargs:
            self.set_size(rcParams["xtick.labelsize"])

    def clear(self):
        self.world = defaultdict(list)
        self.data = defaultdict(list)
        self.angle = defaultdict(list)
        self.text = defaultdict(list)
        self.disp = defaultdict(list)

    def add(
        self,
        axis=None,
        world=None,
        pixel=None,
        angle=None,
        text=None,
        axis_displacement=None,
        data=None,
    ):
        """
        Add a label.

        Parameters
        ----------
        axis : str
            Axis to add label to.
        world : Quantity
            Coordinate value along this axis.
        pixel : [float, float]
            Pixel coordinates of the label. Deprecated and no longer used.
        angle : float
            Angle of the label.
        text : str
            Label text.
        axis_displacement : float
            Displacement from axis.
        data : [float, float]
            Data coordinates of the label.
        """
        required_args = ["axis", "world", "angle", "text", "axis_displacement", "data"]
        if pixel is not None:
            warnings.warn(
                "Setting the pixel coordinates of a label does nothing and is"
                " deprecated, as these can only be accurately calculated when"
                " Matplotlib is drawing a figure. To prevent this warning pass the"
                f" following arguments as keyword arguments: {required_args}",
                AstropyDeprecationWarning,
            )
        if (
            axis is None
            or world is None
            or angle is None
            or text is None
            or axis_displacement is None
            or data is None
        ):
            raise TypeError(
                f"All of the following arguments must be provided: {required_args}"
            )

        self.world[axis].append(world)
        self.data[axis].append(data)
        self.angle[axis].append(angle)
        self.text[axis].append(text)
        self.disp[axis].append(axis_displacement)

        self._stale = True

    def sort(self):
        """
        Sort by axis displacement, which allows us to figure out which parts
        of labels to not repeat.
        """
        for axis in self.world:
            self.world[axis] = sort_using(self.world[axis], self.disp[axis])
            self.data[axis] = sort_using(self.data[axis], self.disp[axis])
            self.angle[axis] = sort_using(self.angle[axis], self.disp[axis])
            self.text[axis] = sort_using(self.text[axis], self.disp[axis])
            self.disp[axis] = sort_using(self.disp[axis], self.disp[axis])
        self._stale = True

    def simplify_labels(self):
        """
        Figure out which parts of labels can be dropped to avoid repetition.
        """
        self.sort()
        for axis in self.world:
            t1 = self.text[axis][0]
            for i in range(1, len(self.world[axis])):
                t2 = self.text[axis][i]
                if len(t1) != len(t2):
                    t1 = self.text[axis][i]
                    continue
                start = 0
                # In the following loop, we need to ignore the last character,
                # hence the len(t1) - 1. This is because if we have two strings
                # like 13d14m15s we want to make sure that we keep the last
                # part (15s) even if the two labels are identical.
                for j in range(len(t1) - 1):
                    if t1[j] != t2[j]:
                        break
                    if t1[j] not in "-0123456789.":
                        start = j + 1
                t1 = self.text[axis][i]
                if start != 0:
                    starts_dollar = self.text[axis][i].startswith("$")
                    self.text[axis][i] = self.text[axis][i][start:]
                    if starts_dollar:
                        self.text[axis][i] = "$" + self.text[axis][i]
                # Remove any empty LaTeX inline math mode string
                if self.text[axis][i] == "$$":
                    self.text[axis][i] = ""

        self._stale = True

    def set_pad(self, value):
        self._pad = value
        self._stale = True

    def get_pad(self):
        return self._pad

    def set_visible_axes(self, visible_axes):
        self._visible_axes = visible_axes
        self._stale = True

    def get_visible_axes(self):
        if self._visible_axes == "all":
            return self.world.keys()
        else:
            return [x for x in self._visible_axes if x in self.world]

    def set_exclude_overlapping(self, exclude_overlapping):
        self._exclude_overlapping = exclude_overlapping

    def _set_xy_alignments(self, renderer):
        """
        Compute and set the x, y positions and the horizontal/vertical alignment of
        each label.
        """
        if not self._stale:
            return

        self.simplify_labels()
        text_size = renderer.points_to_pixels(self.get_size())

        visible_axes = self.get_visible_axes()
        self.xy = {axis: {} for axis in visible_axes}
        self.ha = {axis: {} for axis in visible_axes}
        self.va = {axis: {} for axis in visible_axes}

        for axis in visible_axes:
            for i in range(len(self.world[axis])):
                # In the event that the label is empty (which is not expected
                # but could happen in unforeseen corner cases), we should just
                # skip to the next label.
                if self.text[axis][i] == "":
                    continue

                x, y = self._frame.parent_axes.transData.transform(self.data[axis][i])
                pad = renderer.points_to_pixels(self.get_pad() + self._tick_out_size)

                if isinstance(self._frame, RectangularFrame):
                    # This is just to preserve the current results, but can be
                    # removed next time the reference images are re-generated.
                    if np.abs(self.angle[axis][i]) < 45.0:
                        ha = "right"
                        va = "bottom"
                        dx = -pad
                        dy = -text_size * 0.5
                    elif np.abs(self.angle[axis][i] - 90.0) < 45:
                        ha = "center"
                        va = "bottom"
                        dx = 0
                        dy = -text_size - pad
                    elif np.abs(self.angle[axis][i] - 180.0) < 45:
                        ha = "left"
                        va = "bottom"
                        dx = pad
                        dy = -text_size * 0.5
                    else:
                        ha = "center"
                        va = "bottom"
                        dx = 0
                        dy = pad

                    x = x + dx
                    y = y + dy

                else:
                    # This is the more general code for arbitrarily oriented
                    # axes

                    # Set initial position and find bounding box
                    self.set_text(self.text[axis][i])
                    self.set_position((x, y))
                    bb = super().get_window_extent(renderer)

                    # Find width and height, as well as angle at which we
                    # transition which side of the label we use to anchor the
                    # label.
                    width = bb.width
                    height = bb.height

                    # Project axis angle onto bounding box
                    ax = np.cos(np.radians(self.angle[axis][i]))
                    ay = np.sin(np.radians(self.angle[axis][i]))

                    # Set anchor point for label
                    if np.abs(self.angle[axis][i]) < 45.0:
                        dx = width
                        dy = ay * height
                    elif np.abs(self.angle[axis][i] - 90.0) < 45:
                        dx = ax * width
                        dy = height
                    elif np.abs(self.angle[axis][i] - 180.0) < 45:
                        dx = -width
                        dy = ay * height
                    else:
                        dx = ax * width
                        dy = -height

                    dx *= 0.5
                    dy *= 0.5

                    # Find normalized vector along axis normal, so as to be
                    # able to nudge the label away by a constant padding factor

                    dist = np.hypot(dx, dy)

                    ddx = dx / dist
                    ddy = dy / dist

                    dx += ddx * pad
                    dy += ddy * pad

                    x = x - dx
                    y = y - dy

                    ha = "center"
                    va = "center"

                self.xy[axis][i] = (x, y)
                self.ha[axis][i] = ha
                self.va[axis][i] = va

        self._stale = False

    def _get_bb(self, axis, i, renderer):
        """
        Get the bounding box of an individual label. n.b. _set_xy_alignment()
        must be called before this method.
        """
        if self.text[axis][i] == "":
            return

        self.set_text(self.text[axis][i])
        self.set_position(self.xy[axis][i])
        self.set_ha(self.ha[axis][i])
        self.set_va(self.va[axis][i])
        return super().get_window_extent(renderer)

    @property
    def _all_bboxes(self):
        # List of all tick label bounding boxes
        ret = []
        for axis in self._axis_bboxes:
            ret += self._axis_bboxes[axis]
        return ret

    def _set_existing_bboxes(self, bboxes):
        self._existing_bboxes = bboxes

    @allow_rasterization
    @deprecated_renamed_argument(old_name="bboxes", new_name=None, since="6.0")
    @deprecated_renamed_argument(old_name="ticklabels_bbox", new_name=None, since="6.0")
    @deprecated_renamed_argument(old_name="tick_out_size", new_name=None, since="6.0")
    def draw(self, renderer, bboxes=None, ticklabels_bbox=None, tick_out_size=None):
        # Reset bounding boxes
        self._axis_bboxes = defaultdict(list)

        if not self.get_visible():
            return

        self._set_xy_alignments(renderer)

        for axis in self.get_visible_axes():
            for i in range(len(self.world[axis])):
                # This implicitly sets the label text, position, alignment
                bb = self._get_bb(axis, i, renderer)
                if bb is None:
                    continue

                # TODO: the problem here is that we might get rid of a label
                # that has a key starting bit such as -0:30 where the -0
                # might be dropped from all other labels.
                if (
                    not self._exclude_overlapping
                    or bb.count_overlaps(self._all_bboxes + self._existing_bboxes) == 0
                ):
                    super().draw(renderer)
                    self._axis_bboxes[axis].append(bb)

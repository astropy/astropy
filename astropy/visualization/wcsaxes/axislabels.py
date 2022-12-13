# Licensed under a 3-clause BSD style license - see LICENSE.rst


import matplotlib.transforms as mtransforms
import numpy as np
from matplotlib import rcParams
from matplotlib.text import Text

from .frame import RectangularFrame


class AxisLabels(Text):
    def __init__(self, frame, minpad=1, *args, **kwargs):
        # Use rcParams if the following parameters were not specified explicitly
        if "weight" not in kwargs:
            kwargs["weight"] = rcParams["axes.labelweight"]
        if "size" not in kwargs:
            kwargs["size"] = rcParams["axes.labelsize"]
        if "color" not in kwargs:
            kwargs["color"] = rcParams["axes.labelcolor"]

        self._frame = frame
        super().__init__(*args, **kwargs)
        self.set_clip_on(True)
        self.set_visible_axes("all")
        self.set_ha("center")
        self.set_va("center")
        self._minpad = minpad
        self._visibility_rule = "labels"

    def get_minpad(self, axis):
        try:
            return self._minpad[axis]
        except TypeError:
            return self._minpad

    def set_visible_axes(self, visible_axes):
        self._visible_axes = visible_axes

    def get_visible_axes(self):
        if self._visible_axes == "all":
            return self._frame.keys()
        else:
            return [x for x in self._visible_axes if x in self._frame]

    def set_minpad(self, minpad):
        self._minpad = minpad

    def set_visibility_rule(self, value):
        allowed = ["always", "labels", "ticks"]
        if value not in allowed:
            raise ValueError(
                f"Axis label visibility rule must be one of{' / '.join(allowed)}"
            )

        self._visibility_rule = value

    def get_visibility_rule(self):
        return self._visibility_rule

    def draw(
        self,
        renderer,
        bboxes,
        ticklabels_bbox,
        coord_ticklabels_bbox,
        ticks_locs,
        visible_ticks,
    ):
        if not self.get_visible():
            return

        text_size = renderer.points_to_pixels(self.get_size())
        # Flatten the bboxes for all coords and all axes
        ticklabels_bbox_list = []
        for bbcoord in ticklabels_bbox.values():
            for bbaxis in bbcoord.values():
                ticklabels_bbox_list += bbaxis

        for axis in self.get_visible_axes():
            if self.get_visibility_rule() == "ticks":
                if not ticks_locs[axis]:
                    continue
            elif self.get_visibility_rule() == "labels":
                if not coord_ticklabels_bbox:
                    continue

            padding = text_size * self.get_minpad(axis)

            # Find position of the axis label. For now we pick the mid-point
            # along the path but in future we could allow this to be a
            # parameter.
            x, y, normal_angle = self._frame[axis]._halfway_x_y_angle()

            label_angle = (normal_angle - 90.0) % 360.0
            if 135 < label_angle < 225:
                label_angle += 180
            self.set_rotation(label_angle)

            # Find label position by looking at the bounding box of ticks'
            # labels and the image. It sets the default padding at 1 times the
            # axis label font size which can also be changed by setting
            # the minpad parameter.

            if isinstance(self._frame, RectangularFrame):
                if (
                    len(ticklabels_bbox_list) > 0
                    and ticklabels_bbox_list[0] is not None
                ):
                    coord_ticklabels_bbox[axis] = [
                        mtransforms.Bbox.union(ticklabels_bbox_list)
                    ]
                else:
                    coord_ticklabels_bbox[axis] = [None]

                visible = (
                    axis in visible_ticks and coord_ticklabels_bbox[axis][0] is not None
                )

                if axis == "l":
                    if visible:
                        x = coord_ticklabels_bbox[axis][0].xmin
                    x = x - padding

                elif axis == "r":
                    if visible:
                        x = coord_ticklabels_bbox[axis][0].x1
                    x = x + padding

                elif axis == "b":
                    if visible:
                        y = coord_ticklabels_bbox[axis][0].ymin
                    y = y - padding

                elif axis == "t":
                    if visible:
                        y = coord_ticklabels_bbox[axis][0].y1
                    y = y + padding

            else:  # arbitrary axis
                x = x + np.cos(np.radians(normal_angle)) * (padding + text_size * 1.5)
                y = y + np.sin(np.radians(normal_angle)) * (padding + text_size * 1.5)

            self.set_position((x, y))
            super().draw(renderer)

            bb = super().get_window_extent(renderer)
            bboxes.append(bb)

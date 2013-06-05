from matplotlib.axes import Axes


class ParasiteAxes(Axes):

    def __init__(self, fig, rect, parent, adjustable='datalim'):

        self.fig = fig
        self.rect = rect
        self._parent = parent

        Axes.__init__(self, fig, rect, adjustable=adjustable)

        self.xaxis.set_ticks_position('top')
        self.yaxis.set_ticks_position('right')
        self.set_frame_on(False)

        fig.add_axes(self)

    def draw(self, renderer):

        self.axes.viewLim.set(self._parent.viewLim)
        self.set_position(self._parent.get_position())

        Axes.draw(self, renderer)


class HostAxes(Axes):

    def __init__(self, fig, rect, adjustable='datalim'):

        self.fig = fig
        self.rect = rect

        Axes.__init__(self, fig, rect, adjustable=adjustable)

        self.xaxis.set_ticks_position('bottom')
        self.yaxis.set_ticks_position('left')

        # Not quite API-conforming, but not sure how to make sure both the
        # Parasite and Host axes get added when the user does:
        #
        # >>> fig.add_axes(ax)
        #
        # so we do it here:

        fig.add_axes(self)

        self.twin = ParasiteAxes(fig, rect, self)

    def draw(self, renderer):

        Axes.draw(self, renderer)

from matplotlib.pyplot import Annotation
from kivy.clock import Clock


class Selector(object):

    def __init__(self, canvas):

        self.canvas = canvas

    def connect(self):

        print("here")

        self.cidmotion = self.canvas.mpl_connect("motion_notify_event", self.draw_selection)
        self.cidpress = self.canvas.mpl_connect("button_press_event", self.click)
        self.cidrelease = self.canvas.mpl_connect("button_release_event", self.release)

    def draw_selection(self, event):

        # print(event.x, event.y)
        pass

    def click(self, event):

        print("clicked")
        print(event.x, event.y)

    def release(self, event):

        print("released")
        print(event.x, event.y)


class ShowTooltip(object):

    bbox_style = dict(boxstyle="round",
                      fc="#37abc8",
                      ec="white")

    def __init__(self, plot_obj, plt_idx):

        ax = plot_obj.gca()
        self.padding = ax.get_ylim()[1] * .05

        self.annotations = []

        self.tooltip_info = {}

        self.mp = None


    def _reset(self):

        for an in self.annotations:
            an.remove()
            self.annotations.remove(an)

        self.tooltip_info = {}

    def show_tooltip(self, event):

        for x in event.inaxes.patches:

            if x.contains(event)[0]:

                info = {"text": "{}".format(x.get_height),
                        "pos": (x.get_x() + (x.get_width() / 2),
                                x.get_height() + self.padding)}

                if info == self.tooltip_info:
                    return

                self._reset()

                self.tooltip_info = info

                if self.annotations:
                    map(lambda an: an.remove(), self.annotations)

                annot = event.inaxes.annotate("{}".format(x.get_height()),
                                         xy=info["pos"],
                                         xytext=info["pos"],
                                         bbox=self.bbox_style)

                self.annotations.append(annot)

                event.inaxes.get_figure().tight_layout()
                event.inaxes.get_figure().canvas.draw()
                return

        self._reset()
        event.inaxes.get_figure().canvas.draw()

    def _check_mp(self, mp, event):

        if mp == self.mp:
            self.show_tooltip(event)

    def __call__(self, event):

        # Ignore when the even is outise axes
        if not event.inaxes:
            return

        self.mp = (event.x, event.y)

        mp = (event.x, event.y)
        # Ensure that the mouse is idle for more than 0.5s
        Clock.schedule_once(lambda x: self._check_mp(mp, event), 0.1)


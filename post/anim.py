from post.plots import *
import matplotlib.animation as animation
from matplotlib.animation import FFMpegWriter

class pAnimation():
    def __init__(self,
                 P,
                 name="animation",
                 interval=2,
                 frames=100,
                 xlim=False,
                 ylim=False):

        self.P = P # The plot object
        self.interval = interval
        self.frames = frames
        self.animate()

    def animate(self):
        self.ani = animation.FuncAnimation( self.P.fig,
                                            self.update,
                                            frames=self.frames,
                                            interval=self.interval)

    def update(self, frame):
        self.P.update(frame)

    def save(self, name, fps=15):
        writer = animation.FFMpegWriter(fps=fps,
                                        metadata=dict(artist='Me'),
                                        bitrate=1800)
        self.ani.save(name + '.mp4', writer=writer)




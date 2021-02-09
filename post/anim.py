
from post.plots import *
import matplotlib.animation as animation
from matplotlib.animation import FFMpegWriter

class pAnimation():
    def __init__(self,
                 P,
                 tfilePath,
                 name="animation",
                 interval=10,
                 frames='all',
                 xlim=False,
                 ylim=False):

        self.P = P # The plot object
        self.time = np.loadtxt(tfilePath)
        self.interval = interval
        if frames == 'all':
            self.frames = len(self.time)
        else:
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




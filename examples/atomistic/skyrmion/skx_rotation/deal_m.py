import matplotlib.pyplot as plt
import numpy as np
from moviepy.video.io.bindings import mplfig_to_npimage
from moviepy.editor import VideoClip
import moviepy.editor as mpy
import matplotlib.cm as cm

duration = 10

def read_data(i):
    m = np.load('dyn_npys/m_%d.npy'%i)
    m.shape = (-1,3)
    mx = m[:,2]
    mx.shape=(140,140)
    return mx

fig_mpl, ax = plt.subplots(1,figsize=(5,5), facecolor='white')

mx = read_data(0)
myobj = ax.imshow(np.transpose(mx),origin='lower',cmap=cm.rainbow)
ax.set_title('t= 0 ns')

 
# ANIMATE WITH MOVIEPY (UPDATE THE CURVE FOR EACH t). MAKE A GIF.

def make_frame_mpl(t):
    myobj.set_data(read_data(int(t*20)))
    ax.set_title('t= %0.2f ns'%(t*40*0.01))
    #fig_mpl.savefig('test.png')
    #line.set_ydata()  # <= Update the curve
    return mplfig_to_npimage(fig_mpl) # RGB image of the figure
 
animation = mpy.VideoClip(make_frame_mpl, duration=duration)
animation.write_gif("skx.gif", fps=20)

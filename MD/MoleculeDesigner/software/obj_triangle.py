import matplotlib.transforms
import matplotlib.path
import numpy as np

def triangle(v_0,v_f,ax, start, end, opacity, lw=3):
    
    saturacion_superior=1
    saturacion_inferior=0
    
    if (abs(v_0)>saturacion_superior):
        v_0=saturacion_superior*(v_0/abs(v_0))
        
    if (abs(v_0)<saturacion_inferior):
        v_0=saturacion_inferior*(v_0/abs(v_0))
        
    if (abs(v_f)>saturacion_superior):
        v_f=saturacion_superior*(v_f/abs(v_f))
        
    if (abs(v_f)<saturacion_inferior):
        v_f=saturacion_inferior*(v_f/abs(v_f))
    
    v_0=1-v_0
    v_f=1-v_f
    
    if (v_f<v_0):
        v_0=1-v_0
        v_f=1-v_f
        
    tricoords = [(0,-0.3),(1,0),(0,0.3),(0,-0.6)]
    angle = np.arctan2(end[1]-start[1],end[0]-start[0])
    rot = matplotlib.transforms.Affine2D().rotate(angle)
    tricoords2 = rot.transform(tricoords)
    tri = matplotlib.path.Path(tricoords2, closed=True)
    ax.scatter(end[0],end[1], c=1, s=(2*lw)**2, marker=tri, cmap='gist_gray_r',zorder=10, alpha=opacity)
    ax.autoscale_view()
        
    # Arrow head: Triangle
    tricoords = [(0,-0.06),(0.25,0),(0,0.06),(0,-0.6)]
    angle = np.arctan2(end[1]-start[1],end[0]-start[0])
    rot = matplotlib.transforms.Affine2D().rotate(angle)
    tricoords2 = rot.transform(tricoords)
    tri = matplotlib.path.Path(tricoords2, closed=True)
    ax.scatter(end[0],end[1], c=1, s=(3.3*lw)**2, marker=tri, cmap='gist_gray',zorder=10, alpha=opacity)
    ax.autoscale_view()


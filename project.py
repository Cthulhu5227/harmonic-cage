# harmonic cage coords l = 
# see README for details

import igl
import numpy as np
import pyvista as pv
from scipy.spatial import ConvexHull
from scipy.sparse.linalg import spsolve
from vtkmodules.vtkInteractionStyle import ( vtkInteractorStyleTrackballCamera, vtkInteractorStyleUser,) # should be in when install pyvista

#making the harmonic laplace solver
def harmonic (V, F, b, bc, l):
    laplacian = igl.cotmatrix(V,F).tocsr()# since k only every equals 1
    n = V.shape[0]
    
    free = np.setdiff1d(np.arange(n), b)
    
    free_rows = laplacian[free]                 # rows = the equations to solve
    Qff = free_rows[:, free].tocsc()    # free-vs-free, the matrix we invert
    Qfb = free_rows[:, b]               # free-vs-fixed 
    
    harmony = np.zeros((n, bc.shape[1]))
    sol = spsolve(Qff, -(Qfb @ bc))
    harmony[free] = sol.reshape(free.size, bc.shape[1])
    
    return harmony
    



models = ["models/decimated-knight.off"]  # maybe impliment switching later

# Part 0, load the mesh lol
V, F = igl.read_triangle_mesh(models[0])

# Part 1, Build a cage from the mesh vertices
hull = ConvexHull(V) # there was no specifc way to build a cage, so i just used a convex hull for an easy and quick one
b = hull.vertices
V_cage_unique = V[b].copy()
# cage floats off the mesh surface
centr = V_cage_unique.mean(axis=0)
dirs = V_cage_unique - centr
dirs/= np.linalg.norm(dirs, axis=1, keepdims=True)
diag =np.linalg.norm(V.max(axis=0) - V.min(axis=0))
V_cage_unique += dirs*(diag * 0.03)
faces_pv = np.hstack([np.full((F.shape[0], 1), 3), F])
cage_pts = V_cage_unique.copy() # mutable cage state


bc = np.eye(len(b))
# Part 2, Smooth functions over the mesh / solutions to Laplace equation 
W  = harmonic(V, F, b, bc, 1) # technically this isn't the exact way we solve in the paper but equivalent, was unable to find solver from og paper
#print(f"Weights computed Shape: {W.shape}")
 
# Part 4, Use coordinate functions to interpolate a cage deformation
def deformed_verts():
    return V + W @ (cage_pts - V_cage_unique)

# Finally!, Plot the thing
plotter = pv.Plotter()
plotter.add_text("D to toggle drag mode or R to reset cage",
                 position="lower_left", font_size=14, color="red")

mesh_data = pv.PolyData(deformed_verts(), faces_pv)
cage_data = pv.PolyData(cage_pts.copy())

plotter.add_mesh(mesh_data, color="lightblue", show_edges=False)
plotter.add_points(cage_data, color="red", point_size=20, render_points_as_spheres=True)


# below is implimentation of real time dragging, honestly this took longer then the project itself and was stupid to attempt
# state machine  
drag={"mode": False, "active": False, "selected": None, "last": None}
mode_actor = [None]

trackball= vtkInteractorStyleTrackballCamera()
no_cam = vtkInteractorStyleUser()

# Helpers 
def find_nearest_cage_pt(sx, sy):
    renderer = plotter.renderer
    def d2w(depth): #  # depth to world cords
        renderer.SetDisplayPoint(sx, sy, depth)
        renderer.DisplayToWorld()
        wp = renderer.GetWorldPoint()
        w= wp[3] if wp[3] else 1.0
        return np.array([wp[0] / w,wp[1]/ w, wp[2]/ w])
    near= d2w(0.0)
    ray= d2w(1.0) - near #pew
    rlen=np.linalg.norm(ray)
    if rlen < 1e-10:
        return None
    ray/= rlen
    best,best_d = None, float("inf")
    for i,pt in enumerate(cage_pts):
        v =pt - near
        closest = near + np.dot(v, ray) * ray
        d= np.linalg.norm(pt -closest)
        if d <best_d:
            best_d, best = d,i
    diag=np.linalg.norm(V.max(axis=0)-V.min(axis=0))
    threshold= diag * 0.05
    print(f"  [pick] best pt={best}  dist={best_d:.4f}  threshold={threshold:.4f}")
    return best if best_d < threshold else None


def unproject(sx, sy, ref): # like the opposite of projecting to a screen lol
    renderer = plotter.renderer
    cam = renderer.GetActiveCamera()
    view_dir =np.array(cam.GetFocalPoint()) - np.array(cam.GetPosition())
    view_dir/= np.linalg.norm(view_dir)
    def d2w(depth): # depth to world cords - copied from above
        renderer.SetDisplayPoint(sx, sy, depth)
        renderer.DisplayToWorld()
        wp= renderer.GetWorldPoint()
        w= wp[3] if wp[3] else 1.0
        return np.array([wp[0] / w, wp[1]/w, wp[2]/ w])
    near=d2w(0.0)
    ray=d2w(1.0) - near
    denom=np.dot(ray, view_dir)
    if abs(denom) < 0.000001: # avoid dividing by zero, Its bad fun fact
        return ref.copy()
    t=np.dot(ref-near,view_dir)/denom
    return near+t*ray

def on_left_press(caller, event):
    if not drag["mode"]:
        return
    sx, sy = caller.GetEventPosition()
    print(f"[press] screen=({sx},{sy})")
    idx = find_nearest_cage_pt(sx, sy)
    if idx is not None:
        drag["selected"] =idx
        drag["active"] =True
        drag["last"] =(sx, sy)
       

def on_mouse_move(caller, event):
    if not drag["active"]: # this im sure is sub optimal but first thing that came to mind
        return
    sx, sy = caller.GetEventPosition()
    lx, ly = drag["last"]
    ref =cage_pts[drag["selected"]].copy()
    delta= unproject(sx, sy, ref) - unproject(lx, ly, ref)
    cage_pts[drag["selected"]] += delta
    drag["last"] = (sx, sy)
    mesh_data.points = deformed_verts()
    cage_data.points = cage_pts.copy()

    plotter.render()

def on_left_release(caller, event):
    drag["active"]= False
    drag["selected"] =None

#Wire up interact 
iren = plotter.iren.interactor
iren.AddObserver("LeftButtonPressEvent",on_left_press,1.0)
iren.AddObserver("MouseMoveEvent",on_mouse_move,1.0)
iren.AddObserver("LeftButtonReleaseEvent",on_left_release,1.0)
iren.SetInteractorStyle(trackball)

def toggle_drag():
    drag["mode"] = not drag["mode"]
    if drag["mode"]:
        iren.SetInteractorStyle(no_cam)
        mode_actor[0] = plotter.add_text("DRAG MODE — click+drag red points",
            position="upper_right", font_size=11, color="yellow",)
        print("DRAG MODE")
    else:
        iren.SetInteractorStyle(trackball)
        if mode_actor[0]:
            plotter.remove_actor(mode_actor[0])
            mode_actor[0] =None
        drag["active"]= False
        drag["selected"]= None
        print("VIEW MODE")
    plotter.render()

def reset_cage(): # r key to reset 
    cage_pts[:]= V_cage_unique
    mesh_data.points= deformed_verts()
    cage_data.points =cage_pts.copy()
    cage_edges_data.points= cage_pts.copy()
    plotter.render()

# setup the key events. I forgot to do this and was so confused why this wouldnt work for a few hours. Feel like I need to note this. 
plotter.add_key_event("d", toggle_drag)
plotter.add_key_event("r", reset_cage)
plotter.show()


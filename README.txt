Harmonic Coordinates for Character Articulation Paper Implementation 


An implementation of the paper "Harmonic Coordinates for Character Articulation"
(Joshi et al.). A cage is built around a 3-D mesh, harmonic weights
are solved once using libigl, and the mesh deforms in real time as you drag cage
control points.

Reference papers are included:
  - Harmonic_paper.pdf

Setup (Conda was used for me)

1. Create and activate a conda environment:
     conda create -n harmonic python=3.11
     conda activate harmonic

2.Install core dependencies via conda:
     conda install numpy scipy matplotlib

3.install libigl and PyVista via pip:
     pip install libel
     pip install pyvista


Running the project.....


in the environment run
python project.py

The viewer will open showing the decimated-knight mesh with cage control
points around it.

Controls!

  D Toggle drag mode / view mode
     - View mode (default): rotate, pan, and zoom the camera freely
     - Drag mode: click and drag the red cage points to deform the mesh

  R  Reset the cage to its original position (undoes all deformation)


---- File overview ----

  project.py          Main script — loads mesh, builds cage, solves weights,
                      runs the interactive viewer
  models/             OFF mesh files used as input
  Harmonic_paper.pdf  Original paper this project implements

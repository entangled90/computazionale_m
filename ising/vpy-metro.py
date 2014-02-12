# vising.py
# Simulates the two-dimensional Ising model using the Metropolis algorithm
# This version uses Visual Python, www.vpython.org
# By Dan Schroeder, Weber State University, January 2013
# Modified for VPython 6 compatibility (with help from Bruce Sherwood), July 2013

import visual, visual.controls, numpy, random, math

size = 50                            # number of sites in a lattice row (change if desired)
s = numpy.ones((size, size), int)    # 2D array of dipoles (1=up, -1=down)
T = 2.27                             # temperature in natural units (adjusted by GUI slider)
running = False                      # will be true when simulation is running

# Set up the main graphics window (see www.vpython.org/contents/docs/visual/display.html):
windowSize = 500                           # size in screen pixels
visual.scene.title = "Ising Model"
visual.scene.x = 50                        # move window away from left edge of screen
visual.scene.width = windowSize
visual.scene.height = windowSize + 20      # height in pixels includes title bar; 20 is a guess
visual.scene.fov = 0.01                    # small field of view eliminates 3D perspective
visual.scene.userzoom = visual.scene.userspin = False    # no changing the view!

# Set up the GUI controls (see www.vpython.org/contents/docs/visual/controls.html):
cWindow = visual.controls.controls(title="Control", x=windowSize+60, y=0, width=150, height=visual.scene.height)
pButton = visual.controls.button(pos=(0,65), width=40, height=23, text="Run", action=lambda: runPause())
visual.label(display=cWindow.display, text="Temperature:", box=False, opacity=0, pos=(0,35))
tLabel = visual.label(display=cWindow.display, text="2.27", box=False, opacity=0, pos=(0,25))
tSlider = visual.controls.slider(pos=(0,-80), axis=(0,1,0), length=90, min=0.01, max=10.0, value=2.27, action=lambda: tSliderAdjust())

# Function to run or pause the simulation:
def runPause():
    global running
    running = not running
    pButton.text = "Pause" if running else "Resume"

# Function to update temperature readout when slider is moved:
def tSliderAdjust():
    global T
    T = tSlider.value
    tLabel.text = "%1.2f" % T        # display temperature to two decimal places

# We represent the square 2D array of lattice sites using a vpython "points"
# object, which is basically a 1D array of squares (or circles) placed in the
# 3D space.  The code below arranges them on a grid in the xy plane.
thePoints = visual.points(size=1, size_units="world", shape="square")    # can also try "round"
for i in range(size):
    for j in range(size):
        thePoints.append(pos=(i-size/2+0.5, j-size/2+0.5, 0))

# Function to color the "point" representing site (i,j):
def colorSquare(i, j):
    thePoints.color[i*size + j] = (0.5,0,1) if s[i,j]==1 else (1,1,1)    # purple and white

# Initialize the lattice to a random array, and draw it as we go:
for i in range(size):
    for j in range(size):
        s[i,j] = 1 if random.random()<0.5 else -1
        colorSquare(i,j)

# Function to calculate energy change upon hypothetical flip (with pbc):
def deltaE(i,j):
    leftS = s[size-1,j] if i==0 else s[i-1,j]
    rightS = s[0,j] if i==size-1 else s[i+1,j]
    topS = s[i,size-1] if j==0 else s[i,j-1]
    bottomS = s[i,0] if j==size-1 else s[i,j+1]
    return 2.0 * s[i,j] * (leftS + rightS + topS + bottomS)

# Main simulation loop:
while True:
    if running:
        visual.rate(1000000)             # handle events and run really fast
        i = int(random.random()*size)    # choose a random row and column
        j = int(random.random()*size)
        eDiff = deltaE(i,j)
        if eDiff <= 0 or random.random() < math.exp(-eDiff/T):    # Metropolis!
            s[i,j] = -s[i,j]
            colorSquare(i, j)
    else:
        visual.rate(20)        # don't eat processor time when paused

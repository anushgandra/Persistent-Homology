import numpy as np
from ripser import ripser
from persim import plot_diagrams
import matplotlib.pyplot as plt
import time
import matplotlib.animation as animation

def drawLineColored(X, C):
    for i in range(X.shape[0]-1):
        plt.plot(X[i:i+2, 0], X[i:i+2, 1], c=C[i, :], lineWidth = 3)

def plotCocycle2D(D, X, cocycle, thresh):
    """
    Given a 2D point cloud X, display a cocycle projected
    onto edges under a given threshold "thresh"
    """
    #Plot all edges under the threshold
    N = X.shape[0]
    t = np.linspace(0, 1, 10)
    c = plt.get_cmap('Greys')
    C = c(np.array(np.round(np.linspace(0, 255, len(t))), dtype=np.int32))
    C = C[:, 0:3]

    for i in range(N):
        for j in range(N):
            if D[i, j] <= thresh:
                Y = np.zeros((len(t), 2))
                Y[:, 0] = X[i, 0] + t*(X[j, 0] - X[i, 0])
                Y[:, 1] = X[i, 1] + t*(X[j, 1] - X[i, 1])
                drawLineColored(Y, C)
    #Plot cocycle projected to edges under the chosen threshold
    for k in range(cocycle.shape[0]):
        [i, j, val] = cocycle[k, :]
        if D[i, j] <= thresh:
            [i, j] = [min(i, j), max(i, j)]
            a = 0.5*(X[i, :] + X[j, :])
            plt.text(a[0], a[1], '%g'%val, color='b')
    #Plot vertex labels
    for i in range(N):
        plt.text(X[i, 0], X[i, 1], '%i'%i, color='r')
    #plt.axis('equal')

def pulldata(datafile):

    file = open(datafile,"r")

    lines = file.readlines()
    dictionary = {}

    xmin = np.Inf
    ymin = np.Inf
    xmax = -np.Inf
    ymax = -np.Inf

    for line in lines:
        curr = line.split()
        frame_id = int(float(curr[0]))
        x_pos = float(curr[2])
        y_pos = float(curr[4])

        if frame_id in dictionary:
            dictionary[frame_id].append([x_pos, y_pos])
        else:
            dictionary[frame_id] = [[x_pos,y_pos]]

        if x_pos > xmax:
            xmax = x_pos

        if y_pos > ymax:
            ymax = y_pos

        if x_pos < xmin:
            xmin = x_pos

        if y_pos < ymin:
            ymin = y_pos

    bounds = np.array([xmin, xmax, ymin, ymax])

    file.close()
    return dictionary, bounds

def plotPedestrianCocycles(state, bounds):

    xmin = bounds[0]
    xmax = bounds[1]
    ymin = bounds[2]
    ymax = bounds[3]

    diagrams = ripser(state)['dgms']
    result = ripser(state, coeff=17, do_cocycles=True)
    diagrams = result['dgms']
    cocycles = result['cocycles']
    D = result['dperm2all']

    dgm1 = diagrams[1]
    if len(dgm1):
        idx = np.argmax(dgm1[:, 1] - dgm1[:, 0])

        cocycle = cocycles[1][idx]
        thresh = dgm1[idx, 0]  #Project cocycle onto edges less than or equal to death time
        plotCocycle2D(D, state, cocycle, thresh)        

    plt.scatter(state[:,0], state[:,1], c='c', s=16)  
    plt.ylim([ymin, ymax])
    plt.xlim([xmin, xmax])      
    fig.canvas.draw()


datafile = "seq_eth/obsmat.txt"
dictionary, bounds = pulldata(datafile)

plt.figure()
fig = plt.gcf()
fig.show()
fig.canvas.draw()

for key in dictionary:
    state = np.array(dictionary[key])
    plotPedestrianCocycles(state, bounds)
    plt.pause(0.1)
    plt.clf()









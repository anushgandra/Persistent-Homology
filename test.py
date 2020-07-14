import numpy as np
from ripser import ripser
from persim import plot_diagrams
import matplotlib.pyplot as plt
import tadasets
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
    plt.axis('equal')

file = open("seq_eth/obsmat.txt","r")

lines = file.readlines()
dictionary = {}
# [frame_number pedestrian_ID pos_x pos_z pos_y v_x v_z v_y ]
ped_ids = {}

for line in lines:
    curr = line.split()
    frame_id = float(curr[0])
    x_pos = float(curr[2])
    y_pos = float(curr[4])
    pedestrian_id = float(curr[1])

    if frame_id in dictionary:
        dictionary[frame_id].append([x_pos, y_pos])
    else:
        dictionary[frame_id] = [[x_pos,y_pos]]

    if pedestrian_id in ped_ids:
        ped_ids[pedestrian_id].append([x_pos, y_pos])
    else:
        ped_ids[pedestrian_id] = [[x_pos,y_pos]]

file.close()

fig = plt.figure()
ims = []

for key in dictionary:
    data = np.array(dictionary[key])
    ims.append([plt.scatter(data[:,0], data[:,1], c='k', s=16)])

anim = animation.ArtistAnimation(fig,ims)
ax = plt.gca()
ax.set_aspect(1)
plt.show()
       
maxlen = 0
key = 0

for x in dictionary:
   if(len(dictionary[x]) > maxlen):
       maxlen = len(dictionary[x])
       key = x

data = np.array(dictionary[key])
plt.figure()
plt.scatter(data[:,0],data[:,1])
plt.xlabel("x coordinate")
plt.ylabel("y coordinate")
plt.title("Pedestrian positions")

diagrams = ripser(data)['dgms']
plt.figure()
plot_diagrams(diagrams, show=True)

result = ripser(data, coeff=17, do_cocycles=True)
diagrams = result['dgms']
cocycles = result['cocycles']
D = result['dperm2all']

dgm1 = diagrams[1]
idx = np.argmax(dgm1[:, 1] - dgm1[:, 0])
plt.figure()
plot_diagrams(diagrams, show = False)
plt.scatter(dgm1[idx, 0], dgm1[idx, 1], 20, 'k', 'x')
plt.title("Max 1D birth = %.3g, death = %.3g"%(dgm1[idx, 0], dgm1[idx, 1]))
plt.show()

cocycle = cocycles[1][idx]
thresh = dgm1[idx, 0]  #Project cocycle onto edges less than or equal to death time
plt.figure()
plotCocycle2D(D, data, cocycle, thresh)
plt.title("1-Form Thresh=%g"%thresh)
plt.show()

    
    
    

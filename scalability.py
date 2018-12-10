from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from matplotlib import colors as mcolors
import numpy as np
import sys

def cc(arg):
    return mcolors.to_rgba(arg, alpha=0.6)


if ( __name__ == '__main__' ):
    name = sys.argv[1]
    print("load data from ", name)

    p_list = []
    n_list = []
    t_list = [] 
    with open(name) as f:
        for line in f:
            p = int(line.split()[0])
            p_list.append(p)
            n = int(line.split()[1])
            n_list.append(n)
            t = float(line.split()[2])
            t_list.append(t)
    

    N = len(t_list)
    print("Total size ", N)

    time1 = [t_list[i] for i in range(N) if p_list[i] == 1]
    N1 = len(time1)
    print("Szie of 1p ", N1)

    sca_list = t_list.copy()
    perf_list = t_list.copy()
    print("Matrix & Cores & Time (sec) & SpeedUp & Perfomance (GFLOPS) \\\\")
    for k in range(N1):
        n = n_list[k]
        t = t_list[k]
        top = 2.0/3.0 * n**3
        #print("Matrix %i, total time %f, total operations %f" % (n, t, top))
        for i in range(N):
            if n_list[i] == n:
                sca_list[i] = t / t_list[i]
                perf_list[i] = top / t_list[i] / 10**9
                #print("Proc %i time %f, speedUp %f, perfomance %f" % (p_list[i],  t_list[i], sca_list[i], perf_list[i]))
                # Latex Output
                print("%i & %i & %f & %f & %f \\\\" %  (n, p_list[i], t_list[i], sca_list[i], perf_list[i]))
                print(" ")
        #print("---------------------")

    '''
    # SpeedUp
    fig = plt.figure()
    ax = fig.gca(projection='3d')
    plt.hold(True)
    
    X = np.array(p_list)
    Y = np.array(n_list)
    Z = np.array(sca_list)
    ax.scatter(X, Y, Z)

    ax.plot_trisurf(X, Y, Z)

    ax.set_xlabel('Cores')
    ax.set_ylabel('Matrix size')
    ax.set_zlabel('Speedup')
    ax.set_ylim(4096, 32768)

    #X, Y = np.meshgrid(X, Y)
    #ax.plot_surface(X,Y,Z)
   
    #list1, list2 = zip(*sorted(zip(list1, list2)))
    #Y, X, Z = zip(*sorted(zip(n_list, p_list, sca_list)))
    #verts = [list(zip(X, Y, Z))]
    #ax.add_collection3d(Poly3DCollection(verts), zs=Z)
    
    plt.show()
    '''

    # Perfomance
    fig = plt.figure()
    ax = fig.gca(projection='3d')
    plt.hold(True)
    
    X = np.array(p_list)
    Y = np.array(n_list)
    Z = np.array(perf_list)
    ax.scatter(X, Y, Z)

    ax.plot_trisurf(X, Y, Z, color="purple")

    ax.set_xlabel('Cores')
    ax.set_ylabel('Matrix size')
    ax.set_zlabel('Perfomance, GFLOPS')
    ax.set_ylim(4096, 32768)

    #X, Y = np.meshgrid(X, Y)
    #ax.plot_surface(X,Y,Z)
   
    #list1, list2 = zip(*sorted(zip(list1, list2)))
    #Y, X, Z = zip(*sorted(zip(n_list, p_list, sca_list)))
    #verts = [list(zip(X, Y, Z))]
    #ax.add_collection3d(Poly3DCollection(verts), zs=Z)
    
    plt.show()


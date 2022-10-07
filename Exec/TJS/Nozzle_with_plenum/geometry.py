import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

def main():

    # File coords name
    cfile = './noz_cont1.csv'

    # Read File
    coords = pd.read_csv(cfile, header=None).to_numpy()
    x = coords[:, 0]
    y = coords[: ,1]

    # Get top/bottom indices
    ypos = y[(y > 0)]
    yneg = y[(y < 0)]

    xpos = x[(y > 0)]
    xneg = x[(y < 0)]

    # Search for indices of triangle
    gradind = []
    ygrad = np.gradient(ypos)
    sign = 1
    for i, grad in enumerate(ygrad):
        if grad != 0 and np.sign(grad) != sign:
            gradind.append(i)
            sign*=-1
        elif grad == 0 and np.sign(grad) != np.sign(ygrad[i-1]):
            gradind.append(i)
            sign*=-1            
    
    print(xpos[gradind])
    print(ypos[gradind])

    # Change middle to A* location
    choke = np.where(ypos == min(ypos))
    choke = int(choke[0])
    gradind[1] = choke

    # File Echo
    print('Echo Coordiate Information')
    print('--> Dataset dims = [%i, %i]' % (len(coords[:,0]), len(coords[0,:])))
    print('--> x_lo = [%.2f, %.2f]' % (min(x), min(y)))
    print('--> x_hi = [%.2f, %.2f]\n\n' % (max(x), max(y)))

    # Plot profile
    fig, ax = plt.subplots(1, 1)
    ax.plot(xpos, ypos, c='k', lw=1.0)
    ax.plot(xneg, yneg, c='k', lw=1.0)
    ax.scatter(xpos[gradind], ypos[gradind], c='r')
    ax.plot(xpos[gradind], ypos[gradind], c='r')

    # Basic format
    ax.set_xlim([-5, 10])
    plt.show()

def main2():

    # File coords name
    cfile = './noz_cont1.csv'

    # Read File
    coords = pd.read_csv(cfile, header=None).to_numpy()
    x = coords[:, 0]
    y = coords[: ,1]

    # Convert to meters
    x /= 1000
    y /= 1000

    # Get top/bottom indices
    ypos = y[(y > 0)]
    yneg = y[(y < 0)]

    xpos = x[(y > 0)]
    xneg = x[(y < 0)]

    # Find A*
    ymin = min(ypos)
    Astar = 3.14159 * ymin**2.0
    print('A* = %.4e m^2' % Astar)

    # Find inlet area
    y_inlet = ypos[0]
    Ain = 3.14159 * y_inlet**2.0
    print('A_in = %.4e m^2' % Ain)
    print('Area Ratio = %.4f' % (Ain/Astar))

    # Inlet Mach computation
    gamma = 1.289
    

if __name__ == "__main__": 
    
    main2()
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np

def epsilon(T):
	return((-0.8292e-6 *(T**3)) + (0.1417e-2 * (T**2)) - (0.9297*T) + (233.76 + 5321/T))
	
X = np.linspace(200,400, num=1000)
Y = []
for x in X:
	Y.append(epsilon(x))

plt.plot(X, Y)
plt.show()

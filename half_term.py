import matplotlib.pyplot as plt
import numpy as np
from math import sin, cos

def E(t, a, b):
    return 3*t**2 - a*t**3 - b*cos(t)

def F(t, c ,d):
    return -1 * sin(t) - d*cos(t) - c*t**3

tau = 0.0001
T = 1
t = np.arange(0, T, tau)

x_t1 = [0]*round((T)/tau)
y_t1 = [1]*round((T)/tau)

epsx1 = 0
epsy1 = 0

a = 1
b = 1
c = 0
d = 1


x_t2 = [0]*round((T)/tau)
y_t2 = [1]*round((T)/tau)

epsx2 = 0
epsy2 = 0

x_t3 = [0]*round((T)/tau)
y_t3 = [1]*round((T)/tau)

epsx3 = 0
epsy3 = 0

# (x[n+1] - x[n])/TAU = f(x[n], t[n])
for i in range(1, round((T)/tau)):
    x_t1[i] = x_t1[i-1] + tau*(a * x_t1[i-1] + b * y_t1[i-1] + E(t[i-1], a, b))
    y_t1[i] = y_t1[i-1] + tau*(c * x_t1[i-1] + d * y_t1[i-1] + F(t[i-1], c, d))
    epsx1 = max(epsx1, abs(x_t1[i] - t[i]**3))
    epsy1 = max(epsy1, abs(y_t1[i] - cos(t[i])))

# (x[n+1] - x[n])/TAU = 1/2(f(x[n], t[n]) + f(x[n+1], t[n+1]))
for i in range(1, round((T)/tau)):
    y_t2[i] = (tau * (y_t2[i-1] + F(t[i-1], c, d) + F(t[i], c, d)) + 2 * y_t2[i-1]) / (2 - tau)
    x_t2[i] = (tau * (x_t2[i-1] + y_t2[i-1] + E(t[i-1], a, b) + y_t2[i] + E(t[i], a, b)) + 2 * x_t2[i-1]) / (2 - tau)
    epsx2 = max(epsx2, abs(x_t2[i] - t[i]**3))
    epsy2 = max(epsy2, abs(y_t2[i] - cos(t[i])))
        
# (x[n+1] - x[n])/TAU = 1/2(f(x[n], t[n]) + f(x*[n+1], t[n]))
# (x*[n+1] - x[n])/TAU = f(x[n], t[n])
for i in range(1, round((T)/tau)):
    xh = x_t3[i-1] + tau*(a * x_t3[i-1] + b * y_t3[i-1] + E(t[i-1], a, b))
    yh = y_t3[i-1] + tau*(c * x_t3[i-1] + d * y_t3[i-1] + F(t[i-1], c, d))
    x_t3[i] = x_t3[i-1] + tau / 2 * (a * x_t3[i-1] + b * y_t3[i-1] + E(t[i-1], a, b) + a * xh + b * yh + E(t[i], a, b))
    y_t3[i] = y_t3[i-1] + tau / 2 * (c * x_t3[i-1] + d * y_t3[i-1] + F(t[i-1], c, d) + c * xh + d * yh + F(t[i], c, d))
    epsx3 = max(epsx3, abs(x_t3[i] - t[i]**3))
    epsy3 = max(epsy3, abs(y_t3[i] - cos(t[i])))

# Creating window and charts
fig, ax = plt.subplots(2, 2, figsize=(10, 5.4))
plt.setp(ax, xlim=(0, T), ylim=(-4, 8))

# True functions
ax[0, 0].plot(t, t**3, label="true x(t)")
ax[0, 0].plot(t, np.cos(t), label="true y(t)")
ax[0, 0].set_title("True functions")


# Method 1
ax[0, 1].plot(t, x_t1, label="x(t)")
ax[0, 1].plot(t, y_t1, label="y(t)")
ax[0, 1].set_title("Method 1")
ax[0, 1].text(0.25, 4, "epsx = "+str(epsx1))
ax[0, 1].text(0.25, 4.5, "epsy = "+str(epsy1))

# Method 2
ax[1, 0].plot(t, x_t2, label="x(t)")
ax[1, 0].plot(t, y_t2, label="y(t)")
ax[1, 0].set_title("Method 2")
ax[1, 0].text(0.25, 4, "epsx = "+str(epsx2))
ax[1, 0].text(0.25, 4.5, "epsy = "+str(epsy2))

# Method 3
ax[1, 1].plot(t, x_t3, label="x(t)")
ax[1, 1].plot(t, y_t3, label="y(t)")
ax[1, 1].set_title("Method 3")
ax[1, 1].text(0.25, 4, "epsx = "+str(epsx3))
ax[1, 1].text(0.25, 4.5, "epsy = "+str(epsy3))

for ax in ax.flat:
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.legend()
fig.tight_layout()
plt.show()

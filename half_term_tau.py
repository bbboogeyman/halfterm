import matplotlib.pyplot as plt
import numpy as np
from math import cos, sin

def E(t, a, b):
    return 3*t**2 - a*t**3 - b*cos(t)

def F(t, c ,d):
    return -1 * sin(t) - d*cos(t) - c*t**3

a = 1
b = 1
c = 0
d = 1


TAU_mass = np.arange(0.0001, 0.002, 0.0001)
T_right = 1
T_left = -1

epsx1 = []
epsy1 = []
epsx2 = []
epsy2 = []
epsx3 = []
epsy3 = []
epsx4 = []
epsy4 = []


for tau in TAU_mass:
    t = np.arange(T_left, T_right, tau)
    
    epsx1t = 0
    epsy1t = 0
    epsx2t = 0
    epsy2t = 0
    epsx3t = 0
    epsy3t = 0
    
    x_t1 = [T_left ** 3] * round((T_right - T_left) / tau)
    y_t1 = [cos(T_left)] * round((T_right - T_left) / tau)

    x_t2 = [T_left ** 3] * round((T_right - T_left) / tau)
    y_t2 = [cos(T_left)] * round((T_right - T_left) / tau)


    x_t3 = [T_left ** 3] * round((T_right - T_left) / tau)
    y_t3 = [cos(T_left)] * round((T_right - T_left) / tau)

    # (x[n+1] - x[n])/TAU = f(x[n], t[n])
    for i in range(1, round((T_right - T_left)/tau)):
        x_t1[i] = x_t1[i-1] + tau*(a * x_t1[i-1] + b * y_t1[i-1] + E(t[i-1], a, b))
        y_t1[i] = y_t1[i-1] + tau*(c * x_t1[i-1] + d * y_t1[i-1] + F(t[i-1], c, d))
        epsx1t = max(epsx1t, abs(x_t1[i] - t[i]**3))
        epsy1t = max(epsy1t, abs(y_t1[i] - cos(t[i])))

    # (x[n+1] - x[n])/TAU = 1/2(f(x[n], t[n]) + f(x[n+1], t[n+1]))
    for i in range(1, round((T_right - T_left)/tau)):
        y_t2[i] = (tau * (y_t2[i-1] + F(t[i-1], c, d) + F(t[i], c, d)) + 2 * y_t2[i-1]) / (2 - tau)
        x_t2[i] = (tau * (x_t2[i-1] + y_t2[i-1] + E(t[i-1], a, b) + y_t2[i] + E(t[i], a, b)) + 2 * x_t2[i-1]) / (2 - tau)
        epsx2t = max(epsx2t, abs(x_t2[i] - t[i]**3))
        epsy2t = max(epsy2t, abs(y_t2[i] - cos(t[i])))
        
    # (x[n+1] - x[n])/TAU = 1/2(f(x[n], t[n]) + f(x*[n+1], t[n]))
    # (x*[n+1] - x[n])/TAU = f(x[n], t[n])
    for i in range(1, round((T_right - T_left)/tau)):
        xh = x_t3[i-1] + tau*(a * x_t3[i-1] + b * y_t3[i-1] + E(t[i-1], a, b))
        yh = y_t3[i-1] + tau*(c * x_t3[i-1] + d * y_t3[i-1] + F(t[i-1], c, d))
        x_t3[i] = x_t3[i-1] + tau / 2 * (a * x_t3[i-1] + b * y_t3[i-1] + E(t[i-1], a, b) + a * xh + b * yh + E(t[i], a, b))
        y_t3[i] = y_t3[i-1] + tau / 2 * (c * x_t3[i-1] + d * y_t3[i-1] + F(t[i-1], c, d) + c * xh + d * yh + F(t[i], c, d))
        epsx3t = max(epsx3t, abs(x_t3[i] - t[i]**3))
        epsy3t = max(epsy3t, abs(y_t3[i] - cos(t[i])))
    
    epsx1.append(np.log(epsx1t))
    epsy1.append(np.log(epsy1t))
    epsx2.append(np.log(epsx2t))
    epsy2.append(np.log(epsy2t))
    epsx3.append(np.log(epsx3t))
    epsy3.append(np.log(epsy3t))

# Creating window and chartsы
fig, ax = plt.subplots(2, 2, figsize=(10, 5.4))
plt.setp(ax)

# True functions
ax[0, 0].plot(t, t*t*t, label="true x(t)")
ax[0, 0].plot(t, np.cos(t), label="true y(t)")
ax[0, 0].set_title("True functions")
ax[0, 0].spines['left'].set_position('center')
ax[0, 0].spines['bottom'].set_position('center')

# Method 1
ax[0, 1].plot(np.log(TAU_mass), epsx1, label="log(eps1x)")
ax[0, 1].plot(np.log(TAU_mass), epsy1, label="log(eps1y)")
ax[0, 1].set_title("Method 1")
ax[0, 1].text(-5, -7.0, "kx = "+str(round(np.polyfit(np.log(TAU_mass), epsx1, 1)[0], 3)))
ax[0, 1].text(-5, -7.5, "ky = "+str(round(np.polyfit(np.log(TAU_mass), epsy1, 1)[0], 3)))

# Method 2
ax[1, 0].plot(np.log(TAU_mass), epsx2, label="log(eps2x)")
ax[1, 0].plot(np.log(TAU_mass), epsy2, label="log(eps2y)")
ax[1, 0].set_title("Method 2")
ax[1, 0].text(-5, -16, "kx = "+str(round(np.polyfit(np.log(TAU_mass), epsx2, 1)[0], 3)))
ax[1, 0].text(-5, -16.5, "ky = "+str(round(np.polyfit(np.log(TAU_mass), epsy2, 1)[0], 3)))
# Method 3
ax[1, 1].plot(np.log(TAU_mass), epsx3, label="log(eps3x)")
ax[1, 1].plot(np.log(TAU_mass), epsy3, label="log(eps3y)")
ax[1, 1].set_title("Method 3")
ax[1, 1].text(-5, -15, "kx = "+str(round(np.polyfit(np.log(TAU_mass), epsx3, 1)[0], 3)))
ax[1, 1].text(-5, -15.5, "ky = "+str(round(np.polyfit(np.log(TAU_mass), epsy3, 1)[0], 3)))

# Deleting frames, creating legend and painting charts
for ax in ax.flat:
    ax.set_aspect("equal")
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.legend()
fig.tight_layout()
plt.show()
print('Coefficient values for Method 1:')
print(np.polyfit(np.log(TAU_mass), epsx1, 1)[0], '\n')
print('Coefficient values for Method 2:')
print(np.polyfit(np.log(TAU_mass), epsx2, 1)[0], '\n')
print('Coefficient values for Method 3:')
print(np.polyfit(np.log(TAU_mass), epsx3, 1)[0], '\n')
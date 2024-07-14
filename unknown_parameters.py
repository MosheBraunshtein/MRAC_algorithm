import numpy as np
from scipy.integrate import odeint 
import matplotlib.pyplot as plt 

# define constants 

a1 = 1
a2 = 1
b1 = 1
b2 = 1
h = 2
k = 1 
lambda0 = 20
lambda1 = 10
gamma = 3
a = 1
b = 5
w = 0.01
x_dot_0 = 1
x_0 = 0.01
x_0_model = 0.3
x_m_0 = 0.15





def differential_equation_for_process_response(u,t):
    x_model, x_dot, x, h_estimate, a1_estimate, a2_estimate, b1_estiamate, b2_estimate = u
    x_model_dot =  b*np.sin(w*t) -a*x_model
    x_model_dot_dot = -a*x_model_dot + b*w*np.cos(w*t)
    tmp = -gamma*((x_dot-x_model_dot)+lambda0*(x-x_model))
    h_dot_estimate = tmp*(x_model_dot_dot-lambda0*(x_dot-x_model_dot))
    a1_dot_estimate = tmp*x_dot
    a2_dot_estimate = tmp*x
    b1_dot_estimate = tmp*((x_dot)**3)
    b2_dot_estimate = tmp*np.cos(x)
    h_error = h_estimate-h 
    a1_error = a1_estimate-a1 
    a2_error = a2_estimate-a2 
    b1_error = b1_estiamate-b1 
    b2_error = b2_estimate-b2 
    x_dot_dot = h_error*(x_model_dot_dot-lambda0*(x_dot-x_model_dot)) \
                + a1_error*(x_dot) \
                + a2_error*(x) \
                + b1_error*((x_dot)**3) \
                + b2_error*np.cos(x) \
                - k*((x_dot-x_model_dot) + lambda0*(x_model-x)) \
                + (x_model_dot_dot)  \
                + lambda1*(x_dot-x_model_dot)
    return [x_model_dot, x_dot_dot, x_dot,h_dot_estimate, a1_dot_estimate, a2_dot_estimate, b1_dot_estimate, b2_dot_estimate]



t = np.linspace(0, 100, 1000)
h_0_estimate , a1_0_estimate, a2_0_estimate, b1_0_estimate,b2_0_estimate = 5, 5, 5, 5, 5
initial_conditions = [x_0_model, x_dot_0, x_0,h_0_estimate, a1_0_estimate, a2_0_estimate, b1_0_estimate,b2_0_estimate]
results = odeint(differential_equation_for_process_response, initial_conditions,t)
x_model = results[:, 0]
x = results[:, 2]


h_estimate = results[:, 3]
h = [h]*1000
plt.figure(figsize=(10, 6))
plt.plot(t, h_estimate, label='h_estimate', color='red',linestyle='--')
plt.plot(t, h, label='h', color='red')

a1_estimate = results[:, 4]
a1 = [a1]*1000
plt.plot(t, a1_estimate, label='a1_estimate', color='green',linestyle='--')
plt.plot(t, a1, label='a1', color='green')

plt.xlabel('Time')
plt.ylabel('Values')
plt.title('parameter estimation')
plt.legend()
plt.grid(True)
plt.show()

plt.figure(figsize=(10, 6))
plt.plot(t, x_model, label='x_model')
plt.plot(t, x, label='x')
plt.xlabel('Time')
plt.ylabel('Values')
plt.title('x_model and x over time')
plt.legend()
plt.grid(True)
plt.show()


# Plot error : reference model output - process model response

plt.figure(figsize=(10, 6))
plt.plot(t, x_model-x, label='x_model')
plt.xlabel('Time')
plt.ylabel('error')
plt.title('x_model - x')
plt.legend()
plt.grid(True)
plt.show()


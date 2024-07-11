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
lambda0 = 1
gamma = 3
a = 1
b = 10
w = 0.2
x_dot_0 = 1
x_0 = 0.01
x_0_model = 0.3
procces_model_start_conditions = [0.1,1]
x_m_0 = 0.15
model_start_conditions = [0.15]




def differential_equation_for_process_response(u,t):
    x_model, x_dot, x = u
    x_model_dot =  b*np.sin(w*t) -a*x_model
    x_model_dot_dot = -a*x_model_dot + b*w*np.cos(w*t) 
    x_dot_dot = x_model_dot_dot - ((k/h)+lambda0)*x_dot + (lambda0+(k/h))*x_model_dot +(k*lambda0/h)*(x_model-x)
    return [x_model_dot, x_dot_dot, x_dot]

# function for 
# def compute_process_model_response():


# def compute_reference_model_output():


# plot reference model output & process model response

t = np.linspace(0, 100, 1000)
initial_conditions = [x_0_model, x_dot_0, x_0]
results = odeint(differential_equation_for_process_response, initial_conditions,t)
x_model = results[:, 0]
x = results[:, 2]

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


import matplotlib.pyplot as plt
import numpy as np 
import os

def create_field_plot(x, y, z, field_name, axe):
    axe.contour(x, y, z)
    axe.set(aspect = 1, xlabel = 'x', ylabel = 'y', title = field_name)
    axe.invert_yaxis()

def generate_case_name():
    return 'test'

def generate_case_dir():
    dir_results = './results'
    case_name = generate_case_name()
    path = os.path.join(dir_results, case_name)
    return path 

def generate_time_plot(file_name, field_name, axe):
    data = np.loadtxt(file_name)
    axe.plot(data[:,0], data[:,1])
    axe.set(xlabel = 'time', ylabel = field_name) 
    axe.grid()

def generate_field_plot(file_name, field_name, axe, size):
    x,y,z = np.loadtxt(file_name, unpack=True)
    x = x.reshape(size)
    y = y.reshape(size)
    z = z.reshape(size)
    create_field_plot(x, y, z, field_name, axe)

def savefig_case():
    size = (51, 51) 
    path = generate_case_dir()
    # Generate figure with axis
    fig, axes = plt.subplots(figsize=(10,10), nrows=2, ncols=2)

    time_evolution_file = os.path.join(path, 'out_psimax_t.dat')
    if os.path.exists(time_evolution_file):
        generate_time_plot(time_evolution_file, 'psi_max', axes.flat[0])

    field_psi_file = os.path.join(path, 'out_psi.dat')
    if os.path.exists(field_psi_file):
        generate_field_plot(field_psi_file, 'psi', axes.flat[1], size)

    field_T_file = os.path.join(path, 'out_T.dat')
    if os.path.exists(field_T_file):
        generate_field_plot(field_T_file, 'T', axes.flat[2], size)

    field_C_file = os.path.join(path, 'out_C.dat')
    if os.path.exists(field_C_file):
        generate_field_plot(field_C_file, 'C', axes.flat[3], size)

    #plt.show()
    fig.savefig(os.path.join(path, generate_case_name() + '_pic.png'))

if __name__ == '__main__':
    savefig_case()
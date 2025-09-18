# Basisc square wave: SX-->T1;  flat band --> T2;  T1 + T2 = T
from signal import signal
from scipy.special import jn_zeros
from itertools import combinations
import numpy as np
from qutip import *


def drive(t, args):
    w = args['omega']
    g = args['g']

    # Example: sin²(wt) function
    signl = g * np.sin(w/2  * t)**2
    return signl


# HAMILTONIAN   

def j_ij(Jvalue, i,j, beta):
    return Jvalue/(np.abs(i-j))**beta
 
def position_hamiltonian(args):
    N  =  args['N']
    er = args['er']
    Jvalue = args['J']
    beta = args['beta']
    sx,sy,sz = sigmax(), sigmay(), sigmaz()
    empt = qzero(2**N) + 1j * qzero(2**N)    
    H01, H02, H11 = empt,  empt, empt
    
    # for i in range(N-1):
    #     id = qeye(2**i)    
    #     dim11 = N-2-i
    #     id1 = qeye(2**dim11)
    #     H01 = H01 + Qobj(tensor(id,tensor(sz,tensor(sz,id1))).full())
        
    comb = combinations(np.arange(N), 2)
    for nm in list(comb):
        i,j= np.array(nm)
        id = qeye(2**i)
        dim11 = j-i-1
        id1 = qeye(2**dim11)
        dim12 = N-1-j
        id2 = qeye(2**dim12)
        H01 = H01 + Qobj(tensor(id, tensor(sz, tensor(id1, tensor(sz,id2)))).full())\
            * j_ij(Jvalue, i,j, beta)

        
    for i in range(N):
        id = qeye(2**i)    
        dim11 = N-1-i
        id1 = qeye(2**dim11)
        H02 = H02 + Qobj(tensor(id,tensor(sx,id1)).full()) * (1-er)
        
        
    for i in range(N):
        id = qeye(2**i)    
        dim11 = N-1-i
        id1 = qeye(2**dim11)
        H11 = H11 + Qobj(tensor(id,tensor(sz,id1)).full())

    return H01, H02, H11

# FLOQUET ANALYSIS

def floquet_return_position(args):
    N  =  args['N']
    er = args['er']
    H01, H02, H11 = position_hamiltonian(args)
    H = [H01,[H02, drive]]

    T = 2 * np.pi/args['omega']
    f_modes_0, f_energies = floquet_modes(H, T, args)
    return f_energies

# MAGNETIZATION
def magnetization_position(args):      
    N  =  args['N']
    er = args['er']
    H01, H02, H11 = position_hamiltonian(args)
    H = [H01,[H02, drive]]

    grket = basis(2**N,0)
    times = args['times']
    data = mesolve(H, grket, times, [], [H11/N], args = args)
    return data.expect

def magnetization_stroboscopic(args):      
    N  =  args['N']
    er = args['er']
    maxT = args['maxT']
    w = args['omega']
    T = 2 * np.pi/w
    times = np.arange(0, maxT + 1) * T
    H01, H02, H11 = position_hamiltonian(args)
    H = [H01,[H02, drive]]

    grket = basis(2**N,0)
    data = mesolve(H, grket, times, [], [H11/N], args = args)
    return data.expect

print("Functions defined")

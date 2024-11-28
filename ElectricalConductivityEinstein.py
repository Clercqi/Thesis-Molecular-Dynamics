#Script 5: Berekenen van de elektrische conductiviteit via Einstein (exacte vergelijking)
import numpy as np
import multiprocessing as mp
import matplotlib
matplotlib.use('Agg')

import matplotlib.pyplot as plt

from ase.io import read
from time import time

#Functie voor het bepalen van de volume van de eenheidscel
def get_volume():
    with open('log.lammps', 'r') as fv:
        # Only start reading after 'NPT run at operational temperature'
        v_data = []
        equilibration_run, read = False, False
        for line in fv:
            if read:
                line = line.split()
                if len(line) != 8: break
                v_data.append(float(line[4]))
            if 'NPT run at operational temperature' in line: equilibration_run = True
            if equilibration_run and ' Step ' in line: read = True
        v_data = np.array(v_data)
    print('Average volume: %.2f A^3'%np.mean(v_data))
    return np.mean(v_data)

#Functie voor het oproepen van de atoomladingen
def get_charges():
    charges = []
    with open('inHPC_final.lammps', 'r') as fc:
        for line in fc:
            if 'set type' in line:
                charges.append(float(line.split()[4]))
    charges = np.array(charges)
    print('Charges for (Si, Al, Ca, Mg, Ox):', charges)
    return charges

#Functie voor het dubbelproduct van de vergelijking
def dotp(d, c):
    return np.einsum('i,j,ix,jx', c, c, d, d)

#Functie voor de berekening van de elektrische conductiviteit
def get_EC(atomic_charges):
    traj = read('NVT_op.xyz', format='extxyz', index=':')
    num = traj[0].get_atomic_numbers()
    charges = [atomic_charges[n-1] for n in num]
    charges = np.array(charges)
    print('%d charges'%len(charges))

    pos = np.array([atoms.get_positions() for atoms in traj])
    N_frames = pos.shape[0]
    N_atoms = pos.shape[1]
    print('%d frames, %d atoms'%(N_frames, N_atoms))

    step = 10

    # shape (N_frames/2/step, N_frames/2/step, N_atoms, 3) !
    displace = []
    for deltat in np.arange(stop = N_frames // 2, step = step):
        print('Calc. displace for deltat = %d'%deltat)
        displace_deltat = [pos[t0+deltat] - pos[t0] for t0 in np.arange(stop = N_frames // 2, step = step)]
        displace.append(displace_deltat)
    displace = np.array(displace)
    print(displace.shape)


    dot_product = []
    for i, d_deltat in enumerate(displace):
        dot_products_t0 = []

        t = time()
        pool = mp.Pool(processes=32)
        results = [pool.apply_async(dotp, args=(d_deltat_t0, charges)) for d_deltat_t0 in d_deltat]
        dot_products_t0 = [p.get() for p in results]
        print('Calcs for 1 deltat took %.1f s'%(time()-t))

        dot_product.append(np.mean(dot_products_t0))
    dot_product = np.array(dot_product)
    print(dot_product.shape)

    t = np.array([step * 1000 * i for i in range(len(dot_product))])
    plt.plot(t, dot_product, color='blue', linestyle='-')
    plt.xlabel('timeframes (fs)', fontsize=14)
    plt.ylabel('dot product (A²*C²)', fontsize=14)
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    plt.legend(fontsize=14)
    plt.savefig('dot_product.png')
    
    #Selectie van het middelste bereik (5 - 50 %) waar de elektrische conductiviteit een lineare curve voorstelt
    bereik_timesteps = []
    bereik_gemiddelde_dot_product = []
    for i, element in enumerate(t):
        if element >= 250000 and element <= 2500000:
            bereik_timesteps.append(element) 
            bereik_gemiddelde_dot_product.append(dot_product[i])
    
    #Lineaire regressie (rechte fitten) = het nemen van de afgeleide naar de tijd
    coefficients_gemiddelde_dot_product = np.polyfit(bereik_timesteps, bereik_gemiddelde_dot_product, 1)
    slope_gemiddelde_dot_product = coefficients_gemiddelde_dot_product[0]
    intercept_gemiddelde_dot_product = coefficients_gemiddelde_dot_product[1]
    print('De helling is ', slope_gemiddelde_dot_product)
    print('De intercept is ', intercept_gemiddelde_dot_product)

    #Het dot product en afgeleide samen plotten
    plt.plot(bereik_timesteps, bereik_gemiddelde_dot_product, color='red', linestyle='solid', label='Dot product')
    regressie_lijn_dot_product = slope_gemiddelde_dot_product * np.array(bereik_timesteps) + intercept_gemiddelde_dot_product
    plt.plot(bereik_timesteps, regressie_lijn_dot_product, color ='black', linestyle='solid', label='Linear regression')
    plt.xlabel('timeframes (fs)', fontsize=14)
    plt.ylabel('dot product (A²*C²)', fontsize=14)
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    plt.legend(fontsize=14)
    plt.savefig('Dot_product_regressie.png')
    
    #Berekening van elektrische conductiviteit met de exacte vergelijking
    temperatuur = 1823 #uitgedrukt in Kelvin - Lammps input file
    boltzmann_constante = 1.38064852 * (10) ** (-23) #uitgedrukt in J/K - bron internet
    elementary_charge = 1.602176634 * (10) ** (-19) #uitgedrukt in coulomb -  bron internet
    elektrische_conductiviteit_exact = ((elementary_charge ** 2) * slope_gemiddelde_dot_product) / (6 * boltzmann_constante * temperatuur * get_volume() * 10 ** (-25))
    print('De elektrische conductiviteit berekent via de exacte vergelijking bedraagt ', elektrische_conductiviteit_exact)

get_volume()
atomic_charges = get_charges()
get_EC(atomic_charges)

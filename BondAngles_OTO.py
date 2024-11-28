#Script 6(b): Classificeren van de soorten zuurstof(bridging, nonbridging, free) en bepalen van de bond angle distributie van Si-O-Si, O-Si-O, Al-O-Si, O-Al-O, Al-O-Al
#Files: data.file en NVT_op.trj
import numpy as np
import os
from ase import Atoms
import matplotlib.pyplot as plt

plt.switch_backend('agg') #Stelt agg in als backend

#Deel1: Bepalen van de type's zuurstof en berekenen van de bond angle distributie
#File met de boundary conditions van de box
with open('NVT_op.trj', 'r') as boundary:
    regels = boundary.readlines()
    
#Selectie van de lijnen met de gegevens van de boundary conditions
start_line = 5
end_line = 8
selected_lines = regels[start_line:end_line]
minimum = [] #Een lijst met minimum boundary conditions in volgorde [x, y, z]
maximum = [] #Een lijst met maximum boundary conditions in volgorde [x, y, z]

with open('boundary_conditions.txt', 'w') as conditions:
    for rij in selected_lines:
        rij = rij.split()
        conditions.write(' '.join(rij) + '\n')

data_boundary = np.loadtxt('boundary_conditions.txt', dtype=float)        
for rij in data_boundary:
    minimum.append(float(rij[0]))
    maximum.append(float(rij[1]))
minimum = np.array(minimum)
maximum = np.array(maximum)
print('De minimum boundary conditions zijn ', minimum)
print('De maximum boundary conditions zijn ', maximum) 

#File NVT_op.trj bewerken voor de bepaling van de type's zuurstof
with open('NVT_op.trj', 'r') as bestand:
    lines = bestand.readlines()
    
#Verwijder de lijnen die geen lengte van 9 hebben en voor de eerste 100 snapshots (2144 atomen per snapshot)
lines = [regel for regel in lines if len(regel.split(' ')) == 12]
geselecteerde_lijnen = lines[:214400]

#Schrijf de lijnen naar een nieuw .trj bestand
with open('NVT_op_bewerkt.trj', 'w') as f:
    f.writelines(geselecteerde_lijnen)

#Inlezen van de file NVT_op_bewerkt.trj
data = np.loadtxt('NVT_op_bewerkt.trj', dtype=float, usecols=[0, 1, 5, 6, 7]) #Kolommen in de volgorde [atoom, ID, x, y, z]

#Definieren van de voorwaarden - bond lengths
CO_Si_Ox = 2.2825 #De waarden zijn gehaald uit de berekende bond lengths van de RDFs
CO_Al_Ox = 2.5025
timestep = 1000
total_time = 100000
frequentie = total_time / timestep

#Bepalen van de Al en Si tetraeders (polymerisatiegraad) en de bindingshoeken voor de bidningen O-Si-O en O-Al-O
aantal_tetra_Si = 0
aantal_tetra_Al = 0
bindingen_Ox_Si = 0
bindingen_Ox_Al = 0
coord_Ox_Si = []
coord_Ox_Al = []
hoek_Ox_Si_Ox = []
hoek_Ox_Al_Ox = []
k = 0
l = 1
for regel_n in data:
    atoomsoort = regel_n[1]
    if atoomsoort == 1 and l <= 2144:
        atoomreferentie_Si = [regel_n[2], regel_n[3], regel_n[4]]
        regelnummer_Si = int(regel_n[0])
        for regel_Si in data[regelnummer_Si + k:2144 + k]:
            atoomcontroletype = regel_Si[1]
            atoomcontrole_Si = [regel_Si[2], regel_Si[3], regel_Si[4]]
            positie_atomen_Si = [atoomreferentie_Si, atoomcontrole_Si]
            atoms_Si = Atoms(positions=positie_atomen_Si)
            atoms_Si.set_cell((float(maximum[0]), float(maximum[1]), float(maximum[2])))
            atoms_Si.set_pbc(True)
            afstand_Si = atoms_Si.get_distance(0, 1, mic=True)
            if atoomcontroletype == 5 and afstand_Si <= CO_Si_Ox:
                bindingen_Ox_Si += 1
                if bindingen_Ox_Si == 4: #Berekenen van de hoeken O-Si-O
                    aantal_tetra_Si += 1
                    coord_Ox_Si.append(atoomcontrole_Si)
                    coord_Ox_Si = np.array(coord_Ox_Si)
                    vector1_Ox_Si = coord_Ox_Si[0] - atoomreferentie_Si
                    vector2_Ox_Si = coord_Ox_Si[1] - atoomreferentie_Si
                    vector3_Ox_Si = coord_Ox_Si[2] - atoomreferentie_Si
                    vector4_Ox_Si = coord_Ox_Si[3] - atoomreferentie_Si
                    cos_theta_Ox_Si_Ox_1 = np.dot(vector1_Ox_Si, vector2_Ox_Si) / (np.linalg.norm(vector1_Ox_Si) * np.linalg.norm(vector2_Ox_Si))
                    theta_Ox_Si_Ox_1 = np.arccos(np.clip(cos_theta_Ox_Si_Ox_1, -1.0, 1.0))
                    theta_Ox_Si_Ox_degree_1 = np.degrees(theta_Ox_Si_Ox_1)
                    hoek_Ox_Si_Ox.append(float(theta_Ox_Si_Ox_degree_1))
                    cos_theta_Ox_Si_Ox_2 = np.dot(vector1_Ox_Si, vector3_Ox_Si) / (np.linalg.norm(vector1_Ox_Si) * np.linalg.norm(vector3_Ox_Si))
                    theta_Ox_Si_Ox_2 = np.arccos(np.clip(cos_theta_Ox_Si_Ox_2, -1.0, 1.0))
                    theta_Ox_Si_Ox_degree_2 = np.degrees(theta_Ox_Si_Ox_2)
                    hoek_Ox_Si_Ox.append(float(theta_Ox_Si_Ox_degree_2))
                    cos_theta_Ox_Si_Ox_3 = np.dot(vector1_Ox_Si, vector4_Ox_Si) / (np.linalg.norm(vector1_Ox_Si) * np.linalg.norm(vector4_Ox_Si))
                    theta_Ox_Si_Ox_3 = np.arccos(np.clip(cos_theta_Ox_Si_Ox_3, -1.0, 1.0))
                    theta_Ox_Si_Ox_degree_3 = np.degrees(theta_Ox_Si_Ox_3)
                    hoek_Ox_Si_Ox.append(float(theta_Ox_Si_Ox_degree_3))
                    cos_theta_Ox_Si_Ox_4 = np.dot(vector2_Ox_Si, vector3_Ox_Si) / (np.linalg.norm(vector2_Ox_Si) * np.linalg.norm(vector3_Ox_Si))
                    theta_Ox_Si_Ox_4 = np.arccos(np.clip(cos_theta_Ox_Si_Ox_4, -1.0, 1.0))
                    theta_Ox_Si_Ox_degree_4 = np.degrees(theta_Ox_Si_Ox_4)
                    hoek_Ox_Si_Ox.append(float(theta_Ox_Si_Ox_degree_4))
                    cos_theta_Ox_Si_Ox_5 = np.dot(vector2_Ox_Si, vector4_Ox_Si) / (np.linalg.norm(vector2_Ox_Si) * np.linalg.norm(vector4_Ox_Si))
                    theta_Ox_Si_Ox_5 = np.arccos(np.clip(cos_theta_Ox_Si_Ox_5, -1.0, 1.0))
                    theta_Ox_Si_Ox_degree_5 = np.degrees(theta_Ox_Si_Ox_5)
                    hoek_Ox_Si_Ox.append(float(theta_Ox_Si_Ox_degree_5))
                    cos_theta_Ox_Si_Ox_6 = np.dot(vector3_Ox_Si, vector4_Ox_Si) / (np.linalg.norm(vector3_Ox_Si) * np.linalg.norm(vector4_Ox_Si))
                    theta_Ox_Si_Ox_6 = np.arccos(np.clip(cos_theta_Ox_Si_Ox_6, -1.0, 1.0))
                    theta_Ox_Si_Ox_degree_6 = np.degrees(theta_Ox_Si_Ox_6)
                    hoek_Ox_Si_Ox.append(float(theta_Ox_Si_Ox_degree_6))
                    coord_Ox_Si = []
                elif bindingen_Ox_Si == 5:
                    aantal_tetra_Si -= 1
                    del hoek_Ox_Si_Ox[-6:]
                    bindingen_Ox_Si = 0
                else:
                    coord_Ox_Si.append(atoomcontrole_Si)
    elif atoomsoort == 2 and l <= 2144:
        atoomreferentie_Al = [regel_n[2], regel_n[3], regel_n[4]]
        regelnummer_Al = int(regel_n[0])
        for regel_Al in data[regelnummer_Al + k:2144 + k]:
            atoomcontroletype = regel_Al[1]
            atoomcontrole_Al = [regel_Al[2], regel_Al[3], regel_Al[4]]
            positie_atomen_Al = [atoomreferentie_Al, atoomcontrole_Al]
            atoms_Al = Atoms(positions=positie_atomen_Al)
            atoms_Al.set_cell((float(maximum[0]), float(maximum[1]), float(maximum[2])))
            atoms_Al.set_pbc(True)
            afstand_Al = atoms_Al.get_distance(0, 1, mic=True)
            if atoomcontroletype == 5 and afstand_Al <= CO_Al_Ox:
                bindingen_Ox_Al += 1
                if bindingen_Ox_Al == 4: #Berekenen van de hoeken O-Si-O
                    aantal_tetra_Al += 1
                    coord_Ox_Al.append(atoomcontrole_Al)
                    coord_Ox_Al = np.array(coord_Ox_Al)
                    vector1_Ox_Al = coord_Ox_Al[0] - atoomreferentie_Al
                    vector2_Ox_Al = coord_Ox_Al[1] - atoomreferentie_Al
                    vector3_Ox_Al = coord_Ox_Al[2] - atoomreferentie_Al
                    vector4_Ox_Al = coord_Ox_Al[3] - atoomreferentie_Al
                    cos_theta_Ox_Al_Ox_1 = np.dot(vector1_Ox_Al, vector2_Ox_Al) / (np.linalg.norm(vector1_Ox_Al) * np.linalg.norm(vector2_Ox_Al))
                    theta_Ox_Al_Ox_1 = np.arccos(np.clip(cos_theta_Ox_Al_Ox_1, -1.0, 1.0))
                    theta_Ox_Al_Ox_degree_1 = np.degrees(theta_Ox_Al_Ox_1)
                    hoek_Ox_Al_Ox.append(float(theta_Ox_Al_Ox_degree_1))
                    cos_theta_Ox_Al_Ox_2 = np.dot(vector1_Ox_Al, vector3_Ox_Al) / (np.linalg.norm(vector1_Ox_Al) * np.linalg.norm(vector3_Ox_Al))
                    theta_Ox_Al_Ox_2 = np.arccos(np.clip(cos_theta_Ox_Al_Ox_2, -1.0, 1.0))
                    theta_Ox_Al_Ox_degree_2 = np.degrees(theta_Ox_Al_Ox_2)
                    hoek_Ox_Al_Ox.append(float(theta_Ox_Al_Ox_degree_2))
                    cos_theta_Ox_Al_Ox_3 = np.dot(vector1_Ox_Al, vector4_Ox_Al) / (np.linalg.norm(vector1_Ox_Al) * np.linalg.norm(vector4_Ox_Al))
                    theta_Ox_Al_Ox_3 = np.arccos(np.clip(cos_theta_Ox_Al_Ox_3, -1.0, 1.0))
                    theta_Ox_Al_Ox_degree_3 = np.degrees(theta_Ox_Al_Ox_3)
                    hoek_Ox_Al_Ox.append(float(theta_Ox_Al_Ox_degree_3))
                    cos_theta_Ox_Al_Ox_4 = np.dot(vector2_Ox_Al, vector3_Ox_Al) / (np.linalg.norm(vector2_Ox_Al) * np.linalg.norm(vector3_Ox_Al))
                    theta_Ox_Al_Ox_4 = np.arccos(np.clip(cos_theta_Ox_Al_Ox_4, -1.0, 1.0))
                    theta_Ox_Al_Ox_degree_4 = np.degrees(theta_Ox_Al_Ox_4)
                    hoek_Ox_Al_Ox.append(float(theta_Ox_Al_Ox_degree_4))
                    cos_theta_Ox_Al_Ox_5 = np.dot(vector2_Ox_Al, vector4_Ox_Al) / (np.linalg.norm(vector2_Ox_Al) * np.linalg.norm(vector4_Ox_Al))
                    theta_Ox_Al_Ox_5 = np.arccos(np.clip(cos_theta_Ox_Al_Ox_5, -1.0, 1.0))
                    theta_Ox_Al_Ox_degree_5 = np.degrees(theta_Ox_Al_Ox_5)
                    hoek_Ox_Al_Ox.append(float(theta_Ox_Al_Ox_degree_5))
                    cos_theta_Ox_Al_Ox_6 = np.dot(vector3_Ox_Al, vector4_Ox_Al) / (np.linalg.norm(vector3_Ox_Al) * np.linalg.norm(vector4_Ox_Al))
                    theta_Ox_Al_Ox_6 = np.arccos(np.clip(cos_theta_Ox_Al_Ox_6, -1.0, 1.0))
                    theta_Ox_Al_Ox_degree_6 = np.degrees(theta_Ox_Al_Ox_6)
                    hoek_Ox_Al_Ox.append(float(theta_Ox_Al_Ox_degree_6))
                    coord_Ox_Al = []
                elif bindingen_Ox_Al == 5:
                    aantal_tetra_Al -= 1
                    del hoek_Ox_Al_Ox[-6:]
                    bindingen_Ox_Al = 0
                else:
                    coord_Ox_Al.append(atoomcontrole_Al)
    elif l > 2144:
        k += 2144
        l = 0
    bindingen_Ox_Si = 0
    bindingen_Ox_Al = 0
    coord_Ox_Al = []
    coord_Ox_Si = []
    l += 1
hoek_Ox_Si_Ox_gemiddeld = np.mean(hoek_Ox_Si_Ox)
print('De gemiddelde bindingshoek voor de binding Ox_Si_Ox bedraagt ', hoek_Ox_Si_Ox_gemiddeld)
hoek_Ox_Si_Ox_STD = np.std(hoek_Ox_Si_Ox)
print('De standaardafwijking voor de binding Ox-Si-Ox bedraagt ', hoek_Ox_Si_Ox_STD)
#print('De lijst met bindingshoeken voor Ox-Al-Ox bedraagt ', hoek_Ox_Al_Ox)
hoek_Ox_Al_Ox_gemiddeld = np.mean(hoek_Ox_Al_Ox)
print('De gemiddelde bindingshoek voor de binding Ox_Al_Ox bedraagt ', hoek_Ox_Al_Ox_gemiddeld)
hoek_Ox_Al_Ox_STD = np.std(hoek_Ox_Al_Ox)
print('De standaardafwijking voor de binding Ox-Al-Ox bedraagt ', hoek_Ox_Al_Ox_STD)

#Bepalen van de Ox-Si-Ox hoekverdeling
lengte_Ox_Si_Ox = len(hoek_Ox_Si_Ox)
print('De lengte van de lijst met Ox-Si-Ox hoekwaarden is ', lengte_Ox_Si_Ox)
afgerond_Ox_Si_Ox_waarden = [int(number_Si) for number_Si in hoek_Ox_Si_Ox]
afgerond_Ox_Si_Ox_waarden.sort()
waarden_Ox_Si_Ox = []
frequenties_Ox_Si_Ox = []
for i_Si in afgerond_Ox_Si_Ox_waarden:
    if i_Si not in waarden_Ox_Si_Ox:
        waarden_Ox_Si_Ox.append(i_Si)
        aantal_hoek_Si = afgerond_Ox_Si_Ox_waarden.count(i_Si)
        frequenties_Ox_Si_Ox.append(aantal_hoek_Si)
waarden_Ox_Si_Ox = np.array(waarden_Ox_Si_Ox)
frequenties_Ox_Si_Ox = np.array(frequenties_Ox_Si_Ox)
print('De lijst met alle waarden voor de hoek Ox-Si-Ox is ', waarden_Ox_Si_Ox)
print('De lengte van de lijst met Ox-Si-Ox hoekwaarden is ', len(waarden_Ox_Si_Ox))
print('De lijst met het percentueel voorkomen van de Ox-Si-Ox hoek is ', frequenties_Ox_Si_Ox)
print('De lengte van de lijst met Ox-Si-Ox hoekwaarden is ', len(frequenties_Ox_Si_Ox))

with open('X_Ox_Si_Ox.txt', 'w') as fx_Si:
    for waardex_Si in waarden_Ox_Si_Ox:
        fx_Si.write(str(waardex_Si) + '\n')
        
with open('Y_Ox_Si_Ox.txt', 'w') as fy_Si:
    for waardey_Si in frequenties_Ox_Si_Ox:
        fy_Si.write(str(waardey_Si) + '\n')
        
#Bepalen van de Ox-Al-Ox hoekverdeling
lengte_Ox_Al_Ox = len(hoek_Ox_Al_Ox)
print('De lengte van de lijst met Ox-Al-Ox hoekwaarden is ', lengte_Ox_Al_Ox)
afgerond_Ox_Al_Ox_waarden = [int(number_Al) for number_Al in hoek_Ox_Al_Ox]
afgerond_Ox_Al_Ox_waarden.sort()
waarden_Ox_Al_Ox = []
frequenties_Ox_Al_Ox = []
for i_Al in afgerond_Ox_Al_Ox_waarden:
    if i_Al not in waarden_Ox_Al_Ox:
        waarden_Ox_Al_Ox.append(i_Al)
        aantal_hoek_Al = afgerond_Ox_Al_Ox_waarden.count(i_Al)
        frequenties_Ox_Al_Ox.append(aantal_hoek_Al)
waarden_Ox_Al_Ox = np.array(waarden_Ox_Al_Ox)
frequenties_Ox_Al_Ox = np.array(frequenties_Ox_Al_Ox)
print('De lijst met alle waarden voor de hoek Ox-Al-Ox is ', waarden_Ox_Al_Ox)
print('De lengte van de lijst met Ox-Al-Ox hoekwaarden is ', len(waarden_Ox_Al_Ox))
print('De lijst met het percentueel voorkomen van de Ox-Al-Ox hoek is ', frequenties_Ox_Al_Ox)
print('De lengte van de lijst met Ox-Al-Ox hoekwaarden is ', len(frequenties_Ox_Al_Ox))

with open('X_Ox_Al_Ox.txt', 'w') as fx_Al:
    for waardex_Al in waarden_Ox_Al_Ox:
        fx_Al.write(str(waardex_Al) + '\n')
        
with open('Y_Ox_Al_Ox.txt', 'w') as fy_Al:
    for waardey_Al in frequenties_Ox_Al_Ox:
        fy_Al.write(str(waardey_Al) + '\n')        

#Verwijderen van de aangemaakte files
os.remove('boundary_conditions.txt')
os.remove('NVT_op_bewerkt.trj')
os.remove('data_file_atomen.txt')
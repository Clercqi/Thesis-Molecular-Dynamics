#Script 6(a): Classificeren van de soorten zuurstof(bridging, nonbridging, free) en bepalen van de bond angle distributie van Si-O-Si, O-Si-O, Al-O-Si, O-Al-O, Al-O-Al
#Files: data.file en NVT_op.trj
import numpy as np
import os
from ase import Atoms
import matplotlib.pyplot as plt
#from scipy.stats import norm
#from scipy.interpolate import make_interp_spline

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
    
#####Bepalen van het aantal ionen per type -> File data.file#######
with open('data.file', 'r') as bestand_ionen:
    lines_ionen = bestand_ionen.readlines()
    
#Verwijder de eerste 19 lijnen
del lines_ionen[:19]

#Schrijf uit naar een nieuwe file data_file_atomen.txt
with open('data_file_atomen.txt', 'w') as f_ionen:
    f_ionen.writelines(lines_ionen)
    
#Inlezen van de data_file_atomen.txt en bepalen van het aantal atomen per type ion
data_ionen = np.loadtxt('data_file_atomen.txt', dtype=float)

aantal_Si = 0
aantal_Al = 0
aantal_Ca = 0
aantal_Mg = 0
aantal_Ox = 0

for regel_ionen in data_ionen:
    if regel_ionen[1] == 1:
        aantal_Si += 1
    elif regel_ionen[1] == 2:
        aantal_Al += 1
    elif regel_ionen[1] == 3:
        aantal_Ca += 1
    elif regel_ionen[1] == 4:
        aantal_Mg += 1
    else:
        aantal_Ox += 1
print('Het aantal silicium ionen is ', aantal_Si)
print('Het aantal aluminium ionen is ', aantal_Al)
print('Het aantal calcium ionen is ', aantal_Ca)
print('Het aantal magnesium ionen is ', aantal_Mg)
print('Het aantal zuurstof ionen is ', aantal_Ox)

#Inlezen van de file NVT_op_bewerkt.trj
data = np.loadtxt('NVT_op_bewerkt.trj', dtype=float, usecols=[0, 1, 5, 6, 7]) #Kolommen in de volgorde [atoom, ID, x, y, z]

#Definieren van de voorwaarden - bond lengths
CO_Si_Ox = 2.28 #De waarden zijn gehaald uit de berekende bond lengths van de RDFs
CO_Al_Ox = 2.5
timestep = 1000
total_time = 100000
frequentie = total_time / timestep

#Bepalen van de soorten zuurstof en de bindingshoeken voor de bindingen Si-O-Si, Al-O-Al en Si-O-Al
nonbridging_Ox = 0
bridging_Ox = 0
fractie_Q0 = 0
fractie_Q1 = 0
fractie_Q2 = 0
fractie_Q3 = 0
fractie_Q4 = 0
fractie_Al_Q0 = 0
fractie_Al_Q1 = 0
fractie_Al_Q2 = 0
fractie_Al_Q3 = 0
fractie_Al_Q4 = 0
BO = 0
i = 0
j = 1
k = 1
l = 0
a = 1
b = 0
bindingen_Si = 0
bindingen_Al = 0
bindingen_single = 0
Lijst_BO_Si = []
Lijst_BO_Al = []
positie_Si = []
positie_Al = []
hoek_Si_Ox_Si = []
hoek_Al_Ox_Al = []
hoek_Si_Ox_Al = []
for rij_n in data:
    atoom_type = rij_n[1]
    if atoom_type == 5 and j <= 2144: #De atomen Ox beschouwen en als referentie gebruiken
        atoom_ref = [rij_n[2], rij_n[3], rij_n[4]]
        rijnummer = int(rij_n[0])
        for rij_nn in data[i:rijnummer + i]: #Itereren over de atomen met als referentie een zuurstofatoom
            atoom_controle_type = rij_nn[1]
            atoom_controle = [rij_nn[2], rij_nn[3], rij_nn[4]]
            positions_atoms = [atoom_ref, atoom_controle] #Een lijst van de posities van de atomen
            atoms = Atoms(positions=positions_atoms) #Maak een Atoms object
            atoms.set_cell((float(maximum[0]), float(maximum[1]), float(maximum[2]))) #Definieren van de periodieke boundary conditions
            atoms.set_pbc(True) #Implementeren van de periodieke boundary conditions
            afstand = atoms.get_distance(0, 1, mic=True) #Bereken de minimale afstand tussen de atomen met behulp van de minimale image conventie (houdt rekening met periodieke grnesvoorwaarden)
            if atoom_controle_type == 1 and afstand <= CO_Si_Ox:
                bindingen_Si += 1
                bindingen_single += 1
                if bindingen_Si == 2: #Berekenen van de hoek tussen Si-O-Si
                    bridging_Ox += 1
                    positie_Si.append(atoom_controle)
                    positie_Si = np.array(positie_Si)
                    vector1_Si = positie_Si[0] - atoom_ref
                    vector2_Si = positie_Si[1] - atoom_ref
                    cos_theta_Si_Ox_Si = np.dot(vector1_Si, vector2_Si) / (np.linalg.norm(vector1_Si) * np.linalg.norm(vector2_Si))
                    theta_Si_Ox_Si = np.arccos(np.clip(cos_theta_Si_Ox_Si, -1.0, 1.0))
                    theta_Si_Ox_Si_degree = np.degrees(theta_Si_Ox_Si)
                    hoek_Si_Ox_Si.append(float(theta_Si_Ox_Si_degree))
                    Lijst_BO_Si.append(rijnummer) #index van de BO toevoegen aan een lijst
                    positie_Si = []
                    bindingen_Si = 0
                    nonbridging_Ox -= 1
                elif bindingen_Si == 1 and bindingen_Al == 1: #Berekenen van de hoek tussen Si-O-Al
                    bridging_Ox += 1
                    positie_Si.append(atoom_controle)
                    positie_Si = np.array(positie_Si)
                    positie_Al = np.array(positie_Al)
                    vector1 = positie_Si - atoom_ref
                    vector2 = positie_Al - atoom_ref
                    cos_theta_Si_Ox_Al = np.dot(vector1.flatten(), vector2.flatten()) / (np.linalg.norm(vector1) * np.linalg.norm(vector2))
                    theta_Si_Ox_Al = np.arccos(np.clip(cos_theta_Si_Ox_Al, -1.0, 1.0))
                    theta_Si_Ox_Al_degree = np.degrees(theta_Si_Ox_Al)
                    hoek_Si_Ox_Al.append(float(theta_Si_Ox_Al_degree))
                    positie_Si = []
                    positie_Al = []
                    bindingen_Si = 0
                    bindingen_Al = 0
                    nonbridging_Ox -= 1
                else:
                    positie_Si.append(atoom_controle)
                    nonbridging_Ox += 1
            elif atoom_controle_type == 2 and afstand <= CO_Al_Ox:
                bindingen_Al += 1
                bindingen_single += 1
                if bindingen_Al == 2: #Berekenen van de hoek tussen Al-O-Al
                    bridging_Ox += 1
                    positie_Al.append(atoom_controle)
                    positie_Al = np.array(positie_Al)
                    vector1_Al = positie_Al[0] - atoom_ref
                    vector2_Al = positie_Al[1] - atoom_ref
                    cos_theta_Al_Ox_Al = np.dot(vector1_Al, vector2_Al) / (np.linalg.norm(vector1_Al) * np.linalg.norm(vector2_Al))
                    theta_Al_Ox_Al = np.arccos(np.clip(cos_theta_Al_Ox_Al, -1.0, 1.0))
                    theta_Al_Ox_Al_degree = np.degrees(theta_Al_Ox_Al)
                    hoek_Al_Ox_Al.append(float(theta_Al_Ox_Al_degree))
                    Lijst_BO_Al.append(rijnummer)
                    positie_Al = []
                    bindingen_Al = 0
                    nonbridging_Ox -= 1
                elif bindingen_Al == 1 and bindingen_Si == 1: #Berekenen van de hoek tussen Si-O-Al
                    bridging_Ox += 1
                    positie_Al.append(atoom_controle)
                    positie_Si = np.array(positie_Si)
                    positie_Al = np.array(positie_Al)
                    vector1 = positie_Si - atoom_ref
                    vector2 = positie_Al - atoom_ref
                    cos_theta_Si_Ox_Al = np.dot(vector1.flatten(), vector2.flatten()) / (np.linalg.norm(vector1) * np.linalg.norm(vector2))
                    theta_Si_Ox_Al = np.arccos(np.clip(cos_theta_Si_Ox_Al, -1.0, 1.0))
                    theta_Si_Ox_Al_degree = np.degrees(theta_Si_Ox_Al)
                    hoek_Si_Ox_Al.append(float(theta_Si_Ox_Al_degree))
                    positie_Si = []
                    positie_Al = []
                    bindingen_Al = 0
                    bindingen_Si = 0
                    nonbridging_Ox -= 1
                else: 
                    positie_Al.append(atoom_controle)
                    nonbridging_Ox += 1
    if k == 2144: #Als de eerste frame doorlopen is dan moeten de Si als referentie atomen beschouwd worden en itereren over de eerste frame
        Lijst_BO_Si = np.array(Lijst_BO_Si)
        for ri in data[l:l + 2144]: #Rekening houden met alle snapshots
            atoom_type_Si = ri[1]
            if atoom_type_Si == 1:
                atoom_ref_Si = [ri[2], ri[3], ri[4]]
                rijnummer_2 = int(ri[0])
                for ri_nn in data[l + rijnummer_2:l + 2144]: 
                    atoom_controle_type2 = ri_nn[1]
                    atoom_controle_2 = [ri_nn[2], ri_nn[3], ri_nn[4]]
                    positions_atoms_2 = [atoom_ref_Si, atoom_controle_2] 
                    atoms_2 = Atoms(positions=positions_atoms_2) 
                    atoms_2.set_cell((float(maximum[0]), float(maximum[1]), float(maximum[2])))
                    atoms_2.set_pbc(True) 
                    afstand_2 = atoms_2.get_distance(0, 1, mic=True)
                    if atoom_controle_type2 == 5 and afstand_2 <= CO_Si_Ox: #nagaan of het atoom een zuurstof is en een binding heeft met Si
                        atoom_index = int(ri_nn[0])
                        if atoom_index in Lijst_BO_Si: #Nagaan of de zuurstof type bridging is
                            BO += 1
                if BO == 0:
                    fractie_Q0 += 1
                elif BO == 1:
                    fractie_Q1 += 1
                    BO = 0
                elif BO == 2:
                    fractie_Q2 += 1
                    BO = 0
                elif BO == 3:
                    fractie_Q3 += 1
                    BO = 0
                else: 
                    fractie_Q4 += 1
                    BO = 0
        Lijst_BO_Si = []
        l += 2144
        k = 0
    if a == 2144: #Als de eerste frame doorlopen is dan moeten de Al als referentie atomen beschouwd worden en itereren over de eerste frame
        Lijst_BO_Al = np.array(Lijst_BO_Al)
        for ri_Al in data[b:b + 2144]: #Rekening houden met alle snapshots
            atoom_type_Al = ri_Al[1]
            if atoom_type_Al == 2:
                atoom_ref_Al = [ri_Al[2], ri_Al[3], ri_Al[4]]
                rijnummer_3 = int(ri_Al[0])
                for ri_nn_Al in data[b + rijnummer_3:b + 2144]: 
                    atoom_controle_type3 = ri_nn_Al[1]
                    atoom_controle_3 = [ri_nn_Al[2], ri_nn_Al[3], ri_nn_Al[4]]
                    positions_atoms_3 = [atoom_ref_Al, atoom_controle_3] 
                    atoms_3 = Atoms(positions=positions_atoms_3) 
                    atoms_3.set_cell((float(maximum[0]), float(maximum[1]), float(maximum[2])))
                    atoms_3.set_pbc(True) 
                    afstand_3 = atoms_3.get_distance(0, 1, mic=True)
                    if atoom_controle_type3 == 5 and afstand_3 <= CO_Al_Ox: #nagaan of het atoom een zuurstof is en een binding heeft met Al
                        atoom_index_Al = int(ri_nn_Al[0])
                        if atoom_index_Al in Lijst_BO_Al: #Nagaan of de zuurstof type bridging is
                            BO += 1
                if BO == 0:
                    fractie_Al_Q0 += 1
                elif BO == 1:
                    fractie_Al_Q1 += 1
                    BO = 0
                elif BO == 2:
                    fractie_Al_Q2 += 1
                    BO = 0
                elif BO == 3:
                    fractie_Al_Q3 += 1
                    BO = 0
                else: 
                    fractie_Al_Q4 += 1
                    BO = 0
        Lijst_BO_Al = []
        b += 2144
        a = 0
    elif j > 2144: #Terug bij het eerste atoom van de volgende frequentie  
        i += 2144  #In de for-lus starten bij het begin van de volgende frequentie
        j = 0
    bindingen_single = 0
    bindingen_Si = 0
    bindingen_Al = 0
    positie_Si = []
    positie_Al = []
    j += 1
    k += 1
    a += 1
            
aantal_fractie_Q0 = (fractie_Q0 * 100) / (frequentie * aantal_Si)
print('Het gemiddeld aantal Q0 Si fractie is ', aantal_fractie_Q0)            
aantal_fractie_Q1 = (fractie_Q1 * 100) / (frequentie * aantal_Si)
print('Het gemiddeld aantal Q1 Si fractie is ', aantal_fractie_Q1) 
aantal_fractie_Q2 = (fractie_Q2 * 100) / (frequentie * aantal_Si)
print('Het gemiddeld aantal Q2 Si fractie is ', aantal_fractie_Q2)         
aantal_fractie_Q3 = (fractie_Q3 * 100) / (frequentie * aantal_Si)
print('Het gemiddeld aantal Q3 Si fractie is ', aantal_fractie_Q3) 
aantal_fractie_Q4 = (fractie_Q4 * 100) / (frequentie * aantal_Si)
print('Het gemiddeld aantal Q4 Si fractie is ', aantal_fractie_Q4) 
controle = aantal_fractie_Q0 + aantal_fractie_Q1 + aantal_fractie_Q2 + aantal_fractie_Q3 + aantal_fractie_Q4
print('Het totaal aantal Si atomen is ', aantal_Si)
print('Het percentage aan Q fracties is ', controle)

aantal_fractie_Al_Q0 = (fractie_Al_Q0 * 100) / (frequentie * aantal_Al)
print('Het gemiddeld aantal Q0 Al fractie is ', aantal_fractie_Al_Q0)            
aantal_fractie_Al_Q1 = (fractie_Al_Q1 * 100) / (frequentie * aantal_Al)
print('Het gemiddeld aantal Q1 Al fractie is ', aantal_fractie_Al_Q1) 
aantal_fractie_Al_Q2 = (fractie_Al_Q2 * 100) / (frequentie * aantal_Al)
print('Het gemiddeld aantal Q2 Al fractie is ', aantal_fractie_Al_Q2)         
aantal_fractie_Al_Q3 = (fractie_Al_Q3 * 100) / (frequentie * aantal_Al)
print('Het gemiddeld aantal Q3 Al fractie is ', aantal_fractie_Al_Q3) 
aantal_fractie_Al_Q4 = (fractie_Al_Q4 * 100) / (frequentie * aantal_Al)
print('Het gemiddeld aantal Q4 Al fractie is ', aantal_fractie_Al_Q4) 
controle_Al = aantal_fractie_Al_Q0 + aantal_fractie_Al_Q1 + aantal_fractie_Al_Q2 + aantal_fractie_Al_Q3 + aantal_fractie_Al_Q4
print('Het totaal aantal Al atomen is ', aantal_Al)
print('Het percentage aan Q fracties is ', controle_Al)
    
aantal_nonbridging_Ox = nonbridging_Ox / frequentie
print('Het gemiddeld aantal NBO is ', aantal_nonbridging_Ox)
aantal_bridging_Ox = bridging_Ox / frequentie
print('Het gemiddeld aantal BO is ', aantal_bridging_Ox)
aantal_free_Ox = (aantal_Ox * frequentie - nonbridging_Ox - bridging_Ox) / frequentie
print('Het gemiddeld aantal free Ox is ', aantal_free_Ox)
print('De lijst met bindingshoeken voor Si-Ox-Si bedraagt ', hoek_Si_Ox_Si)
hoek_Si_Ox_Si_gemiddeld = np.mean(hoek_Si_Ox_Si)
print('De gemiddelde bindingshoek voor de binding Si_Ox_Si bedraagt ', hoek_Si_Ox_Si_gemiddeld)
hoek_Si_Ox_Si_STD = np.std(hoek_Si_Ox_Si)
print('De standaardafwijking voor de binding Si_Ox_Si bedraagt ', hoek_Si_Ox_Si_STD)
#print('De lijst met bindingshoeken voor Al-Ox-Al bedraagt ', hoek_Al_Ox_Al)
hoek_Al_Ox_Al_gemiddeld = np.mean(hoek_Al_Ox_Al)
print('De gemiddelde bindingshoek voor de binding Al_Ox_Al bedraagt ', hoek_Al_Ox_Al_gemiddeld)
hoek_Al_Ox_Al_STD = np.std(hoek_Al_Ox_Al)
print('De standaardafwijking voor de binding Al_Ox_Al bedraagt ', hoek_Al_Ox_Al_STD)
#print('De lijst met bindingshoeken voor Si-Ox-Al bedraagt ', hoek_Si_Ox_Al) 
hoek_Si_Ox_Al_gemiddeld = np.mean(hoek_Si_Ox_Al)
print('De gemiddelde bindingshoek voor de binding Si_Ox_Al bedraagt ', hoek_Si_Ox_Al_gemiddeld)
hoek_Si_Ox_Al_STD = np.std(hoek_Si_Ox_Al)
print('De standaardafwijking voor de binding Si_Ox_Al bedraagt ', hoek_Si_Ox_Al_STD)

#Bepalen van de Si-Ox-Si hoekverdeling
lengte_Si_Ox_Si = len(hoek_Si_Ox_Si)
print('De lengte van de lijst met Si-Ox_Si hoekwaarden is ', lengte_Si_Ox_Si)
afgerond_Si_Ox_Si_waarden = [int(number) for number in hoek_Si_Ox_Si]
afgerond_Si_Ox_Si_waarden.sort()
waarden_Si_Ox_Si = []
frequenties_Si_Ox_Si = []
for i in afgerond_Si_Ox_Si_waarden:
    if i not in waarden_Si_Ox_Si:
        waarden_Si_Ox_Si.append(i)
        aantal_hoek = afgerond_Si_Ox_Si_waarden.count(i)
        frequenties_Si_Ox_Si.append(aantal_hoek)
waarden_Si_Ox_Si = np.array(waarden_Si_Ox_Si)
frequenties_Si_Ox_Si = np.array(frequenties_Si_Ox_Si)
print('De lijst met alle waarden voor de hoek Si-Ox-Si is ', waarden_Si_Ox_Si)
print('De lengte van de lijst met Si-Ox_Si hoekwaarden is ', len(waarden_Si_Ox_Si))
print('De lijst met het percentueel voorkomen van de Si-Ox-Si hoek is ', frequenties_Si_Ox_Si)
print('De lengte van de lijst met Si-Ox_Si hoekwaarden is ', len(frequenties_Si_Ox_Si))

with open('X_Si_Ox_Si.txt', 'w') as fx:
    for waardex in waarden_Si_Ox_Si:
        fx.write(str(waardex) + '\n')
        
with open('Y_Si_Ox_Si.txt', 'w') as fy:
    for waardey in frequenties_Si_Ox_Si:
        fy.write(str(waardey) + '\n') 

#Bepalen van de Al-Ox-Al hoekverdeling
lengte_Al_Ox_Al = len(hoek_Al_Ox_Al)
print('De lengte van de lijst met Al-Ox-Al hoekwaarden is ', lengte_Al_Ox_Al)
afgerond_Al_Ox_Al_waarden = [int(number_Al) for number_Al in hoek_Al_Ox_Al]
afgerond_Al_Ox_Al_waarden.sort()
waarden_Al_Ox_Al = []
frequenties_Al_Ox_Al = []
for i_Al in afgerond_Al_Ox_Al_waarden:
    if i_Al not in waarden_Al_Ox_Al:
        waarden_Al_Ox_Al.append(i_Al)
        aantal_hoek_Al = afgerond_Al_Ox_Al_waarden.count(i_Al)
        frequenties_Al_Ox_Al.append(aantal_hoek_Al)
waarden_Al_Ox_Al = np.array(waarden_Al_Ox_Al)
frequenties_Al_Ox_Al = np.array(frequenties_Al_Ox_Al)
print('De lijst met alle waarden voor de hoek Al-Ox-Al is ', waarden_Al_Ox_Al)
print('De lengte van de lijst met Al-Ox-Al hoekwaarden is ', len(waarden_Al_Ox_Al))
print('De lijst met het percentueel voorkomen van de Al-Ox-Al hoek is ', frequenties_Al_Ox_Al)
print('De lengte van de lijst met Al-Ox-Al hoekwaarden is ', len(frequenties_Al_Ox_Al))

with open('X_Al_Ox_Al.txt', 'w') as fx_Al:
    for waardex_Al in waarden_Al_Ox_Al:
        fx_Al.write(str(waardex_Al) + '\n')
        
with open('Y_Al_Ox_Al.txt', 'w') as fy_Al:
    for waardey_Al in frequenties_Al_Ox_Al:
        fy_Al.write(str(waardey_Al) + '\n')  
        
#Bepalen van de Si-Ox-Al hoekverdeling
lengte_Si_Ox_Al = len(hoek_Si_Ox_Al)
print('De lengte van de lijst met Si-Ox-Al hoekwaarden is ', lengte_Si_Ox_Al)
afgerond_Si_Ox_Al_waarden = [int(number_Si) for number_Si in hoek_Si_Ox_Al]
afgerond_Si_Ox_Al_waarden.sort()
waarden_Si_Ox_Al = []
frequenties_Si_Ox_Al = []
for i_Si in afgerond_Si_Ox_Al_waarden:
    if i_Si not in waarden_Si_Ox_Al:
        waarden_Si_Ox_Al.append(i_Si)
        aantal_hoek_Si = afgerond_Si_Ox_Al_waarden.count(i_Si)
        frequenties_Si_Ox_Al.append(aantal_hoek_Si)
waarden_Si_Ox_Al = np.array(waarden_Si_Ox_Al)
frequenties_Si_Ox_Al = np.array(frequenties_Si_Ox_Al)
print('De lijst met alle waarden voor de hoek Si-Ox-Al is ', waarden_Si_Ox_Al)
print('De lengte van de lijst met Si-Ox-Al hoekwaarden is ', len(waarden_Si_Ox_Al))
print('De lijst met het percentueel voorkomen van de Si-Ox-Al hoek is ', frequenties_Si_Ox_Al)
print('De lengte van de lijst met Si-Ox-Al hoekwaarden is ', len(frequenties_Si_Ox_Al))

with open('X_Si_Ox_Al.txt', 'w') as fx_Si:
    for waardex_Si in waarden_Si_Ox_Al:
        fx_Si.write(str(waardex_Si) + '\n')
        
with open('Y_Si_Ox_Al.txt', 'w') as fy_Si:
    for waardey_Si in frequenties_Si_Ox_Al:
        fy_Si.write(str(waardey_Si) + '\n')

#Verwijderen van de aangemaakte files
os.remove('boundary_conditions.txt')
os.remove('NVT_op_bewerkt.trj')
os.remove('data_file_atomen.txt')
os.remove('data_file_MM.txt')
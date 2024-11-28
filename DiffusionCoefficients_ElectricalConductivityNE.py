#Script 4: Bepalen van de diffusiecoefficienten van Al, Ca, Mg, Si en Ox - Berekenen van de elektrische conductiviteit via Nersnt-Einstein (benaderende vergelijking)
#Files: Al-Al-AlMSD.dat, Ca-Ca-CaMSD.dat, Mg-Mg-MgMSD.dat, Ox-Ox-OxMSD.dat, Si-Si-SiMSD.DAT 
import numpy as np
import os
import matplotlib
matplotlib.use('Qt5Agg') #Deze backend biedt GUI-ondersteuning waardoor de figuren worden weergegeven met plt.show()
import matplotlib.pyplot as plt

#Verzamelen van de variabelen voor de berekening van de elektrische conductiviteit via Nernst-Einstein
#Bepalen van het aantal ionen (= aantal atomen) binnen de eenheidscel en het volume van de eenheidscel - inlezen van de file: data.file

#####Bepalen van het aantal ionen per type -> File data.file#######
with open('data.file', 'r') as bestand_ionen:
    lines_ionen = bestand_ionen.readlines()
    
#Verwijder de eerste 19 lijnen
del lines_ionen[:19]

#Schrijf uit naar een nieuwe file data_file_atomen.txt
with open('data_file_atomen.txt', 'w') as f_ionen:
    f_ionen.writelines(lines_ionen)
    
#Inlezen van de data_file_atomen.txt en bepalen van het aantal atomen per type ion
data_ionen = np.loadtxt('data_file_atomen.txt', dtype=int)

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

#####Bepalen van het volume van de eenheidscel - dit volume is de laatste waarde van de coolingdown step/toegevoegde stap -> File log.lammps######
with open('log.lammps', 'r') as bestand_volume:
    lines_volume = bestand_volume.readlines()
    
#Selectie van de nodige lijnen en uitschrijven naar de file log_lammps_volume.txt
start_index = 35533
end_index = 38033
geselecteerde_lijnen = lines_volume[start_index:end_index]

with open('log_lammps_volume.txt', 'w') as f_volume:
    for regel_volume in geselecteerde_lijnen:
        regel_volume = regel_volume.split()
        f_volume.write(' '.join(regel_volume) + '\n')

#Inlezen van de nieuwe file log_lammps_volume.txt
data_volume = np.loadtxt('log_lammps_volume.txt', dtype=float)

#Bepalen van het volume van de eenheidscel na equilibratie op operationele temperatuur
lijst_volume = [] 
timesteps = []   
for regel_vol in data_volume:
    lijst_volume.append(regel_vol[4])
    timesteps.append(float(regel_vol[0]))
    
coefficients_volume = np.polyfit(timesteps, lijst_volume, 1)
slope_volume = coefficients_volume[0]
intercept_volume = coefficients_volume[1]
volume = np.mean(lijst_volume) * (10 ** (-30)) #Eenheid van A³ naar m³
print('Het volume na equilibratie en bij de start van NVT bedraagt ', volume) #Dit is uitgedrukt in m³
    
#####Bepalen van de ladingen per ion type - FF settings -> File inHPC_final.lammps######
with open('inHPC_final.lammps', 'r') as bestand_ladingen:
    lines_ladingen = bestand_ladingen.readlines()
    
start_line = 51
end_line = 56
selected_lines = lines_ladingen[start_line:end_line]

with open('inHPC_final_ladingen.txt', 'w') as f_ladingen:
    for regel_ladingen in selected_lines:
        regel_ladingen = regel_ladingen.split()
        f_ladingen.write(' '.join(regel_ladingen) + '\n')
        
#Inlezen van de nieuwe file inHPC_final_ladingen.txt
data_ladingen = np.loadtxt('inHPC_final_ladingen.txt', dtype=str) 
        
lijst_ladingen = [] #Deze staan in volgorde [Si, Al, Ca, Mg, Ox]
for charges in data_ladingen:
    lijst_ladingen.append(float(charges[4]))
print('De ladingen van de ionen zijn ', lijst_ladingen)
lading_Si = float(lijst_ladingen[0])
lading_Al = float(lijst_ladingen[1])
lading_Ca = float(lijst_ladingen[2])
lading_Mg = float(lijst_ladingen[3])
lading_Ox = float(lijst_ladingen[4])

#####Bepalen van de diffusiecoefficienten van alle elementen########
#Inlezen van de Al-Al-AlMSD.dat file
data_Al = np.loadtxt('Al-Al-AlMSD.dat', dtype=float, skiprows=2)

#Analyseer de gegevens
timesteps_Al = data_Al[:, 0]
MSD_Al = data_Al[:, 4]

#Visualisatie van de MSD van Al over de tijd
plt.plot(timesteps_Al, MSD_Al, color='red', linestyle='solid')
plt.xlabel('time (fs)', fontsize=18)
plt.ylabel('MSD (Å²)', fontsize=18)
plt.title('Mean square displacement of Al', fontsize=18)
plt.xticks(fontsize=18)
plt.yticks(fontsize=18)
plt.show()

#Selectie van het middelste bereik (5 - 50 %) waar de diffusiecoefficient een lineare curve voorstelt
bereik_timesteps_Al = []
bereik_MSD_Al = []
for i_Al, element_Al in enumerate(timesteps_Al):
    if element_Al >= 250000 and element_Al <= 2500000:
        bereik_timesteps_Al.append(element_Al) 
        bereik_MSD_Al.append(MSD_Al[i_Al])
        
#Visualisatie van de beperkte MSD van Al over de tijd
plt.plot(bereik_timesteps_Al, bereik_MSD_Al, color='red', linestyle='solid')
plt.xlabel('time (fs)', fontsize=18)
plt.ylabel('MSD (A²)', fontsize=18)
plt.title('Mean square displacement of Al in selected area', fontsize=18)
plt.xticks(fontsize=18)
plt.yticks(fontsize=18)
plt.show()

#Lineaire regressie (rechte fitten) door het geselecteerde bereik = het nemen van de afgeleide naar de tijd
coefficients_Al = np.polyfit(bereik_timesteps_Al, bereik_MSD_Al, 1)
slope_Al = coefficients_Al[0]
intercept_Al = coefficients_Al[1]
print('De helling is ', slope_Al)
print('De intercept is ', intercept_Al)

#De MSD en afgeleide samen plotten
plt.plot(bereik_timesteps_Al, bereik_MSD_Al, color='red', linestyle='solid', label='MSD_Al')
regressie_lijn_Al = slope_Al * np.array(bereik_timesteps_Al) + intercept_Al
plt.plot(bereik_timesteps_Al, regressie_lijn_Al, color ='black', linestyle='solid', label='Linear regression')
plt.xlabel('time (fs)', fontsize=18)
plt.ylabel('MSD (A²)', fontsize=18)
plt.title('Mean square displacement of Al and linear regression in selected area', fontsize=18)
plt.xticks(fontsize=18)
plt.yticks(fontsize=18)
plt.legend(fontsize=18)
plt.show()

#Berekenen van de diffusiecoefficient van Al
diffusiecoefficient_Al = (slope_Al / 6) * (10 ** (-5)) #Eenheid in m²/s²
print('De diffusiecoefficient van Al bedraagt ', diffusiecoefficient_Al)

#Inlezen van de Ca-Ca-CaMSD.dat file
data_Ca = np.loadtxt('Ca-Ca-CaMSD.dat', dtype=float, skiprows=2)

#Analyseer de gegevens
timesteps_Ca = data_Ca[:, 0]
MSD_Ca = data_Ca[:, 4]

#Visualisatie van de MSD van Ca over de tijd
plt.plot(timesteps_Ca, MSD_Ca, color='blue', linestyle='solid')
plt.xlabel('time (fs)', fontsize=18)
plt.ylabel('MSD (Å²)', fontsize=18)
plt.title('Mean square displacement of Ca', fontsize=18)
plt.xticks(fontsize=18)
plt.yticks(fontsize=18)
plt.show()

#Selectie van het middelste bereik (5 - 50 %) waar de diffusiecoefficient een lineare curve voorstelt
bereik_timesteps_Ca = []
bereik_MSD_Ca = []
for i_Ca, element_Ca in enumerate(timesteps_Ca):
    if element_Ca >= 250000 and element_Ca <= 2500000:
        bereik_timesteps_Ca.append(element_Ca) 
        bereik_MSD_Ca.append(MSD_Ca[i_Ca])
        
#Visualisatie van de beperkte MSD van Ca over de tijd
plt.plot(bereik_timesteps_Ca, bereik_MSD_Ca, color='blue', linestyle='solid')
plt.xlabel('time (fs)', fontsize=18)
plt.ylabel('MSD (A²)', fontsize=18)
plt.title('Mean square displacement of Ca in selected area', fontsize=18)
plt.xticks(fontsize=18)
plt.yticks(fontsize=18)
plt.show()

#Lineaire regressie (rechte fitten) door het geselecteerde bereik = het nemen van de afgeleide naar de tijd
coefficients_Ca = np.polyfit(bereik_timesteps_Ca, bereik_MSD_Ca, 1)
slope_Ca = coefficients_Ca[0]
intercept_Ca = coefficients_Ca[1]
print('De helling is ', slope_Ca)
print('De intercept is ', intercept_Ca)

#De MSD en afgeleide samen plotten
plt.plot(bereik_timesteps_Ca, bereik_MSD_Ca, color='blue', linestyle='solid', label='MSD_Ca')
regressie_lijn_Ca = slope_Ca * np.array(bereik_timesteps_Ca) + intercept_Ca
plt.plot(bereik_timesteps_Ca, regressie_lijn_Ca, color ='black', linestyle='solid', label='Linear regression')
plt.xlabel('time (fs)', fontsize=18)
plt.ylabel('MSD (A²)', fontsize=18)
plt.title('Mean square displacement of Ca and linear regression in selected area', fontsize=18)
plt.xticks(fontsize=18)
plt.yticks(fontsize=18)
plt.legend(fontsize=18)
plt.show()

#Berekenen van de diffusiecoefficient van Ca
diffusiecoefficient_Ca = (slope_Ca / 6) * (10 ** (-5)) #Eenheid in m²/s²
print('De diffusiecoefficient van Ca bedraagt ', diffusiecoefficient_Ca)

#Inlezen van de Mg-Mg-MgMSD.dat file
data_Mg = np.loadtxt('Mg-Mg-MgMSD.dat', dtype=float, skiprows=2)

#Analyseer de gegevens
timesteps_Mg = data_Mg[:, 0]
MSD_Mg = data_Mg[:, 4]

#Visualisatie van de MSD van Mg over de tijd
plt.plot(timesteps_Mg, MSD_Mg, color='green', linestyle='solid')
plt.xlabel('time (fs)', fontsize=18)
plt.ylabel('MSD (Å²)', fontsize=18)
plt.title('Mean square displacement of Mg', fontsize=18)
plt.xticks(fontsize=18)
plt.yticks(fontsize=18)
plt.show()

#Selectie van het middelste bereik (5 - 50 %) waar de diffusiecoefficient een lineare curve voorstelt
bereik_timesteps_Mg = []
bereik_MSD_Mg = []
for i_Mg, element_Mg in enumerate(timesteps_Mg):
    if element_Mg >= 250000 and element_Mg <= 2500000:
        bereik_timesteps_Mg.append(element_Mg) 
        bereik_MSD_Mg.append(MSD_Mg[i_Mg])
        
#Visualisatie van de beperkte MSD van Mg over de tijd
plt.plot(bereik_timesteps_Mg, bereik_MSD_Mg, color='green', linestyle='solid')
plt.xlabel('time (fs)', fontsize=18)
plt.ylabel('MSD (A²)', fontsize=18)
plt.title('Mean square displacement of Mg in selected area', fontsize=18)
plt.xticks(fontsize=18)
plt.yticks(fontsize=18)
plt.show()

#Lineaire regressie (rechte fitten) door het geselecteerde bereik = het nemen van de afgeleide naar de tijd
coefficients_Mg = np.polyfit(bereik_timesteps_Mg, bereik_MSD_Mg, 1)
slope_Mg = coefficients_Mg[0]
intercept_Mg = coefficients_Mg[1]
print('De helling is ', slope_Mg)
print('De intercept is ', intercept_Mg)

#De MSD en afgeleide samen plotten
plt.plot(bereik_timesteps_Mg, bereik_MSD_Mg, color='green', linestyle='solid', label='MSD_Mg')
regressie_lijn_Mg = slope_Mg * np.array(bereik_timesteps_Mg) + intercept_Mg
plt.plot(bereik_timesteps_Mg, regressie_lijn_Mg, color ='black', linestyle='solid', label='Linear regression')
plt.xlabel('time (fs)', fontsize=18)
plt.ylabel('MSD (A²)', fontsize=18)
plt.title('Mean square displacement of Mg and linear regression in selected area', fontsize=18)
plt.xticks(fontsize=18)
plt.yticks(fontsize=18)
plt.legend(fontsize=18)
plt.show()

#Berekenen van de diffusiecoefficient van Mg
diffusiecoefficient_Mg = (slope_Mg / 6) * (10 ** (-5)) #Eenheid in m²/s²
print('De diffusiecoefficient van Mg bedraagt ', diffusiecoefficient_Mg)

#Inlezen van de Si-Si-SiMSD.dat file
data_Si = np.loadtxt('Si-Si-SiMSD.dat', dtype=float, skiprows=2)

#Analyseer de gegevens
timesteps_Si = data_Si[:, 0]
MSD_Si = data_Si[:, 4]

#Visualisatie van de MSD van Si over de tijd
plt.plot(timesteps_Si, MSD_Si, color='black', linestyle='solid')
plt.xlabel('time (fs)', fontsize=18)
plt.ylabel('MSD (Å²)', fontsize=18)
plt.title('Mean square displacement of Si', fontsize=18)
plt.xticks(fontsize=18)
plt.yticks(fontsize=18)
plt.show()

#Selectie van het middelste bereik (5 - 50 %) waar de diffusiecoefficient een lineare curve voorstelt
bereik_timesteps_Si = []
bereik_MSD_Si = []
for i_Si, element_Si in enumerate(timesteps_Si):
    if element_Si >= 250000 and element_Si <= 2500000:
        bereik_timesteps_Si.append(element_Si) 
        bereik_MSD_Si.append(MSD_Si[i_Si])
        
#Visualisatie van de beperkte MSD van Si over de tijd
plt.plot(bereik_timesteps_Si, bereik_MSD_Si, color='magenta', linestyle='solid')
plt.xlabel('time (fs)', fontsize=18)
plt.ylabel('MSD (A²)', fontsize=18)
plt.title('Mean square displacement of Si in selected area', fontsize=18)
plt.xticks(fontsize=18)
plt.yticks(fontsize=18)
plt.show()

#Lineaire regressie (rechte fitten) door het geselecteerde bereik = het nemen van de afgeleide naar de tijd
coefficients_Si = np.polyfit(bereik_timesteps_Si, bereik_MSD_Si, 1)
slope_Si = coefficients_Si[0]
intercept_Si = coefficients_Si[1]
print('De helling is ', slope_Si)
print('De intercept is ', intercept_Si)

#De MSD en afgeleide samen plotten
plt.plot(bereik_timesteps_Si, bereik_MSD_Si, color='magenta', linestyle='solid', label='MSD_Si')
regressie_lijn_Si = slope_Si * np.array(bereik_timesteps_Si) + intercept_Si
plt.plot(bereik_timesteps_Si, regressie_lijn_Si, color ='black', linestyle='solid', label='Linear regression')
plt.xlabel('time (fs)', fontsize=18)
plt.ylabel('MSD (A²)', fontsize=18)
plt.title('Mean square displacement of Si and linear regression in selected area', fontsize=18)
plt.xticks(fontsize=18)
plt.yticks(fontsize=18)
plt.legend(fontsize=18)
plt.show()

#Berekenen van de diffusiecoefficient van Si
diffusiecoefficient_Si = (slope_Si / 6) * (10 ** (-5)) #Eenheid in m²/s²
print('De diffusiecoefficient van Si bedraagt ', diffusiecoefficient_Si)

#Inlezen van de Ox-Ox-OxMSD.dat file
data_Ox = np.loadtxt('Ox-Ox-OxMSD.dat', dtype=float, skiprows=2)

#Analyseer de gegevens
timesteps_Ox = data_Ox[:, 0]
MSD_Ox = data_Ox[:, 4]

#Visualisatie van de MSD van Ox over de tijd
plt.plot(timesteps_Ox, MSD_Ox, color='purple', linestyle='solid')
plt.xlabel('time (fs)', fontsize=18)
plt.ylabel('MSD (Å²)', fontsize=18)
plt.title('Mean square displacement of Ox', fontsize=18)
plt.xticks(fontsize=18)
plt.yticks(fontsize=18)
#plt.show()

#Selectie van het middelste bereik (5 - 50 %) waar de diffusiecoefficient een lineare curve voorstelt
bereik_timesteps_Ox = []
bereik_MSD_Ox = []
for i_Ox, element_Ox in enumerate(timesteps_Ox):
    if element_Ox >= 250000 and element_Ox <= 2500000:
        bereik_timesteps_Ox.append(element_Ox) 
        bereik_MSD_Ox.append(MSD_Ox[i_Ox])
        
#Visualisatie van de beperkte MSD van Ox over de tijd
plt.plot(bereik_timesteps_Ox, bereik_MSD_Ox, color='purple', linestyle='solid')
plt.xlabel('time (fs)', fontsize=18)
plt.ylabel('MSD (A²)', fontsize=18)
plt.title('Mean square displacement of Ox in selected area', fontsize=18)
plt.xticks(fontsize=18)
plt.yticks(fontsize=18)
plt.show()

#Lineaire regressie (rechte fitten) door het geselecteerde bereik = het nemen van de afgeleide naar de tijd
coefficients_Ox = np.polyfit(bereik_timesteps_Ox, bereik_MSD_Ox, 1)
slope_Ox = coefficients_Ox[0]
intercept_Ox = coefficients_Ox[1]
print('De helling is ', slope_Ox)
print('De intercept is ', intercept_Ox)

#De MSD en afgeleide samen plotten
plt.plot(bereik_timesteps_Ox, bereik_MSD_Ox, color='purple', linestyle='solid', label='MSD_Ox')
regressie_lijn_Ox = slope_Ox * np.array(bereik_timesteps_Ox) + intercept_Ox
plt.plot(bereik_timesteps_Ox, regressie_lijn_Ox, color ='black', linestyle='solid', label='Linear regression')
plt.xlabel('time (fs)', fontsize=18)
plt.ylabel('MSD (A²)', fontsize=18)
plt.title('Mean square displacement of Ox and linear regression in selected area', fontsize=18)
plt.xticks(fontsize=18)
plt.yticks(fontsize=18)
plt.legend(fontsize=18)
plt.show()

#Berekenen van de diffusiecoefficient van Ox
diffusiecoefficient_Ox = (slope_Ox / 6) * (10 ** (-5)) #Eenheid in m²/s²
print('De diffusiecoefficient van Ox bedraagt ', diffusiecoefficient_Ox)

#Theoretische berekening van de diffusiecoefficient van Ox - formule
#Bepalen van de molaire massa's van de metaaloxiden - Inlezen van de file: data.file
with open('data.file', 'r') as bestand:
    lijn = bestand.readlines()
    
#Selectie van de nodige lijnen en uitschrijven naar een data_file_MM.txt
start = 11
end = 16
geselecteerde_regels = lijn[start:end + 1]
    
with open('data_file_MM.txt', 'w') as z:
    for regels in geselecteerde_regels:
        regels = regels.split()
        z.write(' '.join(regels) + '\n')

#Inlezen van de nieuwe file data_file_MM.txt
data_MM = np.loadtxt('data_file_MM.txt', dtype=float)

#Definieren van de molaire massa's
Molaire_massa = [] #lijst van de molaire massa's in volgorde [Si, Al, Ca, Mg, O]
for i in data_MM:
    Molaire_massa.append(float(i[1]))
print('De molaire massas zijn:', Molaire_massa)
MM_Si = float(Molaire_massa[0])
MM_Al = float(Molaire_massa[1])
MM_Ca = float(Molaire_massa[2])
MM_Mg = float(Molaire_massa[3])
MM_Ox = float(Molaire_massa[4])
MM_SiO2 = MM_Si + 2 * MM_Ox
MM_Al2O3 = 2 * MM_Al + 3 * MM_Ox
MM_CaO = MM_Ca + MM_Ox
MM_MgO = MM_Mg + MM_Ox

#Omzetten van gewichtsfracties naar molfracties voor het gebruik van de formules
gewichtsfracties = [0.05, 0.4, 0.06, 0.49] #Gewichtsfracties in volgorde [Al2O3, CaO, MgO, SiO2]  
molaire_massa = [MM_Al2O3, MM_CaO, MM_MgO, MM_SiO2] #Molaire massa in volgorde [Al2O3, CaO, MgO, SiO2] - bron internet
print('De molaire massas van de metaaloxiden zijn ', molaire_massa)
aantal_mol = []
for masfractie, molaire in zip(gewichtsfracties, molaire_massa):
    aantal_mol.append(float((masfractie / molaire)))
print('Het aantal mol van de metaaloxide bedraagt ', aantal_mol)
totaal_aantal_mol = sum(aantal_mol)
print('Het totaal aantal mol bedraagt ', totaal_aantal_mol)
molfracties = []
for fractie in aantal_mol:
    molfracties.append(float(fractie / totaal_aantal_mol))
print('De molfracties zijn ', molfracties) #Tabel bestaande uit 4 molfracties in volgorde [Al2O3, CaO, MgO, SiO2]

#Controle berekenening van de diffusie coefficient van zuurstof - formule
diffusie_coefficient_zuurstof = molfracties[3] / (molfracties[3] + molfracties[0]) * diffusiecoefficient_Si + molfracties[0] / (molfracties[0] + molfracties[3]) * diffusiecoefficient_Al #diffusie coefficienten van Si en Al bepalen en als input ingeven voor de berekening van de diffusie coefficient van zuurstof
print('De theoretisch berekende diffusie coefficient van zuurstof bedraagt ', diffusie_coefficient_zuurstof)

########Berekening van elektrische conductiviteit met Nernst-Einstein (benaderende vergelijking)##########
product_Al = diffusiecoefficient_Al * aantal_Al * (lading_Al ** 2)
product_Ca = diffusiecoefficient_Ca * aantal_Ca * (lading_Ca ** 2)
product_Mg = diffusiecoefficient_Mg * aantal_Mg * (lading_Mg ** 2)
product_Si = diffusiecoefficient_Si * aantal_Si * (lading_Si ** 2)
product_Ox = diffusiecoefficient_Ox * aantal_Ox * (lading_Ox ** 2)
totale_som = product_Al + product_Ca + product_Mg + product_Si + product_Ox
temperatuur = 1823 #uitgedrukt in Kelvin - Lammps input file
boltzmann_constante = 1.38064852 * (10) ** (-23) #uitgedrukt in J/K - bron internet
elementary_charge = 1.602176634 * (10) ** (-19) #uitgedrukt in coulomb -  bron internet
elektrische_conductiviteit_NE = ((elementary_charge ** 2) * totale_som) / (temperatuur * volume * boltzmann_constante)
print('De elektrische conductiviteit via de Nernst-Einstein (benaderend) vergelijking bedraagt ', elektrische_conductiviteit_NE) #Uitgedrukt in Siemens/meter (S/m)

#Aangemaakte bestanden terug verwijderen
os.remove('data_file_atomen.txt')
os.remove('log_lammps_volume.txt')
os.remove('inHPC_final_ladingen.txt')
os.remove('data_file_MM.txt')
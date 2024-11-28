#Script 2: Dichtheid, volume na coolingsstap en totale energy van de coolingdown stap naar de operationele temperatuur plotten in functie van de tijdsstappen
#File: log.lammps
import numpy as np
import os
import matplotlib
matplotlib.use('Qt5Agg') #Deze backend biedt GUI-ondersteuning waardoor de figuren worden weergegeven met plt.show()
import matplotlib.pyplot as plt

#Inlezen van de file: log.lammps
with open('log.lammps', 'r') as file:
    lines = file.readlines()
        
#Selectie van de nodige lijnen en uitschrijven naar een log_lammps.txt - NPT run at operational temperature
start_index = 35533
end_index = 38033
geselecteerde_lijnen = lines[start_index:end_index]
    
with open('log_lammps.txt', 'w') as f:
    for regel in geselecteerde_lijnen:
        regel = regel.split()
        f.write(' '.join(regel) + '\n')        
    
#Inlezen van de nieuwe file log_lammps.txt
data = np.loadtxt('log_lammps.txt', dtype=float, usecols=[0, 1, 3, 4, 7]) #Selecteert de kolommen in volgorde [timestep, temperatuur, density, volume, total_energy]
     
#Lijsten opstellen voor de variabelen timesteps, temperatuur, density, volume en total_energy definieren om gegevens in op te slaan
timesteps = []
temperatuur = []
density = []
volume = []
total_energy = []

#Regels in de data overlopen en toevoegen aan de juiste lijst
for rij in data:
    timesteps.append(float(rij[0]) * 0.001)
    temperatuur.append(float(rij[1]))
    density.append(float(rij[2]))
    volume.append(float(rij[3]))
    total_energy.append(float(rij[4]))
    
#Converteer de lijsten naar arrays
timesteps = np.array(timesteps)
temperatuur = np.array(temperatuur)
density = np.array(density)
volume = np.array(volume)
total_energy = np.array(total_energy) 

gemiddelde_temperature = np.mean(temperatuur)
SA_temperature = np.std(temperatuur)
print('De gemiddelde temperatuur is ', gemiddelde_temperature)
print('De standaardafwijking is ', SA_temperature)

#Plot de temperatuur ten opzichte van de tijdsstappen - T in kelvin
plt.plot(timesteps, temperatuur, color='red', linestyle='--')
#plt.xlim(900, 1000)
#plt.ylim(1650, 2200)
plt.xlabel('Time (ps)', fontsize=18)
plt.ylabel('T (K)', fontsize=18)
plt.title('Relaxation time of the temperature during the cooling-to-operational temperature step', fontsize=18)
plt.xticks(fontsize=18)
plt.yticks(fontsize=18)
plt.show()

#Lineaire regressie (rechte fitten) voor het bepalen van de gemiddelde waarde voor de dichtheid
coefficients_dichtheid = np.polyfit(timesteps, density, 1)
slope_dichtheid = coefficients_dichtheid[0]
intercept_dichtheid = coefficients_dichtheid[1]
gemiddelde_dichtheid = np.mean(density)
SA_dichtheid = np.std(density)
print('De gemiddelde dichtheid is ', gemiddelde_dichtheid)
print('De standaardafwijking is ', SA_dichtheid)

#Plot van de dichtheid ten opzichte van de tijdsstappen - dichteheid in g/cm³
plt.plot(timesteps, density, color='blue', linestyle='--')
regressie_lijn_dichtheid = slope_dichtheid * np.array(timesteps) + intercept_dichtheid
plt.plot(timesteps, regressie_lijn_dichtheid, color ='black', linestyle='solid', label='Linear regression')
#plt.xlim(900, 1000)
#plt.ylim(2.5, 2.7)
plt.xlabel('Time (ps)', fontsize=18) 
plt.ylabel('Density (g/cm³)', fontsize=18)
plt.title('Relaxation time of the density during the cooling-to-operational temperature step', fontsize=18)
plt.xticks(fontsize=18)
plt.yticks(fontsize=18)
plt.legend(fontsize=18)
plt.show()

#Lineaire regressie (rechte fitten) voor het bepalen van de gemiddelde waarde voor het volume
coefficients_volume = np.polyfit(timesteps, volume, 1)
slope_volume = coefficients_volume[0]
intercept_volume = coefficients_volume[1]
gemiddelde_volume = np.mean(volume)
SA_volume = np.std(volume)
print('Het gemiddelde volume is ', gemiddelde_volume)
print('De standaardafwijking is ', SA_volume)

#Plot van het volume ten opzichte van de tijdsstappen - volume in A³
plt.plot(timesteps, volume, color='green', linestyle='--')
regressie_lijn_volume = slope_volume * np.array(timesteps) + intercept_volume
plt.plot(timesteps, regressie_lijn_volume, color ='black', linestyle='solid', label='Linear regression')
#plt.xlim(900, 1000)
#plt.ylim(30000, 33000)
plt.xlabel('Time (ps)', fontsize=18)
plt.ylabel('Volume (A³)', fontsize=18)
plt.title('Relaxation time of the volume during the cooling-to-operational temperature step', fontsize=18)
plt.xticks(fontsize=18)
plt.yticks(fontsize=18)
plt.legend(fontsize=18)
plt.show()

#Plot van de totale energy ten opzichte van de tijdsstappen - totale energy in kcal/mol
plt.plot(timesteps, total_energy, color='black', linestyle='--')
#plt.xlim(900, 1000)
plt.xlabel('Time (ps)', fontsize=18) 
plt.ylabel('Total energy (kcal/mol)', fontsize=18)
plt.title('Relaxation time of the total energy during the cooling-to-operational temperature step', fontsize=18)
plt.xticks(fontsize=18)
plt.yticks(fontsize=18)
plt.show()

#Aangemaakte bestanden terug verwijderen
os.remove('log_lammps.txt')
#Script 3: Bepalen van de coordiantie getallen van Al, Ca, Mg en Si
#Files: Cord_N_Al.txt, Cord_N_Ca.txt, Cord_N_Mg.txt, Cord_N_Si.txt
import numpy as np
import os
import matplotlib
matplotlib.use('Qt5Agg') #Deze backend biedt GUI-ondersteuning waardoor de figuren worden weergegeven met plt.show()
import matplotlib.pyplot as plt
from scipy.interpolate import make_interp_spline

#inlezen van de file Cord_N_Al.txt
with open('Cord_N_Al.txt', 'r') as bestand:
    lines_Al = bestand.readlines()
    
#Verwijder de eerste 4 lijnen
del lines_Al[:4]

#Verwijder de lijnen die een lengte van 6 hebben
lines_Al = [regel_Al for regel_Al in lines_Al if len(regel_Al.split(' ')) != 6]

#Schrijf de lijnen naar een nieuw .txt bestand
with open('Cord_N_Al_bewerkt.txt', 'w') as f:
    f.writelines(lines_Al)

#Inlezen van de file Cord_N_Al_bewerkt.txt
data_Al = np.loadtxt('Cord_N_Al_bewerkt.txt', dtype=float)

#Lijsten opstellen voor het bepalen van de average coordinaat getallen
coordinaten_Al = []
atomen_coordinaatgetal_Al = []
percentueel_coordinaatgetal_Al = []
i = 1 
j = 0 
for rij in data_Al:
    if i == 61: #Dit is coordinaatgetal 3
        coordinaten_Al.append(float(rij[1]))
        atomen_coordinaatgetal_Al.append(float(rij[2]))
        percentueel_coordinaatgetal_Al.append(float(rij[3]))
        i += 1
    elif i == 81: #Dit is coordinaatgetal 4
        coordinaten_Al.append(float(rij[1]))
        atomen_coordinaatgetal_Al.append(float(rij[2]))
        percentueel_coordinaatgetal_Al.append(float(rij[3]))
        i += 1
    elif i == 101: #Dit is coordinaatgetal 5
        coordinaten_Al.append(float(rij[1]))
        atomen_coordinaatgetal_Al.append(float(rij[2]))
        percentueel_coordinaatgetal_Al.append(float(rij[3]))
        index_coordinaatgetal_3 = i - 1
        i += 1
    elif i == 121: #Dit is coordinaatgetal 6
        coordinaten_Al.append(float(rij[1]))
        atomen_coordinaatgetal_Al.append(float(rij[2]))
        percentueel_coordinaatgetal_Al.append(float(rij[3]))
        index_coordinaatgetal_4 = i - 1
        i += 1
    else:
        if i > 200:
            j += 1
            if j == 61: #Waarden toevoegen aan de gegevens van cordinaatgetal 3, dus index 0
                atomen_coordinaatgetal_Al[0] += float(rij[2])
                percentueel_coordinaatgetal_Al[0] += float(rij[3])
            elif j == 81: #Waarden toevoegen aan de gegevens van cordinaatgetal 4, dus index 1
                atomen_coordinaatgetal_Al[1] += float(rij[2])
                percentueel_coordinaatgetal_Al[1] += float(rij[3])
            elif j == 101: #Waarden toevoegen aan de gegevens van cordinaatgetal 5, dus index 2
                atomen_coordinaatgetal_Al[2] += float(rij[2])
                percentueel_coordinaatgetal_Al[2] += float(rij[3])
            elif j == 121: #Waarden toevoegen aan de gegevens van cordinaatgetal 6, dus index 3
                atomen_coordinaatgetal_Al[3] += float(rij[2])
                percentueel_coordinaatgetal_Al[3] += float(rij[3])
            elif j == 200:
                j = 0
        i += 1
#Converteer de lijsten atomen_cordinaatgetal_Al en percentueel_coordinaatgetal_Al naar arrays
atomen_coordinaatgetal_Al = np.array(atomen_coordinaatgetal_Al)
coordinaten_Al = np.array(coordinaten_Al)
percentueel_coordinaatgetal_Al = np.array(percentueel_coordinaatgetal_Al)

#Het gemiddelde van het aantal atomen en percentueel voorkomen voor de coordinaatgetallen berekenen
tijdsstap = 100
tijdsduur = 5000000
aantal_frequenties = tijdsduur / tijdsstap
atomen_coordinaatgetal_Al_gemiddeld = atomen_coordinaatgetal_Al / aantal_frequenties
percentueel_coordinaatgetal_Al_gemiddeld = np.array((percentueel_coordinaatgetal_Al / aantal_frequenties) * 100)
afgerond_percentueel_coordinaatgetal_Al_gemiddeld = [round(getal) for getal in percentueel_coordinaatgetal_Al_gemiddeld]
afgerond_coordinaten_Al = [round(num) for num in coordinaten_Al]
print('Het percentueel voorkomen van coordinaatgetal {} bedraagt {} %'.format(afgerond_coordinaten_Al[0], afgerond_percentueel_coordinaatgetal_Al_gemiddeld[0]))
print('Het percentueel voorkomen van coordinaatgetal {} bedraagt {} %'.format(afgerond_coordinaten_Al[1], afgerond_percentueel_coordinaatgetal_Al_gemiddeld[1]))
print('Het percentueel voorkomen van coordinaatgetal {} bedraagt {} %'.format(afgerond_coordinaten_Al[2], afgerond_percentueel_coordinaatgetal_Al_gemiddeld[2]))
print('Het percentueel voorkomen van coordinaatgetal {} bedraagt {} %'.format(afgerond_coordinaten_Al[3], afgerond_percentueel_coordinaatgetal_Al_gemiddeld[3]))

#Plotten van het percentueel voorkomen van de coordinatie getallen als een continue curve
vloeiende_xwaarden = np.linspace(min(coordinaten_Al), max(coordinaten_Al))
spl = make_interp_spline(coordinaten_Al, percentueel_coordinaatgetal_Al_gemiddeld, k=3)
vloeiende_ywaarden = spl(vloeiende_xwaarden)
plt.plot(vloeiende_xwaarden, vloeiende_ywaarden, color='black', label='Al-Ox')
plt.xticks(range(0, 13), fontsize=18)
plt.yticks(np.arange(0, 101, 10), fontsize=18)
plt.xlabel('Coordination number', fontsize=18)
plt.ylabel('% occurence', fontsize=18)
plt.title('Probability of occurence for the coordination numbers of Al-Ox', fontsize=18)
plt.legend(fontsize=18)
plt.show()

#Aangemaakte bestanden terug verwijderen
os.remove('Cord_N_Al_bewerkt.txt')
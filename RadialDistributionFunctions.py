#Script 1: Plotten van RDFs Al-O, Si-O, Mg-O, Ca-O
#File: RDF_op.txt
import numpy as np
import os
import matplotlib
matplotlib.use('Qt5Agg') #Deze backend biedt GUI-ondersteuning waardoor de figuren worden weergegeven met plt.show()
import matplotlib.pyplot as plt

#Inlezen van de file RDF_op
with open('RDF_op.txt', 'r') as bestand:
    lijnen = bestand.readlines()
    
#Verwijder de eerste 4 lijnen
del lijnen[:4]

#Verwijder de lijnen die een lengte van twee hebben
lijnen = [regel for regel in lijnen if len(regel.split(' ')) != 2]
        
#Schrijf de overgebleven lijnen naar een nieuw .txt bestand
with open('RDF_op_bewerkt.txt', 'w') as f:
    f.writelines(lijnen)

#Inlezen van de file RDF_op_bewerkt.txt
data = np.loadtxt('RDF_op_bewerkt.txt', dtype=float, usecols=[0, 1, 10, 18, 24, 28]) #usecols in volgorde [bins, afstand x-as, RDF Si-Ox, RDF Al-Ox, RDF Ca-Ox, RDF Mg-Ox]

#Lijsten opstellen met de waarden van de RDF om het gemiddelde te berekenen
RDF_Si = []
RDF_Al = []
RDF_Ca = []
RDF_Mg = []
afstanden_xas = []
j = 0
i = 0
for rij in data:
        if j < 199:
            RDF_Si.append(float(rij[2]))
            RDF_Al.append(float(rij[3]))
            RDF_Ca.append(float(rij[4]))
            RDF_Mg.append(float(rij[5]))
            afstanden_xas.append(float(rij[1]))
            j += 1
        elif j == 199:
            RDF_Si.append(float(rij[2]))
            RDF_Al.append(float(rij[3]))
            RDF_Ca.append(float(rij[4]))
            RDF_Mg.append(float(rij[5]))
            afstanden_xas.append(float(rij[1]))
            j += 1
        else: #index i moet resetten wanneer het een bepaalde waarde 199 bereikt heeft
            RDF_Si[i] += float(rij[2])
            RDF_Al[i] += float(rij[3])
            RDF_Ca[i] += float(rij[4])
            RDF_Mg[i] += float(rij[5])
            afstanden_xas[i] += float(rij[1])
            if i == 199:
                RDF_Si[i] += float(rij[2])
                RDF_Al[i] += float(rij[3])
                RDF_Ca[i] += float(rij[4])
                RDF_Mg[i] += float(rij[5])
                afstanden_xas[i] += float(rij[1])
                i = -1
            i += 1
            j += 1

#Converteer de lijsten naar arrays
RDF_Si = np.array(RDF_Si)
RDF_Al = np.array(RDF_Al)
RDF_Ca = np.array(RDF_Ca)
RDF_Mg = np.array(RDF_Mg)
afstanden_xas = np.array(afstanden_xas)

#Het gemiddelde van de RDFs berekenen over alle frequenties
tijdsstap = 1000
tijdsduur = 5000000
aantal_frequenties = tijdsduur / tijdsstap
RDF_Si_gemiddeld = RDF_Si / aantal_frequenties
RDF_Al_gemiddeld = RDF_Al / aantal_frequenties
RDF_Ca_gemiddeld = RDF_Ca / aantal_frequenties
RDF_Mg_gemiddeld = RDF_Mg / aantal_frequenties
afstanden_xas_gemiddeld = afstanden_xas / aantal_frequenties

#Plot van de RDF Si-O
plt.plot(afstanden_xas_gemiddeld, RDF_Si_gemiddeld, color='black', linestyle='solid')
plt.xlabel('r(Ä)', fontsize=18)
plt.ylabel('g(r)', fontsize=18)
plt.title('Radial distribution function of Si-Ox', fontsize=18)
plt.xticks(fontsize=18)
plt.yticks(fontsize=18)
plt.show()

#Bepalen van de bond length tussen Si-O -> het maximum van de radial distribution curve
max_value_Si = np.max(RDF_Si_gemiddeld)
max_index_Si = np.argmax(RDF_Si_gemiddeld)
bond_length_Si = afstanden_xas_gemiddeld[max_index_Si]
print('De bond lengte voor de binding tussen Si-O bedraagt:', bond_length_Si)

#Bepalen van de cut off waardetussen Si-O -> het minimum na de eerste piek van de radial distribution curve
min_value_Si = np.min(RDF_Si_gemiddeld[max_index_Si:])
min_index_Si = max_index_Si + np.argmin(RDF_Si_gemiddeld[max_index_Si:])
cut_off_Si = afstanden_xas_gemiddeld[min_index_Si]
print('De cut off afstand voor de binding Si-O bedraagt:', cut_off_Si)

#Plot van de RDF Al-O
plt.plot(afstanden_xas_gemiddeld, RDF_Al_gemiddeld, color='red', linestyle='solid')
plt.xlabel('r(Ä)', fontsize=18)
plt.ylabel('g(r)', fontsize=18)
plt.title('Radial distribution function of Al-Ox', fontsize=18)
plt.xticks(fontsize=18)
plt.yticks(fontsize=18)
plt.show()

#Bepalen van de bond length tussen Al-O -> het maximum van de radial distribution curve
max_value_Al = np.max(RDF_Al_gemiddeld)
max_index_Al = np.argmax(RDF_Al_gemiddeld)
bond_length_Al = afstanden_xas_gemiddeld[max_index_Al]
print('De bond lengte voor de binding tussen Al-O bedraagt:', bond_length_Al)

#Bepalen van de cut off waardetussen Al-O -> het minimum na de eerste piek van de radial distribution curve
min_value_Al = np.min(RDF_Al_gemiddeld[max_index_Al:])
min_index_Al = max_index_Al + np.argmin(RDF_Al_gemiddeld[max_index_Al:])
cut_off_Al = afstanden_xas_gemiddeld[min_index_Al]
print('De cut off afstand voor de binding Al-O bedraagt:', cut_off_Al)

#Plot van de RDF Ca-O
plt.plot(afstanden_xas_gemiddeld, RDF_Ca_gemiddeld, color='blue', linestyle='solid')
plt.xlabel('r(Ä)', fontsize=18)
plt.ylabel('g(r)', fontsize=18)
plt.title('Radial distribution function of Ca-Ox', fontsize=18)
plt.xticks(fontsize=18)
plt.yticks(fontsize=18)
plt.show()

#Bepalen van de bond length tussen Ca-O -> het maximum van de radial distribution curve
max_value_Ca = np.max(RDF_Ca_gemiddeld)
max_index_Ca = np.argmax(RDF_Ca_gemiddeld)
bond_length_Ca = afstanden_xas_gemiddeld[max_index_Ca]
print('De bond lengte voor de binding tussen Ca-O bedraagt:', bond_length_Ca)

#Bepalen van de cut off waardetussen Ca-O -> het minimum na de eerste piek van de radial distribution curve
min_value_Ca = np.min(RDF_Ca_gemiddeld[max_index_Ca:])
min_index_Ca = max_index_Ca + np.argmin(RDF_Ca_gemiddeld[max_index_Ca:])
cut_off_Ca = afstanden_xas_gemiddeld[min_index_Ca]
print('De cut off afstand voor de binding Ca-O bedraagt:', cut_off_Ca)

#Plot van de RDF Mg-O
plt.plot(afstanden_xas_gemiddeld, RDF_Mg_gemiddeld, color='green', linestyle='solid')
plt.xlabel('r(Ä)', fontsize=18)
plt.ylabel('g(r)', fontsize=18)
plt.title('Radial distribution function of Mg-Ox', fontsize=18)
plt.xticks(fontsize=18)
plt.yticks(fontsize=18)
plt.show()

#Bepalen van de bond length tussen Mg-O -> het maximum van de radial distribution curve
max_value_Mg = np.max(RDF_Mg_gemiddeld)
max_index_Mg = np.argmax(RDF_Mg_gemiddeld)
bond_length_Mg = afstanden_xas_gemiddeld[max_index_Mg]
print('De bond lengte voor de binding tussen Mg-O bedraagt:', bond_length_Mg)

#Bepalen van de cut off waardetussen Mg-O -> het minimum na de eerste piek van de radial distribution curve
min_value_Mg = np.min(RDF_Mg_gemiddeld[max_index_Mg:])
min_index_Mg = max_index_Mg + np.argmin(RDF_Mg_gemiddeld[max_index_Mg:])
cut_off_Mg = afstanden_xas_gemiddeld[min_index_Mg]
print('De cut off afstand voor de binding Mg-O bedraagt:', cut_off_Mg)

#Plotten van alle metaal-oxygen paren samen
plt.plot(afstanden_xas_gemiddeld, RDF_Si_gemiddeld, color='black', linestyle='solid', label='Si-Ox')
plt.plot(afstanden_xas_gemiddeld, RDF_Al_gemiddeld, color='red', linestyle='solid', label='Al-Ox')
plt.plot(afstanden_xas_gemiddeld, RDF_Ca_gemiddeld, color='blue', linestyle='solid', label='Ca-Ox')
plt.plot(afstanden_xas_gemiddeld, RDF_Mg_gemiddeld, color='green', linestyle='solid', label='Mg-Ox')
plt.xlim(0, 8)
plt.xlabel('r(Ä)', fontsize=18)
plt.ylabel('g(r)', fontsize=18)
plt.xticks(fontsize=18)
plt.yticks(fontsize=18)
plt.legend(fontsize=18)
plt.show()

#Aangemaakte bestanden terug verwijderen
os.remove('RDF_op_bewerkt.txt')
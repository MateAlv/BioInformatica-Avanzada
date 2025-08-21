from Bio.Seq import Seq
import random
import numpy as np
import matplotlib.pyplot as plt

# 2A) Obtenga la reversa complementaria y compara las estadísticas de la misma con la original.
# 2B) Traduzca la secuencia (obteniendo una secuencia proteica) y obtenga la distribución de los largos de los Marcos 
# Abiertos de Lectura (ORF de Open Reading Frames). Compare los largos de los ORF con el largo de las proteinas analizadas 
# previamente

LARGO_SEQ = 999999
DECIMALES_FREQ = 4

NUCLEOTIDOS = ['A', 'T', 'C', 'G']

AMINOACIDOS = {
    "A": 0.0656,    "C": 0.0328,    "D": 0.0328,
    "E": 0.0328,    "F": 0.0328,    "G": 0.0656,
    "H": 0.0328,    "I": 0.0492,    "K": 0.0328,
    "L": 0.0984,    "M": 0.0164,    "N": 0.0328,
    "P": 0.0656,    "Q": 0.0328,    "R": 0.0984,
    "S": 0.0984,    "T": 0.0656,    "V": 0.0656,
    "W": 0.0164,    "Y": 0.0328
}

CODONES = {
    "ATA":"I", "ATC":"I", "ATT":"I", "ATG":"M",
    "ACA":"T", "ACC":"T", "ACG":"T", "ACT":"T",
    "AAC":"N", "AAT":"N", "AAA":"K", "AAG":"K",
    "AGC":"S", "AGT":"S", "AGA":"R", "AGG":"R",
    "CTA":"L", "CTC":"L", "CTG":"L", "CTT":"L",
    "CCA":"P", "CCC":"P", "CCG":"P", "CCT":"P",
    "CAC":"H", "CAT":"H", "CAA":"Q", "CAG":"Q",
    "CGA":"R", "CGC":"R", "CGG":"R", "CGT":"R",
    "GTA":"V", "GTC":"V", "GTG":"V", "GTT":"V",
    "GCA":"A", "GCC":"A", "GCG":"A", "GCT":"A",
    "GAC":"D", "GAT":"D", "GAA":"E", "GAG":"E",
    "GGA":"G", "GGC":"G", "GGG":"G", "GGT":"G",
    "TCA":"S", "TCC":"S", "TCG":"S", "TCT":"S",
    "TTC":"F", "TTT":"F", "TTA":"L", "TTG":"L",
    "TAC":"Y", "TAT":"Y", "TAA":"*", "TAG":"*",
    "TGC":"C", "TGT":"C", "TGA":"*", "TGG":"W"
}


def seq_adn_random():
    seq_random = ''
    for i in range(LARGO_SEQ):
        seq_random += random.choice(NUCLEOTIDOS)

    return seq_random

def determinar_frecuencia(secuencia, simbolos):
    total = len(secuencia)
    frecuencias = {}
    
    for simbolo in simbolos:
        conteo = secuencia.count(simbolo)
        frecuencia = round(conteo / total, DECIMALES_FREQ)  
        frecuencias[simbolo] = frecuencia
    
    return frecuencias

def main_2a():
    seq_nucleotidos = seq_adn_random()
    seq_nuc_lista = list(seq_nucleotidos)
    seq_prot = ''
    
    for i in range(0, len(seq_nuc_lista) - 2, 3):
        codon = ''
        for j in range(3):
            codon += seq_nuc_lista[i + j]
        seq_prot += CODONES.get(codon, 'X')
    
    freq_adn = determinar_frecuencia(seq_nucleotidos, NUCLEOTIDOS)
    freq_amino = determinar_frecuencia(seq_prot, list(AMINOACIDOS.keys()))

    print(freq_adn)
    print(freq_amino)


    # Gráfico de frecuencias de nucleótidos
    plt.figure(figsize=(8,4))
    plt.bar(freq_adn.keys(), freq_adn.values(), color='skyblue')
    plt.xlabel('Nucleótidos')
    plt.ylabel('Frecuencia')
    plt.title('Frecuencia de nucleótidos en ADN aleatorio')
    plt.ylim(0, 1)
    plt.show()

    # Gráfico de frecuencias de aminoácidos
    plt.figure(figsize=(12,5))
    plt.bar(freq_amino.keys(), freq_amino.values(), color='salmon')
    plt.xlabel('Aminoácidos')
    plt.ylabel('Frecuencia')
    plt.title('Frecuencia de aminoácidos en proteína traducida')
    plt.ylim(0, max(list(freq_amino.values()))* 1.1)
    plt.show()
    
main_2a()

# El gráfico nos muestra como a partir de secuencias de nucleótidos aleatorios, obviamente la frequencia de cada base
# está perfectamente distribuída. Esto solamente evidencia el correcto funcionamiento de la librería random de python,
# ya que la instrucción era general bases totalmente al azar y se distribuyeron de forma 100% pareja.
# Con la distribución de aminoácidos pasa otra historia. Al haber muchos codones que producen el mismo aminoácido en
# la cadena, esto provoca que esos aminoácidos que se "repiten" aparezcan sobre-representados, lo que se evidencia
# claramente en el gráfico. Es fácil ver con la distribución generada, cuales AA son generados por más de un codón.
# Pensando más por el lado biológico, el análisis se queda corto ya que la evolución te sesga totalmente los genomas.
# Sin embargo, sí se puede extraer la conclusión de que los AA que surgen de varios codones distintos, van a producirse
# más en los casos de mutaciones y deberían estar mayormente representados.
# - Mate




def obtener_largos_orfs(proteina):
    orfs = []
    i = 0
    while i < len(proteina):
        if proteina[i] == 'M': 
            largo = 1
            i += 1
            while i < len(proteina) and proteina[i] != '*':
                largo += 1
                i += 1
            orfs.append(largo)
        else:
            i += 1
    return orfs

def main_2b():
    seq_nucleotidos = seq_adn_random()
    seq_nuc_lista = list(seq_nucleotidos)
    
    marcos = []
    for offset in range(3):
        seq_prot = ''
        for i in range(offset, len(seq_nuc_lista) - 2, 3):
            codon = ''.join(seq_nuc_lista[i:i+3])
            seq_prot += CODONES.get(codon, 'X')
        marcos.append(seq_prot)
    
    orfs_todos = []
    for proteina in marcos:
        orfs = obtener_largos_orfs(proteina)
        orfs_todos.extend(orfs)
    
    plt.figure(figsize=(12,5))
    plt.hist(orfs_todos, bins=50, color='lightgreen', edgecolor='black')
    plt.xlabel('Largo del ORF (número de aminoácidos)')
    plt.ylabel('Cantidad de ORFs')
    plt.title('Distribución de largos de ORFs en proteína traducida')
    plt.show()
    
    print(f"Número total de ORFs encontrados: {len(orfs_todos)}")
    print(f"Largo promedio de los ORFs: {np.mean(orfs_todos):.2f}")

main_2b()

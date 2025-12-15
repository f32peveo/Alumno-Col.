import loki as lk
import matplotlib.pyplot as plt
import numpy as np
import os

def cargar_energias(path="databaseStateEnergyHe.txt"):
    base = os.path.dirname(os.path.abspath(__file__))
    path = os.path.abspath(os.path.join(base, path))
    
    energy_dict = {}
    with open(path, "r") as f:
        for linea in f:
            if not linea.startswith('%') and linea.strip():
                parts = linea.strip().split()
                species = ' '.join(parts[:-1])
                energy = float(parts[-1])
                energy_dict[species] = energy
    return energy_dict

def mostrar_matriz(matriz, speciesList, titulo, cmap='Blues'):
    base = os.path.dirname(os.path.abspath(__file__))
    plt.figure(figsize=(10, 6))
    plt.imshow(matriz, cmap=cmap, interpolation='nearest', aspect='auto')
    plt.title(titulo)
    plt.xlabel("Especies")
    plt.ylabel("Reacciones")
    plt.xticks(range(len(speciesList)), speciesList, rotation=90)
    plt.tight_layout()
    plt.savefig(os.path.join(base,f"{titulo.replace(' ', '_').lower()}.png"))
    plt.show()
    plt.close()

def calcular_grados(matriz):
    return [sum(col) for col in zip(*matriz)]

def mostrar_histogramas_por_especie(speciesList, reactantsDegree, productsDegree, ponderado=False):
    fig, axes = plt.subplots(1, 2, figsize=(12, 6))
    label = "Grado Ponderado" if ponderado else "Grado"

    axes[0].bar(speciesList, reactantsDegree, color='blue', alpha=0.7)
    axes[0].set_title(f"{label} de Nodos - Reactivos")
    axes[0].set_xlabel("Especies")
    axes[0].set_ylabel(label)
    axes[0].tick_params(axis='x', rotation=90)

    axes[1].bar(speciesList, productsDegree, color='orange', alpha=0.7)
    axes[1].set_title(f"{label} de Nodos - Productos")
    axes[1].set_xlabel("Especies")
    axes[1].set_ylabel(label)
    axes[1].tick_params(axis='x', rotation=90)

    plt.tight_layout()
    plt.show()
    plt.close()

def mostrar_histograma_global(reactantsDegree, productsDegree, rango=(0, 100)):
    plt.hist(reactantsDegree, bins=50, color='blue', alpha=0.7, label='Reactivos', range=rango)
    plt.hist(productsDegree, bins=50, color='orange', alpha=0.7, label='Productos', range=rango)
    plt.legend()
    plt.xlabel("Grado Ponderado")
    plt.ylabel("Número de especies")
    plt.title("Grado de conectividad")
    plt.show()
    plt.close()


def main(file):
    base = os.path.dirname(os.path.abspath(__file__))
    file = os.path.abspath(os.path.join(base, file))

    uniqueSpecies, reactions = lk.parseChemFile(file)

    reactantsMatrix = []
    productsMatrix = []

    for reaction in reactions:
        reactantsRow = {s: 0 for s in uniqueSpecies}
        productsRow = {s: 0 for s in uniqueSpecies}

        for s, c in zip(reaction['lhsSpecies'], reaction['lhsStoichiometricCoeffs']):
            reactantsRow[s] = c
        for s, c in zip(reaction['rhsSpecies'], reaction['rhsStoichiometricCoeffs']):
            productsRow[s] = c

        reactantsMatrix.append(list(reactantsRow.values()))
        productsMatrix.append(list(productsRow.values()))

    np.save(os.path.join(base, 'reactantsMatrix_w.npy'), reactantsMatrix)
    np.save(os.path.join(base, 'productsMatrix_w.npy'), productsMatrix)

    reactantsDegree = calcular_grados(reactantsMatrix)
    productsDegree = calcular_grados(productsMatrix)

    np.save(os.path.join(base, 'reactantsDegree_w.npy'), reactantsDegree)
    np.save(os.path.join(base, 'productsDegree_w.npy'), productsDegree)

    energy_dict = cargar_energias()
    ordered_species = [s for s, _ in sorted(energy_dict.items(), key=lambda x: x[1]) if s in uniqueSpecies]
    np.save(os.path.join(base, 'ordered_species.npy'), np.array(ordered_species))

    reactantsDegree_sorted = [reactantsDegree[uniqueSpecies.index(s)] for s in ordered_species]
    productsDegree_sorted = [productsDegree[uniqueSpecies.index(s)] for s in ordered_species]

    np.save(os.path.join(base, 'reactantsDegree_sorted.npy'), np.array(reactantsDegree_sorted))
    np.save(os.path.join(base, 'productsDegree_sorted.npy'), np.array(productsDegree_sorted))

    mostrar_matriz(reactantsMatrix, ordered_species, "Matriz Ponderada Reactivos", cmap='Blues')
    mostrar_matriz(productsMatrix, ordered_species, "Matriz Ponderada Productos", cmap='Blues')
    plt.show()
    plt.close()

    mostrar_histogramas_por_especie(ordered_species, reactantsDegree_sorted, productsDegree_sorted, ponderado=True)
    mostrar_histograma_global(reactantsDegree, productsDegree)

if __name__ == "__main__":
    main('helium.chem')

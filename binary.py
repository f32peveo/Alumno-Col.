import loki as lk
import matplotlib.pyplot as plt
import numpy as np

def cargar_energias(path="databaseStateEnergyHe.txt"):
    energy_dict = {}
    with open(path, "r") as f:
        for linea in f:
            if not linea.startswith('%') and linea.strip():
                parts = linea.strip().split()
                species = ' '.join(parts[:-1])
                energy = float(parts[-1])
                energy_dict[species] = energy
    return energy_dict

def main(file):
    uniqueSpecies, reactions = lk.parseChemFile(file)

    reactantsMatrix = []
    productsMatrix = []

    for reaction in reactions:
        reactantsRow = {species: 0 for species in uniqueSpecies}
        productsRow = {species: 0 for species in uniqueSpecies}

        for species in reaction['lhsSpecies']:
            reactantsRow[species] = 1
        for species in reaction['rhsSpecies']:
            productsRow[species] = 1

        reactantsMatrix.append(list(reactantsRow.values()))
        productsMatrix.append(list(productsRow.values()))

    np.save('reactantsMatrix.npy', reactantsMatrix)
    np.save('productsMatrix.npy', productsMatrix)

    reactantsDegree = [sum(col) for col in zip(*reactantsMatrix)]
    productsDegree = [sum(col) for col in zip(*productsMatrix)]

    np.save('reactantsDegree.npy', reactantsDegree)
    np.save('productsDegree.npy', productsDegree)

    # Ordenar por energía
    energy_dict = cargar_energias()
    ordered_species = [s for s, _ in sorted(energy_dict.items(), key=lambda x: x[1]) if s in uniqueSpecies]
    np.save('ordered_species.npy', np.array(ordered_species))

    # Reordenar los grados
    reactantsDegree_sorted = [reactantsDegree[uniqueSpecies.index(s)] for s in ordered_species]
    productsDegree_sorted = [productsDegree[uniqueSpecies.index(s)] for s in ordered_species]
    np.save(r"C:\Users\isafe\OneDrive\Desktop\Escritorio\FÍSICA\TFG\TFG\reactantsDegree_sorted.npy", np.array(reactantsDegree_sorted))
    np.save(r"C:\Users\isafe\OneDrive\Desktop\Escritorio\FÍSICA\TFG\TFG\rproductsDegree_sorted.npy", np.array(productsDegree_sorted))

    # Gráficas
    mostrar_matriz_binaria(reactantsMatrix, ordered_species, "Matriz binaria de Reactivos")
    mostrar_matriz_binaria(productsMatrix, ordered_species, "Matriz binaria de Productos")
    plt.show()

    fig, axes = plt.subplots(1, 2, figsize=(12, 6))

    axes[0].bar(ordered_species, reactantsDegree_sorted, color='blue', alpha=0.7)
    axes[0].set_title("Grado de Nodos - Reactivos")
    axes[0].set_xlabel("Especies")
    axes[0].set_ylabel("Grado")
    axes[0].tick_params(axis='x', rotation=90)

    axes[1].bar(ordered_species, productsDegree_sorted, color='orange', alpha=0.7)
    axes[1].set_title("Grado de Nodos - Productos")
    axes[1].set_xlabel("Especies")
    axes[1].set_ylabel("Grado")
    axes[1].tick_params(axis='x', rotation=90)

    plt.tight_layout()
    plt.show()

    plt.hist(reactantsDegree, bins=50, color='blue', alpha=0.7, label='Reactivos', range=(0, 100))
    plt.hist(productsDegree, bins=50, color='orange', alpha=0.7, label='Productos', range=(0, 100))
    plt.legend()
    plt.show()

def mostrar_matriz_binaria(matriz, speciesList, titulo):
    plt.figure(figsize=(10, 6))
    plt.imshow(matriz, cmap='binary', interpolation='nearest', aspect='auto')
    plt.title(titulo)
    plt.xlabel("Especies")
    plt.ylabel("Reacciones")
    plt.xticks(range(len(speciesList)), speciesList, rotation=90)
    plt.tight_layout()

if __name__ == "__main__":
    main('helium.chem')

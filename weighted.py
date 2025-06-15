import loki as lk
import matplotlib.pyplot as plt
import numpy as np

# === FUNCIONES AUXILIARES ===

def mostrar_matriz(matriz, speciesList, titulo, cmap='Blues'):
    plt.figure(figsize=(10, 6))
    plt.imshow(matriz, cmap=cmap, interpolation='nearest', aspect='auto')
    plt.title(titulo)
    plt.xlabel("Especies")
    plt.ylabel("Reacciones")
    plt.xticks(range(len(speciesList)), speciesList, rotation=90)
    plt.tight_layout()

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

def mostrar_histograma_global(reactantsDegree, productsDegree, rango=(0, 100)):
    plt.hist(reactantsDegree, bins=50, color='blue', alpha=0.7, label='Reactivos', range=rango)
    plt.hist(productsDegree, bins=50, color='orange', alpha=0.7, label='Productos', range=rango)
    plt.legend()
    plt.xlabel("Grado Ponderado")
    plt.ylabel("Frecuencia")
    plt.title("Distribución de Grado Ponderado")
    plt.show()

# === PROGRAMA PRINCIPAL ===

def main(file):
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

    # Guardar matrices con sufijo _w
    np.save('reactantsMatrix_w.npy', reactantsMatrix)
    np.save('productsMatrix_w.npy', productsMatrix)

    mostrar_matriz(reactantsMatrix, uniqueSpecies, "Matriz Ponderada Reactivos", cmap='Blues')
    mostrar_matriz(productsMatrix, uniqueSpecies, "Matriz Ponderada Productos", cmap='Blues')
    plt.show()

    reactantsDegree = calcular_grados(reactantsMatrix)
    productsDegree = calcular_grados(productsMatrix)

    # Guardar grados con sufijo _w
    np.save('reactantsDegree_w.npy', reactantsDegree)
    np.save('productsDegree_w.npy', productsDegree)

    mostrar_histogramas_por_especie(uniqueSpecies, reactantsDegree, productsDegree, ponderado=True)
    mostrar_histograma_global(reactantsDegree, productsDegree)

if __name__ == "__main__":
    main('helium.chem')

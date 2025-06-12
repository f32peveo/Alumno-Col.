import numpy as np
import matplotlib.pyplot as plt
import os
from loki import parseChemFile

def load_energies():
    path = os.path.join(os.path.dirname(__file__), "databaseStateEnergyHe.txt")
    energy_dict = {}
    with open(path, 'r') as file:
        for line in file:
            line = line.strip()
            if not line or line.startswith('%'):
                continue  # ignorar encabezados y líneas vacías
            parts = line.split()
            if len(parts) == 2:
                species, energy = parts
                energy_dict[species] = float(energy)
            else:
                print(f"⚠️ Línea ignorada por formato inesperado: {line}")
    return energy_dict

def plot_connectivity_vs_energy(energy_dict, connectivity, species_names, title):
    x = []
    y = []
    labels = []

    for i, species in enumerate(species_names):
        if species in energy_dict:
            x.append(energy_dict[species])
            y.append(connectivity[i])
            labels.append(species)
        else:
            print(f"⚠️ Warning: '{species}' not found in energy data.")

    plt.figure(figsize=(10, 6))
    plt.scatter(x, y, color='purple', alpha=0.7)

    for xi, yi, label in zip(x, y, labels):
        if yi > np.percentile(y, 95) or xi == min(x) or xi == max(x):
            plt.text(xi, yi, label, fontsize=8, ha='right')

    max_idx = np.argmax(y)
    plt.text(x[max_idx], y[max_idx], labels[max_idx],
             fontsize=9, fontweight='bold', color='red')

    plt.title(title)
    plt.xlabel("Energy (eV)")
    plt.ylabel("Connectivity")
    plt.grid(True)
    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    # Cargar conectividad desde archivo local
    connectivity = np.load("reactantsDegree.npy")  # o "productsDegree.npy"

    # Cargar especies desde archivo .chem con ruta completa
    chem_path = r"C:\Users\isafe\OneDrive\Desktop\Escritorio\FÍSICA\TFG\helium.chem"
    species, _ = parseChemFile(chem_path)

    # Cargar energías desde archivo local
    energies = load_energies()

    # Graficar
    plot_connectivity_vs_energy(energies, connectivity, species, title="Reactants Connectivity vs Energy")

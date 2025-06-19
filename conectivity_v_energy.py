import numpy as np
import matplotlib.pyplot as plt
import os

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
    # Cargar arrays alineados y ordenados
    ordered_species = np.load("ordered_species.npy", allow_pickle=True)
    reactantsDegree_sorted = np.load("reactantsDegree_sorted.npy")
    productsDegree_sorted = np.load("productsDegree_sorted.npy")

    # Cargar energías desde archivo local
    energies = load_energies()

    # Plot para reactivos
    plot_connectivity_vs_energy(
        energies, reactantsDegree_sorted, ordered_species, title="Reactants Connectivity vs Energy"
    )

    # Plot para productos
    plot_connectivity_vs_energy(
        energies, productsDegree_sorted, ordered_species, title="Products Connectivity vs Energy"
    )
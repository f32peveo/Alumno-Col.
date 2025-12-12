import numpy as np
import matplotlib.pyplot as plt
import os
from matplotlib import transforms


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

def plot_connectivity_vs_energy(energy_dict, connectivity, species_names, title, xscale='log', yscale='log'):
    x = []
    y = []
    labels = []

    for i, species in enumerate(species_names):
        print(f"{species} - {energy_dict.get(species, 'No data')}")
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

    # plt.title(title)
    plt.xlabel("Energy (eV)")
    plt.ylabel("Connectivity")
    plt.xscale(xscale)
    plt.yscale(yscale)
    plt.grid(True)
    plt.tight_layout()
    plt.show()

def plot_connectivity_vs_energy_xspace(energy_dict, connectivity, species_names, title, xscale='log', yscale='log'):
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

    # Construimos el eje x xspace
    x_axis = np.linspace(min(x), max(x), len(x))

    plt.figure(figsize=(10, 6))
    plt.scatter(x_axis, y, color='purple', alpha=0.7)

    # for xi, yi, label in zip(x_axis, y, labels):
    #     if yi > np.percentile(y, 95) or xi == min(x_axis) or xi == max(x_axis):
    #         plt.text(xi, yi, label, fontsize=8, ha='right')

    # max_idx = np.argmax(y)
    # plt.text(x[max_idx], y[max_idx], labels[max_idx],
    #          fontsize=9, fontweight='bold', color='red')
    
    # for xi, yi, label in zip(x_axis, y, labels):
        # plt.text(xi, yi, label, fontsize=4, ha='right')
# Ajustar posición de las etiquetas para evitar solapamiento

    ax = plt.gca()
    # for xi, yi, label in zip(x_axis, y, labels):
        # Offset en x y y para separar etiquetas
        # offset_x = 2
        # offset_y = 2
        # trans_offset = transforms.offset_copy(ax.transData, fig=plt.gcf(),
                                            #  x=offset_x, y=offset_y, units='points')
        # ax.text(xi, yi, label, fontsize=6, ha='left', va='bottom', transform=trans_offset)
    
    for i, (xi, yi, label) in enumerate(zip(x_axis, y, labels)):
        offset_x = 2
        offset_y = 4 if i % 2 == 0 else -4

        trans_offset = transforms.offset_copy(
            ax.transData, fig=plt.gcf(),
            x=offset_x, y=offset_y, units='points'
        )
        ax.text(xi, yi, label, fontsize=6, ha='left', va='bottom', transform=trans_offset)
             

    # plt.title(title)
    plt.xlabel("Energía (eV)")
    plt.ylabel("Conectividad")
    plt.xscale(xscale)
    plt.yscale(yscale)
    plt.grid(True)
    plt.tight_layout()
    # plt.show()
    # Guardar grafica
    plt.savefig(f"{title.replace(' ', '_')}.png", bbox_inches='tight', dpi=300)

if __name__ == "__main__":
    base_path = os.path.dirname(os.path.abspath(__file__))

    # Cargar arrays alineados y ordenados
    ordered_species = np.load(os.path.join(base_path, "ordered_species.npy"), allow_pickle=True)
    reactantsDegree_sorted = np.load(os.path.join(base_path, "reactantsDegree_sorted.npy"))
    productsDegree_sorted = np.load(os.path.join(base_path, "productsDegree_sorted.npy"))

    # Cargar energías desde archivo local
    energies = load_energies()

    # # Plot para reactivos
    # plot_connectivity_vs_energy(
    #     energies, reactantsDegree_sorted, ordered_species, title="Reactants Connectivity vs Energy (Lin -Log)", xscale='log', yscale='linear'
    # )

    # # Plot para productos
    # plot_connectivity_vs_energy(
    #     energies, productsDegree_sorted, ordered_species, title="Products Connectivity vs Energy (Lin - Log)", xscale='log', yscale='linear'
    # )

    # Plot para reactivos log - log
    # plot_connectivity_vs_energy(
        # energies, reactantsDegree_sorted, ordered_species, title="Reactants Connectivity vs Energy (Log-Log)", xscale='log', yscale='log'
    # )

    # Plot para productos log - log
    # plot_connectivity_vs_energy(
       # energies, productsDegree_sorted, ordered_species, title="Products Connectivity vs Energy (Log-Log)", xscale='log', yscale='log'
   # )

    # Plot para reactivos con eje x xspace
    plot_connectivity_vs_energy_xspace(
        energies, reactantsDegree_sorted, ordered_species, title="Reactants Connectivity vs Energy (Xspace)", xscale='linear', yscale='linear'
    )

    # Plot para productos con eje x xspace
    plot_connectivity_vs_energy_xspace(
        energies, productsDegree_sorted, ordered_species, title="Products Connectivity vs Energy (Xspace)", xscale='linear', yscale='linear'
    )
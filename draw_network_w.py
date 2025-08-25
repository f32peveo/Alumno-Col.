import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
import loki as lk
import community as community_louvain

# Cargar energías y ordenar especies
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

def construir_red_bipartita(uniqueSpecies, reactions):
    B = nx.Graph()

    B.add_nodes_from(uniqueSpecies, bipartite='species', tipo='especie')

    for i, reaction in enumerate(reactions):
        reaction_id = f"R{i}"
        B.add_node(reaction_id, bipartite='reaction', tipo='reaccion')

        for s, c in zip(reaction['lhsSpecies'], reaction['lhsStoichiometricCoeffs']):
            if s in uniqueSpecies:
                B.add_edge(reaction_id, s, weight=c)
        for s, c in zip(reaction['rhsSpecies'], reaction['rhsStoichiometricCoeffs']):
            if s in uniqueSpecies:
                B.add_edge(reaction_id, s, weight=c)

    return B

def dibujar_red_bipartita(B, species_nodes, reaction_nodes, titulo="Red Bipartita", height=40):
    plt.figure(figsize=(18, height))
    pos = nx.bipartite_layout(B, species_nodes)
    nx.draw(B, pos, with_labels=True, node_size=50,
            node_color=['lightblue' if n in species_nodes else 'lightgreen' for n in B.nodes()],
            edge_color='gray', font_size=8, font_color='black')
    plt.title(titulo)
    plt.axis('off')
    plt.tight_layout()
    # plt.show()
    plt.gcf().set_size_inches(6, 4)           # tamaño físico razonable
    plt.savefig("red_bipartita.pdf", bbox_inches="tight")
    plt.close()


def main():
    uniqueSpecies, reactions = lk.parseChemFile("helium.chem")
    species_nodes = [s for s in uniqueSpecies if s != 'e']
    reaction_nodes = [f"R{i}" for i in range(len(reactions))]
    B = construir_red_bipartita(uniqueSpecies, reactions)

    # Dibujar toda la red bipartita en una sola imagen
    altura = max(10, len(species_nodes) * 0.5 + len(reaction_nodes) * 0.2)
    dibujar_red_bipartita(B, species_nodes, reaction_nodes, titulo="Red Bipartita Completa", height=altura)

    # Proyectar la red bipartita en una red de especies
    G = nx.bipartite.projected_graph(B, species_nodes, multigraph=False)
    plt.figure(figsize=(12, 8))
    pos = nx.spring_layout(G, seed=42)
    nx.draw(G, pos, with_labels=True, node_color='lightblue', edge_color='gray', font_size=8, font_color='black')
    # plt.title("Red Proyectada de Especies")
    plt.show()
    # plt.savefig("red_proyectada_especies.png", bbox_inches='tight', dpi=300)

    # Aplicar la detección de comunidades
    partition = community_louvain.best_partition(G, weight='weight', resolution=1.0)
    # Dibujar la red con las comunidades
    plt.figure(figsize=(12, 8))
    pos = nx.spring_layout(G, seed=42)
    cmap = plt.cm.get_cmap('viridis', max(partition.values()) + 1)
    node_colors = [cmap(partition[n]) if n in partition else (0.5, 0.5, 0.5, 1.0) for n in G.nodes()]  # Gray for missing

    nx.draw(G, pos, with_labels=True, node_size=50, node_color=node_colors, edge_color='gray', font_size=8)

    # nx.draw(G, pos, with_labels=True, node_size=50, node_color=[cmap(partition[n]) for n in B.nodes()], edge_color='gray', font_size=8)
    # plt.title("Red Bipartita con Comunidades")
    plt.axis('off')
    plt.tight_layout()
    # plt.show()
    # Guardar imagen
    plt.savefig("red_bipartita_comunidades.png")

    # Calcular modularidad
    modularidad = community_louvain.modularity(partition, G)
    print(f"Modularidad de la red bipartita: {modularidad:.4f}")


if __name__ == "__main__":
    main()
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
import loki as lk

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

def dibujar_red_bipartita(B, species_nodes, reaction_nodes, titulo="Red Bipartita"):
    plt.figure(figsize=(18, 10))
    pos = {}
    pos.update(nx.bipartite_layout(B, species_nodes, align='vertical', scale=100))

    sizes_species = [max(10, sum(B[u][v]['weight'] for u in B.neighbors(v)) * 5) for v in species_nodes]
    sizes_reactions = [20 for _ in reaction_nodes]

    nx.draw_networkx_nodes(B, pos, nodelist=species_nodes, node_size=sizes_species, node_color='skyblue', label='Especies')
    nx.draw_networkx_nodes(B, pos, nodelist=reaction_nodes, node_size=sizes_reactions, node_color='salmon', label='Reacciones')
    nx.draw_networkx_edges(B, pos, width=1, alpha=0.6)
    nx.draw_networkx_labels(B, pos, font_size=7)

    plt.title(titulo)
    plt.axis('off')
    plt.legend()
    plt.tight_layout()
    plt.show()

def main():
    uniqueSpecies, reactions = lk.parseChemFile("helium.chem")
    species_nodes = [s for s in uniqueSpecies if s != 'e']
    reaction_nodes = [f"R{i}" for i in range(len(reactions))]
    B = construir_red_bipartita(uniqueSpecies, reactions)

    # Mostrar por bloques de 100 reacciones
    bloque = 100
    for start in range(0, len(reaction_nodes), bloque):
        end = min(start + bloque, len(reaction_nodes))
        sub_reactions = reaction_nodes[start:end]
        # Especies conectadas a estas reacciones
        sub_species = set()
        for r in sub_reactions:
            sub_species.update(B.neighbors(r))
        sub_species = [s for s in sub_species if s != 'e'] 
        # Subgrafo
        subB = B.subgraph(sub_species + sub_reactions)
        dibujar_red_bipartita(subB, sub_species, sub_reactions, titulo=f"Reacciones {start+1}-{end}")

if __name__ == "__main__":
    main()
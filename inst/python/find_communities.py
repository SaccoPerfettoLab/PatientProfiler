import pandas as pd
import networkx as nx
import os

def check_sign(src, target):
    contains_inhibited_1 = 'INHIBITED' in src
    contains_inhibited_2 = 'INHIBITED' in target
    if contains_inhibited_1 == contains_inhibited_2:
        return 1
    else:
        return -1

def check_reltype(sign):
    if sign==-1:
        return "inhibition"
    else:
        return "activation"

def check_node_sign(node):
    if "INHIBITED" in node:
        return -1
    else:
        return 1

def extract_edges(patients_list,nodes_dir,edges_dir):
    edges_hist={}
    edges_list=[]
    mechanism_map={}
    mf_map={}

    print("Processing patients data...")
    for i in range(len(patients_list)):
        nodes_df = pd.read_excel(os.path.join(nodes_dir, f'nodes_{patients_list[i]}.xlsx'), engine='openpyxl')
        edges_df = pd.read_excel(os.path.join(edges_dir,  f'edges_{patients_list[i]}.xlsx'), engine='openpyxl')
        active_genes = list(nodes_df[nodes_df['carnival_activity'] == 100]["gene_name"])
        filtered_edges = edges_df[edges_df['source'].isin(active_genes)]

        for index, row in nodes_df.iterrows():
            gene_name = row['gene_name']
            sign = int(row["carnival_activity"])
            if sign <=0:
                gene_name+="_INHIBITED"
            mf=row["mf"]
            if not pd.isna(mf):
                mf=mf.replace("*","")
                if gene_name not in mf_map.keys() or len(mf_map[gene_name])<= len(mf):
                    mf_map[gene_name] = mf

        for index, row in filtered_edges.iterrows():
            source_gene = row['source']
            target_gene = row['target']
            if source_gene not in active_genes:
                source_gene = source_gene + "_INHIBITED"
            if target_gene not in active_genes:
                target_gene = target_gene + "_INHIBITED"
            e = (source_gene, target_gene)
            mechanism=row["mechanism"]

            if not pd.isna(mechanism):
                if e not in mechanism_map.keys():
                    mechanism_map[e] = mechanism

            edges_list.append((e, i))
            if e not in edges_hist.keys():
                edges_hist[e] = 1
            else:
                edges_hist[e] += 1

    return edges_hist, edges_list, mechanism_map, mf_map

def write_to_txt(communities, output_dir, patients_list, id_2_edge, mechanism_map, mf_map):
    for i, community in enumerate(communities):

        patients_community = [patients_list[node] for node in community if node < len(patients_list)]

        interactome = [
            (id_2_edge[node][0], id_2_edge[node][1], check_sign(id_2_edge[node][0], id_2_edge[node][1]),
             "NaN" if id_2_edge[node] not in mechanism_map.keys() else mechanism_map[id_2_edge[node]])
            for node in community if node >= len(patients_list)
        ]

        # Save community to directory
        community_dir = os.path.join(output_dir, f"community_{i + 1}")
        os.makedirs(community_dir, exist_ok=True)

        # Save patients to file
        with open(os.path.join(community_dir, "patients.txt"), "w") as f:
            f.write("\n".join(patients_community))

        # Save interactome edges to file
        with open(os.path.join(community_dir, "interactome.txt"), "w") as f:
            for edge in interactome:

                f.write(f"{edge[0]}\t{edge[1]}\t{edge[2]}\t{edge[3]}\t{edge[4]}\n")

def write_output(communities,output_dir,patients_list,id_2_edge,mechanism_map,mf_map):

    for i, community in enumerate(communities):
        patients_community = [patients_list[node] for node in community if node < len(patients_list)]

        interactome = [
            (id_2_edge[node][0], id_2_edge[node][1], check_sign(id_2_edge[node][0], id_2_edge[node][1]),
             "NaN" if id_2_edge[node] not in mechanism_map.keys() else mechanism_map[id_2_edge[node]])
            for node in community if node >= len(patients_list)
        ]

        community_dir = os.path.join(output_dir, f"community_{i + 1}")
        os.makedirs(community_dir, exist_ok=True)

        with open(os.path.join(community_dir, "patients.txt"), "w") as f:
            f.write("\n".join(patients_community))

        node_df = pd.DataFrame({
            "gene_name": [edge[0] for edge in interactome] + [edge[1] for edge in interactome],
            "sign": [check_node_sign(edge[0]) for edge in interactome] + [check_node_sign(edge[1]) for edge in interactome],
            "mf": ["NaN" if edge[0] not in mf_map.keys() else mf_map[edge[0]] for edge in interactome] + ["NaN" if edge[1] not in mf_map.keys() else mf_map[edge[1]] for edge in interactome]
        }).drop_duplicates()
        node_df.to_csv(os.path.join(community_dir, f"nodes_{i + 1}.csv"), index=False)

        edge_df = pd.DataFrame({
            "source": [edge[0] for edge in interactome],
            "target": [edge[1] for edge in interactome],
            "sign": [edge[2] for edge in interactome],
            "mechanism": [edge[3] for edge in interactome]
        })
        edge_df.to_csv(os.path.join(community_dir, f"edges_{i + 1}.csv"), index=False)

        '''
        # Salva file in formato .sif
        sif_path = os.path.join(community_dir, f"interactome_{i + 1}.sif")
        with open(sif_path, "w") as sif_file:
            for edge in interactome:
                sif_file.write(f"{edge[0]}\t{check_reltype(edge[2])}\t{edge[1]}\n")
                

        '''
        #pyreadr.write_rds(os.path.join(community_dir,f"interactome_{i+1}.rds"), edge_df)
        

def find_communities(path,t_lower=4,t_upper=30,output_dir="output_communities"):
    """
       Identifies communities in a bipartite graph created from biological data.

       This function processes patient data, extracts relevant edges and nodes, constructs a bipartite graph,
       and identifies patient communities and the corresponding representative interactions using the Louvain algorithm.
       The resulting communities are saved in various formats for further analysis.

       Args:
           path (str): The root directory containing the patient data, nodes, and edges files.
               The directory should contain subdirectories named 'nodes' and 'edges', and a patient list file.
           t_lower (int, optional): The lower threshold for filtering edges based on their frequency across patients.
               Default is 4.
           t_upper (int, optional): The upper threshold for filtering edges based on their frequency across patients.
               If not provided, no upper limit is applied, and all edges above `t_lower` will be included.
           output_dir (str, optional): The directory where the results (e.g., community files, graphs) will be saved.
               Default is 'output_communities'.

       Raises:
           FileNotFoundError: If the required files or directories are missing in the specified `path`.
           ValueError: If the edge or node extraction process encounters unexpected data formats or missing values.

       Returns:
           None: All outputs are saved directly in the specified `output_dir`.

       Side Effects:
           - Creates output subdirectories for each identified community.
           - Writes patient lists, interactome edges, and node/edge details to CSV and .sif files.
           - Saves interactome edge data in .rds format.

       Notes:
           - The function assumes the Louvain algorithm is used for community detection.
           - Requires pandas, networkx, pyreadr, and os libraries for proper functionality.

       Example:
           >>> find_communities("/path/to/data", t_lower=5, output_dir="results")
           Processing patient data...
           Extracting edges and constructing bipartite graph...
           Detected 15 communities.
           Results saved in 'results' directory.
       """
    nodes_dir = None
    edges_dir = None
    patients_file = None

    for root, dirs, files in os.walk(path):
        if len(dirs) == 0:
            continue
        for dir_name in dirs:
            if 'nodes' in dir_name:
                nodes_dir = os.path.join(root, dir_name)
            else:
                edges_dir = os.path.join(root, dir_name)

        patients_file = os.path.join(root, files[0])

    patients_list=pd.read_table(patients_file).iloc[:,0].tolist()
    edge_hist, edges_list, mechanism_map, mf_map = extract_edges(patients_list, nodes_dir, edges_dir)
    print("Extracting edges and constructing bipartite graph")
    edges_hist = {k: v for k, v in edge_hist.items() if v >= t_lower and v <= t_upper}

    edge_2_id = {}
    biadjacency_list = []
    counter=len(patients_list)

    for k, v in edges_hist.items():
        edge_2_id[k] = counter
        counter += 1

    for e in edges_list:
        if e[0] in edges_hist.keys():
            biadjacency_list.append((edge_2_id[e[0]], e[1]))



    bipartite_graph = nx.Graph()
    bipartite_graph.add_nodes_from(range(122, counter), bipartite=0)
    bipartite_graph.add_nodes_from(range(122), bipartite=1)
    bipartite_graph.add_edges_from(biadjacency_list)

    communities = nx.community.louvain_communities(bipartite_graph, seed=277, resolution=0.99)
    print(f"Detected {len(communities)} communities")

    # Map IDs back to edges
    id_2_edge = {v: k for k, v in edge_2_id.items()}

    # Prepare output directory
    os.makedirs(output_dir, exist_ok=True)

    write_output(communities, output_dir, patients_list, id_2_edge, mechanism_map, mf_map)

    print(f"Results saved in {output_dir} directory.")



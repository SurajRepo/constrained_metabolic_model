import numpy as np

def _map_expression_reactions(reactions, grRules, genes, values):
    ## Helper Function to map gene expression values to reactions based on GPR rules
    gene_reactions = {}

    gene_indexes = {gene: index for index, gene in enumerate(genes)}

    for reaction in reactions:
        gene_reactions[reaction] = np.ones(values.shape[1])

    for rxn, associated_genes in enumerate(grRules):
        if len(associated_genes) > 0:
            # print(f"Mapping expression for reaction: {reactions[rxn]} with GPR: {associated_genes}")
            gene_sets = associated_genes.split(' or ')
            max_gene_values = np.zeros(values.shape[1])

            for gene_item in gene_sets:
                if ' and ' in gene_item:
                    gene_and = gene_item.split(' and ')
                    gene_and_values = [values[gene_indexes[g.strip()], :] for g in gene_and if g.strip() in gene_indexes]
                    if gene_and_values:
                        max_gene_values = np.maximum(max_gene_values, np.min(gene_and_values))
                else:
                    if gene_item.strip() in gene_indexes:
                        max_gene_values = np.maximum(max_gene_values, values[gene_indexes[gene_item.strip()], :])

            gene_reactions[reactions[rxn]] = max_gene_values
        else:
            gene_reactions[reactions[rxn]] = np.ones(values.shape[1])
    return gene_reactions


def update_reaction_bounds(default_model, df, gamma=2):     
    model = default_model.copy()   
    genes_in_model = [gene.id for gene in model.genes]
    metabolic_genes = df[df.index.isin(genes_in_model)]
    metabolic_genes = metabolic_genes[metabolic_genes.values != 0]
    if metabolic_genes.shape[0] == 0:
        print("No metabolic genes found in the expression data.")
        return model
                    
    values = metabolic_genes.values.reshape(-1, 1)
    genes = metabolic_genes.index.tolist()
    for rxn in model.reactions:
        expr = _map_expression_reactions([rxn.id], [str(rxn.gpr)], genes, values)
        if expr[rxn.id][0] != 1.0 and expr[rxn.id][0] != 0.0:
            rxn.lower_bound = rxn.lower_bound * gamma * expr[rxn.id][0]
            rxn.upper_bound = rxn.upper_bound * gamma * expr[rxn.id][0]
            # print(f"Updated bounds for reaction {rxn.id}: Lower Bound = {rxn.lower_bound}, Upper Bound = {rxn.upper_bound}")
    return model
    
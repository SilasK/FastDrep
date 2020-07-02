import scipy.spatial as sp
import scipy.cluster.hierarchy as hc
import ete3


def linkage_to_newick(Z, labels):
    """
    Input :  Z = linkage matrix, labels = leaf labels
    Output:  Newick formatted tree string
    """
    tree = hc.to_tree(Z, False)

    def buildNewick(node, newick, parentdist, leaf_names):
        if node.is_leaf():
            return f"{leaf_names[node.id]}:{(parentdist - node.dist)/2}{newick}"
        else:
            if len(newick) > 0:
                newick = f"):{(parentdist - node.dist)/2}{newick}"
            else:
                newick = ");"
            newick = buildNewick(node.get_left(), newick, node.dist, leaf_names)
            newick = buildNewick(node.get_right(), f",{newick}", node.dist, leaf_names)
            newick = f"({newick}"
            return newick

    return buildNewick(tree, "", tree.dist, labels)

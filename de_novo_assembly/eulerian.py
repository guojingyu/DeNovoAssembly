"""


Author Jingyu Guo
"""
import networkx as nx
import logging
import datetime
import sys
import operator

# Configure logging
timestamp = datetime.datetime.now().strftime("%Y-%m-%d")
logging.basicConfig(filename=timestamp+"_denovo_assembly.log",
                    stream=sys.stdout,
                    level=logging.DEBUG)
stderrLogger=logging.StreamHandler()
stderrLogger.setFormatter(logging.Formatter(logging.BASIC_FORMAT))
logger = logging.getLogger(__name__)
logger.addHandler(stderrLogger)


def is_euler(graph):
    """
    Return true or false. Every node in graph must have equal in degree and out
    :param graph: graph as nx graph obj
    :return: true or false for graph is Eulerian
    """
    return has_euler_circuit(graph) or has_euler_path(graph)

def has_euler_circuit(graph):
    """
     Check if the graph contains odd number of degree nodes to has euler
     path
     :param graph: networkx graph obj
     :return: true or false
     """
    if not nx.is_strongly_connected(graph):
        return False
    for node in graph.nodes_iter():
        if graph.in_degree(node) != graph.out_degree(node):
            return False
    return True

def has_euler_path(graph):
    """
    Check if the graph contains odd number of degree nodes to has euler
    path. If true, also return start and end node; Otherwise, return false
    and none for start and end node
    :param graph: networkx graph obj
    :return: true or false, and start node and end node (depending on if
    true or false)
    """
    flag = False
    start = None
    end = None
    for node in graph.nodes_iter():
        node['end'] = False
        node['start'] = False
        if graph.degree(node) % 2 != 0:
            if graph.in_degree(node) - 1 == graph.out_degree(node):
                node['end'] = True
                end = node
            elif graph.in_degree(node) + 1 == graph.out_degree(node):
                node['start'] = True
                start = node
            flag = True
    return flag, start, end


def eulerian_random_walk(DBG):
    """
    the method to reconstruct the genome

    Euler walk methods, while left out from implementation, can be
    described as a restarted random walking -- as an analogy here,
    Euler walks through the graph and try to pass through every
    edge, and when he is stuck, restart on an explored node with
    unexplored edges and follow the unexplored edges (repeat the explored
    edges first if needed)... until there is no more unexplored edges.

    :param DBG: DeBruijnGraph obj
    :return: dictionary containing a list of edges that represents the
    assembly of the genome, each is a contig
    """

    def find_eulerian_path(graph, start):
        copy_graph = graph.__class__(graph)  # copy graph structure
        degree = copy_graph.in_degree
        edges = copy_graph.in_edges_iter
        get_node = operator.itemgetter(0)

        stack = [start]
        last_node = None
        while stack:
            current_node = stack[-1]
            # if there is no more edge to explore
            if degree(current_node) == 0:
                if last_node is not None:
                    yield (last_node, current_node)
                last_node = current_node
                stack.pop() # remove [-1]
            else: # move on and remove the edge through
                random_edge = next(edges(current_node))
                stack.append(get_node(random_edge))
                copy_graph.remove_edge(*random_edge)

    assembly = []
    for subg in DBG.subgraphs:
        # nx eulerian assumes strong connection for directed graph
        if nx.is_eulerian(subg):
            subg['eulerian'] = True
            logger.info("Eulerian : Eulerian circuit found for the "
                        "subgraph of De Bruijn Graph built from input "
                        "sequence : ")
            logger.info(subg.nodes())
            a_euler_circuit = list(nx.eulerian_circuit(subg))
            # assembly removed the last edge to source
            assembly.append(make_contig_from_path(a_euler_circuit[:-1]))
        elif (not nx.is_eulerian(subg)) and subg.euler_path_exists:
                # TODO: add a method to find the Euler path
                subg['eulerian'] = True
                logger.info("Eulerian : Eulerian path found for the "
                            "subgraph of De Bruijn Graph built from input "
                            "sequence : ")
                logger.info(subg.nodes())

                a_euler_path = list(find_eulerian_path(subg,
                                                       subg.euler_path_start))
                assert a_euler_path[-1][1] in subg.euler_path_end
                assembly.append(make_contig_from_path(a_euler_path))
        else:
            subg['eulerian'] = False
            logger.warn("Eulerian : no Eulerian feature found for the "
                        "subgraph of De Bruijn Graph built from input "
                        "sequence : ")
            logger.warn(subg.nodes())
            # TODO consider to export the graph here
    return assembly

def make_contig_from_path(path):
    """
    Assemble the sequence from list of (l_node, r_node)
    :param path: a list of kmers
    :return: a string of sequence
    """
    return reduce(lambda x, y: x + y[-1], [l + r[-1] for l, r in path])
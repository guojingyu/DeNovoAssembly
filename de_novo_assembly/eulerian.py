#!/usr/bin/env python

"""
Eulerian features and functions of De Bruijn Graph

Author Jingyu Guo
"""
import networkx as nx
import matplotlib.pyplot as plt
import logging
import datetime
import sys
import operator
from utils import clock_now

# Configure logging
timestamp = datetime.datetime.now().strftime("%Y-%m-%d")
logging.basicConfig(filename="../" + timestamp+"_denovo_assembly.log",
                    stream=sys.stdout,
                    level=logging.DEBUG)
stderrLogger=logging.StreamHandler()
stderrLogger.setFormatter(logging.Formatter(logging.BASIC_FORMAT))
logger = logging.getLogger(__name__)
logger.addHandler(stderrLogger)


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
    for n in graph.nodes_iter():
        graph.node[n]['end'] = False
        graph.node[n]['start'] = False
        if graph.degree(n) % 2 != 0:
            if graph.out_degree(n) == graph.in_degree(n) - 1:
                graph.node[n]['end'] = True
                end = n
            elif graph.in_degree(n) == graph.out_degree(n) - 1:
                graph.node[n]['start'] = True
                start = n
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
    assembly = []
    for subg in DBG.subgraphs:
        # nx eulerian assumes strong connection for directed graph
        if subg.graph['euler_circuit']:
            a_euler_circuit = list(nx.eulerian_circuit(subg))
            # assembly removed the last edge to source
            assembly.append(make_contig_from_path(a_euler_circuit[:-1]))
            logger.info(clock_now() + " Eulerian : Eulerian circuit found for "
                                     "the subgraph of De Bruijn Graph built "
                                      "from the input sequence : ")
            logger.info(a_euler_circuit)
        elif subg.graph['euler_path']:
            # add a temp edge here
            subg.add_edge(subg.graph['euler_path_end'], subg.graph[
                'euler_path_start'], seq="[temp_edge]")
            a_euler_path = list(find_eulerian_path(subg,subg.graph[
                'euler_path_start']))[:-1]
            print a_euler_path
            assembly.append(make_contig_from_path(a_euler_path))
            logger.info(clock_now() + " Eulerian : Eulerian path found "
                                          "for the subgraph of De Bruijn Graph built from the input "
                            "sequence : ")
            logger.info(a_euler_path)
        else:
            subg.graph['eulerian'] = False
            logger.warn(clock_now() + " Eulerian : no Eulerian feature found "
                                      "for the subgraph of De Bruijn Graph built from input sequence : ")
            logger.warn(subg.nodes())
    return assembly


def find_eulerian_path(graph, start):
    """
    a method to find the Eulerian path with a starting node
    :param graph: potentially a subgraph that satisify the
    :param start: starting node
    :return:
    """
    if not graph.graph['euler_path']:
        raise RuntimeError("Eulerian : call Eulerian path method on "
                           "graph/subgraph having no Eulerian path: " +
                           graph.nodes())

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
            stack.pop()  # remove [-1]
        else:  # move on and remove the edge through
            random_edge = next(edges(current_node))
            stack.append(get_node(random_edge))
            copy_graph.remove_edge(*random_edge)


def make_contig_from_path(path):
    """
    Assemble the sequence from list of (l_node, r_node)
    :param path: a list of kmers
    :return: a string of sequence
    """
    return reduce(lambda x,y: x+y[-1],[l + r[-1] for l, r in path])
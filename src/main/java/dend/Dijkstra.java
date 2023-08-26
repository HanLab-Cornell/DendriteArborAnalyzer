package dend;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Set;

/**
 * This class implements Dijkstra's algorithm to find the shortest path from
 * each node to the
 * source vertex.
 *
 * @author Diana Bank
 * @author Inle Bush
 * @author Elena Zhong
 * @author Chun Han Lab, Cornell University
 */

public class Dijkstra {

	/* Vertices marked as visited */
	private Set<Vertex> settledNodes = new HashSet<>();
	/* Vertices marked as unvisited */
	private Set<Vertex> unSettledNodes = new HashSet<>();
	/* List of previous vertices */
	private Map<Vertex, Vertex> predecessors = new HashMap<>();
	/* Map of distance from each Vertex to a specified "source" Vertex */
	private Map<Vertex, Double> distance = new HashMap<>();

	/**
	 * Initializes and runs the main algorithm to find shortest paths between all
	 * nodes
	 * and the source Vertex, saves result.
	 *
	 * @param source the Vertex to use as the source
	 */
	public Dijkstra(Vertex source) {
		/* Vertices marked as visited */
		settledNodes = new HashSet<>();
		/* Vertices marked as unvisited */
		unSettledNodes = new HashSet<>();
		/* List of previous vertices */
		predecessors = new HashMap<>();
		/* Map of distance from each Vertex to a specified "source" Vertex */
		distance = new HashMap<>();

		distance.put(source, 0.0);
		unSettledNodes.add(source);

		while (unSettledNodes.size() > 0) {
			Vertex node = getMinimum(unSettledNodes);
			settledNodes.add(node);
			unSettledNodes.remove(node);
			findMinimalDistances2(node);
		}
	}

	/**
	 * Getter method to return the map of distances
	 *
	 * @return distances
	 */
	public Map<Vertex, Double> getDistances() {
		return distance;
	}

	/**
	 * Given a Vertex node, this method calculates the shortest distances between
	 * the node Vertex and
	 * its neighbors
	 *
	 * @param node
	 */
	private void findMinimalDistances2(Vertex node) {
		List<Vertex> adjacentNodes = getNeighbors(node);
		for (Vertex target : adjacentNodes) {
			if (target != null) {
				Double dw = getShortestDistance(node) +
						getDistance(node, target);
				if (!settledNodes.contains(target) &&
						!unSettledNodes.contains(target)) {
					distance.put(target, dw);
					predecessors.put(target, node);
					unSettledNodes.add(target);
				} else if (distance.containsKey(target) &&
						distance.get(target) > dw) {
					distance.put(target, dw);
					predecessors.put(target, node);
				}
			}
		}
	}

	/**
	 * Get the distance between two vertices on opposite ends of an edge
	 *
	 * @param node
	 * @param target
	 * @return edge length
	 */
	private double getDistance(Vertex node, Vertex target) {
		for (Edge e : node.getBranches()) {
			if (e.getV2().equals(target) || e.getV1().equals(target)) {
				return e.getLength();
			}
		}
		throw new RuntimeException("Should not happen");
	}

	/**
	 * Gets a list of Vertices immediately connected to the node Vertex by one edge
	 *
	 * @param node
	 * @return neighbors - List of neighbor vertices
	 */
	private ArrayList<Vertex> getNeighbors(Vertex node) {

		ArrayList<Vertex> neighbors = new ArrayList<>();
		ArrayList<Edge> edges = node.getBranches();
		for (Edge e : edges) {
			if (e.getV1().equals(node))
				neighbors.add(e.getV2());
			else if (e.getV2().equals(node))
				neighbors.add(e.getV1());
		}
		return neighbors;
	}

	/**
	 * Returns the Vertex with the shortest distance from a set of vertices
	 *
	 * @param vertexes
	 * @return
	 */
	private Vertex getMinimum(Set<Vertex> vertices) {
		Vertex minimum = null;
		for (Vertex vertex : vertices) {
			if (minimum == null) {
				minimum = vertex;
			} else {
				if (getShortestDistance(vertex) < getShortestDistance(
						minimum)) {
					minimum = vertex;
				}
			}
		}
		return minimum;
	}

	/**
	 * Returns the calculated shortest distance for a destination vertex from the
	 * map of distances
	 *
	 * @param destination destination vertex
	 * @return shortest distance
	 */
	public double getShortestDistance(Vertex destination) {
		Double d = distance.get(destination);
		if (d == null) {
			return Integer.MAX_VALUE;
		} else {
			return d;
		}
	}

	/**
	 * This method returns the path from the source node to the selected target node
	 * and returns NULL
	 * if no path exists
	 *
	 * @param target target node
	 * @return path
	 */
	public LinkedList<Vertex> getPath(Vertex target) {
		LinkedList<Vertex> path = new LinkedList<>();
		Vertex step = target;
		// check if a path exists
		if (predecessors.get(step) == null) {
			return null;
		}
		path.add(step);
		while (predecessors.get(step) != null) {
			step = predecessors.get(step);
			path.add(step);
		}
		// Put it into the correct order
		Collections.reverse(path);
		return path;
	}

}

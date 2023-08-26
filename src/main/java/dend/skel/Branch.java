package dend.skel;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.Set;

import dend.Edge;
import dend.Point;
import dend.Vertex;

/**
 * This class represents a branch, which is essentially an edge with specified
 * path distance,
 * Euclidean distance, and order.
 *
 * @author Diana Bank
 * @author Inle Bush
 * @author Elena Zhong
 * @author Chun Han Lab, Cornell University
 */

public class Branch implements Comparable<Branch> {

	private Edge edge = null;
	/** vertex on graph */
	private Vertex vertex = null;
	/** path distances */
	private double pathDist = 0.0;
	/** Euclidean distance */
	private double dist = 0.0;
	/** Strahler order */
	private int order = 0;
	/** branch length */
	private double length = 0;
	/** visited */
	private boolean visited = false;
	private int predorder = 10000;
	private Set<Point> points = null;

	/**
	 * Create a branch from an edge.
	 *
	 * @param e edge
	 */
	public Branch(Edge e) {
		edge = e;
		vertex = e.getV1();
		length = e.getLength();
		points = new HashSet<>(e.getSlabs());
	}

	/**
	 * Create a branch from an edge with specified path distance, Euclidean
	 * distance, and order.
	 *
	 * @param e        edge
	 * @param pathDist path distance from a source vertex
	 * @param dist     Euclidean distance from a source vertex
	 * @param order    Strahler order
	 */
	public Branch(Edge e, double pathDist, double dist, int order) {
		edge = e;
		vertex = e.getV1();
		length = e.getLength();
		this.pathDist = pathDist;
		this.dist = dist;
		this.order = order;
		points = new HashSet<>(e.getSlabs());
	}

	/**
	 * Get the vertex belonging to the branch.
	 *
	 * @return vertex of the branch
	 */
	public Vertex getVertex() {
		return vertex;
	}

	/**
	 * Get the path distance to the branch.
	 *
	 * @return path distance to the branch
	 */
	public double getPathDist() {
		return pathDist;
	}

	/**
	 * Get the Euclidean distance to the branch.
	 *
	 * @return Euclidean distance to the branch
	 */
	public double getDist() {
		return dist;
	}

	/**
	 * Get order of the branch.
	 *
	 * @return Strahler order of the branch
	 */
	public int getOrder() {
		return order;
	}

	/**
	 * Get length of the branch.
	 *
	 * @return length of the branch
	 */
	public double getLength() {
		return length;
	}

	/**
	 * Get the edge belonging to the branch.
	 *
	 * @return edge of the branch
	 */
	public Edge getEdge() {
		return edge;
	}

	/**
	 * Set vertex.
	 *
	 * @param vertex vertex belonging to the branch
	 */
	public void setVertex(Vertex vertex) {
		this.vertex = vertex;
	}

	/**
	 * Set path distance.
	 *
	 * @param pathDist path distance from a source vertex
	 */
	public void setPathDist(double pathDist) {
		this.pathDist = pathDist;
	}

	/**
	 * Set Euclidean distance.
	 *
	 * @param dist Euclidean distance from a source vertex
	 */
	public void setDist(double dist) {
		this.dist = dist;
	}

	/**
	 * Set Strahler order.
	 *
	 * @param order Strahler order
	 */
	public void setOrder(int order) {
		this.order = order;
	}

	/**
	 * Set length of branch.
	 *
	 * @param length length of branch
	 */
	public void setLength(int length) {
		this.length = length;
	}

	/**
	 * Set edge.
	 *
	 * @param edge edge belonging to the branch
	 */
	public void setEdge(Edge edge) {
		this.edge = edge;
	}

	/**
	 * Returns a string representation of the branch.
	 *
	 * @return a string representation of the branch
	 */
	@Override
	public String toString() {
		return "Order: " + order + " Branch Length: " + length +
				" Path Dist: " + pathDist + " Euc Dist: " + dist;
	}

	/**
	 * Compares this branch with another branch for order. Returns a negative
	 * integer, zero, or a
	 * positive integer as the order of this branch is less than, equal to, or
	 * greater than that of the
	 * specified branch.
	 *
	 * @param other the branch to be compared
	 * @return a negative integer, zero, or a positive integer as the order of this
	 *         branch less than,
	 *         equal to, or greater than that of the specified branch.
	 */
	@Override
	public int compareTo(Branch other) {
		Integer thisOrder = order;
		return thisOrder.compareTo(other.order);
	}

	/**
	 * Check visit status.
	 *
	 * @return true if the branch was already visited
	 */
	public boolean isVisited() {
		return visited;
	}

	/**
	 * Set branch as visited or not.
	 *
	 * @param visited boolean
	 */
	public void setVisited(boolean visited) {
		this.visited = visited;
	}

	/** @return the predorder */
	public int getPredorder() {
		return predorder;
	}

	/** @param predorder the predorder to set */
	public void setPredorder(int predorder) {
		this.predorder = predorder;
	}

	/** @return the points */
	public Set<Point> getPoints() {
		return points;
	}

	/** @param points the points to set */
	public void setPoints(Set<Point> points) {
		this.points = points;
	}

	public void addPoint(ArrayList<Point> ps) {
		points.addAll(ps);
	}

	public boolean isNeighbor(Branch bb, int i) {
		for (Point pp : points) {
			for (Point pp2 : bb.getPoints()) {
				if (pp2.isNeighbor(pp, i)) {
					return true;
				}
			}
		}
		return false;
	}

}

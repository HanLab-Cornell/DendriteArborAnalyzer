/*
 * #%L
 * Modified based on Graph.java of the ImageJ AnalyzeSkeleton_ plugin by Ignacio Arganda-Carreras.
 *
 * Major Modifications:
 * Added functions: clone(), cutEdge(), getLargeEdges().
 *
 * Diana Bank, Inle Bush, and Elena Zhong, Chun Han Lab, Cornell University
 * #L%
 */
package dend;

import static java.util.function.Function.identity;
import static java.util.stream.Collectors.toMap;

import java.util.ArrayList;
import java.util.Map;
import java.util.Stack;
import java.util.function.Function;

/**
 * This class represents an undirected graph to allow visiting the skeleton in
 * an efficient way.
 * Modified based on {@code Graph.java} of the ImageJ {@code AnalyzeSkeleton_}
 * plugin by Ignacio
 * Arganda-Carreras.
 * <p>
 * For more information about the original plugin, visit the AnalyzeSkeleton
 * home page:
 * <A target="_blank" href=
 * "http://fiji.sc/AnalyzeSkeleton">http://fiji.sc/AnalyzeSkeleton</A>
 *
 * @author Diana Bank
 * @author Inle Bush
 * @author Elena Zhong
 * @author Chun Han Lab, Cornell University
 */

public class Graph {
	/** list of edges */
	private ArrayList<Edge> edges = null;
	/** list of vertices */
	private ArrayList<Vertex> vertices = null;
	/** root vertex */
	private Vertex root = null;

	// --------------------------------------------------------------------------
	/** Empty constructor. */
	public Graph() {
		edges = new ArrayList<>();
		vertices = new ArrayList<>();
	}

	// --------------------------------------------------------------------------
	/**
	 * Add edge to the graph.
	 *
	 * @param e edge to be added
	 * @return false if the edge could not be added, true otherwise
	 */
	public boolean addEdge(Edge e) {
		if (edges.contains(e))
			return false;
		else {
			// Set vertices from e as neighbors (undirected graph)
			e.getV1().setBranch(e);
			if (!e.getV1().equals(e.getV2()))
				e.getV2().setBranch(e);
			// Add edge to the list of edges in the graph
			edges.add(e);
			return true;
		}
	}

	// --------------------------------------------------------------------------
	/**
	 * Add vertex to the graph.
	 *
	 * @param v vertex to be added
	 * @return false if the vertex could not be added, true otherwise
	 */
	public boolean addVertex(Vertex v) {
		if (vertices.contains(v))
			return false;
		else {
			vertices.add(v);
			return true;
		}
	}
	// --------------------------------------------------------------------------

	/**
	 * Get list of vertices in the graph.
	 *
	 * @return list of vertices in the graph
	 */
	public ArrayList<Vertex> getVertices() {
		return vertices;
	}

	// --------------------------------------------------------------------------
	/**
	 * Get list of edges in the graph.
	 *
	 * @return list of edges in the graph
	 */
	public ArrayList<Edge> getEdges() {
		return edges;
	}

	// --------------------------------------------------------------------------
	/** Set root vertex. */
	void setRoot(Vertex v) {
		root = v;
	}

	// --------------------------------------------------------------------------
	/**
	 * Get root vertex.
	 *
	 * @return root vertex of the graph
	 */
	public Vertex getRoot() {
		return root;
	}

	// --------------------------------------------------------------------------
	/**
	 * Depth first search method to detect cycles in the graph.
	 *
	 * @return list of BACK edges
	 */
	ArrayList<Edge> depthFirstSearch() {
		ArrayList<Edge> backEdges = new ArrayList<>();

		// Create empty stack
		Stack<Vertex> stack = new Stack<>();

		// Mark all vertices as non-visited
		for (final Vertex v : vertices)
			v.setVisited(false);

		// Push the root into the stack
		stack.push(root);

		int visitOrder = 0;

		while (!stack.empty()) {
			Vertex u = stack.pop();

			if (!u.isVisited()) {
				// IJ.log(" Visiting vertex " + u.getPoints().get(0));

				// If the vertex has not been visited yet, then
				// the edge from the predecessor to this vertex
				// is mark as TREE
				if (u.getPredecessor() != null)
					u.getPredecessor().setType(Edge.TREE);

				// mark as visited
				u.setVisited(true, visitOrder++);

				for (final Edge e : u.getBranches()) {
					// For the undefined branches:
					// We push the unvisited vertices in the stack,
					// and mark the edge to the others as BACK
					if (e.getType() == Edge.UNDEFINED) {
						final Vertex ov = e.getOppositeVertex(u);
						if (!ov.isVisited()) {
							stack.push(ov);
							ov.setPredecessor(e);
						} else {
							e.setType(Edge.BACK);
							backEdges.add(e);
						}

					}
				}
			}
		}

		return backEdges;

	} // end method depthFirstSearch

	// --------------------------------------------------------------------------
	/**
	 * Clone original graph.
	 *
	 * @return cloned graph
	 */
	@Override
	public Graph clone() {
		final Graph clone = new Graph();
		final Map<Vertex, Vertex> vertexMap = vertices.stream().collect(toMap(
				identity(), Vertex::cloneUnconnected));
		final Function<Edge, Edge> cloner = e -> e.clone(
				vertexMap.get(e.getV1()),
				vertexMap.get(e.getV2()));
		final Map<Edge, Edge> edgeMap = edges.stream().collect(toMap(identity(),
				cloner));
		vertices.forEach(v -> vertexMap.get(v).setPredecessor(edgeMap.get(v
				.getPredecessor())));
		// Iterate in the order of the originals to preserve order in the cloned
		// graph (makes testing easier)
		edges.stream().map(edgeMap::get).forEach(clone::addEdge);
		vertices.stream().map(vertexMap::get).forEach(clone::addVertex);
		clone.setRoot(vertexMap.get(root));
		return clone;
	} // end method clone

	// --------------------------------------------------------------------------
	/**
	 * Cut edge e at point p.
	 *
	 * @param e edge to be cut
	 * @param p point belongs to slabs of edge e
	 * @return new vertex at point p
	 */
	public Vertex cutEdge(Edge e, Point p) {
		edges.remove(e);

		// Truncates list of slabs at point p
		ArrayList<Point> slabs1 = new ArrayList<>(e.getSlabs().subList(0,
				e.getSlabs().indexOf(p)));
		ArrayList<Point> slabs2 = new ArrayList<>(e.getSlabs()
				.subList(e.getSlabs().indexOf(p) + 1, e.getSlabs().size()));

		// Finds vertices connected to new lists of slabs
		Vertex v1 = e.getV2();
		for (Point p0 : slabs1) {
			for (Point p1 : e.getV1().getPoints()) {
				if (p1.infDist(p0) < 1) {
					v1 = e.getV1();
					break;
				}
			}
		}

		// Constructs new vertex at point p
		Vertex v0 = new Vertex();
		v0.addPoint(p);

		double lz = e.getLength() / e.getSlabs().size();

		// Constructs new edges
		Edge e1 = new Edge(v0, v1, slabs1, slabs1.size() * lz);
		Edge e2 = new Edge(v0, e.getOppositeVertex(v1), slabs2,
				slabs2.size() * lz);

		addVertex(v0);
		addEdge(e1);
		addEdge(e2);

		return v0;
	}// end method cutEdge

	// --------------------------------------------------------------------------
	/**
	 * Get new list of edges greater than l in the graph.
	 * 
	 * @param l length with which to threshold
	 * @return new list of edges in the graph with length greater than l
	 */
	public ArrayList<Edge> getLargeEdges(double l) {
		ArrayList<Edge> largeEdges = new ArrayList<>(edges);
		largeEdges.removeIf((e) -> e.getLength() <= l);
		return largeEdges;
	}// end method getLargeEdges

}// end class Graph

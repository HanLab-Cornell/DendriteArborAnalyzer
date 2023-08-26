/*
 * #%L
 * Modified based on Vertex.java of the ImageJ AnalyzeSkeleton_ plugin by Ignacio Arganda-Carreras.
 *
 * Major Modifications:
 * Added fields: marked;
 * Added functions: cloneUnconnected(), aboutEqual(), equals(), isMarked(), setMarked().
 *
 * Diana Bank, Inle Bush, and Elena Zhong, Chun Han Lab, Cornell University
 * #L%
 */
package dend;

import static java.util.stream.Collectors.toList;

import java.util.ArrayList;
import java.util.List;

/** This class represents a vertex or node in a graph. Modified based on {@code Vertex.java} of the
 * ImageJ {@code AnalyzeSkeleton_} plugin by Ignacio Arganda-Carreras.
 * <p>
 * For more information about the original plugin, visit the AnalyzeSkeleton home page:
 * <A target="_blank" href="http://fiji.sc/AnalyzeSkeleton">http://fiji.sc/AnalyzeSkeleton</A>
 *
 * @author Diana Bank
 * @author Inle Bush
 * @author Elena Zhong
 * @author Chun Han Lab, Cornell University */

public class Vertex {
	/** list of points belonging to the vertex */
	private ArrayList<Point> points= null;
	/** list of projecting edges from this vertex */
	private ArrayList<Edge> branches= null;
	/** visit status (for depth-first search, DFS) */
	private boolean visited= false;
	/** previously visited edge in DFS */
	private Edge precedessor= null;
	/** DFS visit order */
	private int visitOrder= -1;
	private boolean marked= false;

	// --------------------------------------------------------------------------
	/** Create empty vertex. */
	public Vertex() {
		points= new ArrayList<>();
		branches= new ArrayList<>();
	}

	// --------------------------------------------------------------------------
	/** Add point to the vertex.
	 *
	 * @param p input point */
	public void addPoint(Point p) {
		points.add(p);
	}

	// --------------------------------------------------------------------------
	/** Check if a point belongs to the vertex list of points.
	 *
	 * @param p input points
	 * @return true if the point is in the vertex point list */
	public boolean isVertexPoint(Point p) {
		return points != null && points.contains(p);
	}

	// --------------------------------------------------------------------------
	/** Convert list of points to String.
	 *
	 * @return printable version of the list of points */
	public String pointsToString() {
		StringBuilder sb= new StringBuilder();
		for (final Point p : points)
			sb.append(p.toString()).append(" ");

		return sb.toString();
	}

	// --------------------------------------------------------------------------
	/** Get list of points.
	 *
	 * @return list of points */
	public ArrayList<Point> getPoints() {
		return points;
	}

	// --------------------------------------------------------------------------
	/** Add a new branch to the vertex.
	 *
	 * @param e neighbor edge */
	public void setBranch(Edge e) {
		branches.add(e);
	}

	// --------------------------------------------------------------------------
	/** Get branch list.
	 *
	 * @return list of branch vertices */
	public ArrayList<Edge> getBranches() {
		return branches;
	}

	// --------------------------------------------------------------------------
	/** Set vertex as visited or not.
	 *
	 * @param b boolean flag */
	public void setVisited(boolean b) {
		visited= b;
	}

	// --------------------------------------------------------------------------
	/** Set vertex as visited or not.
	 *
	 * @param b          boolean flag
	 * @param visitOrder visit order */
	public void setVisited(boolean b, int visitOrder) {
		visited= b;
		this.visitOrder= visitOrder;
	}

	/** Check visit status.
	 *
	 * @return true if the vertex was already visited (DFS) */
	public boolean isVisited() {
		return visited;
	}

	/** Set predecessor (for DFS).
	 *
	 * @param pred predecessor edge in DFS visit. */
	public void setPredecessor(Edge pred) {
		precedessor= pred;
	}

	/** Get predecessor edge.
	 *
	 * @return predecessor edge. */
	public Edge getPredecessor() {
		return precedessor;
	}

	/** Get DFS visit order.
	 *
	 * @return visit order */
	public int getVisitOrder() {
		return visitOrder;
	}

	/** Clones the Vertex disconnected from its {@link Graph}
	 *
	 * @return The Vertex without branches or a predecessor */
	public Vertex cloneUnconnected() {
		final Vertex clone= new Vertex();
		clone.setVisited(visited, visitOrder);
		final List<Point> clonedPoints= points.stream().map(Point::clone)
			.collect(
				toList());
		clone.points.addAll(clonedPoints);
		return clone;
	}

	/** Compare two vertices with loose criterion.
	 *
	 * @param v2 vertex to be compared with
	 * @return true if infinity norm between two points from two vertices each is less than 1 */
	public boolean aboutEqual(Vertex v2) {
		if (v2.getPoints().size() == 0) { throw new NullPointerException(); }
		for (Point p : v2.getPoints()) {
			for (Point p2 : getPoints()) {
				if (p2.infDist(p) <= 1) { return true; }
			}
		}
		return false;
	}

	/** Compare two vertices with strict criterion.
	 *
	 * @param v2 vertex to be compared with
	 * @return true if two vertices contain exactly the same points */
	@Override
	public boolean equals(Object v2) {
		if (v2 == null) { return false; }
		Vertex vt2= (Vertex) v2;
		return vt2.getPoints().equals(getPoints());
	}

	/** Check marking status (for branch ID)
	 *
	 * @return true if the vertex is marked */
	public boolean isMarked() {
		return marked;
	}

	/** Set vertex as marked or not (for branch ID)
	 *
	 * @param marked boolean */
	public void setMarked(boolean marked) {
		this.marked= marked;
	}

}// end class Vertex

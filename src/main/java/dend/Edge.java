/*
 * #%L
 * Modified based on Edge.java of the ImageJ AnalyzeSkeleton_ plugin by Ignacio Arganda-Carreras.
 *
 * Major Modifications:
 * Added functions: clone().
 *
 * Diana Bank, Inle Bush, and Elena Zhong, Chun Han Lab, Cornell University
 * #L%
 *
 */package dend;

import java.util.ArrayList;
import java.util.stream.Collectors;

/** This class represents the edge between two vertices in an undirected graph. Modified based on
 * {@code Edge.java} of the ImageJ {@code AnalyzeSkeleton_} plugin by Ignacio Arganda-Carreras.
 * <p>
 * For more information about the original plugin, visit the AnalyzeSkeleton home page:
 * <A target="_blank" href="http://fiji.sc/AnalyzeSkeleton">http://fiji.sc/AnalyzeSkeleton</A>
 *
 * @author Diana Bank
 * @author Inle Bush
 * @author Elena Zhong
 * @author Chun Han Lab, Cornell University */

public class Edge {
	/** "tree" edge classification constant for Depth-first search (DFS) */
	public final static int TREE= 0;
	/** "back" edge classification constant for Depth-first search (DFS) */
	public final static int BACK= 1;
	/** not yet defined edge classification constant for Depth-first search (DFS) */
	public final static int UNDEFINED= -1;

	/** DFS classification */
	private int type= Edge.UNDEFINED;

	/** vertex at one extreme of the edge */
	private Vertex v1= null;
	/** vertex at the other extreme of the edge */
	private Vertex v2= null;
	/** list of slab voxels belonging to this edge */
	private ArrayList<Point> slabs= null;
	/** length of the edge */
	private double length= 0;
	/** average color of edge */
	private double color= 0;
	/** average color of inner third of edge */
	private double color3rd= 0;
	/** length calculated by running average ovr 5 Pixel */
	private double length_ra= 0;

	/** Create an edge of specific vertices and list of slab voxels.
	 *
	 * @param v1     first vertex
	 * @param v2     second vertex
	 * @param slabs  list of slab voxels
	 * @param length calibrated edge length */
	public Edge(
		Vertex v1,
		Vertex v2,
		ArrayList<Point> slabs,
		double length) {
		this.v1= v1;
		this.v2= v2;
		this.slabs= slabs;
		this.length= length;
	}

	/** Create an edge of specific vertices and list of slab voxels.
	 *
	 * @param v1        first vertex
	 * @param v2        second vertex
	 * @param slabs     list of slab voxels
	 * @param length    calibrated edge length
	 * @param color3rd  average color value of the inner third
	 * @param color     average color value
	 * @param length_ra calibrated edge length calculated with running average ofer 5 Pixel */
	public Edge(
		Vertex v1,
		Vertex v2,
		ArrayList<Point> slabs,
		double length,
		double color3rd,
		double color,
		double length_ra) {
		this.v1= v1;
		this.v2= v2;
		this.slabs= slabs;
		this.length= length;
		this.color= color;
		this.color3rd= color3rd;
		this.length_ra= length_ra;

	}

	/** Get first vertex.
	 *
	 * @return first vertex of the edge */
	public Vertex getV1() {
		return v1;
	}

	/** Get second vertex.
	 *
	 * @return second vertex of the edge */
	public Vertex getV2() {
		return v2;
	}

	/** Get list of slab voxels belonging to the edge.
	 *
	 * @return list of slab voxels */
	public ArrayList<Point> getSlabs() {
		return slabs;
	}

	/** Set DFS type (BACK or TREE)
	 *
	 * @param type DFS classification (BACK or TREE) */
	public void setType(int type) {
		this.type= type;
	}

	/** Get DFS edge type
	 *
	 * @return DFS classification type */
	public int getType() {
		return type;
	}

	/** Get opposite vertex from a given one.
	 *
	 * @param v input vertex
	 * @return opposite vertex in the edge */
	public Vertex getOppositeVertex(Vertex v) {
		if (v1.equals(v))
			return v2;
		else if (v2.equals(v))
			return v1;
		else return null;
	}

	/** Set edge length
	 *
	 * @param length calibrated edge length */
	public void setLength(double length) {
		this.length= length;
	}

	/** Get edge length
	 *
	 * @return calibrated edge length */
	public double getLength() {
		return length;
	}

	/** Get edge length_ra (running average)
	 *
	 * @return calibrated edge length (running average) */
	public double getLength_ra() {
		return length_ra;
	}

	/** Set color
	 *
	 * @param color color of vertex */
	public void setColor(double color) {
		this.color= color;
	}

	/** Get color
	 *
	 * @return color */
	public double getColor() {
		return color;
	}

	/** Set color3rd
	 *
	 * @param color color3rd of vertex */
	public void setColor3rd(double color) {
		color3rd= color;
	}

	/** Get color3rd
	 *
	 * @return color */
	public double getColor3rd() {
		return color3rd;
	}

	/** Clones the Edge with all its properties
	 * <p>
	 * NB Does not clone the vertices!
	 * </p>
	 *
	 * @param v1 One endpoint of the edge, can be null
	 * @param v2 The other endpoint of the edge, can be null
	 * @return new edge */
	public Edge clone(final Vertex v1, final Vertex v2) {
		final ArrayList<Point> clonedSlabs;
		if (slabs != null) {
			clonedSlabs= slabs.stream().map(Point::clone).collect(Collectors
				.toCollection(ArrayList::new));
		} else {
			clonedSlabs= null;
		}
		final Edge clonedEdge= new Edge(v1, v2, clonedSlabs, length, color3rd,
			color, length_ra);
		clonedEdge.setType(type);
		return clonedEdge;
	}

}// end class Edge

/*
 * #%L
 * Modified based on Point.java of the ImageJ AnalyzeSkeleton_ plugin by Ignacio Arganda-Carreras.
 *
 * Major Modifications:
 * Added functions: clone() and infDist().
 *
 * Diana Bank, Inle Bush, and Elena Zhong, Chun Han Lab, Cornell University
 * #L%
 */
package dend;

/** This class represents a 3D point or position on a 3D image. Modified based on {@code Point.java}
 * of the ImageJ {@code AnalyzeSkeleton_} plugin by Ignacio Arganda-Carreras.
 * <p>
 * For more information about the original plugin, visit the AnalyzeSkeleton home page:
 * <A target="_blank" href="http://fiji.sc/AnalyzeSkeleton">http://fiji.sc/AnalyzeSkeleton</A>
 *
 * @author Diana Bank
 * @author Inle Bush
 * @author Elena Zhong
 * @author Chun Han Lab, Cornell University
 * */

public class Point {
	/** x- coordinate */
	public int x= 0;
	/** y- coordinate */
	public int y= 0;
	/** z- coordinate */
	public int z= 0;

	/** Create point from integer coordinates.
	 *
	 * @param x x- coordinate
	 * @param y y- coordinate
	 * @param z z- coordinate */
	public Point(int x, int y, int z) {
		this.x= x;
		this.y= y;
		this.z= z;
	}

	/** Convert point to string. */
	@Override
	public String toString() {
		return "(" + x + ", " + y + ", " + z + ")";
	}

	/** Override equals method to compare points.
	 *
	 * @param o input object
	 * @return true if the input object is equal to this Point */
	@Override
	public boolean equals(Object o) {
		if (this == o) return true;

		if (o == null || getClass() != o.getClass()) return false;

		final Point p= (Point) o;
		return p.x == x && p.y == y && p.z == z;
	}

	@Override
	public Point clone() {
		return new Point(x, y, z);
	}

	/** Calculate L-infinity norm between source point and point p.
	 *
	 * @param p point to be compared with
	 * @return double L-infinity norm */
	public double infDist(Point p) {
		return Math.max(Math.abs(x - p.x), Math.abs(y - p.y));
	}

	public boolean isNeighbor(Point p, int d) {
		return infDist(p) == d;
	}

}// end class point

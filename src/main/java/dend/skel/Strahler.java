/*
 * #%L
 * Modified based on Strahler.java of the ImageJ Strahler Analysis plugin by Tiago Ferreira.
 *
 * Major Modifications:
 * Modified root selection to always select rootpoint in the ROI;
 * Added eliminating endpoints at each pruning cycle;
 * Added displaying Euclidean distance and shortest path distance from rootpoint to each branch;
 * Added displaying branch ID;
 * Added displaying standard deviation of branch length at each Strahler Order;
 * Deleted unused functions.
 *
 * Diana Bank, Inle Bush, and Elena Zhong, Chun Han Lab, Cornell University
 * #L%
 */
package dend.skel;

import java.awt.Choice;
import java.awt.Color;
import java.awt.Font;
import java.awt.Rectangle;
import java.awt.image.ColorModel;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.ListIterator;
import java.util.Map;
import java.util.Objects;

import org.apache.commons.math3.stat.descriptive.SummaryStatistics;

import dend.*;
import dend.processing.Binary;
import dend.sholl.gui.EnhancedGenericDialog;
import ij.IJ;
import ij.ImageJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.WindowManager;
import ij.gui.DialogListener;
import ij.gui.GenericDialog;
import ij.gui.ImageCanvas;
import ij.gui.Overlay;
import ij.gui.Roi;
import ij.gui.TextRoi;
import ij.measure.Calibration;
import ij.measure.ResultsTable;
import ij.plugin.PlugIn;
import ij.plugin.ZProjector;
import ij.process.ImageProcessor;
import sc.fiji.analyzeSkeleton.AnalyzeSkeleton_;
import sc.fiji.skeletonize3D.Skeletonize3D_;

/**
 * This class implements the {@code Dendrite Arbor Analyzer} plugin for ImageJ.
 * Modified based on
 * {@code Strahler.java} of the ImageJ {@code Strahler Analysis} plugin by Tiago
 * Ferreira.
 * <p>
 * For more information about the original plugin, visit the hIPNAT repository
 * {@literal https://github.com/tferr/hIPNAT} and the the plugin's documentation
 * page:
 * {@literal http://imagej.net/Strahler_Analysis}
 *
 * @author Diana Bank
 * @author Inle Bush
 * @author Elena Zhong
 * @author Chun Han Lab, Cornell University
 */
public class Strahler implements PlugIn, DialogListener {

	protected static final String URL = "http://imagej.net/Strahler_Analysis";

	/** Default value for max. number of pruning cycles */
	int maxOrder = 30;

	/** Default option for loop detection */
	private int pruneChoice = AnalyzeSkeleton_.SHORTEST_BRANCH;

	/** Default option for verbose mode */
	private boolean verbose = true;

	/** Default option for tabular option */
	private boolean tabular = true;

	/** Default option for branch ID */
	private boolean branchid = false;

	/** Default option for reverse order */
	private boolean reverseOrder = false;

	/** Default option for frequency table */
	private boolean freqtbl = true;

	/** Default option for branch length frequency table bin size */
	private int binSize = 10;

	/** Default option for path distance frequency table bin size */
	private int binSize2 = 10;

	/** Default option for Euclidean distance frequency table bin size */
	private int binSize3 = 10;

	/** iternum */
	private int iternum = 3;

	/** Maximum length to be consider Higher Order */
	private double largeEdgeThreshold = 60;

	/**
	 * Number of Orders considered Higher Order. Non-positive indicates no length
	 * thresholding
	 */
	private int higherOrderMaxOrder = 3;

	/** Title of main results window */
	private static final String STRAHLER_TABLE = "Strahler Table";

	/** Title of detailed euclidean branch distances to root table */
	private static final String BRANCH_DIST = "Euclidean Branch Distances";

	/** Title of accurate branch distances to root table */
	private static final String BRANCH_PATHS = "Accurate Branch Distances";

	/** Title of comparison results table */
	private static final String COMPARE_TABLE = "Branch Table";

	/** Title of frequency results table for branch length */
	private static final String FREQUENCY_TABLE_LENGTH = "Branch Length Frequency Table";

	/** Title of frequency results table for path distance */
	private static final String FREQUENCY_TABLE_PATH = "Path Distance Frequency Table";

	/** Title of frequency results table for distance */
	private static final String FREQUENCY_TABLE_DISTANCE = "Euclidean Distance Frequency Table";

	/** Title of frequency results table for parent order */
	private static final String FREQUENCY_TABLE_PARENT = "Parent Order Frequency Table";

	/*
	 * Grayscale image for intensity-based pruning of skeleton loops. While it is
	 * unlikely that the iterative pruning of terminal branches will cause new
	 * loops on pre-existing skeletons, offering the option to resolve loops
	 * with intensity based methods remains useful specially when analyzing
	 * non-thinned grayscale images.
	 */
	ImagePlus grayscaleImp = null;
	int grayscaleImpChoice;

	ImagePlus srcImp; // Image to be analyzed (we'll be working on a copy)
	boolean validRootRoi; // Flag assessing validity of 'root-protective' ROI
	String title; // Title of active image
	Roi rootRoi; // Reference to the "root-protecting" ROI
	ImageProcessor ip; // Image Processor for image being analyzed

	/**
	 * Calls {@link fiji.Debug#run(String, String) fiji.Debug.run()} so that the
	 * plugin can be debugged
	 * from an IDE
	 *
	 * @param args the arguments as specified in {@code plugins.config}
	 */
	public static void main(final String[] args) {

		new ImageJ(); // start ImageJ

		ImagePlus impb = IJ.openImage();
		ImagePlus imp = IJ.openImage();
		// imp.setRoi(517, 517, 6, 6);
		imp.setRoi(500, 500, 30, 30); // set ROIs at different locations for debugging purposes
		// imp.setRoi(518, 522, 4, 6);
		// imp.setRoi(486, 463, 25, 23);
		// imp.setRoi(871, 81, 69, 75);

		// throw an exception if the user did not define an ROI around the root
		if (imp.getRoi() == null) {
			throw new RuntimeException("No ROI found");
		}
		impb.show();
		imp.show();
		IJ.runPlugIn(imp, "dend.skel.Strahler", null);
		WindowManager.addWindow(imp.getWindow());
	}

	/**
	 * This method is called when the plugin is loaded.
	 *
	 * @param arg the arguments as specified in {@code plugins.config}
	 */
	@Override
	public void run(final String arg) {

		// Retrieve analysis image and its ROI
		srcImp = WindowManager.getCurrentImage();
		if (!validRequirements(srcImp))
			return;

		title = srcImp.getTitle();
		rootRoi = srcImp.getRoi();
		validRootRoi = rootRoi != null && rootRoi.getType() == Roi.RECTANGLE;

		if (srcImp.getNSlices() > 1) {
			final String warning = "3D images are currently supported with the following limitations:\n" +
					"    - 'Root-protecting' ROIs are not yet supported\n" +
					"    - Lengths are estimated from Z-projections\n \n" +
					"These issues will be addressed in future releases.";
			if (IJ.macroRunning())
				IJ.log(warning);
			else
				IJ.showMessage("Warning", warning);
			validRootRoi = false;
		}
		if (!getSettings())
			return;

		// Work on a skeletonized copy since we'll be modifying the image
		if (rootRoi != null)
			srcImp.killRoi();
		final ImagePlus imp = srcImp.duplicate();
		if (rootRoi != null)
			srcImp.setRoi(rootRoi);
		ip = imp.getProcessor();
		skeletonizeWithoutHermits(imp);

		// Initialize ResultsTable: main and detailed info
		final ResultsTable rt = Utils.getTable(STRAHLER_TABLE);
		ResultsTable tt = Utils.getTable(COMPARE_TABLE);
		ResultsTable ft = Utils.getTable(FREQUENCY_TABLE_LENGTH);
		ResultsTable ft2 = Utils.getTable(FREQUENCY_TABLE_PATH);
		ResultsTable ft3 = Utils.getTable(FREQUENCY_TABLE_DISTANCE);
		ResultsTable ft4 = Utils.getTable(FREQUENCY_TABLE_PARENT);

		// Analyze root
		ImagePlus rootImp;
		ImageProcessor rootIp = null;
		dend.SkeletonResult rootResult = null;
		ArrayList<Point> rootEndpointsList = null;
		int nRootEndpoints = 0, nRootJunctions = 0;
		Vertex rootVertex = null;
		ArrayList<Point> juncts = new ArrayList<>();
		AnalyzeSkeleton2_ root = null;

		// The code for finding distances to the root hinges on
		// having a valid root ROI
		// if none is found, the distances calculated will be inaccurate
		if (validRootRoi) {

			// Duplicate entire canvas. Ignore tree(s) outside ROI
			rootImp = imp.duplicate();
			rootIp = rootImp.getProcessor();
			rootIp.setValue(0.0);
			rootIp.fillOutside(rootRoi);

			// Get root properties
			root = new AnalyzeSkeleton2_();
			root.setup("", rootImp);

			rootResult = root.run(pruneChoice, false, false, grayscaleImp, true,
					false);
			rootImp.flush();

			// Our array juncts will hold all Points inside the root ROI
			juncts.addAll(rootResult.getListOfJunctionVoxels());
			juncts.addAll(rootResult.getListOfEndPoints());
			nRootJunctions = sum(rootResult.getJunctions());

			// Remove end-points at ROI boundaries
			rootEndpointsList = rootResult.getListOfEndPoints();
			final ListIterator<Point> it = rootEndpointsList.listIterator();
			final Rectangle r = rootRoi.getBounds();
			while (it.hasNext()) {
				final Point p = it.next();
				if (p.x == r.x || p.y == r.y ||
						p.x == (int) (r.x + r.getWidth() - 1) ||
						p.y == (int) (r.y + r.getHeight() - 1))
					it.remove();
			}
			rootResult.setListOfEndPoints(rootEndpointsList);
			nRootEndpoints = rootEndpointsList.size();
		}

		// Initialize display images. Use Z-projections to populate
		// iteration stack when dealing with 3D skeletons
		final int nSlices = imp.getNSlices();
		ZProjector zp = null;

		ImageStack imgStack = srcImp.getImageStack();

		final ImageStack iterationStack = new ImageStack(imp.getWidth(),
				imp.getHeight());
		if (nSlices > 1) {
			zp = new ZProjector(imp);
			zp.setMethod(ZProjector.MAX_METHOD);
			zp.setStartSlice(1);
			zp.setStopSlice(nSlices);
		}

		// Initialize AnalyzeSkeleton2_
		final AnalyzeSkeleton2_ as = new AnalyzeSkeleton2_();
		as.setup("", imp);

		// Perform the iterative pruning
		int order = 1, nEndpoints = 0, nJunctions = 0, nJunctions2 = 0;
		int maxBranches = 0;

		dend.SkeletonResult ogSr = as.run(0, false, false, grayscaleImp, true,
				false);
		Graph graphOG = null;

		// Array to hold list of Vertices that are end points of branches
		ArrayList<Point> ends = ogSr.getListOfEndPoints();
		ArrayList<Vertex> endVertices = new ArrayList<>();

		// Find the graph with the most number of branches which should
		// represent the bulk of the neuron skeleton structure
		for (Graph g : ogSr.getGraph()) {
			if (g.getVertices().size() > maxBranches) {
				maxBranches = g.getVertices().size();
				graphOG = g;
			}
		}

		// this loop adds vertices to our list that lie on the endpoints of the graph
		// i.e., the outermost branches that we want to track
		if (graphOG != null) {
			for (Vertex v : graphOG.getVertices()) {
				for (Point p : v.getPoints()) {
					if (ends.contains(p)) {
						endVertices.add(v);
						break;
					}
				}
			}
		}

		// Finds the vertex on the graph associated to the correct root points
		// taken from junction and endpoint voxels in the root ROI
		Vertex[] verts = new Vertex[graphOG.getVertices().size()];
		verts = graphOG.getVertices().toArray(verts);

		rootVertex = null;
		ArrayList<Edge> rootBranches = new ArrayList<>();
		for (Point j : juncts) {
			rootVertex = as.findPointVertex(verts, j);
			if (rootVertex != null) {
				rootBranches = rootVertex.getBranches();
				break;
			}
		}
		Edge[] rootEdges = new Edge[rootBranches.size()];
		rootEdges = rootBranches.toArray(rootEdges);

		Point rootPoint = null;

		// If a root vertex / junction does not exist, we have to create one
		if (rootVertex == null) {
			Rectangle rect = rootRoi.getBounds();
			int x = rect.x + rect.width / 2;
			int y = rect.y + rect.height / 2;

			Point center = new Point(x, y, 0);

			double minDist = Double.MAX_VALUE;
			Point minPoint = null;

			// Finds the point closest to centroid of the ROI
			for (Point p : rootResult.getListOfSlabVoxels()) {
				double dist = calculateDistance(p, center, srcImp);
				if (dist < minDist) {
					minDist = dist;
					minPoint = p;
				}
			}

			// Cuts the edge at that point
			for (Edge e : graphOG.getEdges()) {
				if (e.getSlabs().contains(minPoint)) {
					rootVertex = graphOG.cutEdge(e, minPoint);
					rootPoint = minPoint;
					break;
				}
			}

			rootBranches = rootVertex.getBranches();
			rootEdges = new Edge[rootBranches.size()];
			rootEdges = rootBranches.toArray(rootEdges);
		} else {
			rootPoint = rootVertex.getPoints().get(0);
		}

		// Copy input image
		ColorModel cm = ColorModel.getRGBdefault();
		imgStack = new ImageStack(srcImp.getWidth(), srcImp.getHeight(), cm);
		for (int i = 1; i <= srcImp.getStack().getSize(); i++)
			imgStack.addSlice(srcImp.getImageStack().getSliceLabel(i),
					srcImp.getImageStack().getProcessor(i).duplicate());

		ArrayList<Edge> largeAndRootBranches = null;
		if (higherOrderMaxOrder > 0) {
			largeAndRootBranches = graphOG.getLargeEdges(largeEdgeThreshold);
			largeAndRootBranches.addAll(rootBranches);
		}

		// This do-while loop performs the function of pruning end point branches
		// at each order, re-skeletonizing the image and classifying each branch
		// by Strahler order
		do {

			IJ.showStatus("Retrieving measurements for order " + order + "...");
			IJ.showProgress(order, maxOrder);

			// (Re)skeletonize image
			if (order > 1)
				skeletonizeWithoutHermits(imp);

			// Get properties of loop-resolved tree(s)
			SkeletonResult sr = as.run(pruneChoice, false, false, grayscaleImp,
					true, false);

			ArrayList<Point> totalPoints = new ArrayList<>();
			totalPoints.addAll(sr.getListOfEndPoints());
			totalPoints.addAll(sr.getListOfSlabVoxels());
			totalPoints.addAll(sr.getListOfJunctionVoxels());
			totalPoints.addAll(sr.getListOfStartingSlabVoxels());

			nEndpoints = sum(sr.getEndPoints());
			nJunctions = sum(sr.getJunctions());

			if (order == 1) {
				// Do not include root in 1st order calculations
				nEndpoints -= nRootEndpoints;
				nJunctions -= nRootJunctions;
			}

			// Is it worth proceeding?
			if (nEndpoints == 0 || nJunctions2 == nJunctions) {
				break;
			}

			// Add current tree(s) to iteration stack
			ImageProcessor ipd;
			if (nSlices > 1 && zp != null) {
				zp.doProjection();
				ipd = zp.getProjection().getProcessor();
			} else {
				ipd = ip.duplicate();
			}

			iterationStack.addSlice("Order " + IJ.pad(order, 2), ipd);

			// Remember main results
			nJunctions2 = nJunctions;

			// Eliminate end-points (i.e., single pixel branches). If setTBThreshold,
			// protects large branches from pruning in first higherOrderMaxOrder rounds
			if (order <= higherOrderMaxOrder && !Objects.isNull(largeAndRootBranches)) {
				as.run(pruneChoice, true, false, grayscaleImp, true, false,
						largeAndRootBranches);
			} else {
				as.run(pruneChoice, true, false, grayscaleImp, true, false,
						rootBranches);
			}
		} while (order++ <= maxOrder && nJunctions > 0);

		// Set counter to the de facto order
		order -= 1;

		// Safety check
		if (iterationStack.getSize() < 1) {
			error("Enable \"detailed\" mode.");
			return;
		}

		// Create iteration stack
		final Calibration cal = srcImp.getCalibration();
		final ImagePlus imp2 = new ImagePlus("StrahlerIteration_" + title,
				iterationStack);
		imp2.setCalibration(cal);

		// Generate Strahler mask
		zp = new ZProjector(imp2);
		zp.setMethod(ZProjector.SUM_METHOD);
		zp.setStartSlice(1);
		zp.setStopSlice(order);
		zp.doProjection();
		final ImageProcessor ip3 = zp.getProjection().getProcessor()
				.convertToShortProcessor(false);
		ip3.multiply(1 / 255.0); // map intensities to Strahler orders
		final ImagePlus imp3 = new ImagePlus("StrahlerMask_" + title, ip3);
		imp3.setCalibration(cal);

		// Segment branches by order
		ImagePlus maskImp = imp3.duplicate(); // Calibration is retained
		IJ.setThreshold(maskImp, 1, 6);
		IJ.run(maskImp, "Convert to Mask", "");

		// Analyze segmented order
		AnalyzeSkeleton2_ maskAs = new AnalyzeSkeleton2_();
		maskAs.setup("", maskImp);
		maskAs.run(0, false, false, grayscaleImp, true, false);
		SkeletonResult maskSr = maskAs.run(pruneChoice, false, false,
				grayscaleImp, true, false);
		maskImp.flush();

		// Trees if the user requested it
		int nBranches = sum(maskSr.getBranches());
		maxBranches = nBranches;

		// List to hold all branches
		List<Branch> branchesList2 = new ArrayList<>();

		// Measure segmented orders
		double prevNbranches = Double.NaN;

		for (int i = order; i >= 1; i--) {
			List<Branch> branchesList3 = new ArrayList<>(branchesList2);

			// Segment branches by order
			maskImp = imp3.duplicate(); // Calibration is retained
			IJ.setThreshold(maskImp, i, i);
			IJ.run(maskImp, "Convert to Mask", "");

			// Analyze segmented order
			maskAs = new AnalyzeSkeleton2_();
			maskAs.setup("", maskImp);
			maskSr = maskAs.run(pruneChoice, false, false, grayscaleImp, true,
					false);

			maskImp.flush();

			nBranches = sum(maskSr.getBranches());

			// Run Dijkstra's algorithm to find the shortest path from each vertex to the
			// root
			Dijkstra dij1 = new Dijkstra(rootVertex);
			Map<Vertex, Double> branchDist1 = dij1.getDistances();

			double minxy = 10000;
			for (double xy : branchDist1.values()) {
				if (xy < minxy) {
					minxy = xy;
				}
			}

			// Cut any branch that contains root point as a slab point
			if (i == order) {
				ArrayList<Edge> ess = new ArrayList<>(
						maskSr.getGraph()[0].getEdges());
				for (Edge es : ess) {
					if (es.getSlabs().contains(rootPoint)) {
						maskSr.getGraph()[0].cutEdge(es, rootPoint);
						nBranches += 1;
					}
				}
			}

			// Go through all edges, construct branch, and assign order, distance, and path
			// distance
			for (int k = 0; k < maskSr.getNumOfTrees(); k++) {
				List<Edge> edgesg = maskSr.getGraph()[k].getEdges();
				if (edgesg.toArray().length != 0) {
					for (int k1 = 0; k1 < edgesg.toArray().length; k1++) {
						Edge newe = edgesg.get(k1);
						Branch newb = new Branch(newe);
						if (reverseOrder) {
							newb.setOrder(i);
						} else {
							newb.setOrder(order - i + 1);
						}

						// Calculate Euclidean distance
						Point p1 = newe.getV1().getPoints().get(0);
						Point p2 = newe.getV2().getPoints().get(0);
						double dist1 = calculateDistance(p1, rootPoint, srcImp);
						double dist2 = calculateDistance(p2, rootPoint, srcImp);
						newb.setDist(Math.min(dist1, dist2));

						// Find path distance for this branch from Dijkstra's results
						Double pathdist1 = null;
						Double pathdist2 = null;
						for (Map.Entry<Vertex, Double> entry : branchDist1
								.entrySet()) {
							if (entry.getKey().aboutEqual(newe.getV1())) {
								pathdist1 = entry.getValue();
							}
							if (entry.getKey()
									.aboutEqual(newe.getV2())) {
								pathdist2 = entry.getValue();
							}
						}

						if (pathdist1 != null && pathdist2 != null) {
							newb.addPoint(newe.getV1().getPoints());
							newb.addPoint(newe.getV2().getPoints());
							if (pathdist1 < pathdist2) {
								newb.setPathDist(pathdist1);

							} else {
								newb.setPathDist(pathdist2);
							}
						} else if (pathdist1 != null) {
							newb.setPathDist(pathdist1);
							newb.addPoint(newe.getV1().getPoints());
						} else if (pathdist2 != null) {
							newb.setPathDist(pathdist2);
							newb.addPoint(newe.getV2().getPoints());
						} else {
							newb.setPathDist(-1);
						}
						if (reverseOrder) {
							newb.setPredorder(-10000);
						}

						if (i == order) {
							newb.setPredorder(0);
						} else {
							for (int iter = 1; iter <= iternum && Math
									.abs(newb.getPredorder()) == 10000; iter++) {
								for (Branch bb : branchesList3) {
									if (!reverseOrder &&
											bb.isNeighbor(newb, iter)) {
										if (bb.getOrder() < newb
												.getPredorder()) {
											newb.setPredorder(
													bb.getOrder());
										}
									} else if (reverseOrder &&
											bb.isNeighbor(newb, iter)) {
										if (bb.getOrder() > newb
												.getPredorder()) {
											newb.setPredorder(
													bb.getOrder());
										}
									}
								}
							}
						}

						branchesList2.add(newb);
					}
				}
			}

			// Add edges to summary statistics calculation
			SummaryStatistics ss = new SummaryStatistics();
			for (Graph gg : maskSr.getGraph()) {
				for (Edge ee : gg.getEdges()) {
					ss.addValue(ee.getLength());
				}
			}

			// Log measurements
			rt.incrementCounter();
			rt.addValue("Image", title);
			if (reverseOrder) {
				rt.addValue("Strahler Order", i);
			} else {
				rt.addValue("Strahler Order", order - i + 1);
			}
			rt.addValue("# Branches", nBranches);
			rt.addValue("Ramification ratios", prevNbranches / nBranches);
			rt.addValue("Average branch length",
					average(maskSr.getAverageBranchLength()));
			rt.addValue("Standard deviation",
					ss.getStandardDeviation());
			rt.addValue("Unit", cal.getUnit());

			int co = order - i + 1;
			if (co != 1) {
				ft4.incrementCounter();
				int[] pOrder = new int[order - 1];

				if (reverseOrder) {
					ft4.addValue("Strahler Order", i);
				} else {
					ft4.addValue("Strahler Order", co);
				}

				for (Branch bb : branchesList2) {
					if (!reverseOrder && bb.getOrder() == co) {
						pOrder[bb.getPredorder() - 1] += 1;
					} else if (reverseOrder && bb.getOrder() == i) {
						pOrder[order - bb.getPredorder()] += 1;
					}
				}
				for (int ii = 0; ii < order - 1; ii += 1) {
					if (!reverseOrder) {
						ft4.addValue("Order" + (ii + 1), pOrder[ii]);
					} else {
						ft4.addValue("Order" + (order - ii), pOrder[ii]);
					}
				}
			}

			// Remember results for previous order
			prevNbranches = nBranches;
		}

		// Display Strahler Mask
		if (tabular) {
			ip3.setMinAndMax(0, order + 1);
			ColorMaps.applyMagmaColorMap(imp3, 1, false);
			if (validRootRoi)
				imp3.setRoi(rootRoi);
			imp3.show();
			addCalibrationBar(imp3, 0, "White");
		}

		// Display results table
		if (verbose) {
			rt.show(STRAHLER_TABLE);
			displayAllResults(branchesList2, tt);
		}

		// Display frequency table
		if (freqtbl) {
			ResultsTable tfs = displayFrequencyTable(branchesList2, ft, order,
					binSize, 1);
			tfs.show(FREQUENCY_TABLE_LENGTH);
			ResultsTable tfs2 = displayFrequencyTable(branchesList2, ft2, order,
					binSize2, 2);
			tfs2.show(FREQUENCY_TABLE_PATH);
			ResultsTable tfs3 = displayFrequencyTable(branchesList2, ft3, order,
					binSize3, 3);
			tfs3.show(FREQUENCY_TABLE_DISTANCE);
			ft4.show(FREQUENCY_TABLE_PARENT);
		}

		if (branchid) {
			ImageProcessor ip4 = ip3.duplicate();
			Font font = new Font("SansSerif", Font.PLAIN, 5);
			Overlay overlay = new Overlay();

			// Assign branch ID labels to vertices
			int counter = 0;
			for (Branch b : branchesList2) {
				counter++;
				String TheText = Integer.toString(counter);
				Point pp1 = null;
				if (b.getEdge().getV1().isMarked() == false) {
					pp1 = b.getEdge().getV1().getPoints().get(0);
					b.getEdge().getV1().setMarked(true);
				} else {
					pp1 = b.getEdge().getV2().getPoints().get(0);
					b.getEdge().getV2().setMarked(true);
				}
				TextRoi roi = new TextRoi(pp1.x, pp1.y, TheText, font);
				roi.setStrokeColor(Color.white);
				overlay.add(roi);

				for (Point pp : b.getEdge().getSlabs()) {
					ip4.putPixel(pp.x, pp.y, counter);
				}
				for (Point pp : b.getEdge().getV1().getPoints()) {
					ip4.putPixel(pp.x, pp.y, counter);
				}
				for (Point pp : b.getEdge().getV2().getPoints()) {
					ip4.putPixel(pp.x, pp.y, counter);
				}
			}

			// Display image with branch ID
			ip4.setMinAndMax(-1, counter + 1);

			final ImagePlus imp4 = new ImagePlus("BranchID_" + title, ip4);
			imp4.setOverlay(overlay);
			Roi[] ovlArray = imp4.getOverlay().toArray();
			for (Roi roi : ovlArray) {
				imp4.setRoi(roi);
				IJ.run(imp4, "Draw", "slice");
			}
			if (validRootRoi)
				imp4.setRoi(rootRoi);

			ColorMaps.applyMagmaColorMap(imp4, 1, true);
			imp4.show();
		}

		IJ.showProgress(0, 0);
		IJ.showTime(imp, imp.getStartTime(), "Strahler Analysis concluded... ");
		imp.flush();

	}

	/**
	 * Displays all results for Branches in the skeleton Info includes: Strahler
	 * order, branch length,
	 * path distance, distance
	 *
	 * @param branches List of Branch objects
	 * @param tt       Results table to write to
	 */
	public void displayAllResults(List<Branch> branches, ResultsTable tt) {
		int counter = 0;
		for (Branch b : branches) {
			counter++;
			tt.incrementCounter();
			tt.addValue("Branch #", counter);
			tt.addValue("Strahler Order", b.getOrder());
			tt.addValue("Branch Length", b.getLength());
			tt.addValue("Path Distance", b.getPathDist());
			tt.addValue("Distance", b.getDist());
			tt.addValue("Parent Order", b.getPredorder());
		}
		tt.show(COMPARE_TABLE);
	}

	/**
	 * Displays Frequency Table for selected field
	 *
	 * @param branches List of Branch objects
	 * @param tf       Results table to write to
	 * @param order    Maximum branch order
	 * @param binSize  Histogram bin size
	 * @param field    1:branch length; 2:path distance; 3:Euclidean distance
	 */
	public ResultsTable displayFrequencyTable(List<Branch> branches,
			ResultsTable tf,
			int order, int binSize, int field) {
		ArrayList<Double>[] a0 = new ArrayList[order];
		double maxVal = 0;
		for (Branch b : branches) {
			int bo = b.getOrder();
			double bl;
			if (field == 1) {
				bl = b.getLength();
			} else if (field == 2) {
				bl = b.getPathDist();
			} else {
				bl = b.getDist();
			}
			maxVal = Math.max(maxVal, bl);
			if (a0[bo - 1] == null)
				a0[bo - 1] = new ArrayList<>();
			a0[bo - 1].add(bl);
		}
		// int nbins= Math.round((float) (Math.ceil(maxVal) / binSize));
		for (double bin = 0; bin <= maxVal; bin += binSize) {
			tf.incrementCounter();
			tf.addValue("Bins", bin);
		}

		for (int k = 0; k < order; k++) {
			Double[] data = a0[k].toArray(new Double[a0[k].size()]);
			Arrays.parallelSort(data);
			int i = 0;
			int j = 0;
			for (double bin = 0; bin <= maxVal; bin += binSize) {
				int count = 0;
				while (j < data.length && data[j] < bin + binSize) {
					j += 1;
					count += 1;
				}
				String str = "order" + (k + 1);
				tf.setValue(str, i, count);
				i++;
			}

		}

		return tf;
	}

	/**
	 * Display the Euclidean distance from each branch endpoint to the root of the
	 * skeleton in a
	 * ResultsTable called "Branch Distances"
	 *
	 * @param order       the Strahler order of the branches
	 * @param maxBranches the maximum number of branches in one Strahler order
	 * @param distances   the distances of each branch to the root
	 * @param bt          the results table to populate
	 */
	public void displayBranchDistances(int order, int maxBranches,
			double[][] distances, ResultsTable bt) {

		// Set up the heading of the results table
		bt.incrementCounter();
		bt.addValue("Branch #", "Average Distance");

		// calculate and display the average distance in the results table
		for (int i = 0; i < order; i++) {
			double average = 0;
			int counter = 0;
			for (int j = 0; j < maxBranches; j++) {
				if (distances[j][i] > 0.0) {
					counter++;
					average += distances[j][i];
				}
			}
			average /= counter;
			bt.addValue("Strahler Order " + (order - i), average);
		}
		bt.incrementCounter();

		// Loop through and display branch lengths in the results table
		for (int i = 0; i < maxBranches; i++) {
			bt.incrementCounter();
			bt.addValue("Branch #", i);
			for (int j = 0; j < order; j++) {
				bt.addValue("Strahler Order " + (order - j), distances[i][j]);
			}
		}

		bt.show(BRANCH_DIST);
	}

	/**
	 * Displays the exact distances from each vertex in the graph to the root in a
	 * results table
	 *
	 * @param distances map of vertices to distances
	 * @param endVerts  list of endpoint vertices
	 * @param dt        our exact distances ResultsTable
	 */
	public void displayExactDistances(Map<Vertex, Double> distances,
			ArrayList<Vertex> endVerts, ResultsTable dt) {
		dt.incrementCounter();

		int counter = 0;

		for (Map.Entry<Vertex, Double> entry : distances.entrySet()) {
			counter++;
			dt.incrementCounter();

			Vertex key = entry.getKey();
			Double val = entry.getValue();

			dt.addValue("Vertex #", counter);
			dt.addValue("Distance", val);

			if (endVerts.contains(key))
				dt.addValue("EndPoint Vertex", "Yes");
			else
				dt.addValue("EndPoint Vertex", "No");
		}
		dt.show(BRANCH_PATHS);
	}

	/**
	 * Outputs the comparison values for the exact and estimated branch distances
	 *
	 * @param g          graph
	 * @param branchDist map of vertices to distances
	 * @param ct         results table
	 * @param rootPoint  root point
	 */
	public void displayComparisonTable(Graph g, Map<Vertex, Double> branchDist,
			ResultsTable ct, Point rootPoint) {

		ct.incrementCounter();
		int count = 1;

		for (Vertex v : g.getVertices()) {
			if (v.getPoints().size() >= 1) {
				Point branchPoint = v.getPoints().get(0);
				double dist = calculateDistance(branchPoint, rootPoint, srcImp);
				if (branchDist.get(v) != null) {
					count++;
					ct.incrementCounter();
					ct.addValue("Branch #", count);
					ct.addValue("Exact Distance", branchDist.get(v));
					ct.addValue("Estimated Distance", dist);
				}
			}
		}
		ct.show(COMPARE_TABLE);

	}

	/**
	 * Calculate Euclidean distance between two points in 3D.
	 *
	 * @param point1 first point coordinates
	 * @param point2 second point coordinates
	 * @return distance (in the corresponding units)
	 */
	private double calculateDistance(Point point1, Point point2,
			ImagePlus image) {
		return Math.sqrt(Math
				.pow((point1.x - point2.x) * image.getCalibration().pixelWidth, 2) +
				Math.pow((point1.y - point2.y) * image.getCalibration().pixelHeight,
						2)
				+
				Math.pow((point1.z - point2.z) * image.getCalibration().pixelDepth,
						2));
	}

	/**
	 * Checks if image to be analyzed fulfills analysis requirements and warns the
	 * user if required
	 * dependencies are present (i.e,, if all the required update sites have been
	 * subscribed).
	 *
	 * @param imp the image to be analyzed
	 * @return {@code true}, if assessment was successful. If {@code false} a macro
	 *         friendly
	 *         {@link Utils#error} is displayed.
	 */
	boolean validRequirements(final ImagePlus imp) {
		boolean validImp = imp != null && imp.getBitDepth() == 8;
		final boolean validSetup = Utils.validSkelDependencies();
		if (!validImp) {
			final String msg = imp == null ? "An 8-bit image is required but none was found."
					: imp.getTitle() + " is not an 8-bit image.";
			if (IJ.macroRunning()) {
				Utils.error("Invalid image", msg, imp);
			} else {
				final GenericDialog gd = new GenericDialog("Invalid Image");
				gd.addMessage(msg);
				gd.enableYesNoCancel("OK", "Analyze Sample Image");
				gd.hideCancelButton();
				gd.showDialog();
				if (!gd.wasOKed() && !gd.wasCanceled()) {
					final LSystemsTree lst = new LSystemsTree();
					srcImp = lst.sampleTree();
					srcImp.setRoi(58, 130, 25, 35);
					srcImp.show();
					new ij.plugin.Zoom().run("in");
					validImp = true;
				}
			}
		}
		return validSetup && validImp;
	}

	/**
	 * Displays an error message that will not disrupt macro calls. This is useful
	 * for batch processing
	 * of images: Even if the analysis of a particular image fails, remaining images
	 * can still be
	 * analyzed by the same macro
	 *
	 * @param errorMsg the error message
	 */
	private void error(final String errorMsg) {
		Utils.error("Strahler Analysis", errorMsg, srcImp);
	}

	/**
	 * Gets the analysis parameters from the user.
	 *
	 * @return {@code true} if the dialog input is valid and dialog was not
	 *         dismissed.
	 */
	private boolean getSettings() {

		final EnhancedGenericDialog gd = new EnhancedGenericDialog(
				"Dendrite Arbor Analyzer");
		final Font headerFont = new Font("SansSerif", Font.BOLD, 12);
		gd.setSmartRecording(true);

		// Part 2: Loop elimination
		gd.setInsets(0, 0, 0);
		gd.addMessage("Elimination of Skeleton Loops:", headerFont);
		gd.addChoice("Method:", AnalyzeSkeleton_.pruneCyclesModes,
				AnalyzeSkeleton_.pruneCyclesModes[pruneChoice]);

		// 8-bit grayscale is the only image type recognized by
		// AnalyzeSkeleton2_,
		// so we'll provide the user with a pre-filtered list of valid choices
		final ArrayList<Integer> validIds = new ArrayList<>();
		final ArrayList<String> validTitles = new ArrayList<>();
		final int[] ids = WindowManager.getIDList();
		for (int i = 0; i < ids.length; ++i) {
			final ImagePlus imp = WindowManager.getImage(ids[i]);
			if (imp.getBitDepth() == 8) {
				validIds.add(ids[i]);
				validTitles.add(imp.getTitle());
			}
		}
		gd.addChoice("8-bit grayscale image:",
				validTitles.toArray(new String[validTitles.size()]), title);

		// Part 4: Output
		gd.setInsets(0, 0, 0);
		gd.addMessage("Output Options:", headerFont);
		gd.addCheckbox("Display Tables", verbose);
		gd.addCheckbox("Display Strahler Mask", tabular);
		gd.addCheckbox("Display Branch ID", branchid);
		gd.addCheckbox("Reverse Branch Order", reverseOrder);
		gd.addCheckbox("Display Frequency Table", freqtbl);
		gd.addNumericField("Branch Length Bin Size", binSize, 0);
		gd.addNumericField("Path Distance Bin Size", binSize2, 0);
		gd.addNumericField("Euclidean Distance Bin Size",
				binSize3, 0);
		gd.addNumericField("Iteration Number for Parent Order", iternum, 0);
		gd.addNumericField("Maximum length to be Higher Order",
				largeEdgeThreshold, 0);
		gd.addNumericField("# Orders to be Higher Order (0 for no threshold)",
				higherOrderMaxOrder, 0);

		gd.addDialogListener(this);
		dialogItemChanged(gd, null);

		gd.showDialog();

		// Set grayscale image for intensity-based pruning of skeleton loops
		if (pruneChoice == AnalyzeSkeleton_.LOWEST_INTENSITY_VOXEL ||
				pruneChoice == AnalyzeSkeleton_.LOWEST_INTENSITY_BRANCH) {
			grayscaleImp = WindowManager
					.getImage(validIds.get(grayscaleImpChoice));
		} else {
			grayscaleImp = null;
		}

		return gd.wasOKed();
	}

	/* Retrieve dialog options using the DialogListener interface */
	@Override
	public boolean dialogItemChanged(final GenericDialog gd,
			final java.awt.AWTEvent e) {

		pruneChoice = gd.getNextChoiceIndex();
		grayscaleImpChoice = gd.getNextChoiceIndex();
		verbose = gd.getNextBoolean();
		tabular = gd.getNextBoolean();
		branchid = gd.getNextBoolean();
		reverseOrder = gd.getNextBoolean();
		freqtbl = gd.getNextBoolean();
		binSize = (int) gd.getNextNumber();
		binSize2 = (int) gd.getNextNumber();
		binSize3 = (int) gd.getNextNumber();
		iternum = (int) gd.getNextNumber();
		largeEdgeThreshold = (double) gd.getNextNumber();
		higherOrderMaxOrder = (int) gd.getNextNumber();

		// Enable/Disable key components of GenericDialog
		if (!IJ.macroRunning()) {
			final Choice cImgChoice = (Choice) gd.getChoices().elementAt(0);

			cImgChoice.setEnabled(
					pruneChoice == AnalyzeSkeleton_.LOWEST_INTENSITY_VOXEL ||
							pruneChoice == AnalyzeSkeleton_.LOWEST_INTENSITY_BRANCH ||
							pruneChoice == AnalyzeSkeleton_.SHORTEST_BRANCH ||
							pruneChoice == AnalyzeSkeleton_.NONE);
		}

		return !gd.wasCanceled();
	}

	/**
	 * Returns the sum of the values in the input array, or zero if the array is
	 * empty or {@code null}.
	 *
	 * @param array
	 *              array of values to be summed
	 * @return the sum of elements in the array. Returns zero if array is
	 *         {@code null} or empty.
	 */
	int sum(final int[] array) {
		if (array != null && array.length > 0)
			return Arrays.stream(array).sum();
		return 0;
	}

	/**
	 * Returns the sum of the values in the input array, or zero if the array is
	 * empty or {@code null}.
	 *
	 * @param array
	 *              array of values to be summed
	 * @return the sum of elements in the array. Returns zero if array is
	 *         {@code null} or empty.
	 */
	double sum(final double[] array) {
		if (array != null && array.length > 0)
			return Arrays.stream(array).sum();
		return 0;
	}

	/**
	 * Returns the arithmetic mean of the values in the input array, or
	 * {@code Double.NaN} if the array is empty or {@code null}.
	 *
	 * @param array
	 *              array of values to be averaged
	 * @return the arithmetic mean of the array. Returns {@code Double.NaN} if
	 *         array is {@code null} or empty.
	 */
	double average(final double[] array) {
		if (array != null && array.length > 0)
			return sum(array) / array.length;
		// org.apache.commons.math3.stat.StatUtils?
		return Double.NaN;
	}

	/* Paints point positions. */
	void paintPoints(final ImageStack stack, final ArrayList<Point> points,
			final int value, final String sliceLabel) {
		if (points != null) {
			final ImageProcessor ipp = stack.getProcessor(1)
					.createProcessor(stack.getWidth(), stack.getHeight());
			// ip.createProcessor(stack.getWidth(), stack.getHeight());
			for (int j = 0; j < points.size(); j++) {
				final Point point = points.get(j);
				ipp.putPixel(point.x, point.y, value);
			}
			stack.addSlice(sliceLabel, ipp);
		}
	}

	/*
	 * Skeletonization method that erodes the thinned structure in order to
	 * eliminate isolated pixels. Thinning and pruning may give rise to single
	 * point arbors. These 'debris' trees have 1 end-point but no branches or
	 * junctions. If present they overestimate the total number of end-points
	 */
	private void skeletonizeWithoutHermits(final ImagePlus imp) {
		final Skeletonize3D_ thin = new Skeletonize3D_();
		thin.setup("", imp);
		thin.run(null);
		Binary.removeIsolatedPixels(imp);
	}

	/**
	 * Runs {@link ij.plugin.CalibrationBar} on the specified image using sensible
	 * settings.
	 *
	 * @param imp     the image to processed
	 * @param nLabels the n. of labels in the calibration bar
	 * @param color   Labels' foreground color as per {@link ij.plugin.Colors}
	 **/
	private void addCalibrationBar(final ImagePlus imp, final int nLabels,
			final String color) {
		final ImageCanvas ic = imp.getCanvas();
		double zoom = imp.getHeight() > 200 ? 1.0 : 0.8;
		final double mag = ic != null ? ic.getMagnification() : 1.0;
		if (zoom <= 1 && mag < 1)
			zoom = 1.0 / mag;
		IJ.run(imp, "Calibration Bar...",
				"fill=None label=" + color + " number=" + nLabels + " zoom=" +
						zoom + " overlay");
	}

}

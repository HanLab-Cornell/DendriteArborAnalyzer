/*
 * #%L
 * Modified based on AnalyzeSkeleton_.java of the ImageJ AnalyzeSkeleton_ plugin by Ignacio Arganda-Carreras.
 *
 * Major Modifications:
 * Set silent run to true;
 * Added root protection.
 *
 * Diana Bank, Inle Bush, and Elena Zhong, Chun Han Lab, Cornell University
 * #L%
 */
package dend;

import java.awt.AWTEvent;
import java.awt.Checkbox;
import java.awt.Font;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.Iterator;
import java.util.List;
import java.util.ListIterator;

import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.Prefs;
import ij.WindowManager;
import ij.gui.DialogListener;
import ij.gui.GenericDialog;
import ij.gui.Roi;
import ij.measure.Calibration;
import ij.measure.ResultsTable;
import ij.plugin.filter.PlugInFilter;
import ij.plugin.frame.Recorder;
import ij.process.ByteProcessor;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;

/**
 * This class performs pruning and skeleton analysis. Modified based on
 * {@code AnalyzeSkeleton_.java} of the ImageJ {@code AnalyzeSkeleton_} plugin
 * by Ignacio
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

public class AnalyzeSkeleton2_ implements PlugInFilter, DialogListener {
	private static final String PRUNE_MODE_INDEX_KEY = "sc.fiji.analyzeSkeleton.pruneModeIndex";
	private static final String PRUNE_ENDS_KEY = "sc.fiji.analyzeSkeleton.pruneEnds";
	private static final String CALCULATE_PATH_KEY = "sc.fiji.analyzeSkeleton.shortestPath";
	private static final String VERBOSE_KEY = "sc.fiji.analyzeSkeleton.showDetailedInfo";
	private static final String DISPLAY_SKELETONS_KEY = "sc.fiji.analyzeSkeleton.displayLabeledSkeletons";
	private static final int DEFAULT_PRUNE_MODE_INDEX = AnalyzeSkeleton2_.SHORTEST_BRANCH;
	private static final boolean DEFAULT_PRUNE_ENDS = true;
	private static final boolean DEFAULT_CALCULATE_SHORTEST_PATH = false;
	private static final boolean DEFAULT_VERBOSE = false;
	private static final boolean DEFAULT_PROTECT_ROI = false;
	private static final boolean DEFAULT_DISPLAY_SKELETONS = false;
	private static final String HELP_URL = "http://fiji.sc/wiki/index.php/AnalyzeSkeleton";

	/** end point flag */
	public static byte END_POINT = 30;
	/** junction flag */
	public static byte JUNCTION = 70;
	/** slab flag */
	public static byte SLAB = 127;
	/** shortest path flag */
	public static byte SHORTEST_PATH = 96;

	private GenericDialog settingsDialog;

	/** working image plus */
	private ImagePlus imRef = null;

	/** working image width */
	private int width = 0;
	/** working image height */
	private int height = 0;
	/** working image depth */
	private int depth = 0;
	/** working image stack */
	private ImageStack inputImage = null;

	/** visit flags */
	private boolean[][][] visited = null;

	// Measures
	/** total number of end points voxels */
	private int totalNumberOfEndPoints = 0;
	/** total number of junctions voxels */
	public int totalNumberOfJunctionVoxels = 0;
	/** total number of slab voxels */
	private int totalNumberOfSlabs = 0;
	/**
	 * auxiliary variable to store the length of longest shortest path in a graph
	 */
	private double shortestPath = 0;

	// Shortest path variables
	/** list of longest shortest paths from the skeletons in the image */
	private ArrayList<Double> shortestPathList;
	/** list containing longest shortest path points (one per tree) */
	private ArrayList<Point>[] shortestPathPoints = null;
	/** shortest path x start position */
	private int spx = 0;
	/** shortest path y start position */
	private int spy = 0;
	/** shortest path z start position */
	private int spz = 0;
	/** shortest path start position array */
	private double[][] spStartPosition;
	/** shortest path output stack */
	private ImageStack shortPathImage = null;

	// Tree fields
	/** image stack containing all skeletons marked with their corresponding id */
	ImageStack labeledSkeletons = null;

	/** number of branches for every specific tree */
	private int[] numberOfBranches = null;
	/** number of end points voxels of every tree */
	private int[] numberOfEndPoints = null;
	/** number of junctions voxels of every tree */
	private int[] numberOfJunctionVoxels = null;
	/** number of slab voxels of every specific tree */
	private int[] numberOfSlabs = null;
	/** number of junctions of every specific tree */
	private int[] numberOfJunctions = null;
	/** number of triple points in every tree */
	private int[] numberOfTriplePoints = null;
	/** number of quadruple points in every tree */
	private int[] numberOfQuadruplePoints = null;
	/** list of end points in every tree */
	private ArrayList<Point> endPointsTree[] = null;
	/** list of junction voxels in every tree */
	private ArrayList<Point> junctionVoxelTree[] = null;
	/** list of special slab coordinates where circular tree starts */
	private ArrayList<Point> startingSlabTree[] = null;

	/** average branch length */
	private double[] averageBranchLength = null;

	/** maximum branch length */
	private double[] maximumBranchLength = null;

	/** list of end point coordinates in the entire image */
	private ArrayList<Point> listOfEndPoints = null;
	/** list of junction coordinates in the entire image */
	public ArrayList<Point> listOfJunctionVoxels = null;
	/** list of slab coordinates in the entire image */
	private ArrayList<Point> listOfSlabVoxels = null;
	/** list of slab coordinates in the entire image */
	private ArrayList<Point> listOfStartingSlabVoxels = null;

	/**
	 * list of groups of junction voxels that belong to the same tree junction (in
	 * every tree)
	 */
	private ArrayList<ArrayList<Point>> listOfSingleJunctions[] = null;
	/** array of junction vertex per tree */
	private Vertex[][] junctionVertex = null;

	/**
	 * stack image containing the corresponding skeleton tags (end point, junction
	 * or slab)
	 */
	public ImageStack taggedImage = null;

	/** auxiliary temporary point */
	private Point auxPoint = null;

	/** number of trees (skeletons) in the image */
	private int numOfTrees = 0;

	/** pruning option */
	private boolean bPruneCycles = true;

	/** dead-end pruning option */
	public static boolean pruneEnds = DEFAULT_PRUNE_ENDS;

	/** protective-ROI option (branches inside ROI are spared from pruning) */
	public static boolean protectRoi = DEFAULT_PROTECT_ROI;

	/** calculate largest shortest path option */
	public static boolean calculateShortestPath = DEFAULT_CALCULATE_SHORTEST_PATH;

	/** array of graphs (one per tree) */
	private Graph[] graph = null;

	/** auxiliary list of slabs */
	private ArrayList<Point> slabList = null;
	/** auxiliary final vertex */
	private Vertex auxFinalVertex = null;

	/** prune cycle options */
	public static final String[] pruneCyclesModes = { "none",
			"shortest branch",
			"lowest intensity voxel",
			"lowest intensity branch" };
	/** no pruning mode index */
	public static final int NONE = 0;
	/** shortest branch pruning mode index */
	public static final int SHORTEST_BRANCH = 1;
	/** lowest pixel intensity pruning mode index */
	public static final int LOWEST_INTENSITY_VOXEL = 2;
	/** lowest intensity branch pruning mode index */
	public static final int LOWEST_INTENSITY_BRANCH = 3;

	/** original grayscale image (for lowest pixel intensity pruning mode) */
	private ImageStack originalImage = null;

	/** prune cycle options index */
	public static int pruneIndex = DEFAULT_PRUNE_MODE_INDEX;

	/** x- neighborhood offset */
	private int x_offset = 1;
	/** y- neighborhood offset */
	private int y_offset = 1;
	/** z- neighborhood offset */
	private int z_offset = 1;

	/** boolean flag to display extra information in result tables */
	public static boolean verbose = DEFAULT_VERBOSE;

	/** silent run flag, to distinguish between GUI and plugin calls */
	protected boolean silent = true;

	/** debugging flag */
	private static final boolean debug = false;

	/** flag to output the labeled skeletons in a new image */
	public static boolean displaySkeletons = DEFAULT_DISPLAY_SKELETONS;

	/**
	 * returns this tagged image
	 *
	 * @return tagged image
	 */
	public ImageStack getTaggedImage() {
		return taggedImage;
	}

	/**
	 * sets this tagged image
	 *
	 * @param img image
	 */
	public void setTaggedImage(ImageStack img) {
		taggedImage = img;
	}

	/* ----------------------------------------------------------------------- */
	/**
	 * This method is called once when the filter is loaded.
	 *
	 * @param arg argument specified for this plugin
	 * @param imp currently active image
	 * @return flag word that specifies the filters capabilities
	 */
	@Override
	public int setup(String arg, ImagePlus imp) {
		imRef = imp;

		if (arg.equals("about")) {
			showAbout();
			return DONE;
		}

		return DOES_8G;
	} // end method setup

	/* ----------------------------------------------------------------------- */
	/**
	 * Process the image: tag skeleton and show results.
	 *
	 * @see ij.plugin.filter.PlugInFilter#run(ij.process.ImageProcessor)
	 */
	@Override
	public void run(ImageProcessor ip) {

		loadDialogSettings();
		createSettingsDialog();
		settingsDialog.showDialog();
		if (settingsDialog.wasCanceled()) {
			return;
		}
		setSettingsFromDialog();
		saveDialogSettings();

		pruneIndex = 1;
		// pre-checking if another image is needed and also setting bPruneCycles
		ImagePlus origIP = null;
		switch (pruneIndex) {

			// No pruning
			case AnalyzeSkeleton2_.NONE:
				bPruneCycles = false;
				break;
			// Pruning cycles by shortest branch
			case AnalyzeSkeleton2_.SHORTEST_BRANCH:
				bPruneCycles = true;
				break;
			// Pruning cycles by lowest pixel intensity
			case AnalyzeSkeleton2_.LOWEST_INTENSITY_VOXEL:
			case AnalyzeSkeleton2_.LOWEST_INTENSITY_BRANCH:
				// Select original image between the open images
				int[] ids = WindowManager.getIDList();
				if (ids == null || ids.length < 1) {
					IJ.showMessage("You should have at least one image open.");
					return;
				}

				String[] titles = new String[ids.length];
				for (int i = 0; i < ids.length; ++i) {
					titles[i] = WindowManager.getImage(ids[i]).getTitle();
				}

				final GenericDialog gd2 = new GenericDialog("Image selection");

				gd2.addMessage("Select original grayscale image:");
				final String current = WindowManager.getCurrentImage().getTitle();
				gd2.addChoice("original_image", titles, current);

				gd2.showDialog();

				if (gd2.wasCanceled())
					return;

				// Get original stack
				origIP = WindowManager.getImage(ids[gd2.getNextChoiceIndex()]);

				bPruneCycles = true;
				break;
			default:
		}

		// now we have all the information that's needed for running the plugin
		// as if it was called from somewhere else

		run(pruneIndex, 50.0, false, origIP, false, verbose);

		if (debug)
			IJ.log("num of skeletons = " + numOfTrees);

		// Show labeled skeletons
		if (AnalyzeSkeleton2_.displaySkeletons) {
			ImagePlus labeledSkeletons = new ImagePlus(
					imRef.getShortTitle() + "-labeled-skeletons",
					this.labeledSkeletons.duplicate());
			IJ.run(labeledSkeletons, "Fire", null);
			labeledSkeletons.show();
		}

		// Show results table
		showResults();

	} // end run method

	/** Disables dialog components that are irrelevant to GUI-based analysis. */
	@Override
	public boolean dialogItemChanged(GenericDialog gd, AWTEvent e) {
		if (imRef.getRoi() == null && null != gd &&
				null != gd.getCheckboxes()) {
			Checkbox roiOption = (Checkbox) gd.getCheckboxes().elementAt(1);
			roiOption.setEnabled(false);
			if (Recorder.record)
				roiOption.setState(false);
		}
		return true;
	}

	/**
	 * This method is intended for non-interactively using this plugin.
	 * <p>
	 *
	 * @param pruneIndex The pruneIndex, as asked by the initial gui dialog.
	 * @param pruneEnds  flag to prune end-point-ending branches
	 * @param shortPath  flag to calculate the longest shortest path
	 * @param origIP     original grayscale input image (for lowest pixel intensity
	 *                   pruning mode)
	 * @param silent     silent run
	 * @param verbose    flag to display running information
	 */
	public SkeletonResult run(
			int pruneIndex,
			boolean pruneEnds,
			boolean shortPath,
			ImagePlus origIP,
			boolean silent,
			boolean verbose) {
		Roi roi = null;
		return run(pruneIndex, pruneEnds, shortPath, origIP, silent, verbose,
				roi);
	}

	/**
	 * This method is intended for non-interactively using this plugin.
	 * <p>
	 *
	 * @param pruneIndex The pruneIndex, as asked by the initial gui dialog.
	 * @param pruneEnds  flag to prune end-point-ending branches
	 * @param shortPath  flag to calculate the longest shortest path
	 * @param origIP     original grayscale input image (for lowest pixel intensity
	 *                   pruning mode)
	 * @param silent     silent run
	 * @param verbose    flag to display running information
	 * @param roi        points inside this region are spared from elimination when
	 *                   pruning end
	 *                   branches. ROI can be associated to a single image in the
	 *                   stack or all images as
	 *                   per {@link ij.gui.Roi#getPosition ij.gui.Roi.getPosition()}
	 */
	public SkeletonResult run(
			int pruneIndex,
			boolean pruneEnds,
			boolean shortPath,
			ImagePlus origIP,
			boolean silent,
			boolean verbose,
			Roi roi) {
		AnalyzeSkeleton2_.pruneIndex = pruneIndex;
		this.silent = silent;
		AnalyzeSkeleton2_.pruneEnds = pruneEnds;
		AnalyzeSkeleton2_.calculateShortestPath = shortPath;
		AnalyzeSkeleton2_.verbose = verbose;

		switch (pruneIndex) {
			// No pruning
			case AnalyzeSkeleton2_.NONE:
				bPruneCycles = false;
				break;
			// Pruning cycles by shortest branch
			case AnalyzeSkeleton2_.SHORTEST_BRANCH:
				bPruneCycles = true;
				break;
			// Pruning cycles by lowest pixel intensity
			case AnalyzeSkeleton2_.LOWEST_INTENSITY_VOXEL:
			case AnalyzeSkeleton2_.LOWEST_INTENSITY_BRANCH:
				// calculate neighborhood size given the calibration
				calculateNeighborhoodOffsets(origIP.getCalibration());
				originalImage = origIP.getStack();
				bPruneCycles = true;
				break;
			default:
		}

		width = imRef.getWidth();
		height = imRef.getHeight();
		depth = imRef.getStackSize();
		inputImage = imRef.getStack();

		// initialize visit flags
		resetVisited();

		// Tag skeleton, differentiate trees and visit them
		processSkeleton(inputImage, null);

		// prune ends
		if (pruneEnds) {
			pruneEndBranches(inputImage, taggedImage, roi);
		}

		// Prune cycles if necessary
		if (bPruneCycles) {
			if (pruneCycles(inputImage, originalImage,
					AnalyzeSkeleton2_.pruneIndex)) {
				// initialize visit flags
				resetVisited();
				// Recalculate analysis over the new image
				bPruneCycles = false;
				processSkeleton(inputImage, null);
			}
		}

		// Calculate triple points (junctions with exactly 3 branches)
		calculateTripleAndQuadruplePoints();

		if (shortPath) {
			if (debug)
				IJ.log("Calculating longest shortest paths...");

			// Copy input image
			shortPathImage = new ImageStack(width, height,
					inputImage.getColorModel());
			for (int i = 1; i <= inputImage.getSize(); i++)
				shortPathImage.addSlice(inputImage.getSliceLabel(i),
						inputImage.getProcessor(i).duplicate());

			shortestPathList = new ArrayList<>();

			shortestPathPoints = new ArrayList[numOfTrees];

			// Visit skeleton and measure distances.
			// and apply Warshall algorithm
			spStartPosition = new double[numOfTrees][3];
			for (int i = 0; i < numOfTrees; i++) {
				shortestPathPoints[i] = new ArrayList<>();
				// Warshall algorithm including tag positions
				shortestPath = warshallAlgorithm(graph[i],
						shortestPathPoints[i]);
				shortestPathList.add(shortestPath);
				spStartPosition[i][0] = spx * imRef.getCalibration().pixelWidth;
				spStartPosition[i][1] = spy * imRef.getCalibration().pixelHeight;
				spStartPosition[i][2] = spz * imRef.getCalibration().pixelDepth;
			}

			if (!silent) {
				// Display short paths in a new stack
				ImagePlus shortIP = new ImagePlus("Longest shortest paths",
						shortPathImage);
				shortIP.show();

				// Set same calibration as the input image
				shortIP.setCalibration(imRef.getCalibration());

				// We apply the Fire LUT and reset the min and max to be between 0-255.
				IJ.run(shortIP, "Fire", null);

				// IJ.resetMinAndMax();
				shortIP.resetDisplayRange();
				shortIP.updateAndDraw();
			}

		}

		// Return the analysis results
		return assembleResults();
	}

	/**
	 * This method is intended for non-interactively using this plugin.
	 * <p>
	 *
	 * @param pruneIndex The pruneIndex, as asked by the initial gui dialog.
	 * @param pruneEnds  flag to prune end-point-ending branches
	 * @param shortPath  flag to calculate the longest shortest path
	 * @param origIP     original grayscale input image (for lowest pixel intensity
	 *                   pruning mode)
	 * @param silent     silent run
	 * @param verbose    flag to display running information
	 * @param branches   points inside these branches are spared from elimination
	 *                   when pruning end
	 *                   branches.
	 */
	public SkeletonResult run(
			int pruneIndex,
			boolean pruneEnds,
			boolean shortPath,
			ImagePlus origIP,
			boolean silent,
			boolean verbose,
			ArrayList<Edge> branches) {
		AnalyzeSkeleton2_.pruneIndex = pruneIndex;
		this.silent = silent;
		AnalyzeSkeleton2_.pruneEnds = pruneEnds;
		AnalyzeSkeleton2_.calculateShortestPath = shortPath;
		AnalyzeSkeleton2_.verbose = verbose;

		switch (pruneIndex) {
			// No pruning
			case AnalyzeSkeleton2_.NONE:
				bPruneCycles = false;
				break;
			// Pruning cycles by shortest branch
			case AnalyzeSkeleton2_.SHORTEST_BRANCH:
				bPruneCycles = true;
				break;
			// Pruning cycles by lowest pixel intensity
			case AnalyzeSkeleton2_.LOWEST_INTENSITY_VOXEL:
			case AnalyzeSkeleton2_.LOWEST_INTENSITY_BRANCH:
				// calculate neighborhood size given the calibration
				calculateNeighborhoodOffsets(origIP.getCalibration());
				originalImage = origIP.getStack();
				bPruneCycles = true;
				break;
			default:
		}

		width = imRef.getWidth();
		height = imRef.getHeight();
		depth = imRef.getStackSize();
		inputImage = imRef.getStack();

		// initialize visit flags
		resetVisited();

		// Tag skeleton, differentiate trees and visit them
		processSkeleton(inputImage, null);

		// prune ends
		if (pruneEnds) {
			pruneEndBranches(inputImage, taggedImage, branches);
		}

		// Prune cycles if necessary
		if (bPruneCycles) {
			if (pruneCycles(inputImage, originalImage,
					AnalyzeSkeleton2_.pruneIndex)) {
				// initialize visit flags
				resetVisited();
				// Recalculate analysis over the new image
				bPruneCycles = false;
				processSkeleton(inputImage, null);
			}
		}

		// Calculate triple points (junctions with exactly 3 branches)
		calculateTripleAndQuadruplePoints();

		if (shortPath) {
			if (debug)
				IJ.log("Calculating longest shortest paths...");

			// Copy input image
			shortPathImage = new ImageStack(width, height,
					inputImage.getColorModel());
			for (int i = 1; i <= inputImage.getSize(); i++)
				shortPathImage.addSlice(inputImage.getSliceLabel(i),
						inputImage.getProcessor(i).duplicate());

			shortestPathList = new ArrayList<>();

			shortestPathPoints = new ArrayList[numOfTrees];

			// Visit skeleton and measure distances.
			// and apply Warshall algorithm
			spStartPosition = new double[numOfTrees][3];
			for (int i = 0; i < numOfTrees; i++) {
				shortestPathPoints[i] = new ArrayList<>();
				// Warshall algorithm including tag positions
				shortestPath = warshallAlgorithm(graph[i],
						shortestPathPoints[i]);
				shortestPathList.add(shortestPath);
				spStartPosition[i][0] = spx * imRef.getCalibration().pixelWidth;
				spStartPosition[i][1] = spy * imRef.getCalibration().pixelHeight;
				spStartPosition[i][2] = spz * imRef.getCalibration().pixelDepth;
			}

			if (!silent) {
				// Display short paths in a new stack
				ImagePlus shortIP = new ImagePlus("Longest shortest paths",
						shortPathImage);
				shortIP.show();

				// Set same calibration as the input image
				shortIP.setCalibration(imRef.getCalibration());

				// We apply the Fire LUT and reset the min and max to be between 0-255.
				IJ.run(shortIP, "Fire", null);

				// IJ.resetMinAndMax();
				shortIP.resetDisplayRange();
				shortIP.updateAndDraw();
			}

		}

		// Return the analysis results
		return assembleResults();
	}

	/**
	 * This method is intended for non-interactively using this plugin.
	 * <p>
	 *
	 * @param pruneIndex      The pruneIndex, as asked by the initial gui dialog.
	 * @param thresholdLength maximum length of the branches to prune (all branches
	 *                        below that value are
	 *                        removed)
	 * @param shortPath       flag to calculate the longest shortest path
	 * @param origIP          original input image
	 * @param silent          silent run
	 * @param verbose         flag to display running information
	 */
	public SkeletonResult run(
			int pruneIndex,
			double thresholdLength,
			boolean shortPath,
			ImagePlus origIP,
			boolean silent,
			boolean verbose) {

		AnalyzeSkeleton2_.pruneIndex = pruneIndex;
		this.silent = silent;
		AnalyzeSkeleton2_.pruneEnds = true;
		AnalyzeSkeleton2_.calculateShortestPath = shortPath;
		AnalyzeSkeleton2_.verbose = verbose;

		switch (pruneIndex) {
			// No pruning
			case AnalyzeSkeleton2_.NONE:
				bPruneCycles = false;
				break;
			// Pruning cycles by shortest branch
			case AnalyzeSkeleton2_.SHORTEST_BRANCH:
				bPruneCycles = true;
				break;
			// Pruning cycles by lowest pixel intensity
			case AnalyzeSkeleton2_.LOWEST_INTENSITY_VOXEL:
			case AnalyzeSkeleton2_.LOWEST_INTENSITY_BRANCH:
				// calculate neighborhood size given the calibration
				calculateNeighborhoodOffsets(origIP.getCalibration());
				originalImage = origIP.getStack();
				bPruneCycles = true;
				break;
			default:
		}

		width = imRef.getWidth();
		height = imRef.getHeight();
		depth = imRef.getStackSize();
		inputImage = imRef.getStack();

		// initialize visit flags
		resetVisited();

		// Tag skeleton, differentiate trees and visit them
		processSkeleton(inputImage, null);

		// prune ends
		pruneEndBranches(inputImage, taggedImage, thresholdLength);

		// Prune cycles if necessary
		if (bPruneCycles) {
			if (pruneCycles(inputImage, originalImage,
					AnalyzeSkeleton2_.pruneIndex)) {
				// initialize visit flags
				resetVisited();
				// Recalculate analysis over the new image
				bPruneCycles = false;
				processSkeleton(inputImage, null);
			}
		}

		// Calculate triple points (junctions with exactly 3 branches)
		calculateTripleAndQuadruplePoints();

		if (shortPath) {
			if (debug)
				IJ.log("Calculating longest shortest paths...");

			// Copy input image
			shortPathImage = new ImageStack(width, height,
					inputImage.getColorModel());
			for (int i = 1; i <= inputImage.getSize(); i++)
				shortPathImage.addSlice(inputImage.getSliceLabel(i),
						inputImage.getProcessor(i).duplicate());

			shortestPathList = new ArrayList<>();

			shortestPathPoints = new ArrayList[numOfTrees];

			// Visit skeleton and measure distances.
			// and apply Warshall algorithm
			spStartPosition = new double[numOfTrees][3];
			for (int i = 0; i < numOfTrees; i++) {
				shortestPathPoints[i] = new ArrayList<>();
				// Warshall algorithm including tag positions
				shortestPath = warshallAlgorithm(graph[i],
						shortestPathPoints[i]);
				shortestPathList.add(shortestPath);
				spStartPosition[i][0] = spx * imRef.getCalibration().pixelWidth;
				spStartPosition[i][1] = spy * imRef.getCalibration().pixelHeight;
				spStartPosition[i][2] = spz * imRef.getCalibration().pixelDepth;
			}

			if (!silent) {
				// Display short paths in a new stack
				ImagePlus shortIP = new ImagePlus("Longest shortest paths",
						shortPathImage.duplicate());
				shortIP.show();

				// Set same calibration as the input image
				shortIP.setCalibration(imRef.getCalibration());

				// We apply the Fire LUT and reset the min and max to be between 0-255.
				IJ.run(shortIP, "Fire", null);

				// IJ.resetMinAndMax();
				shortIP.resetDisplayRange();
				shortIP.updateAndDraw();
			}

		}

		// Return the analysis results
		return assembleResults();
	}

	/**
	 * Get the graphs of the current skeletons
	 *
	 * @return array of graphs (one per tree/skeleton)
	 */
	public Graph[] getGraphs() {
		return graph;
	}

	/**
	 * Get the list of points (including junctions and end points) of the largest
	 * shortest paths in the
	 * skeleton image (one per tree).
	 *
	 * @return array with the lists of points of the shortest paths
	 */
	public ArrayList<Point>[] getShortestPathPoints() {
		return shortestPathPoints;
	}

	/**
	 * A simpler standalone running method, for analysis without pruning or showing
	 * images.
	 * <p>
	 * This one just calls run(AnalyzeSkeleton_.NONE, false, null, true, false)
	 */
	public SkeletonResult run() {
		return run(1, 50, false, null, true, false);
		// return run(1, true, false, null, true, false);
	}

	/**
	 * Prune end branches outside the specified ROI
	 *
	 * @param stack       input skeleton image
	 * @param taggedImage tagged skeleton image
	 * @param roi         'protective' ROI: points inside this region are spared
	 *                    from pruning
	 */
	private void pruneEndBranches(ImageStack stack, ImageStack taggedImage,
			Roi roi) {

		if (debug)
			IJ.log("Pruning end-point branches...");
		for (int t = 0; t < numOfTrees; t++) {
			if (debug)
				IJ.log("Pruning tree #" + t);

			Graph g = graph[t];
			ArrayList<Vertex> vertices = g.getVertices();
			ListIterator<Vertex> vit = vertices.listIterator();

			if (debug)
				IJ.log("Initial number of vertices: " +
						graph[t].getVertices().size());

			while (vit.hasNext()) {
				Vertex v = vit.next();
				// Check if the vertex is an end point
				if (v.getBranches().size() == 1 &&
						isEndPoint(v.getPoints().get(0), roi)) {
					if (debug)
						IJ.log("Pruning branch starting at " +
								v.getPoints().get(0));

					// Remove end point voxels
					ArrayList<Point> points = v.getPoints();
					final int nPoints = points.size();

					for (int i = 0; i < nPoints; i++) {
						Point p = points.get(i);
						setPixel(stack, p.x, p.y, p.z, (byte) 0);
						setPixel(taggedImage, p.x, p.y, p.z, (byte) 0);
						numberOfEndPoints[t]--;
						totalNumberOfEndPoints--;
						Iterator<Point> pit = listOfEndPoints.listIterator();
						while (pit.hasNext()) {
							Point ep = pit.next();
							if (ep.equals(p)) {
								pit.remove();
								break;
							}
						}
					}

					// Remove branch voxels
					Edge branch = v.getBranches().get(0);
					points = branch.getSlabs();
					final int nSlabs = points.size();
					for (int i = 0; i < nSlabs; i++) {
						Point p = points.get(i);
						setPixel(stack, p.x, p.y, p.z, (byte) 0);
						setPixel(taggedImage, p.x, p.y, p.z, (byte) 0);
						numberOfSlabs[t]--;
						totalNumberOfSlabs--;
						Iterator<Point> pit = listOfSlabVoxels.listIterator();
						while (pit.hasNext()) {
							Point ep = pit.next();
							if (ep.equals(p)) {
								pit.remove();
								break;
							}
						}
					}

					// remove the Edge from the Graph
					ArrayList<Edge> gEdges = graph[t].getEdges();
					Iterator<Edge> git = gEdges.listIterator();
					while (git.hasNext()) {
						Edge e = git.next();
						if (e.equals(branch)) {
							git.remove();
							break;
						}
					}

					// remove the Edge from the opposite Vertex
					Vertex opp = branch.getOppositeVertex(v);
					ArrayList<Edge> oppBranches = opp.getBranches();
					Iterator<Edge> oppIt = oppBranches.listIterator();
					while (oppIt.hasNext()) {
						Edge oppBranch = oppIt.next();
						if (oppBranch.equals(branch)) {
							oppIt.remove();
							break;
						}
					}

					// remove the Edge from the Vertex
					v.getBranches().remove(0);

					// remove the Vertex from the Graph
					vit.remove();
				}
			}

			if (debug)
				IJ.log("Final number of vertices: " +
						graph[t].getVertices().size());
		}

		return;
	}

	/**
	 * Prune end branches that aren't included in the specified "root" branches
	 *
	 * @param stack       input skeleton image
	 * @param taggedImage tagged skeleton image
	 * @param branches
	 */
	private void pruneEndBranches(ImageStack stack, ImageStack taggedImage,
			ArrayList<Edge> branches) {

		ArrayList<Point> branchPoints = new ArrayList<>();
		for (Edge e : branches) {
			branchPoints.addAll(e.getSlabs());
			branchPoints.addAll(e.getV1().getPoints());
			branchPoints.addAll(e.getV2().getPoints());
		}

		if (debug)
			IJ.log("Pruning end-point branches...");
		for (int t = 0; t < numOfTrees; t++) {
			if (debug)
				IJ.log("Pruning tree #" + t);

			Graph g = graph[t];
			ArrayList<Vertex> vertices = g.getVertices();
			ListIterator<Vertex> vit = vertices.listIterator();

			if (debug)
				IJ.log("Initial number of vertices: " +
						graph[t].getVertices().size());

			while (vit.hasNext()) {
				Vertex v = vit.next();
				// Check if the vertex is an end point
				boolean nonRoot = true;

				ArrayList<Point> vertexPoints = new ArrayList<>();
				vertexPoints.addAll(v.getPoints());
				if (v.getBranches().size() == 1) {
					vertexPoints.addAll(v.getBranches().get(0).getSlabs());
					// vertexPoints.addAll(v.getBranches().get(0).getV1().getPoints());
					// vertexPoints.addAll(v.getBranches().get(0).getV2().getPoints());
				}

				for (Point p : vertexPoints) {
					if (branchPoints.contains(p)) {
						nonRoot = false;
						break;
					}
				}

				// If end branch is NOT a root it can be pruned
				if (v.getBranches().size() == 1 && nonRoot) {
					if (debug)
						IJ.log("Pruning branch starting at " +
								v.getPoints().get(0));

					// Remove end point voxels
					ArrayList<Point> points = v.getPoints();
					final int nPoints = points.size();

					for (int i = 0; i < nPoints; i++) {
						Point p = points.get(i);
						setPixel(stack, p.x, p.y, p.z, (byte) 0);
						setPixel(taggedImage, p.x, p.y, p.z, (byte) 0);
						numberOfEndPoints[t]--;
						totalNumberOfEndPoints--;
						Iterator<Point> pit = listOfEndPoints.listIterator();
						while (pit.hasNext()) {
							Point ep = pit.next();
							if (ep.equals(p)) {
								pit.remove();
								break;
							}
						}
					}

					// Remove branch voxels
					Edge branch = v.getBranches().get(0);
					points = branch.getSlabs();
					final int nSlabs = points.size();
					for (int i = 0; i < nSlabs; i++) {
						Point p = points.get(i);
						setPixel(stack, p.x, p.y, p.z, (byte) 0);
						setPixel(taggedImage, p.x, p.y, p.z, (byte) 0);
						numberOfSlabs[t]--;
						totalNumberOfSlabs--;
						Iterator<Point> pit = listOfSlabVoxels.listIterator();
						while (pit.hasNext()) {
							Point ep = pit.next();
							if (ep.equals(p)) {
								pit.remove();
								break;
							}
						}
					}

					// remove the Edge from the Graph
					ArrayList<Edge> gEdges = graph[t].getEdges();
					Iterator<Edge> git = gEdges.listIterator();
					while (git.hasNext()) {
						Edge e = git.next();
						if (e.equals(branch)) {
							git.remove();
							break;
						}
					}

					// remove the Edge from the opposite Vertex
					Vertex opp = branch.getOppositeVertex(v);
					ArrayList<Edge> oppBranches = opp.getBranches();
					Iterator<Edge> oppIt = oppBranches.listIterator();
					while (oppIt.hasNext()) {
						Edge oppBranch = oppIt.next();
						if (oppBranch.equals(branch)) {
							oppIt.remove();
							break;
						}
					}

					// remove the Edge from the Vertex
					v.getBranches().remove(0);

					// remove the Vertex from the Graph
					vit.remove();
				}
			}

			if (debug)
				IJ.log("Final number of vertices: " +
						graph[t].getVertices().size());
		}

		return;
	}

	/**
	 * Prune end branches of a specific length
	 *
	 * @param stack       input skeleton image
	 * @param taggedImage tagged skeleton image
	 * @param length      limit length to prune the branches (in calibrated units)
	 */
	private void pruneEndBranches(
			ImageStack stack,
			ImageStack taggedImage,
			double length) {

		if (debug)
			IJ.log("Pruning end-point branches...");
		for (int t = 0; t < numOfTrees; t++) {
			if (debug)
				IJ.log("Pruning tree #" + t);

			Graph g = graph[t];
			ArrayList<Vertex> vertices = g.getVertices();
			ListIterator<Vertex> vit = vertices.listIterator();

			if (debug) {
				// System.out.println("Initial number of vertices: " +
				// graph[t].getVertices().size());
				IJ.log("Initial number of vertices: " +
						graph[t].getVertices().size());
			}

			while (vit.hasNext()) {
				Vertex v = vit.next();
				//
				if (v.getBranches().size() == 1 &&
						v.getBranches().get(0).getLength() <= length) {

					// System.out.println(v.getBranches() + "");

					if (debug)
						IJ.log("Pruning branch starting at " +
								v.getPoints().get(0));
					// Remove end point voxels
					ArrayList<Point> points = v.getPoints();
					final int nPoints = points.size();

					for (int i = 0; i < nPoints; i++) {
						Point p = points.get(i);
						setPixel(stack, p.x, p.y, p.z, (byte) 3);
						setPixel(taggedImage, p.x, p.y, p.z, (byte) 3);
						numberOfEndPoints[t]--;
						totalNumberOfEndPoints--;
						Iterator<Point> pit = listOfEndPoints.listIterator();
						while (pit.hasNext()) {
							Point ep = pit.next();
							if (ep.equals(p)) {
								// System.out.println("removed endpoint voxel");
								pit.remove();
								break;
							}
						}
					}

					// Remove branch voxels
					Edge branch = v.getBranches().get(0);
					points = branch.getSlabs();
					final int nSlabs = points.size();
					for (int i = 0; i < nSlabs; i++) {
						Point p = points.get(i);
						setPixel(stack, p.x, p.y, p.z, (byte) 3);
						setPixel(taggedImage, p.x, p.y, p.z, (byte) 3);
						numberOfSlabs[t]--;
						totalNumberOfSlabs--;
						Iterator<Point> pit = listOfSlabVoxels.listIterator();
						while (pit.hasNext()) {
							Point ep = pit.next();
							if (ep.equals(p)) {
								// System.out.println("removed branch voxel");
								pit.remove();
								break;
							}
						}
					}

					// remove the Edge from the Graph
					ArrayList<Edge> gEdges = graph[t].getEdges();
					Iterator<Edge> git = gEdges.listIterator();
					while (git.hasNext()) {
						Edge e = git.next();
						if (e.equals(branch)) {
							// System.out.println("removed edge from graph");
							// System.out.println(e);
							git.remove();
							break;
						}
					}

					// remove the Edge from the opposite Vertex
					Vertex opp = branch.getOppositeVertex(v);
					ArrayList<Edge> oppBranches = opp.getBranches();
					Iterator<Edge> oppIt = oppBranches.listIterator();
					while (oppIt.hasNext()) {
						Edge oppBranch = oppIt.next();
						if (oppBranch.equals(branch)) {
							oppIt.remove();
							break;
						}
					}

					// remove the Edge from the Vertex
					if (!v.getBranches().isEmpty())
						v.getBranches().remove(0);

					// remove the Vertex from the Graph
					vit.remove();
				}
			}

			if (debug) {
				// System.out.println("Final number of vertices: " +
				// graph[t].getVertices().size());
				IJ.log("Final number of vertices: " +
						graph[t].getVertices().size());
			}
		}

		return;
	}

	// ---------------------------------------------------------------------------
	/**
	 * Calculate the neighborhood size based on the calibration of the image.
	 *
	 * @param calibration image calibration
	 */
	private void calculateNeighborhoodOffsets(Calibration calibration) {
		double max = calibration.pixelDepth;
		if (calibration.pixelHeight > max)
			max = calibration.pixelHeight;
		if (calibration.pixelWidth > max)
			max = calibration.pixelWidth;

		x_offset = (int) Math.round(max / calibration.pixelWidth) > 1 ? (int) Math.round(max / calibration.pixelWidth)
				: 1;
		y_offset = (int) Math.round(max / calibration.pixelHeight) > 1 ? (int) Math.round(max / calibration.pixelHeight)
				: 1;
		z_offset = (int) Math.round(max / calibration.pixelDepth) > 1 ? (int) Math.round(max / calibration.pixelDepth)
				: 1;

		if (debug) {
			IJ.log("x_offset = " + x_offset);
			IJ.log("y_offset = " + y_offset);
			IJ.log("z_offset = " + z_offset);
		}

	}// end method calculateNeighborhoodOffsets

	// ---------------------------------------------------------------------------
	/**
	 * Process skeleton: tag image, mark trees and visit.
	 *
	 * @param inputImage2 input skeleton image to process
	 * @param rootPoint   root point
	 */
	public void processSkeleton(ImageStack inputImage2, Point rootPoint) {
		// Initialize global lists of points
		listOfEndPoints = new ArrayList<>();
		listOfJunctionVoxels = new ArrayList<>();
		listOfSlabVoxels = new ArrayList<>();
		listOfStartingSlabVoxels = new ArrayList<>();
		totalNumberOfEndPoints = 0;
		totalNumberOfJunctionVoxels = 0;
		totalNumberOfSlabs = 0;

		// Prepare data: classify voxels and tag them.
		if (rootPoint != null)
			taggedImage = tagImage(inputImage2, rootPoint);
		else
			taggedImage = tagImage(inputImage2);

		// Show tags image.
		if (!bPruneCycles && !silent) {
			displayTagImage(taggedImage);
		}

		// Mark trees
		labeledSkeletons = markTrees(taggedImage);

		if (numOfTrees == 0)
			return;

		// Ask memory for every tree
		initializeTrees();

		// Divide groups of end-points and junction voxels
		if (numOfTrees > 1)
			divideVoxelsByTrees(labeledSkeletons);
		if (numOfTrees == 1) {
			if (debug)
				IJ.log("list of end points size = " + listOfEndPoints.size());
			endPointsTree[0] = listOfEndPoints;
			numberOfEndPoints[0] = listOfEndPoints.size();
			junctionVoxelTree[0] = listOfJunctionVoxels;
			numberOfJunctionVoxels[0] = listOfJunctionVoxels.size();
			startingSlabTree[0] = listOfStartingSlabVoxels;
		}

		// Calculate number of junctions (skipping neighbor junction voxels)
		groupJunctions(labeledSkeletons);

		// Mark all unvisited
		resetVisited();

		// Visit skeleton and measure distances.
		for (int i = 0; i < numOfTrees; i++)
			visitSkeleton(taggedImage, labeledSkeletons, i + 1);

	} // end method processSkeleton

	// -----------------------------------------------------------------------
	/**
	 * Prune cycles from tagged image and update it.
	 *
	 * @param inputImage    input skeleton image
	 * @param originalImage original gray-scale image
	 * @param pruningMode   (SHORTEST_BRANCH, LOWEST_INTENSITY_VOXEL,
	 *                      LOWEST_INTENSITY_BRANCH)
	 * @return true if the input image was pruned or false if there were no cycles
	 */
	private boolean pruneCycles(
			ImageStack inputImage,
			final ImageStack originalImage,
			final int pruningMode) {
		boolean pruned = false;

		for (int iTree = 0; iTree < numOfTrees; iTree++) {
			// For circular trees we just remove one slab
			if (startingSlabTree[iTree].size() == 1) {
				setPixel(inputImage, startingSlabTree[iTree].get(0), (byte) 0);
				pruned = true;
			} else // For the rest, we do depth-first search to detect the cycles
			{
				// DFS
				ArrayList<Edge> backEdges = graph[iTree].depthFirstSearch();

				if (debug) {
					IJ.log(" --------------------------- ");
					final String[] s = new String[] { "UNDEFINED", "TREE",
							"BACK" };
					for (final Edge e : graph[iTree].getEdges()) {
						IJ.log(" edge " + e.getV1().getPoints().get(0) + " - " +
								e.getV2().getPoints().get(0) + " : " +
								s[e.getType() + 1]);
					}
				}

				// If DFS returned backEdges, we need to delete the loops
				if (backEdges.size() > 0) {
					// Find all edges of each loop (backtracking the predecessors)
					for (final Edge e : backEdges) {
						ArrayList<Edge> loopEdges = new ArrayList<>();
						loopEdges.add(e);

						Edge minEdge = e;

						// backtracking (starting at the vertex with higher order index
						final Vertex finalLoopVertex = e.getV1()
								.getVisitOrder() < e.getV2().getVisitOrder() ? e.getV1() : e.getV2();

						Vertex backtrackVertex = e.getV1().getVisitOrder() < e
								.getV2().getVisitOrder() ? e.getV2() : e.getV1();

						// backtrack until reaching final loop vertex
						while (!finalLoopVertex.equals(backtrackVertex)) {
							// Extract predecessor
							final Edge pre = backtrackVertex.getPredecessor();
							// Update shortest loop edge if necessary
							if (pruningMode == AnalyzeSkeleton2_.SHORTEST_BRANCH &&
									pre.getSlabs().size() < minEdge.getSlabs()
											.size())
								minEdge = pre;
							// Add to loop edge list
							loopEdges.add(pre);
							// Extract predecessor
							backtrackVertex = pre.getV1().equals(
									backtrackVertex) ? pre.getV2() : pre.getV1();
						}

						// Prune cycle
						if (pruningMode == AnalyzeSkeleton2_.SHORTEST_BRANCH) {
							// Remove middle slab from the shortest loop edge
							Point removeCoords = null;
							if (minEdge.getSlabs().size() > 0)
								removeCoords = minEdge.getSlabs()
										.get(minEdge.getSlabs().size() / 2);
							else
								removeCoords = minEdge.getV1().getPoints()
										.get(0);
							setPixel(inputImage, removeCoords, (byte) 0);
						} else if (pruningMode == AnalyzeSkeleton2_.LOWEST_INTENSITY_VOXEL) {
							removeLowestIntensityVoxel(loopEdges, inputImage,
									originalImage);
						} else if (pruningMode == AnalyzeSkeleton2_.LOWEST_INTENSITY_BRANCH) {
							cutLowestIntensityBranch(loopEdges, inputImage,
									originalImage);
						}
					} // endfor backEdges

					pruned = true;
				}
			}
		}

		return pruned;
	}// end method pruneCycles

	// -----------------------------------------------------------------------
	/**
	 * Cut the a list of edges in the lowest pixel intensity voxel (calculated from
	 * the original
	 * -grayscale- image).
	 *
	 * @param loopEdges         list of edges to be analyzed
	 * @param inputImage2       input skeleton image
	 * @param originalGrayImage original gray image
	 */
	private void removeLowestIntensityVoxel(
			final ArrayList<Edge> loopEdges,
			ImageStack inputImage2,
			ImageStack originalGrayImage) {
		Point lowestIntensityVoxel = null;

		double lowestIntensityValue = Double.MAX_VALUE;

		for (final Edge e : loopEdges) {
			for (final Point p : e.getSlabs()) {
				final double avg = getAverageNeighborhoodValue(originalGrayImage,
						p,
						x_offset, y_offset, z_offset);
				if (avg < lowestIntensityValue) {
					lowestIntensityValue = avg;
					lowestIntensityVoxel = p;
				}
			}
			// Check vertices
			/*
			 * for(final Point p : e.getV1().getPoints())
			 * {
			 * final double avg = getAverageNeighborhoodValue(originalGrayImage, p,
			 * this.x_offset, this.y_offset, this.z_offset);
			 * if(avg < lowestIntensityValue)
			 * {
			 * lowestIntensityValue = avg;
			 * lowestIntensityVoxel = p;
			 * }
			 * }
			 * for(final Point p : e.getV2().getPoints())
			 * {
			 * final double avg = getAverageNeighborhoodValue(originalGrayImage, p,
			 * this.x_offset, this.y_offset, this.z_offset);
			 * if(avg < lowestIntensityValue)
			 * {
			 * lowestIntensityValue = avg;
			 * lowestIntensityVoxel = p;
			 * }
			 * }
			 */
		}

		// Cut loop in the lowest intensity pixel value position
		if (debug)
			IJ.log("Cut loop at coordinates: " + lowestIntensityVoxel);
		setPixel(inputImage2, lowestIntensityVoxel, (byte) 0);
	}// end method removeLowestIntensityVoxel

	// -----------------------------------------------------------------------
	/**
	 * Cut the a list of edges in the lowest pixel intensity branch.
	 *
	 * @param loopEdges   list of edges to be analyzed
	 * @param inputImage2 input skeleton image
	 */
	private void cutLowestIntensityBranch(
			final ArrayList<Edge> loopEdges,
			ImageStack inputImage2,
			ImageStack originalGrayImage) {
		Edge lowestIntensityEdge = null;

		double lowestIntensityValue = Double.MAX_VALUE;

		Point cutPoint = null;

		for (final Edge e : loopEdges) {
			// Calculate average intensity of the edge neighborhood
			double min_val = Double.MAX_VALUE;
			Point darkestPoint = null;

			double edgeIntensity = 0;
			double n_vox = 0;

			// Check slab points
			for (final Point p : e.getSlabs()) {
				final double avg = getAverageNeighborhoodValue(originalGrayImage,
						p,
						x_offset, y_offset, z_offset);
				// Keep track of the darkest slab point of the edge
				if (avg < min_val) {
					min_val = avg;
					darkestPoint = p;
				}
				edgeIntensity += avg;
				n_vox++;
			}
			// Check vertices
			for (final Point p : e.getV1().getPoints()) {
				edgeIntensity += getAverageNeighborhoodValue(originalGrayImage,
						p,
						x_offset, y_offset, z_offset);
				n_vox++;
			}
			for (final Point p : e.getV2().getPoints()) {
				edgeIntensity += getAverageNeighborhoodValue(originalGrayImage,
						p,
						x_offset, y_offset, z_offset);
				n_vox++;
			}

			if (n_vox != 0)
				edgeIntensity /= n_vox;
			if (debug) {
				IJ.log("Loop edge between " + e.getV1().getPoints().get(0) +
						" and " + e.getV2().getPoints().get(0) + ":");
				IJ.log("avg edge intensity = " + edgeIntensity +
						" darkest slab point = " + darkestPoint.toString());
			}
			// Keep track of the lowest intensity edge
			if (edgeIntensity < lowestIntensityValue) {
				lowestIntensityEdge = e;
				lowestIntensityValue = edgeIntensity;
				cutPoint = darkestPoint;
			}
		}

		// Cut loop in the lowest intensity branch medium position
		Point removeCoords = null;
		if (lowestIntensityEdge.getSlabs().size() > 0)
			removeCoords = cutPoint;
		else {
			IJ.error("Lowest intensity branch without slabs?!: vertex " +
					lowestIntensityEdge.getV1().getPoints().get(0));
			removeCoords = lowestIntensityEdge.getV1().getPoints().get(0);
		}

		if (debug)
			IJ.log("Cut loop at coordinates: " + removeCoords);
		setPixel(inputImage2, removeCoords, (byte) 0);

	}// end method cutLowestIntensityBranch

	// -----------------------------------------------------------------------
	/**
	 * Display tag image on a new window.
	 *
	 * @param taggedImage tag image to be displayed
	 */
	void displayTagImage(ImageStack taggedImage) {
		final int slices = imRef.getNSlices();
		final int frames = imRef.getNFrames();
		final int channels = imRef.getNChannels();

		ImagePlus tagIP = IJ.createHyperStack("Tagged skeleton", width, height,
				channels, slices, frames,
				imRef.getBitDepth());
		tagIP.setStack(taggedImage.duplicate(), channels, slices, frames);
		tagIP.show();

		// Set same calibration as the input image
		tagIP.setCalibration(imRef.getCalibration());

		// We apply the Fire LUT and reset the min and max to be between 0-255.
		IJ.run(tagIP, "Fire", null);

		// IJ.resetMinAndMax();
		tagIP.resetDisplayRange();
		tagIP.updateAndDraw();
	} // end method displayTagImage

	// -----------------------------------------------------------------------
	/**
	 * Divide the end point, junction and special (starting) slab voxels in the
	 * corresponding tree
	 * lists.
	 *
	 * @param treeIS tree image
	 */
	private void divideVoxelsByTrees(ImageStack treeIS) {
		// Add end points to the corresponding tree
		for (int i = 0; i < totalNumberOfEndPoints; i++) {
			final Point p = listOfEndPoints.get(i);
			endPointsTree[(int) getFloatPixel(treeIS, p) - 1].add(p);
		}

		// Add junction voxels to the corresponding tree
		for (int i = 0; i < totalNumberOfJunctionVoxels; i++) {
			final Point p = listOfJunctionVoxels.get(i);
			junctionVoxelTree[(int) getFloatPixel(treeIS, p) - 1].add(p);
		}

		// Add special slab voxels to the corresponding tree
		for (int i = 0; i < listOfStartingSlabVoxels.size(); i++) {
			final Point p = listOfStartingSlabVoxels.get(i);
			startingSlabTree[(int) getFloatPixel(treeIS, p) - 1].add(p);
		}

		// Assign number of end points and junction voxels per tree
		for (int iTree = 0; iTree < numOfTrees; iTree++) {
			numberOfEndPoints[iTree] = endPointsTree[iTree].size();
			numberOfJunctionVoxels[iTree] = junctionVoxelTree[iTree].size();
		}

	} // end divideVoxelsByTrees

	// -----------------------------------------------------------------------
	/** Ask memory for trees. */
	private void initializeTrees() {
		numberOfBranches = new int[numOfTrees];
		numberOfEndPoints = new int[numOfTrees];
		numberOfJunctionVoxels = new int[numOfTrees];
		numberOfJunctions = new int[numOfTrees];
		numberOfSlabs = new int[numOfTrees];
		numberOfTriplePoints = new int[numOfTrees];
		numberOfQuadruplePoints = new int[numOfTrees];
		averageBranchLength = new double[numOfTrees];
		maximumBranchLength = new double[numOfTrees];
		endPointsTree = new ArrayList[numOfTrees];
		junctionVoxelTree = new ArrayList[numOfTrees];
		startingSlabTree = new ArrayList[numOfTrees];
		listOfSingleJunctions = new ArrayList[numOfTrees];

		graph = new Graph[numOfTrees];

		for (int i = 0; i < numOfTrees; i++) {
			endPointsTree[i] = new ArrayList<>();
			junctionVoxelTree[i] = new ArrayList<>();
			startingSlabTree[i] = new ArrayList<>();
			listOfSingleJunctions[i] = new ArrayList<>();
		}
		junctionVertex = new Vertex[numOfTrees][];
	}// end method initializeTrees

	// -----------------------------------------------------------------------
	/** Show results table. */
	private void showResults() {
		final ResultsTable rt = new ResultsTable();

		final String[] head = { "Skeleton", "# Branches", "# Junctions",
				"# End-point voxels",
				"# Junction voxels", "# Slab voxels", "Average Branch Length",
				"# Triple points", "# Quadruple points",
				"Maximum Branch Length",
				"Longest Shortest Path", "spx", "spy", "spz" };

		for (int i = 0; i < numOfTrees; i++) {
			rt.incrementCounter();

			rt.addValue(head[1], numberOfBranches[i]);
			rt.addValue(head[2], numberOfJunctions[i]);
			rt.addValue(head[3], numberOfEndPoints[i]);
			rt.addValue(head[4], numberOfJunctionVoxels[i]);
			rt.addValue(head[5], numberOfSlabs[i]);
			rt.addValue(head[6], averageBranchLength[i]);
			rt.addValue(head[7], numberOfTriplePoints[i]);
			rt.addValue(head[8], numberOfQuadruplePoints[i]);
			rt.addValue(head[9], maximumBranchLength[i]);
			if (null != shortestPathList) {
				rt.addValue(head[10], shortestPathList.get(i));
				rt.addValue(head[11], spStartPosition[i][0]);
				rt.addValue(head[12], spStartPosition[i][1]);
				rt.addValue(head[13], spStartPosition[i][2]);
			}

			if (0 == i % 100)
				rt.show("Results");
		}
		rt.show("Results");

		// Extra information
		if (AnalyzeSkeleton2_.verbose) {
			// New results table
			final ResultsTable extra_rt = new ResultsTable();

			final String[] extra_head = { "Branch", "Skeleton ID",
					"Branch length", "V1 x", "V1 y",
					"V1 z", "V2 x", "V2 y", "V2 z", "Euclidean distance",
					"running average length", "average intensity (inner 3rd)",
					"average intensity" };

			// Edge comparator (by branch length)
			Comparator<Edge> comp = new Comparator<Edge>() {
				@Override
				public int compare(Edge o1, Edge o2) {
					final double diff = o1.getLength() - o2.getLength();
					if (diff < 0)
						return 1;
					else if (diff == 0)
						return 0;
					else
						return -1;
				}

				@Override
				public boolean equals(Object o) {
					return false;
				}
			};
			// Display branch information for each tree
			for (int i = 0; i < numOfTrees; i++) {
				final ArrayList<Edge> listEdges = graph[i].getEdges();
				// Sort branches by length
				Collections.sort(listEdges, comp);
				for (final Edge e : listEdges) {
					extra_rt.incrementCounter();
					extra_rt.addValue(extra_head[1], i + 1);
					extra_rt.addValue(extra_head[2], e.getLength());
					extra_rt.addValue(extra_head[3],
							e.getV1().getPoints().get(0).x *
									imRef.getCalibration().pixelWidth);
					extra_rt.addValue(extra_head[4],
							e.getV1().getPoints().get(0).y *
									imRef.getCalibration().pixelHeight);
					extra_rt.addValue(extra_head[5],
							e.getV1().getPoints().get(0).z *
									imRef.getCalibration().pixelDepth);
					extra_rt.addValue(extra_head[6],
							e.getV2().getPoints().get(0).x *
									imRef.getCalibration().pixelWidth);
					extra_rt.addValue(extra_head[7],
							e.getV2().getPoints().get(0).y *
									imRef.getCalibration().pixelHeight);
					extra_rt.addValue(extra_head[8],
							e.getV2().getPoints().get(0).z *
									imRef.getCalibration().pixelDepth);
					extra_rt.addValue(extra_head[9],
							this.calculateDistance(e.getV1().getPoints().get(0),
									e.getV2().getPoints().get(0)));
					extra_rt.addValue(extra_head[10], e.getLength_ra());
					extra_rt.addValue(extra_head[11], e.getColor3rd());
					extra_rt.addValue(extra_head[12], e.getColor());
				}

			}
			extra_rt.show("Branch information");
		}

	}// end method showResults

	/**
	 * Returns one of the two result images in an ImageStack object.
	 *
	 * @param longestShortestPath Get the tagged longest shortest paths instead of
	 *                            the standard tagged
	 *                            image
	 *
	 * @return The results image with a tagged skeleton
	 */
	public ImageStack getResultImage(boolean longestShortestPath) {
		if (longestShortestPath) {
			return shortPathImage;
		}
		return taggedImage;
	}

	/**
	 * Returns the analysis results in a SkeletonResult object.
	 * <p>
	 *
	 * @return The results of the skeleton analysis.
	 */
	protected SkeletonResult assembleResults() {
		SkeletonResult result = new SkeletonResult(numOfTrees);
		result.setBranches(numberOfBranches);
		result.setJunctions(numberOfJunctions);
		result.setEndPoints(numberOfEndPoints);
		result.setJunctionVoxels(numberOfJunctionVoxels);
		result.setSlabs(numberOfSlabs);
		result.setAverageBranchLength(averageBranchLength);
		result.setTriples(numberOfTriplePoints);
		result.setQuadruples(numberOfQuadruplePoints);
		result.setMaximumBranchLength(maximumBranchLength);

		result.setListOfEndPoints(listOfEndPoints);
		result.setListOfJunctionVoxels(listOfJunctionVoxels);
		result.setListOfSlabVoxels(listOfSlabVoxels);
		result.setListOfStartingSlabVoxels(listOfStartingSlabVoxels);

		result.setShortestPathList(shortestPathList);
		result.setSpStartPosition(spStartPosition);

		result.setGraph(graph);

		result.calculateNumberOfVoxels();

		return result;
	}

	// -----------------------------------------------------------------------
	/**
	 * Visit skeleton starting at end-points, junctions and slab of circular
	 * skeletons, and record
	 * measurements.
	 *
	 * @param taggedImage tag skeleton image
	 * @param treeImage   skeleton image with tree classification
	 * @param currentTree number of the tree to be visited
	 */
	private void visitSkeleton(ImageStack taggedImage, ImageStack treeImage,
			int currentTree) {
		// tree index
		final int iTree = currentTree - 1;

		if (debug) {
			// Junction vertices in the tree
			IJ.log("this.junctionVertex[" + iTree + "].length = " +
					junctionVertex[iTree].length);
			for (int i = 0; i < junctionVertex[iTree].length; i++) {
				IJ.log(" vertices points: " + junctionVertex[iTree][i]);
			}
		}

		// Create new graph
		graph[iTree] = new Graph();
		// Add all junction vertices
		for (int i = 0; i < junctionVertex[iTree].length; i++)
			graph[iTree].addVertex(junctionVertex[iTree][i]);

		if (debug)
			IJ.log(" Analyzing tree number " + currentTree);
		// length of branches
		double branchLength = 0;

		maximumBranchLength[iTree] = 0;
		numberOfSlabs[iTree] = 0;

		// Visit branches starting at end points
		for (int i = 0; i < numberOfEndPoints[iTree]; i++) {
			final Point endPointCoord = endPointsTree[iTree].get(i);

			if (debug)
				IJ.log(
						"\n*** visit from end point: " + endPointCoord + " *** ");

			// Skip when visited
			if (isVisited(endPointCoord)) {
				// if(this.initialPoint[iTree] == null)
				// IJ.error("WEIRD:" + " (" + endPointCoord.x + ", " + endPointCoord.y + ", " +
				// endPointCoord.z +
				// ")");
				if (debug)
					IJ.log("visited = (" + endPointCoord.x + ", " +
							endPointCoord.y + ", " + endPointCoord.z + ")");
				continue;
			}

			// Initial vertex
			Vertex v1 = new Vertex();
			v1.addPoint(endPointCoord);
			graph[iTree].addVertex(v1);
			if (i == 0)
				graph[iTree].setRoot(v1);

			// slab list for the edge
			slabList = new ArrayList<>();

			// Otherwise, visit branch until next junction or end point.
			double[] properties = visitBranch(endPointCoord, iTree);
			double length = properties[0];
			double color3rd = properties[1];
			double color = properties[2];
			double length_ra = properties[3];

			// If length is 0, it means the tree is formed by only one voxel.
			if (length == 0) {
				// If there is an adjacent visited junction, count it
				// as a single voxel branch
				final Point aux = getVisitedJunctionNeighbor(endPointCoord, v1);
				if (null != aux) {
					auxFinalVertex = findPointVertex(junctionVertex[iTree], aux);
					length += calculateDistance(endPointCoord, aux);

					// Add the length to the first point of the vertex (to prevent later from having
					// euclidean distances larger than the actual distance)
					length += calculateDistance(
							auxFinalVertex.getPoints().get(0), endPointCoord);
					// Add branch to graph
					if (debug)
						IJ.log("adding branch from " + v1.getPoints().get(0) +
								" to " + auxFinalVertex.getPoints().get(0));
					graph[iTree].addVertex(auxFinalVertex);
					graph[iTree].addEdge(new Edge(v1, auxFinalVertex, slabList,
							length, color3rd, color, length_ra));
					// increase number of branches
					numberOfBranches[iTree]++;

					if (debug)
						IJ.log(
								"increased number of branches, length = " + length);

					branchLength += length;
				} else if (debug)
					IJ.log("set initial point to final point");
				continue;
			}

			// If the final point is a slab, then we add the path to the
			// neighbor junction voxel not belonging to the initial vertex
			// (unless it is a self loop)
			if (isSlab(auxPoint)) {
				final Point aux = auxPoint;
				// IJ.log("Looking for " + this.auxPoint + " in the list of vertices...");
				auxPoint = getVisitedJunctionNeighbor(auxPoint, v1);
				auxFinalVertex = findPointVertex(junctionVertex[iTree],
						auxPoint);
				if (auxPoint == null) {
					// IJ.log("Point "+ aux + " has not neighbor end junction! (inner loop)");
					// Inner loop
					auxFinalVertex = v1;
					auxPoint = aux;
				}
				length += calculateDistance(auxPoint, aux);

				// Add the length to the first point of the vertex (to prevent later from having
				// euclidean distances larger than the actual distance)
				length += calculateDistance(auxFinalVertex.getPoints().get(0),
						auxPoint);
			}

			// Add branch to graph
			if (debug)
				IJ.log("adding branch from " + v1.getPoints().get(0) + " to " +
						auxFinalVertex.getPoints().get(0) + ", aux point = " +
						auxPoint);
			graph[iTree].addVertex(auxFinalVertex);
			graph[iTree].addEdge(new Edge(v1, auxFinalVertex, slabList, length,
					color3rd, color, length_ra));

			// increase number of branches
			numberOfBranches[iTree]++;

			if (debug)
				IJ.log("increased number of branches, length = " + length);

			branchLength += length;

			// update maximum branch length
			if (length > maximumBranchLength[iTree]) {
				maximumBranchLength[iTree] = length;
			}
		}

		// If there is no end points, set the first junction as root.
		if (numberOfEndPoints[iTree] == 0 &&
				junctionVoxelTree[iTree].size() > 0)
			graph[iTree].setRoot(junctionVertex[iTree][0]);

		if (debug)
			IJ.log(" --------------------------- ");

		// Now visit branches starting at junctions
		// 08/26/2009 Changed the loop to visit first the junction voxels that are
		// forming a single junction.
		for (int i = 0; i < junctionVertex[iTree].length; i++) {
			for (int j = 0; j < junctionVertex[iTree][i].getPoints()
					.size(); j++) {
				final Point junctionCoord = junctionVertex[iTree][i].getPoints()
						.get(j);

				if (debug)
					IJ.log(
							"\n*** visit from junction " + junctionCoord + " *** ");

				// Mark junction as visited
				setVisited(junctionCoord, true);

				Point nextPoint = getNextUnvisitedVoxel(junctionCoord);

				while (nextPoint != null) {
					// Do not count adjacent junctions
					if (!isJunction(nextPoint)) {

						if (debug)
							IJ.log("visiting " + nextPoint);

						// Create graph edge
						slabList = new ArrayList<>();
						slabList.add(nextPoint);
						numberOfSlabs[iTree]++;

						// Calculate distance from junction to that point
						double length = calculateDistance(junctionCoord,
								nextPoint);

						// Visit branch
						auxPoint = null;

						double[] properties = visitBranch(nextPoint, iTree);
						length += properties[0];
						double color3rd = properties[1];
						double color = properties[2];
						double length_ra = properties[3];

						// Increase total length of branches
						branchLength += length;

						// Increase number of branches
						if (length != 0) {
							if (auxPoint == null)
								auxPoint = nextPoint;

							numberOfBranches[iTree]++;

							// Initial vertex
							Vertex initialVertex = null;
							for (int k = 0; k < junctionVertex[iTree].length; k++)
								if (junctionVertex[iTree][k]
										.isVertexPoint(junctionCoord)) {
									initialVertex = junctionVertex[iTree][k];
									break;
								}

							// If the final point is a slab, then we add the path to the
							// neighbor junction voxel not belonging to the initial vertex
							// (unless it is a self loop)
							if (isSlab(auxPoint)) {
								final Point aux = auxPoint;
								// IJ.log("Looking for " + this.auxPoint + " in the list of vertices...");
								auxPoint = getVisitedJunctionNeighbor(auxPoint,
										initialVertex);
								auxFinalVertex = findPointVertex(
										junctionVertex[iTree], auxPoint);
								if (auxPoint == null) {
									// IJ.log("Point "+ aux + " has not neighbor end junction! (inner loop)");
									// Inner loop
									auxFinalVertex = initialVertex;
									auxPoint = aux;
								}
								length += calculateDistance(auxPoint, aux);
							}

							if (debug)
								IJ.log(
										"increased number of branches, length = " +
												length + " (last point = " + auxPoint +
												")");
							// update maximum branch length
							if (length > maximumBranchLength[iTree]) {
								maximumBranchLength[iTree] = length;
							}

							// Add the distance between the main vertex of the junction
							// and the initial junction vertex of the branch (this prevents from
							// having branches in the graph larger than the calculated branch length)
							length += calculateDistance(
									initialVertex.getPoints().get(0),
									junctionCoord);

							// Create graph branch
							// Add branch to graph
							if (debug)
								IJ.log("adding branch from " +
										initialVertex.getPoints().get(0) + " to " +
										auxFinalVertex.getPoints().get(0));
							graph[iTree].addEdge(new Edge(initialVertex,
									auxFinalVertex, slabList, length, color3rd,
									color, length_ra));
						}
					} else
						setVisited(nextPoint, true);

					nextPoint = getNextUnvisitedVoxel(junctionCoord);
				}
			}
		}

		if (debug)
			IJ.log(" --------------------------- ");

		// Finally visit branches starting at slabs (special case for circular trees)
		if (startingSlabTree[iTree].size() == 1) {
			if (debug)
				IJ.log("visit from slabs");

			final Point startCoord = startingSlabTree[iTree].get(0);

			// Create circular graph (only one vertex)
			final Vertex v1 = new Vertex();
			v1.addPoint(startCoord);
			graph[iTree].addVertex(v1);

			slabList = new ArrayList<>();
			slabList.add(startCoord);

			numberOfSlabs[iTree]++;

			// visit branch until finding visited voxel.
			double[] properties = visitBranch(startCoord, iTree);
			final double length = properties[0];
			double color3rd = properties[1];
			double color = properties[2];
			double length_ra = properties[3];

			if (length != 0) {
				// increase number of branches
				numberOfBranches[iTree]++;
				branchLength += length;

				// update maximum branch length
				if (length > maximumBranchLength[iTree]) {
					maximumBranchLength[iTree] = length;
				}
			}

			// Create circular edge
			graph[iTree].addEdge(
					new Edge(v1, v1, slabList, length, color3rd, color, length_ra));
		}

		if (debug)
			IJ.log(" --------------------------- ");

		if (numberOfBranches[iTree] == 0)
			return;
		// Average length
		averageBranchLength[iTree] = branchLength / numberOfBranches[iTree];

		if (debug) {
			IJ.log("Num of vertices = " + graph[iTree].getVertices().size() +
					" num of edges = " + graph[iTree].getEdges().size());
			for (int i = 0; i < graph[iTree].getVertices().size(); i++) {
				Vertex v = graph[iTree].getVertices().get(i);
				IJ.log(" vertex " + v.getPoints().get(0) + " has neighbors: ");
				for (int j = 0; j < v.getBranches().size(); j++) {
					final Vertex v1 = v.getBranches().get(j).getV1();
					final Vertex oppositeVertex = v1.equals(v) ? v.getBranches().get(j).getV2() : v1;
					IJ.log(j + ": " + oppositeVertex.getPoints().get(0));
				}

			}

			IJ.log(" --------------------------- ");
			for (int i = 0; i < junctionVertex[iTree].length; i++) {
				IJ.log("Junction #" + i + " is formed by: ");
				for (int j = 0; j < junctionVertex[iTree][i].getPoints()
						.size(); j++)
					IJ.log(
							j + ": " + junctionVertex[iTree][i].getPoints().get(j));
			}
		}

	} // end visitSkeleton

	/* ----------------------------------------------------------------------- */
	/**
	 * Color the different trees in the skeleton.
	 *
	 * @param taggedImage
	 *
	 * @return image with every tree tagged with a different number
	 */
	private ImageStack markTrees(ImageStack taggedImage) {
		if (debug)
			IJ.log("=== Mark Trees ===");
		// Create output image
		ImageStack outputImage = new ImageStack(width, height);
		for (int z = 0; z < depth; z++) {
			outputImage.addSlice(taggedImage.getSliceLabel(z + 1),
					new FloatProcessor(width, height));
		}

		numOfTrees = 0;

		int color = 0;

		// Visit trees starting at end points
		for (int i = 0; i < totalNumberOfEndPoints; i++) {
			Point endPointCoord = listOfEndPoints.get(i);

			if (isVisited(endPointCoord))
				continue;

			color++;

			if (color == Integer.MAX_VALUE) {
				IJ.error("More than " + (Integer.MAX_VALUE - 1) +
						" skeletons in the image. AnalyzeSkeleton can only process up to " +
						(Integer.MAX_VALUE - 1));
				return null;
			}

			if (debug)
				IJ.log("-- Visit tree from end-point:");
			// Visit the entire tree.
			int numOfVoxelsInTree = visitTree(endPointCoord, outputImage, color);

			// increase number of trees
			numOfTrees++;
		}

		// Visit trees starting at junction points
		// (some circular trees do not have end points)
		// Visit trees starting at end points
		for (int i = 0; i < totalNumberOfJunctionVoxels; i++) {
			Point junctionCoord = listOfJunctionVoxels.get(i);
			if (isVisited(junctionCoord))
				continue;

			color++;

			if (color == Short.MAX_VALUE) {
				IJ.error("More than " + (Short.MAX_VALUE - 1) +
						" skeletons in the image. AnalyzeSkeleton can only process up to 255");
				return null;
			}

			if (debug)
				IJ.log("-- Visit tree from junction:");

			// else, visit branch until next junction or end point.
			int length = visitTree(junctionCoord, outputImage, color);

			if (length == 0) {
				color--; // the color was not used
				continue;
			}

			// increase number of trees
			numOfTrees++;
		}

		// Check for unvisited slab voxels
		// (just in case there are circular trees without junctions)
		for (int i = 0; i < listOfSlabVoxels.size(); i++) {
			Point p = listOfSlabVoxels.get(i);
			if (isVisited(p) == false) {
				// Mark that voxel as the start point of the circular skeleton
				listOfStartingSlabVoxels.add(p);

				color++;

				if (color == Short.MAX_VALUE) {
					IJ.error("More than " + (Short.MAX_VALUE - 1) +
							" skeletons in the image. AnalyzeSkeleton can only process up to 255");
					return null;
				}

				if (debug)
					IJ.log("-- Visit tree from slab:");

				// else, visit branch until next junction or end point.
				int length = visitTree(p, outputImage, color);

				if (length == 0) {
					color--; // the color was not used
					continue;
				}

				// increase number of trees
				numOfTrees++;
			}
		}

		// System.out.println("Number of trees = " + this.numOfTrees);

		// Show tree image.
		if (debug) {
			ImagePlus treesIP = new ImagePlus("Trees skeleton", outputImage);
			treesIP.show();

			// Set same calibration as the input image
			treesIP.setCalibration(imRef.getCalibration());

			// We apply the Fire LUT and reset the min and max to be between 0-255.
			IJ.run("Fire");

			// IJ.resetMinAndMax();
			treesIP.resetDisplayRange();
			treesIP.updateAndDraw();
		}

		// Reset visited variable
		resetVisited();

		// IJ.log("Number of trees: " + this.numOfTrees + ", # colors = " + color);

		return outputImage;

	} /* end markTrees */

	// --------------------------------------------------------------
	/**
	 * Visit tree marking the voxels with a reference tree color.
	 *
	 * @param startingPoint starting tree point
	 * @param outputImage   3D image to visit
	 * @param color         reference tree color
	 * @return number of voxels in the tree
	 */
	private int visitTree(Point startingPoint, ImageStack outputImage,
			int color) {
		int numOfVoxels = 0;

		if (debug)
			IJ.log("visiting " + startingPoint + " color = " + color);

		if (isVisited(startingPoint))
			return 0;

		// Set pixel color
		this.setPixel(outputImage, startingPoint.x, startingPoint.y,
				startingPoint.z, color);
		setVisited(startingPoint, true);

		ArrayList<Point> toRevisit = new ArrayList<>();

		// Add starting point to revisit list if it is a junction
		if (isJunction(startingPoint))
			toRevisit.add(startingPoint);

		Point nextPoint = getNextUnvisitedVoxel(startingPoint);

		while (nextPoint != null || toRevisit.size() != 0) {
			if (nextPoint != null) {
				if (!isVisited(nextPoint)) {
					numOfVoxels++;
					if (debug)
						IJ.log("visiting " + nextPoint + " color = " + color);

					// Set color and visit flat
					this.setPixel(outputImage, nextPoint.x, nextPoint.y,
							nextPoint.z, color);
					setVisited(nextPoint, true);

					// If it is a junction, add it to the revisit list
					if (isJunction(nextPoint))
						toRevisit.add(nextPoint);

					// Calculate next point to visit
					nextPoint = getNextUnvisitedVoxel(nextPoint);
				}
			} else // revisit list
			{
				nextPoint = toRevisit.get(0);
				if (debug)
					IJ.log("visiting " + nextPoint + " color = " + color);

				// Calculate next point to visit
				nextPoint = getNextUnvisitedVoxel(nextPoint);
				// Maintain junction in the list until there is no more branches
				if (nextPoint == null)
					toRevisit.remove(0);
			}
		}

		return numOfVoxels;
	} // end method visitTree

	// -----------------------------------------------------------------------
	/**
	 * Visit a branch and calculate length.
	 *
	 * @param startingPoint starting coordinates
	 * @return branch length
	 *
	 * @deprecated
	 */
	@Deprecated
	private double visitBranch(Point startingPoint) {
		double length = 0;

		// mark starting point as visited
		setVisited(startingPoint, true);

		// Get next unvisited voxel
		Point nextPoint = getNextUnvisitedVoxel(startingPoint);

		if (nextPoint == null)
			return 0;

		Point previousPoint = startingPoint;

		// We visit the branch until we find an end point or a junction
		while (nextPoint != null && isSlab(nextPoint)) {
			// Add length
			length += calculateDistance(previousPoint, nextPoint);

			// Mark as visited
			setVisited(nextPoint, true);

			// Move in the graph
			previousPoint = nextPoint;
			nextPoint = getNextUnvisitedVoxel(previousPoint);
		}

		if (nextPoint != null) {
			// Add distance to last point
			length += calculateDistance(previousPoint, nextPoint);

			// Mark last point as visited
			setVisited(nextPoint, true);
		}

		auxPoint = previousPoint;

		return length;
	}// end visitBranch

	// -----------------------------------------------------------------------
	/**
	 * Visit a branch and calculate length in a specific tree
	 *
	 * @param startingPoint starting coordinates
	 * @param iTree         tree index
	 * @return branch length and color
	 */
	private double[] visitBranch(Point startingPoint, int iTree) {
		// IJ.log("startingPoint = (" + startingPoint.x + ", " + startingPoint.y + ", "
		// + startingPoint.z +
		// ")");
		double length = 0;
		double intensity = 0.0;
		double intensity3rd = 0.0;
		double length_ra = 0.0;
		double[] ret = new double[4];

		List<Point> pointHistory = new ArrayList<>(0);
		pointHistory.add(0, startingPoint);

		// mark starting point as visited
		setVisited(startingPoint, true);

		// Get next unvisited voxel
		Point nextPoint = getNextUnvisitedVoxel(startingPoint);

		if (nextPoint == null)
			return ret;

		Point previousPoint = startingPoint;

		// We visit the branch until we find an end point or a junction
		while (nextPoint != null && isSlab(nextPoint)) {
			numberOfSlabs[iTree]++;

			// Add slab voxel to the edge
			slabList.add(nextPoint);

			// Add length
			length += calculateDistance(previousPoint, nextPoint);
			pointHistory.add(0, nextPoint);

			length_ra += calculateDistance(pointHistory);

			// Mark as visited
			setVisited(nextPoint, true);

			// Move in the graph
			previousPoint = nextPoint;
			if (debug)
				IJ.log("visiting " + previousPoint);
			nextPoint = getNextUnvisitedVoxel(previousPoint);
		}

		// If we find an unvisited end-point or junction, we set it
		// as final vertex of the branch
		if (nextPoint != null) {
			// Add distance to last point
			length += calculateDistance(previousPoint, nextPoint);
			pointHistory.add(0, nextPoint);
			length_ra += calculateDistance(pointHistory);

			// Mark last point as visited
			setVisited(nextPoint, true);

			// Mark final vertex
			if (isEndPoint(nextPoint)) {
				if (debug)
					IJ.log("found unvisited end point: " + nextPoint);
				auxFinalVertex = new Vertex();
				auxFinalVertex.addPoint(nextPoint);
			} else if (isJunction(nextPoint)) {
				if (debug)
					IJ.log("found unvisited junction point: " + nextPoint);
				auxFinalVertex = findPointVertex(junctionVertex[iTree],
						nextPoint);
				// Add the length to the first point of the vertex (to prevent later from having
				// euclidean distances larger than the actual distance)
				length += calculateDistance(auxFinalVertex.getPoints().get(0),
						nextPoint);
				length_ra += calculateDistance(auxFinalVertex.getPoints().get(0),
						nextPoint);

			}

			auxPoint = nextPoint;
		} else
			auxPoint = previousPoint;

		// calculate average intensity (thickness) value, but only take the inner third
		// of a branch.
		// at both ends the intensity (thickness) is most likely affected by junctions.

		int size = pointHistory.size();
		int start = (int) (size / 3.0);
		int end = (int) (2 * size / 3.0);

		for (int i = 0; i < size; i++) {
			int value = getPixel(inputImage, pointHistory.get(i));
			if (value < 0) {
				value += 256;
			}
			intensity += value;
			if (i >= start && i < end) {
				intensity3rd += value;
			}
		}

		intensity /= size;
		intensity3rd /= end - start;

		ret[0] = length;
		ret[1] = intensity3rd;
		ret[2] = intensity;
		ret[3] = length_ra;

		return ret;
	} // end visitBranch

	// -----------------------------------------------------------------------
	/**
	 * Find vertex in an array given a specific vertex point.
	 *
	 * @param vertex array of search
	 * @param p      vertex point
	 * @return vertex containing that point
	 */
	public Vertex findPointVertex(Vertex[] vertex, Point p) {
		int j = 0;
		for (j = 0; j < vertex.length; j++)
			if (vertex[j].isVertexPoint(p)) {
				if (debug)
					IJ.log(" " + p + " belongs to junction " +
							vertex[j].getPoints().get(0));
				return vertex[j];
			}
		if (debug)
			IJ.log("point " + p +
					" was not found in vertex list! (vertex.length= " +
					vertex.length + ")");
		return null;
	}

	// -----------------------------------------------------------------------
	/**
	 * Calculate Euclidean distance between two points in 3D.
	 *
	 * @param point1 first point coordinates
	 * @param point2 second point coordinates
	 * @return distance (in the corresponding units)
	 */
	private double calculateDistance(Point point1, Point point2) {
		return Math.sqrt(Math
				.pow((point1.x - point2.x) * imRef.getCalibration().pixelWidth, 2) +
				Math.pow((point1.y - point2.y) * imRef.getCalibration().pixelHeight,
						2)
				+
				Math.pow((point1.z - point2.z) * imRef.getCalibration().pixelDepth,
						2));
	}

	// -----------------------------------------------------------------------
	/**
	 * Calculate linear corrected distance between two points in 3D. Uses the
	 * PointHistory of a branch
	 * and its 5 last points.
	 *
	 * @param Points - the last visited Points (most recent has Index 0)
	 * @return linear corrected distance between the last two Points (in the
	 *         corresponding units)
	 */
	private double calculateDistance(List<Point> Points) {
		int indexOfLast = Points.size() - 1;

		// no Distance to be calculated here...
		if (indexOfLast < 1)
			return 0;

		// Point Of Interest
		int poi = 5;

		// poi is indexOflast if List is shorter than 5
		if (indexOfLast < 5) {
			poi = indexOfLast;
		}
		return Math.sqrt(Math
				.pow((Points.get(poi).x - Points.get(0).x) *
						imRef.getCalibration().pixelWidth, 2)
				+
				Math.pow((Points.get(poi).y - Points.get(0).y) *
						imRef.getCalibration().pixelHeight, 2)
				+
				Math.pow((Points.get(poi).z - Points.get(0).z) *
						imRef.getCalibration().pixelDepth, 2))
				/
				poi;
	}

	// -----------------------------------------------------------------------
	/**
	 * Calculate number of junction skipping neighbor junction voxels
	 *
	 * @param treeIS tree stack
	 */
	private void groupJunctions(ImageStack treeIS) {
		// Mark all unvisited
		resetVisited();

		for (int iTree = 0; iTree < numOfTrees; iTree++) {
			// Visit list of junction voxels
			for (int i = 0; i < numberOfJunctionVoxels[iTree]; i++) {
				Point pi = junctionVoxelTree[iTree].get(i);

				if (!isVisited(pi))
					fusionNeighborJunction(pi, listOfSingleJunctions[iTree]);
			}
		}

		// Count number of single junctions for every tree in the image
		for (int iTree = 0; iTree < numOfTrees; iTree++) {
			if (debug)
				IJ.log("this.listOfSingleJunctions[" + iTree + "].size() = " +
						listOfSingleJunctions[iTree].size());

			numberOfJunctions[iTree] = listOfSingleJunctions[iTree].size();

			// Create array of junction vertices for the graph
			junctionVertex[iTree] = new Vertex[listOfSingleJunctions[iTree]
					.size()];

			for (int j = 0; j < listOfSingleJunctions[iTree].size(); j++) {
				final ArrayList<Point> list = listOfSingleJunctions[iTree]
						.get(j);
				junctionVertex[iTree][j] = new Vertex();
				for (final Point p : list)
					junctionVertex[iTree][j].addPoint(p);

			}
		}

		// Mark all unvisited
		resetVisited();
	}

	// -----------------------------------------------------------------------
	/** Reset visit variable and set it to false. */
	private void resetVisited() {
		// Reset visited variable
		visited = null;
		visited = new boolean[width][height][depth];
		/*
		 * for(int i = 0; i < this.width; i ++)
		 * for(int j = 0; j < this.height; j++)
		 * for(int k = 0; k < this.depth; k++)
		 * this.visited[i][j][k] = false;
		 */
	}

	// -----------------------------------------------------------------------
	/**
	 * Fusion neighbor junctions voxels into the same list.
	 *
	 * @param startingPoint       starting junction voxel
	 * @param singleJunctionsList list of single junctions
	 */
	private void fusionNeighborJunction(Point startingPoint,
			ArrayList<ArrayList<Point>> singleJunctionsList) {
		// Create new group of junctions
		ArrayList<Point> newGroup = new ArrayList<>();
		newGroup.add(startingPoint);

		// Mark the starting junction as visited
		setVisited(startingPoint, true);

		// Look for neighbor junctions and add them to the new group
		ArrayList<Point> toRevisit = new ArrayList<>();
		toRevisit.add(startingPoint);

		Point nextPoint = getNextUnvisitedJunctionVoxel(startingPoint);

		while (nextPoint != null || toRevisit.size() != 0) {
			if (nextPoint != null && !isVisited(nextPoint)) {
				// Add to the group
				newGroup.add(nextPoint);
				// Mark as visited
				setVisited(nextPoint, true);

				// add it to the revisit list
				toRevisit.add(nextPoint);

				// Calculate next junction point to visit
				nextPoint = getNextUnvisitedJunctionVoxel(nextPoint);
			} else // revisit list
			{
				nextPoint = toRevisit.get(0);
				// IJ.log("visiting " + nextPoint + " color = " + color);

				// Calculate next point to visit
				nextPoint = getNextUnvisitedJunctionVoxel(nextPoint);
				// Maintain junction in the list until there is no more branches
				if (nextPoint == null)
					toRevisit.remove(0);
			}
		}

		// Add group to the single junction list
		singleJunctionsList.add(newGroup);

	}// end method fusionNeighborJunction

	// -----------------------------------------------------------------------
	/**
	 * Check if two groups of voxels are neighbors.
	 *
	 * @param g1 first group
	 * @param g2 second group
	 *
	 * @return true if the groups have any neighbor voxel
	 */
	boolean checkNeighborGroups(ArrayList<Point> g1, ArrayList<Point> g2) {
		for (int i = 0; i < g1.size(); i++) {
			Point pi = g1.get(i);
			for (int j = 0; j < g2.size(); j++) {
				Point pj = g2.get(j);
				if (isNeighbor(pi, pj))
					return true;
			}
		}
		return false;
	}

	// -----------------------------------------------------------------------
	/**
	 * Calculate number of triple and quadruple points in the skeleton. Triple and
	 * quadruple points are
	 * junctions with exactly 3 and 4 branches respectively.
	 */
	private void calculateTripleAndQuadruplePoints() {
		for (int iTree = 0; iTree < numOfTrees; iTree++) {
			// Visit the groups of junction voxels
			for (int i = 0; i < numberOfJunctions[iTree]; i++) {

				ArrayList<Point> groupOfJunctions = listOfSingleJunctions[iTree]
						.get(i);

				// Count the number of slab and end-points neighbors of every voxel in the group
				int nBranch = 0;
				for (int j = 0; j < groupOfJunctions.size(); j++) {
					Point pj = groupOfJunctions.get(j);

					// Get neighbors and check the slabs or end-points
					byte[] neighborhood = this.getNeighborhood(taggedImage, pj.x,
							pj.y, pj.z);
					for (int k = 0; k < 27; k++)
						if (neighborhood[k] == AnalyzeSkeleton2_.SLAB ||
								neighborhood[k] == AnalyzeSkeleton2_.END_POINT)
							nBranch++;
				}
				// If the junction has only 3 slab/end-point neighbors, then it is a triple
				// point
				if (nBranch == 3)
					numberOfTriplePoints[iTree]++;
				else if (nBranch == 4) // quadruple point if 4
					numberOfQuadruplePoints[iTree]++;
			}

		}

	}// end calculateTripleAndQuadruplePoints

	/* ----------------------------------------------------------------------- */
	/**
	 * Calculate if two points are neighbors.
	 *
	 * @param point1 first point
	 * @param point2 second point
	 * @return true if the points are neighbors (26-pixel neighborhood)
	 */
	private boolean isNeighbor(Point point1, Point point2) {
		return Math.sqrt(Math.pow(point1.x - point2.x, 2) +
				Math.pow(point1.y - point2.y, 2) +
				Math.pow(point1.z - point2.z, 2)) <= Math.sqrt(3);
	}

	/* ----------------------------------------------------------------------- */
	/**
	 * Check if the point is slab.
	 *
	 * @param point actual point
	 * @return true if the point has slab status
	 */
	private boolean isSlab(Point point) {
		return getPixel(taggedImage, point.x, point.y,
				point.z) == AnalyzeSkeleton2_.SLAB;
	}

	/* ----------------------------------------------------------------------- */
	/**
	 * Check if the point is a junction.
	 *
	 * @param point actual point
	 * @return true if the point has slab status
	 */
	public boolean isJunction(Point point) {
		return getPixel(taggedImage, point.x, point.y,
				point.z) == AnalyzeSkeleton2_.JUNCTION;
	}

	/* ----------------------------------------------------------------------- */
	/**
	 * Check if the point is an end point.
	 *
	 * @param point actual point
	 * @return true if the point has slab status
	 */
	private boolean isEndPoint(Point point) {
		return getPixel(taggedImage, point.x, point.y,
				point.z) == AnalyzeSkeleton2_.END_POINT;
	}

	/* ----------------------------------------------------------------------- */
	/**
	 * Check if the point is an 'unprotected' end point.
	 *
	 * @param point actual point
	 * @param roi   Only points outside this ROI will have end-point status. ROI can
	 *              be associated to a
	 *              single image in the stack (or all) as per
	 *              {@link ij.gui.Roi#getPosition
	 *              ij.gui.Roi.getPosition()}
	 * @return true if the point has end-point status
	 */
	private boolean isEndPoint(Point point, Roi roi) {
		boolean endPoint = isEndPoint(point);
		if (endPoint && roi != null) {
			// Is roi associated with all the images in the ImageStack?
			boolean roiContainsZ = roi.getPosition() == 0 ? true : roi.getPosition() == point.z;

			// Maintain end-point status only if point lies outside roi boundaries
			endPoint = !(roi.contains(point.x, point.y) && roiContainsZ);
		}
		return endPoint;
	}

	/* ----------------------------------------------------------------------- */
	/**
	 * Check if the point is a junction.
	 *
	 * @param x x- voxel coordinate
	 * @param y y- voxel coordinate
	 * @param z z- voxel coordinate
	 * @return true if the point has slab status
	 */
	public boolean isJunction(int x, int y, int z) {
		return getPixel(taggedImage, x, y, z) == AnalyzeSkeleton2_.JUNCTION;
	}

	/* ----------------------------------------------------------------------- */
	/**
	 * Get next unvisited neighbor voxel.
	 *
	 * @param point starting point
	 * @return unvisited neighbor or null if all neighbors are visited
	 */
	private Point getNextUnvisitedVoxel(Point point) {
		Point unvisitedNeighbor = null;

		// Check neighbors status
		for (int x = -1; x < 2; x++)
			for (int y = -1; y < 2; y++)
				for (int z = -1; z < 2; z++) {
					if (x == 0 && y == 0 && z == 0)
						continue;

					if (getPixel(inputImage, point.x + x, point.y + y,
							point.z + z) != 0 &&
							isVisited(point.x + x, point.y + y,
									point.z + z) == false) {
						unvisitedNeighbor = new Point(point.x + x, point.y + y,
								point.z + z);
						break;
					}

				}

		return unvisitedNeighbor;
	}// end getNextUnvisitedVoxel

	/* ----------------------------------------------------------------------- */
	/**
	 * Get next unvisited junction neighbor voxel.
	 *
	 * @param point starting point
	 * @return unvisited neighbor or null if all neighbors are visited
	 */
	private Point getNextUnvisitedJunctionVoxel(Point point) {
		Point unvisitedNeighbor = null;

		// Check neighbors status
		for (int x = -1; x < 2; x++)
			for (int y = -1; y < 2; y++)
				for (int z = -1; z < 2; z++) {
					if (x == 0 && y == 0 && z == 0)
						continue;

					if (getPixel(inputImage, point.x + x, point.y + y,
							point.z + z) != 0 &&
							isVisited(point.x + x, point.y + y,
									point.z + z) == false
							&&
							isJunction(point.x + x, point.y + y, point.z + z)) {
						unvisitedNeighbor = new Point(point.x + x, point.y + y,
								point.z + z);
						break;
					}

				}

		return unvisitedNeighbor;
	}// end getNextUnvisitedJunctionVoxel

	// -----------------------------------------------------------------------
	/**
	 * Get next visited junction neighbor voxel excluding the ones belonging to a
	 * give vertex
	 *
	 * @param point   starting point
	 * @param exclude exclusion vertex
	 * @return unvisited neighbor or null if all neighbors are visited
	 */
	private Point getVisitedJunctionNeighbor(Point point, Vertex exclude) {
		Point finalNeighbor = null;

		// Check neighbors status
		for (int x = -1; x < 2; x++)
			for (int y = -1; y < 2; y++)
				for (int z = -1; z < 2; z++) {
					if (x == 0 && y == 0 && z == 0)
						continue;

					final Point neighbor = new Point(point.x + x, point.y + y,
							point.z + z);

					if (getPixel(inputImage, neighbor) != 0 &&
							isVisited(neighbor) && isJunction(neighbor) &&
							!exclude.getPoints().contains(neighbor)) {
						finalNeighbor = neighbor;
						break;
					}

				}

		return finalNeighbor;
	}// end getNextUnvisitedJunctionVoxel

	// -----------------------------------------------------------------------
	/**
	 * Check if a voxel is visited taking into account the borders. Out of range
	 * voxels are considered
	 * as visited.
	 *
	 * @param point
	 * @return true if the voxel is visited
	 */
	private boolean isVisited(Point point) {
		return isVisited(point.x, point.y, point.z);
	}

	/* ----------------------------------------------------------------------- */
	/**
	 * Check if a voxel is visited taking into account the borders. Out of range
	 * voxels are considered
	 * as visited.
	 *
	 * @param x x- voxel coordinate
	 * @param y y- voxel coordinate
	 * @param z z- voxel coordinate
	 * @return true if the voxel is visited
	 */
	private boolean isVisited(int x, int y, int z) {
		if (x >= 0 && x < width && y >= 0 && y < height && z >= 0 && z < depth)
			return visited[x][y][z];
		return true;
	}

	/* ----------------------------------------------------------------------- */
	/**
	 * Set value in the visited flags matrix.
	 *
	 * @param x x- voxel coordinate
	 * @param y y- voxel coordinate
	 * @param z z- voxel coordinate
	 * @param b
	 */
	private void setVisited(int x, int y, int z, boolean b) {
		if (x >= 0 && x < width && y >= 0 && y < height && z >= 0 && z < depth)
			visited[x][y][z] = b;
	}

	/* ----------------------------------------------------------------------- */
	/**
	 * Set value in the visited flags matrix.
	 *
	 * @param point voxel coordinates
	 * @param b     visited flag value
	 */
	private void setVisited(Point point, boolean b) {
		int x = point.x;
		int y = point.y;
		int z = point.z;

		setVisited(x, y, z, b);
	}

	/* ----------------------------------------------------------------------- */
	/**
	 * Tag skeleton dividing the voxels between end points, junctions and slabs.
	 *
	 * @param inputImage2 skeleton image to be tagged
	 * @return tagged skeleton image
	 */
	private ImageStack tagImage(ImageStack inputImage2) {
		// Create output image
		ImageStack outputImage = new ImageStack(width, height,
				inputImage2.getColorModel());

		// Tag voxels
		for (int z = 0; z < depth; z++) {
			outputImage.addSlice(inputImage2.getSliceLabel(z + 1),
					new ByteProcessor(width, height));
			for (int x = 0; x < width; x++)
				for (int y = 0; y < height; y++) {
					if (getPixel(inputImage2, x, y, z) != 0) {
						int numOfNeighbors = getNumberOfNeighbors(inputImage2, x,
								y, z);
						if (numOfNeighbors < 2) {
							setPixel(outputImage, x, y, z,
									AnalyzeSkeleton2_.END_POINT);
							totalNumberOfEndPoints++;
							Point endPoint = new Point(x, y, z);
							listOfEndPoints.add(endPoint);
						} else if (numOfNeighbors > 2) {
							setPixel(outputImage, x, y, z,
									AnalyzeSkeleton2_.JUNCTION);
							Point junction = new Point(x, y, z);
							listOfJunctionVoxels.add(junction);
							totalNumberOfJunctionVoxels++;
						} else {
							setPixel(outputImage, x, y, z,
									AnalyzeSkeleton2_.SLAB);
							Point slab = new Point(x, y, z);
							listOfSlabVoxels.add(slab);
							totalNumberOfSlabs++;
						}
					}
				}
		}

		return outputImage;
	}// end method tagImage

	/* ----------------------------------------------------------------------- */
	/**
	 * Tag skeleton dividing the voxels between end points, junctions and slabs.
	 *
	 * @param inputImage2 skeleton image to be tagged
	 * @param rootPoint   - point to designate as new root junction if one doesn't
	 *                    exist
	 * @return tagged skeleton image
	 */
	private ImageStack tagImage(ImageStack inputImage2, Point rootPoint) {
		// Create output image
		ImageStack outputImage = new ImageStack(width, height,
				inputImage2.getColorModel());

		// Tag voxels
		for (int z = 0; z < depth; z++) {
			outputImage.addSlice(inputImage2.getSliceLabel(z + 1),
					new ByteProcessor(width, height));
			for (int x = 0; x < width; x++)
				for (int y = 0; y < height; y++) {
					if (getPixel(inputImage2, x, y, z) != 0) {
						int numOfNeighbors = getNumberOfNeighbors(inputImage2, x,
								y, z);
						if (numOfNeighbors < 2) {
							setPixel(outputImage, x, y, z,
									AnalyzeSkeleton2_.END_POINT);
							totalNumberOfEndPoints++;
							Point endPoint = new Point(x, y, z);
							listOfEndPoints.add(endPoint);
						} else if (numOfNeighbors > 2) {
							setPixel(outputImage, x, y, z,
									AnalyzeSkeleton2_.JUNCTION);
							Point junction = new Point(x, y, z);
							listOfJunctionVoxels.add(junction);
							totalNumberOfJunctionVoxels++;
						} else {
							setPixel(outputImage, x, y, z,
									AnalyzeSkeleton2_.SLAB);
							Point slab = new Point(x, y, z);
							listOfSlabVoxels.add(slab);
							totalNumberOfSlabs++;
						}
					}
				}
		}

		setPixel(outputImage, rootPoint.x, rootPoint.y, rootPoint.z,
				AnalyzeSkeleton2_.JUNCTION);

		return outputImage;
	}// end method tagImage

	/* ----------------------------------------------------------------------- */
	/**
	 * Get number of neighbors of a voxel in a 3D image (0 border conditions).
	 *
	 * @param image 3D image (ImageStack)
	 * @param x     x- coordinate
	 * @param y     y- coordinate
	 * @param z     z- coordinate (in image stacks the indexes start at 1)
	 * @return corresponding 27-pixels neighborhood (0 if out of image)
	 */
	private int getNumberOfNeighbors(ImageStack image, int x, int y, int z) {
		int n = 0;
		byte[] neighborhood = getNeighborhood(image, x, y, z);

		for (int i = 0; i < 27; i++)
			if (neighborhood[i] != 0)
				n++;
		// We return n-1 because neighborhood includes the actual voxel.
		return n - 1;
	}// end method getNumberOfNeighbors

	// -----------------------------------------------------------------------
	/**
	 * Get average 3x3x3 neighborhood pixel value of a given point
	 *
	 * @param image input image
	 * @param p     image coordinates
	 */
	private double getAverageNeighborhoodValue(ImageStack image, Point p) {
		byte[] neighborhood = getNeighborhood(image, p);

		double avg = 0;
		for (int i = 0; i < neighborhood.length; i++)
			avg += neighborhood[i] & 0xFF;
		if (neighborhood.length > 0)
			return avg / neighborhood.length;
		else
			return 0;
	}// end method getAverageNeighborhoodValue

	// -----------------------------------------------------------------------
	/**
	 * Get average neighborhood pixel value of a given point.
	 *
	 * @param image    input image
	 * @param p        image coordinates
	 * @param x_offset x- neighborhood offset
	 * @param y_offset y- neighborhood offset
	 * @param z_offset z- neighborhood offset
	 * @return average neighborhood pixel value
	 */
	public static double getAverageNeighborhoodValue(
			final ImageStack image,
			final Point p,
			final int x_offset,
			final int y_offset,
			final int z_offset) {
		byte[] neighborhood = getNeighborhood(image, p, x_offset, y_offset,
				z_offset);

		double avg = 0;
		for (int i = 0; i < neighborhood.length; i++)
			avg += neighborhood[i] & 0xFF;
		if (neighborhood.length > 0)
			return avg / neighborhood.length;
		else
			return 0;
	}// end method getAverageNeighborhoodValue

	// -----------------------------------------------------------------------
	/**
	 * Get neighborhood of a pixel in a 3D image (0 border conditions).
	 *
	 * @param image    3D image (ImageStack)
	 * @param p        point coordinates
	 * @param x_offset x- neighborhood offset
	 * @param y_offset y- neighborhood offset
	 * @param z_offset z- neighborhood offset
	 * @return corresponding neighborhood (0 if out of image)
	 */
	public static byte[] getNeighborhood(
			final ImageStack image,
			final Point p,
			final int x_offset,
			final int y_offset,
			final int z_offset) {
		final byte[] neighborhood = new byte[(2 * x_offset + 1) *
				(2 * y_offset + 1) * (2 * z_offset + 1)];

		for (int l = 0, k = p.z - z_offset; k <= p.z + z_offset; k++)
			for (int j = p.y - y_offset; j <= p.y + y_offset; j++)
				for (int i = p.x - x_offset; i <= p.x + x_offset; i++, l++)
					neighborhood[l] = getPixel(image, i, j, k);
		return neighborhood;
	} // end getNeighborhood

	// -----------------------------------------------------------------------
	/**
	 * Get neighborhood of a pixel in a 3D image (0 border conditions).
	 *
	 * @param image 3D image (ImageStack)
	 * @param p     3D point coordinates
	 * @return corresponding 27-pixels neighborhood (0 if out of image)
	 */
	private byte[] getNeighborhood(ImageStack image, Point p) {
		return getNeighborhood(image, p.x, p.y, p.z);
	}

	// -----------------------------------------------------------------------
	/**
	 * Get neighborhood of a pixel in a 3D image (0 border conditions).
	 *
	 * @param image 3D image (ImageStack)
	 * @param x     x- coordinate
	 * @param y     y- coordinate
	 * @param z     z- coordinate (in image stacks the indexes start at 1)
	 * @return corresponding 27-pixels neighborhood (0 if out of image)
	 */
	private byte[] getNeighborhood(ImageStack image, int x, int y, int z) {
		byte[] neighborhood = new byte[27];

		neighborhood[0] = getPixel(image, x - 1, y - 1, z - 1);
		neighborhood[1] = getPixel(image, x, y - 1, z - 1);
		neighborhood[2] = getPixel(image, x + 1, y - 1, z - 1);

		neighborhood[3] = getPixel(image, x - 1, y, z - 1);
		neighborhood[4] = getPixel(image, x, y, z - 1);
		neighborhood[5] = getPixel(image, x + 1, y, z - 1);

		neighborhood[6] = getPixel(image, x - 1, y + 1, z - 1);
		neighborhood[7] = getPixel(image, x, y + 1, z - 1);
		neighborhood[8] = getPixel(image, x + 1, y + 1, z - 1);

		neighborhood[9] = getPixel(image, x - 1, y - 1, z);
		neighborhood[10] = getPixel(image, x, y - 1, z);
		neighborhood[11] = getPixel(image, x + 1, y - 1, z);

		neighborhood[12] = getPixel(image, x - 1, y, z);
		neighborhood[13] = getPixel(image, x, y, z);
		neighborhood[14] = getPixel(image, x + 1, y, z);

		neighborhood[15] = getPixel(image, x - 1, y + 1, z);
		neighborhood[16] = getPixel(image, x, y + 1, z);
		neighborhood[17] = getPixel(image, x + 1, y + 1, z);

		neighborhood[18] = getPixel(image, x - 1, y - 1, z + 1);
		neighborhood[19] = getPixel(image, x, y - 1, z + 1);
		neighborhood[20] = getPixel(image, x + 1, y - 1, z + 1);

		neighborhood[21] = getPixel(image, x - 1, y, z + 1);
		neighborhood[22] = getPixel(image, x, y, z + 1);
		neighborhood[23] = getPixel(image, x + 1, y, z + 1);

		neighborhood[24] = getPixel(image, x - 1, y + 1, z + 1);
		neighborhood[25] = getPixel(image, x, y + 1, z + 1);
		neighborhood[26] = getPixel(image, x + 1, y + 1, z + 1);

		return neighborhood;
	} // end getNeighborhood

	// -----------------------------------------------------------------------
	/**
	 * Get pixel in 3D image (0 border conditions)
	 *
	 * @param image 3D image
	 * @param x     x- coordinate
	 * @param y     y- coordinate
	 * @param z     z- coordinate (in image stacks the indexes start at 1)
	 * @return corresponding pixel (0 if out of image)
	 */
	public static byte getPixel(final ImageStack image, final int x,
			final int y, final int z) {
		final int width = image.getWidth();
		final int height = image.getHeight();
		final int depth = image.getSize();

		if (x >= 0 && x < width && y >= 0 && y < height && z >= 0 && z < depth)
			return ((byte[]) image.getPixels(z + 1))[x + y * width];
		else
			return 0;
	} // end getPixel

	/* ----------------------------------------------------------------------- */
	/**
	 * Get pixel in 3D image (0 border conditions).
	 *
	 * @param image 3D image
	 * @param x     x- coordinate
	 * @param y     y- coordinate
	 * @param z     z- coordinate (in image stacks the indexes start at 1)
	 * @return corresponding pixel (0 if out of image)
	 */
	private float getFloatPixel(
			ImageStack image,
			int x,
			int y,
			int z) {
		if (x >= 0 && x < width && y >= 0 && y < height && z >= 0 && z < depth)
			return ((float[]) image.getPixels(z + 1))[x + y * width];
		else
			return 0;
	} // end getFloatPixel

	/* ----------------------------------------------------------------------- */
	/**
	 * Get pixel in 3D image (0 border conditions).
	 *
	 * @param image 3D image
	 * @param point point to be evaluated
	 * @return corresponding pixel (0 if out of image)
	 */
	private float getFloatPixel(
			ImageStack image,
			Point point) {
		return getFloatPixel(image, point.x, point.y, point.z);
	} // end getFloatPixel

	/* ----------------------------------------------------------------------- */
	/**
	 * Get pixel in 3D image (0 border conditions).
	 *
	 * @param image 3D image
	 * @param point point to be evaluated
	 * @return corresponding pixel (0 if out of image)
	 */
	private byte getPixel(ImageStack image, Point point) {
		return getPixel(image, point.x, point.y, point.z);
	} // end getPixel

	/* ----------------------------------------------------------------------- */
	/**
	 * Set pixel in 3D image.
	 *
	 * @param image 3D image
	 * @param p     point coordinates
	 * @param value pixel value
	 */
	private void setPixel(ImageStack image, Point p, byte value) {
		if (p.x >= 0 && p.x < width && p.y >= 0 && p.y < height && p.z >= 0 &&
				p.z < depth)
			((byte[]) image.getPixels(p.z + 1))[p.x + p.y * width] = value;
	} // end setPixel

	/* ----------------------------------------------------------------------- */
	/**
	 * Set pixel in 3D image.
	 *
	 * @param image 3D image
	 * @param x     x- coordinate
	 * @param y     y- coordinate
	 * @param z     z- coordinate (in image stacks the indexes start at 1)
	 * @param value pixel value
	 */
	private void setPixel(ImageStack image, int x, int y, int z, byte value) {
		if (x >= 0 && x < width && y >= 0 && y < height && z >= 0 && z < depth)
			((byte[]) image.getPixels(z + 1))[x + y * width] = value;
	} // end setPixel

	/* ----------------------------------------------------------------------- */
	/**
	 * Set pixel in 3D (short) image.
	 *
	 * @param image 3D image
	 * @param x     x- coordinate
	 * @param y     y- coordinate
	 * @param z     z- coordinate (in image stacks the indexes start at 1)
	 * @param value pixel value
	 */
	private void setPixel(ImageStack image, int x, int y, int z, float value) {
		if (x >= 0 && x < width && y >= 0 && y < height && z >= 0 && z < depth)
			((float[]) image.getPixels(z + 1))[x + y * width] = value;
	} // end setPixel

	/* ----------------------------------------------------------------------- */
	/** Show plug-in information. */
	void showAbout() {
		IJ.showMessage(
				"About AnalyzeSkeleton...",
				"This plug-in filter analyzes a 2D/3D image skeleton.\n");
	} // end showAbout

	/**
	 * Determine the longest shortest path using the APSP (all pairs shortest path)
	 * warshall algorithm
	 *
	 * @param graph              the graph of a tree
	 * @param shortestPathPoints list to store the longest shortest path points
	 * @return longest shortest path length
	 * @author Huub Hovens
	 */
	private double warshallAlgorithm(Graph graph,
			ArrayList<Point> shortestPathPoints) {
		// local fields
		/** vertex 1 of an edge */
		Vertex v1 = null;
		/** vertex 2 of an edge */
		Vertex v2 = null;
		/** the equivalent of row in a matrix */
		int row = 0;
		/** the equivalent of column in a matrix */
		int column = 0;
		/** the value of the longest shortest path */
		double maxPath = 0;
		/** row that contains the longest shortest path value */
		int a = 0;
		/** column that contains the longest shortest path value */
		int b = 0;

		ArrayList<Edge> edgeList = graph.getEdges();
		ArrayList<Vertex> vertexList = graph.getVertices();

		// check for paths of only one vertex
		if (vertexList.size() == 1) {
			shortestPathPoints.add(vertexList.get(0).getPoints().get(0));
			spx = vertexList.get(0).getPoints().get(0).x;
			spy = vertexList.get(0).getPoints().get(0).y;
			spz = vertexList.get(0).getPoints().get(0).z;
			return 0;
		}

		// create empty adjacency and predecessor matrix

		/**
		 * the matrix that contains the length of the shortest path from vertex a to
		 * vertex b
		 */
		double[][] adjacencyMatrix = new double[vertexList.size()][vertexList
				.size()];
		/**
		 * the matrix that contains the predecessor vertex of vertex b in the shortest
		 * path from vertex a
		 * to b
		 */
		int[][] predecessorMatrix = new int[vertexList.size()][vertexList
				.size()];

		// initial conditions for both matrices
		/*
		 * create 2D-adjacency array with the distance between the nodes.
		 * distance i --> i = 0
		 * distance i -/> j = infinite (edge does not exist)
		 * distance i --> j = length of branch
		 */

		/*
		 * create 2D-predecessor array using interconnected vertices.
		 * the predecessor matrix contains the predecessor of j in a path from node i to
		 * j.
		 * distance i --> i = NIL (there is no edge between a single vertex)
		 * distance i -/> j = NIL (there is no edge between the vertices)
		 * distance i --> j = i (the initial matrix only contains paths with a single
		 * edge)
		 *
		 * I'm using -1 as NIL since it cannot refer to an index (and thus a vertex)
		 */

		// applying initial conditions
		for (int i = 0; i < vertexList.size(); i++) {
			for (int j = 0; j < vertexList.size(); j++) {
				adjacencyMatrix[i][j] = Double.POSITIVE_INFINITY;
				predecessorMatrix[i][j] = -1;
			}
		}

		for (Edge edge : edgeList) {
			v1 = edge.getV1();
			v2 = edge.getV2();
			// use the index of the vertices as the index in the matrix
			row = vertexList.indexOf(v1);
			if (row == -1) {
				IJ.log("Vertex " + v1.getPoints().get(0) +
						" not found in the list of vertices!");
				continue;
			}

			column = vertexList.indexOf(v2);
			if (column == -1) {
				IJ.log("Vertex " + v2.getPoints().get(0) +
						" not found in the list of vertices!");
				continue;
			}

			/*
			 * the diagonal is 0.
			 *
			 * Because not every vertex is a 'v1 vertex'
			 * the [column][column] statement is needed as well.
			 *
			 * in an undirected graph the adjacencyMatrix is symmetric
			 * thus A = Transpose(A)
			 */
			adjacencyMatrix[row][row] = 0;
			adjacencyMatrix[column][column] = 0;
			adjacencyMatrix[row][column] = edge.getLength();
			adjacencyMatrix[column][row] = edge.getLength();

			/*
			 * the diagonal remains -1.
			 * for the rest I use the index of the vertex so I can later refer to the
			 * vertexList
			 * for the correct information
			 *
			 * Determining what belongs where requires careful consideration of the
			 * definition
			 *
			 * the array contains the predecessor of "column" in a path from "row" to
			 * "column"
			 * therefore in the other statement it is the other way around.
			 */
			predecessorMatrix[row][row] = -1;
			predecessorMatrix[column][column] = -1;
			predecessorMatrix[row][column] = row;
			predecessorMatrix[column][row] = column;
		}
		// matrices now have their initial conditions

		// the warshall algorithm with k as candidate vertex and i and j walk through
		// the adjacencyMatrix
		// the predecessor matrix is updated at the same time.

		for (int k = 0; k < vertexList.size(); k++) {
			for (int i = 0; i < vertexList.size(); i++) {
				for (int j = 0; j < vertexList.size(); j++) {
					if (adjacencyMatrix[i][k] +
							adjacencyMatrix[k][j] < adjacencyMatrix[i][j]) {
						adjacencyMatrix[i][j] = adjacencyMatrix[i][k] +
								adjacencyMatrix[k][j];
						predecessorMatrix[i][j] = predecessorMatrix[k][j];

					}
				}
			}
		}

		// find the maximum of all shortest paths
		for (int i = 0; i < vertexList.size(); i++) {
			for (int j = 0; j < vertexList.size(); j++) {
				// sometimes infinities still remain
				if (adjacencyMatrix[i][j] > maxPath &&
						adjacencyMatrix[i][j] != Double.POSITIVE_INFINITY) {
					maxPath = adjacencyMatrix[i][j];
					a = i;
					b = j;

				}
			}
		}

		// trace back the longest shortest path
		reconstructPath(predecessorMatrix, a, b, edgeList, vertexList,
				shortestPathPoints);

		// !important return maxPath;
		return maxPath;

	}
	// end method warshallAlgorithm

	/**
	 * Reconstruction and visualisation of the longest shortest path found by the
	 * APSP warshall
	 * algorithm
	 *
	 * @param predecessorMatrix  the Matrix which contains the predecessor of vertex
	 *                           b in the shortest
	 *                           path from a to b
	 * @param startIndex         the index of the row which contains the longest
	 *                           shortest path
	 * @param endIndex           the index of the column which contains the longest
	 *                           shortest path
	 * @param edgeList           the list of edges
	 * @param vertexList         the list of vertices
	 * @param shortestPathPoints contains points of the longest shortest path for
	 *                           each graph
	 * @author Huub Hovens
	 */
	private void reconstructPath(
			int[][] predecessorMatrix,
			int startIndex,
			int endIndex,
			ArrayList<Edge> edgeList,
			ArrayList<Vertex> vertexList,
			ArrayList<Point> shortestPathPoints) {
		// We know the first and last vertex of the longest shortest path, namely a and
		// b
		// using the predecessor matrix we can now determine the path that is taken from
		// a to b
		// remember a and b are indices and not the actual vertices.

		int b = endIndex;
		final int a = startIndex;

		while (b != a) {
			Vertex predecessor = vertexList.get(predecessorMatrix[a][b]);
			Vertex endvertex = vertexList.get(b);
			ArrayList<Edge> sp_edgeslist = new ArrayList<>();
			Double lengthtest = Double.POSITIVE_INFINITY;
			Edge shortestedge = null;

			// search all edges for a combination of the two vertices
			for (Edge edge : edgeList) {

				if (edge.getV1() == predecessor && edge.getV2() == endvertex ||
						edge.getV1() == endvertex && edge.getV2() == predecessor) {
					// sometimes there are multiple edges between two vertices so add them to a list
					// for a second test
					sp_edgeslist.add(edge);
				}

			}
			// the second test
			// this test looks which edge has the shortest length in sp_edgeslist
			for (Edge edge : sp_edgeslist) {
				if (edge.getLength() < lengthtest) {
					shortestedge = edge;
					lengthtest = edge.getLength();
				}

			}
			// add vertex 1 points
			Vertex v1 = shortestedge.getV2() != predecessor ? shortestedge.getV2() : shortestedge.getV1();
			for (Point p : v1.getPoints()) {
				if (!shortestPathPoints.contains(p)) {
					shortestPathPoints.add(p);
					// setPixel(this.shortPathImage, p.x, p.y, p.z, SHORTEST_PATH);
				}
			}

			// add slab points of the shortest edge to the list of points
			ArrayList<Point> slabs = shortestedge.getSlabs();
			// reverse order if needed
			if (shortestedge.getV2() != predecessor)
				Collections.reverse(slabs);
			for (Point p : slabs) {
				shortestPathPoints.add(p);
				setPixel(shortPathImage, p.x, p.y, p.z, SHORTEST_PATH);
			}

			// add vertex 2 points too
			Vertex v2 = shortestedge.getV2() != predecessor ? shortestedge.getV1() : shortestedge.getV2();
			for (Point p : v2.getPoints()) {
				if (!shortestPathPoints.contains(p)) {
					shortestPathPoints.add(p);
					// setPixel(this.shortPathImage, p.x, p.y, p.z, SHORTEST_PATH);
				}
			}

			// now make the index of the endvertex the index of the predecessor so that the
			// path now goes from
			// a to predecessor and repeat cycle
			b = predecessorMatrix[a][b];
		}
		if (shortestPathPoints.size() != 0) {
			spx = shortestPathPoints.get(0).x;
			spy = shortestPathPoints.get(0).y;
			spz = shortestPathPoints.get(0).z;
		}

	}
	// end method reconstructPath

	/**
	 * Get the stack containing all the trees labeld with their corresponding
	 * skeleton id.
	 *
	 * @return labeled-skeleton image stack
	 */
	public ImageStack getLabeledSkeletons() {
		return labeledSkeletons;
	}

	private void createSettingsDialog() {
		settingsDialog = new GenericDialog("Analyze Skeleton");
		Font headerFont = new Font("SansSerif", Font.BOLD, 12);

		settingsDialog.setInsets(0, 0, 0);
		settingsDialog.addMessage("Elimination of Loops:", headerFont);
		settingsDialog.addChoice("Prune cycle method: ",
				AnalyzeSkeleton2_.pruneCyclesModes,
				AnalyzeSkeleton2_.pruneCyclesModes[pruneIndex]);

		settingsDialog.setInsets(20, 0, -15); // default top inset for 1st checkbox is 15
		settingsDialog.addMessage("Elimination of End-points:", headerFont);
		settingsDialog.addCheckbox("Prune ends", pruneEnds);
		settingsDialog.addCheckbox("Exclude ROI from pruning", protectRoi);

		settingsDialog.setInsets(20, 0, 0); // default top inset for subsequent checkboxes is 0
		settingsDialog.addMessage("Results and Output:", headerFont);
		settingsDialog.addCheckbox("Calculate largest shortest path",
				calculateShortestPath);
		settingsDialog.addCheckbox("Show detailed info", verbose);
		settingsDialog.addCheckbox("Display labeled skeletons",
				displaySkeletons);

		settingsDialog.addHelp(HELP_URL);
		dialogItemChanged(settingsDialog, null);
	}

	private void loadDialogSettings() {
		String index = Prefs.get(PRUNE_MODE_INDEX_KEY,
				String.valueOf(DEFAULT_PRUNE_MODE_INDEX));
		pruneIndex = Integer.parseInt(index); // fails to find the key
		pruneEnds = Prefs.get(PRUNE_ENDS_KEY, DEFAULT_PRUNE_ENDS);
		calculateShortestPath = Prefs.get(CALCULATE_PATH_KEY,
				DEFAULT_CALCULATE_SHORTEST_PATH);
		verbose = Prefs.get(VERBOSE_KEY, DEFAULT_VERBOSE);
		displaySkeletons = Prefs.get(DISPLAY_SKELETONS_KEY,
				DEFAULT_DISPLAY_SKELETONS);
	}

	private void saveDialogSettings() {
		Prefs.set(PRUNE_MODE_INDEX_KEY, String.valueOf(pruneIndex));
		Prefs.set(PRUNE_ENDS_KEY, pruneEnds);
		Prefs.set(CALCULATE_PATH_KEY, calculateShortestPath);
		Prefs.set(VERBOSE_KEY, verbose);
		Prefs.set(DISPLAY_SKELETONS_KEY, displaySkeletons);
	}

	private void setSettingsFromDialog() {
		pruneIndex = settingsDialog.getNextChoiceIndex();
		pruneEnds = settingsDialog.getNextBoolean();
		protectRoi = settingsDialog.getNextBoolean();
		calculateShortestPath = settingsDialog.getNextBoolean();
		verbose = settingsDialog.getNextBoolean();
		displaySkeletons = settingsDialog.getNextBoolean();
	}
}// end class AnalyzeSkeleton2_

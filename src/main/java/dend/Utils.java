/*
 * #%L
 * Modified based on Utils.java of the ImageJ Strahler Analysis plugin by Tiago Ferreira.
 *
 * Major Modifications:
 * Deleted unused functions.
 * Edited displayed text.
 *
 * Diana Bank, Inle Bush, and Elena Zhong, Chun Han Lab, Cornell University
 * #L%
 */
package dend;

import java.awt.Window;
import java.util.Arrays;
import java.util.List;

import ij.IJ;
import ij.ImagePlus;
import ij.WindowManager;
import ij.measure.ResultsTable;
import ij.text.TextWindow;

/** This class contains utilities. Modified based on {@code Utils.java} of the ImageJ
 * {@code Strahler Analysis} plugin by Tiago Ferreira.
 * <p>
 * For more information about the original plugin, visit the hIPNAT repository
 * {@literal https://github.com/tferr/hIPNAT} and the the plugin's documentation page:
 * {@literal http://imagej.net/Strahler_Analysis}
 *
 * @author Diana Bank
 * @author Inle Bush
 * @author Elena Zhong
 * @author Chun Han Lab, Cornell University */
public class Utils {

	/** Private constructor to prevent class instantiation. */
	private Utils() {}

	/** Returns the ResultsTable of the specified window title (if open) or a new ResultsTable with
	 * appropriated properties (precision of 5 decimal places, no row numbers, "NaN" padding of empty
	 * cells)
	 *
	 * @param title the window title of the table
	 * @return a referenced to the opened ResultsTable or a new one if the window of the specified title
	 *         is not associated to a valid ResultsTable */
	public static ResultsTable getTable(final String title) {
		ResultsTable rt= null;
		final Window window= WindowManager.getWindow(title);
		if (window != null)
			rt= ((TextWindow) window).getTextPanel().getResultsTable();
		if (rt == null)
			rt= new ResultsTable();
		rt.setPrecision(5);
		rt.setNaNEmptyCells(true);
		rt.showRowNumbers(false);
		return rt;
	}

	/** Macro-friendly error message.
	 *
	 * If a macro is running it will not be aborted and the error message is displayed in the "Log"
	 * window (or the Java console if ImageJ is not present). Otherwise displays a regular
	 * {@link ij.IJ#error(String) IJ.error()}.
	 *
	 * @param errorTitle The error title
	 * @param errorMsg   The error message
	 * @param imp        The Image to be mentioned in the message. It is ignored if {@code null} */
	public static void error(final String errorTitle, final String errorMsg,
		final ImagePlus imp) {
		String title= "Dendrite Arbor Analyzer";
		if (errorTitle != null)
			title= errorTitle + " (" + title + ")";
		final String impMsg= imp == null ? "" :
			"Error while processing " + imp.getTitle();
		if (IJ.macroRunning())
			IJ.log("\n>>> " + title + ": " + impMsg + "\n" + errorMsg);
		else IJ.error(title, impMsg + "\n" + errorMsg);
	}

	public static boolean validSkelDependencies() {
		return classExists(
			Arrays.asList("sc.fiji.analyzeSkeleton.AnalyzeSkeleton_",
				"sc.fiji.skeletonize3D.Skeletonize3D_"));
	}

	public static boolean classExists(final List<String> classStringNames) {
		if (classStringNames != null) {
			for (final String cls : classStringNames) {
				try {
					Class.forName(cls);
				} catch (final ClassNotFoundException e) {
					IPNAT.handleException(e);
					return false;
				}
			}
		}
		return true;
	}
}

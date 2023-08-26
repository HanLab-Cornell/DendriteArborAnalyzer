/*
 * #%L
 * Modified based on IPNAT.java of the ImageJ Strahler Analysis plugin by Tiago Ferreira.
 *
 * Major Modifications:
 * Deleted unused functions.
 *
 * Diana Bank, Inle Bush, and Elena Zhong, Chun Han Lab, Cornell University
 * #L%
 */

package dend;

import ij.IJ;

/** This class preserves the exception handler function of the original plugin. Modified based on
 * {@code IPNAT.java} of the ImageJ {@code Strahler Analysis} plugin by Tiago Ferreira.
 * <p>
 * For more information about the original plugin, visit the hIPNAT repository
 * {@literal https://github.com/tferr/hIPNAT} and the the plugin's documentation page:
 * {@literal http://imagej.net/Strahler_Analysis}
 *
 * @author Diana Bank
 * @author Inle Bush
 * @author Elena Zhong
 * @author Chun Han Lab, Cornell University
 * */

public class IPNAT {

	private IPNAT() {}

	public static void handleException(final Exception e) {
		IJ.setExceptionHandler(new dend.ExceptionHandler());
		IJ.handleException(e);
		IJ.setExceptionHandler(null); // Revert to the default behavior
	}

}

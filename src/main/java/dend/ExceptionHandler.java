/*
 * #%L
 * Modified based on ExceptionHandler.java of the ImageJ Strahler Analysis plugin by Tiago Ferreira.
 *
 * Major Modifications:
 * Edited displayed text.
 *
 * Diana Bank, Inle Bush, and Elena Zhong, Chun Han Lab, Cornell University
 * #L%
 */
package dend;

import java.io.CharArrayWriter;
import java.io.PrintWriter;

import ij.IJ;
import ij.text.TextWindow;

/** hIPNAT's ExceptionHandler. Modified based on {@code ExceptionHandler.java} of the ImageJ
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
public class ExceptionHandler implements IJ.ExceptionHandler {

	String CLASS_NOT_FOUND= "A Java file was not found";
	String METHOD_NOT_FOUND= "Your IJ installation is likeley outdated";
	String UNNAMED_ERROR= "An error occured";
	String TIPS= "Troubleshooting tips:\n" +
		"  - Ensure you are subscribed to the Java-8 and Fiji update sites\n" +
		"  - Run the updater to install missing files or update deprecated ones\n" +
		"  - Useful resources:\n" +
		"    - http://imagej.net/Troubleshooting\n" +
		"    - http://imagej.net/Frequently_Asked_Questions\n" +
		"    - http://forum.imagej.net/";

	@Override
	public void handle(final Throwable t) {
		final CharArrayWriter writer= new CharArrayWriter();
		final PrintWriter pw= new PrintWriter(writer);
		t.printStackTrace(pw);
		final String trace= writer.toString();
		String tMsg;
		switch (trace.substring(0, trace.indexOf(":"))) {
		case "java.lang.ClassNotFoundException":
			tMsg= CLASS_NOT_FOUND;
			break;
		case "java.lang.NoSuchMethodException":
			tMsg= METHOD_NOT_FOUND;
			break;
		default:
			tMsg= UNNAMED_ERROR;
		}
		tMsg+= " (details below). " + TIPS + "\n \n" + trace;
		if (IJ.getInstance() != null) {
			new TextWindow("Exception", "Dendrite Arbor Analyzer", tMsg, 500,
				250);
		} else IJ.log(tMsg);
	}

}
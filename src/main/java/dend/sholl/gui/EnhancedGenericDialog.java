/*
 * #%L
 * Sholl_Analysis plugin for ImageJ
 * %%
 * Copyright (C) 2016 Tiago Ferreira
 * %%
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public
 * License along with this program.  If not, see
 * <http://www.gnu.org/licenses/gpl-3.0.html>.
 * #L%
 */
package dend.sholl.gui;

import java.awt.BorderLayout;
import java.awt.Button;
import java.awt.Color;
import java.awt.Component;
import java.awt.Dimension;
import java.awt.Font;
import java.awt.Frame;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.Panel;
import java.awt.ScrollPane;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;
import java.awt.event.WindowEvent;

import javax.swing.JPopupMenu;
import javax.swing.UIManager;

import fiji.util.gui.GenericDialogPlus;
import ij.IJ;
import ij.plugin.frame.Recorder;

/** Enhances GenericDialog with a few additional features, including scrollbars as soon as the
 * Dialog is too large to be displayed, ability to use the "help" button to display a drop-down
 * menu, and labels featuring clickable hyperlinks. Customizations are ignored if running
 * headless. */
public class EnhancedGenericDialog extends GenericDialogPlus {
	private static final long serialVersionUID= 1L;

	private String labelOfHelpActionButton= null;
	private ActionListener helpActionButtonListener= null;
	private MouseAdapter helpActionMouseListener= null;

	/** {@link ij.gui.GenericDialog#GenericDialog(String, Frame) GenericDialog constructor}
	 *
	 * @param title Dialog's title */
	public EnhancedGenericDialog(final String title) {
		super(title);
	}

	/** {@link ij.gui.GenericDialog#GenericDialog(String) GenericDialog constructor}
	 *
	 * @param title  Dialog's title
	 * @param parent Parent frame */
	public EnhancedGenericDialog(final String title, final Frame parent) {
		super(title, parent);
	}

	/** Appends a message to the specified dialog consisting of one or more lines of text, and assigns a
	 * functional URL to it.
	 *
	 * @param text  The contents of the clickable label
	 * @param font  the label font. If {@code null}, the GenericDialog's default font is used
	 * @param color the label color. If {@code null}, the GenericDialog's default foreground color is
	 *              used
	 * @param url   The URL to be opened by the default browser of the OS */
	public void addHyperlinkMessage(final String text, final Font font,
		final Color color, final String url) {
		super.addMessage(text, font, color);
		Utils.addClickableURLtoLabel(super.getMessage(), url, color);
	}

	/** Allows users to visit the manuscript from a dialog prompt */
	public void addCitationMessage() {
		super.addMessage(Utils.citationMsg(), null, Utils.infoColor());
		Utils.addClickableURLtoLabel(super.getMessage(), Utils.citationURL(),
			Utils.infoColor());
	}

	@Override
	public void actionPerformed(final ActionEvent e) {
		final Object source= e.getSource();
		if (!Utils.isHeadless() && source != null &&
			labelOfHelpActionButton != null &&
			source.toString().contains(labelOfHelpActionButton) &&
			helpActionButtonListener != null) {
			helpActionButtonListener.actionPerformed(e);
		} else {
			super.actionPerformed(e);
		}
	}

	/** Adds a "Help" button and attaches the specified listener to it. In an IJ1 GenericDialog, event
	 * listeners triggered by the help button are only notified when the "Cancel" button is hidden. This
	 * method overcomes that limitation. NB: Actions triggered by menu items will not be macro
	 * recordable unless they trigger e.g., {@link ij.IJ#doCommand(String ) IJ.doCommand()} calls
	 *
	 * @param buttonLabel the label of the customized "Help" button
	 * @param listener    the ActionListener monitoring action events */
	public void assignListenerToHelpButton(final String buttonLabel,
		final ActionListener listener) {
		if (!Utils.isHeadless() && buttonLabel != null && listener != null) {
			super.addHelp("");
			super.setHelpLabel(buttonLabel);
			labelOfHelpActionButton= buttonLabel;
			helpActionButtonListener= listener;
		}
	}

	/** Attaches the specified JPopupMenu to the GenericDialog "Help" button renamed "More &#187;". NB:
	 * Actions triggered by menu items will not be macro recordable unless they trigger e.g.,
	 * {@link ij.IJ#doCommand(String ) IJ.doCommand()} calls.
	 *
	 * @param popupmenu the JPopupMenu to be attached to the "More &#187;" button.
	 * @see #assignPopupToHelpButton(String, JPopupMenu) */
	public void assignPopupToHelpButton(final JPopupMenu popupmenu) {
		assignPopupToHelpButton("More \u00bb", popupmenu);
	}

	/** Adds a "Help" button and attaches the specified JPopupMenu to it. NB: Actions triggered by menu
	 * items will not be macro recordable unless they trigger e.g., {@link ij.IJ#doCommand(String)
	 * IJ.doCommand()} calls.
	 *
	 * @param buttonLabel the label of the "Help" button.
	 * @param popupmenu   the JPopupMenu to be attached to the "Help" button. */
	public void assignPopupToHelpButton(final String buttonLabel,
		final JPopupMenu popupmenu) {
		if (!Utils.isHeadless() && buttonLabel != null && popupmenu != null) {
			// Ensure swing component is displayed with a java.awt look and feel
			try {
				UIManager
					.setLookAndFeel(UIManager.getSystemLookAndFeelClassName());
				popupmenu.updateUI();
			} catch (final Exception ignored) {}
			super.addHelp("");
			super.setHelpLabel(buttonLabel);
			labelOfHelpActionButton= buttonLabel;
			helpActionMouseListener= new MouseAdapter() {
				@Override
				public void mousePressed(final MouseEvent e) {
					popupmenu.show((Button) e.getSource(), 0, 0);
				}
			};

			helpActionButtonListener= new ActionListener() {
				@Override
				public void actionPerformed(final ActionEvent e) {
					final MouseEvent me= new MouseEvent(
						(Component) e.getSource(), MouseEvent.MOUSE_CLICKED,
						e.getWhen(), MouseEvent.MOUSE_PRESSED, 0, 0, 0, true);
					helpActionMouseListener.mousePressed(me);
				}
			};
		}
	}

	/** Adds AWT scroll bars to the dialog before displaying it, when not running headless. Scroll bars
	 * are only added if the largest dimension of the Dialog reaches ~90% of the primary display as
	 * reported by {@link IJ#getScreenSize()}. Dialog remains fully recordable. */
	public void showScrollableDialog() {
		if (!Utils.isHeadless())
			addScrollBars();
		super.showDialog();
	}

	/** Closes the dialog without Recording the command. */
	public void disposeWithouRecording() {
		dispatchEvent(new WindowEvent(this, WindowEvent.WINDOW_CLOSING));
		dispose();
		if (Recorder.record)
			Recorder.setCommand(null);
	}

	/** Adds AWT scroll bars to the GenericDialog. From an old version of bio-formats Window.Tools,
	 * licensed under GNU GPLv2 (April 2013)
	 *
	 * @see <a href=
	 *      "http://git.openmicroscopy.org/?p=bioformats.git;a=blob;f=components/loci-plugins/src/loci/plugins/util/WindowTools.java;hb=HEAD">
	 *      git.openmicroscopy</a> */
	private void addScrollBars() {

		final GridBagLayout layout= (GridBagLayout) getLayout();

		// extract components
		final int count= getComponentCount();
		final Component[] c= new Component[count];
		final GridBagConstraints[] gbc= new GridBagConstraints[count];
		for (int i= 0; i < count; i++ ) {
			c[i]= getComponent(i);
			gbc[i]= layout.getConstraints(c[i]);
		}

		// clear components
		removeAll();
		layout.invalidateLayout(this);

		// create new container panel
		final Panel newPane= new Panel();
		final GridBagLayout newLayout= new GridBagLayout();
		newPane.setLayout(newLayout);
		for (int i= 0; i < count; i++ ) {
			newLayout.setConstraints(c[i], gbc[i]);
			newPane.add(c[i]);
		}

		// HACK - get preferred size for container panel
		// NB: don't know a better way:
		// - newPane.getPreferredSize() doesn't work
		// - newLayout.preferredLayoutSize(newPane) doesn't work
		final Frame f= new Frame();
		f.setLayout(new BorderLayout());
		f.add(newPane, BorderLayout.CENTER);
		f.pack();
		final Dimension size= newPane.getSize();
		f.remove(newPane);
		f.dispose();

		// compute best size for scrollable viewport
		size.width+= 35; // initially 25;
		size.height+= 30; // initially 15;
		final Dimension screen= IJ.getScreenSize();
		final int maxWidth= 9 * screen.width / 10; // initially 7/8;
		final int maxHeight= 8 * screen.height / 10; // initially 3/4
		if (size.width > maxWidth)
			size.width= maxWidth;
		if (size.height > maxHeight)
			size.height= maxHeight;

		// Create scroll pane. This has some problems: In some installations
		// horizontal scrollbars will be displayed even if not required. In
		// others, setScrollPosition will not be respected and panel will not
		// be anchored to top
		final ScrollPane scroll= new ScrollPane(
			ScrollPane.SCROLLBARS_AS_NEEDED) {

			private static final long serialVersionUID= 1L;

			@Override
			public Dimension getPreferredSize() {
				return size;
			}
		};

		scroll.add(newPane);
		scroll.validate();
		scroll.setScrollPosition(0, 0);

		// ensure constraints
		final GridBagConstraints constraints= new GridBagConstraints();
		constraints.gridwidth= GridBagConstraints.REMAINDER;
		constraints.fill= GridBagConstraints.BOTH;
		constraints.weightx= 1.0;
		constraints.weighty= 1.0;
		layout.setConstraints(scroll, constraints);

		// add scroll pane to original container
		this.add(scroll);

	}

}
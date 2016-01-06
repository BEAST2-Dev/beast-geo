package beast.app.beauti;


import java.awt.Color;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;

import javax.swing.Box;
import javax.swing.JButton;
import javax.swing.JColorChooser;
import javax.swing.JDialog;
import javax.swing.JOptionPane;
import javax.swing.border.EmptyBorder;

import beast.app.beauti.BeautiDoc;
import beast.core.BEASTInterface;
import beast.core.Input;
import beast.app.draw.InputEditor;

public class ColorInputEditor2 extends InputEditor.Base {
	private static final long serialVersionUID = 1L;
	
	JButton button;
	
	public ColorInputEditor2(BeautiDoc doc) {
		super(doc);
	}

	@Override
	public Class<?> type() {
		return Color.class;
	}

	@Override
	public void init(Input<?> input, BEASTInterface plugin, int itemNr, ExpandOption bExpandOption, boolean bAddButtons) {
        m_bAddButtons = bAddButtons;
        m_input = input;
        m_beastObject = plugin;
        this.itemNr= itemNr;
        
        addInputLabel();

        
        Color color = (Color) m_input.get();
        button = new JButton("color");
        if (color != null) {
        	button.setBackground(color);
        }
        button.addActionListener(new ActionListener() {
			
			@Override
			public void actionPerformed(ActionEvent e) {
				editColor();
			}

		});
        add(button);
        
        add(Box.createHorizontalGlue());
        addValidationLabel();
    } // init
	
	private void editColor() {
		Color color = (Color) m_input.get();
        JColorChooser chooser = new JColorChooser();
        if (color != null) {
        	chooser.setColor(color);
        }
	
        JOptionPane optionPane = new JOptionPane(chooser,
                JOptionPane.QUESTION_MESSAGE,
                JOptionPane.OK_CANCEL_OPTION,
                null,
                null,
                null);
        optionPane.setBorder(new EmptyBorder(12, 12, 12, 12));

        //Frame frame = (doc != null ? doc.getFrame(): Frame.getFrames()[0]);
        final JDialog dialog = optionPane.createDialog(null, "Choose colour for " + m_input.getName());
        dialog.pack();

        dialog.setVisible(true);

        int result = JOptionPane.CANCEL_OPTION;
        Integer value = (Integer) optionPane.getValue();
        if (value != null && value != -1) {
            result = value;
        }
        if (result != JOptionPane.CANCEL_OPTION) {
        	color = chooser.getColor();
        	m_input.setValue(color, m_beastObject);
        }
        button.setBackground(color);
	}
}

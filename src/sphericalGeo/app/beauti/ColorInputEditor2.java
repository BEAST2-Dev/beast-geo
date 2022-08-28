package sphericalGeo.app.beauti;




import java.util.Optional;

import beastfx.app.inputeditor.BeautiDoc;
import beast.base.core.BEASTInterface;
import beast.base.core.Input;
import beastfx.app.inputeditor.InputEditor;
import javafx.scene.control.Button;
import javafx.scene.control.ButtonType;
import javafx.scene.control.ColorPicker;
import javafx.scene.control.Dialog;
import javafx.scene.paint.Color;
import javafx.scene.text.Font;
import javafx.scene.text.Text;

public class ColorInputEditor2 extends InputEditor.Base {
	
	Button button;
	
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
        button = new Button("color");
        if (color != null) {
        	button.setStyle("-fx-background-color: " + color + ";");
        }
        button.setOnAction(e->editColor());
        getChildren().add(button);
        
        //add(Box.createHorizontalGlue());
        addValidationLabel();
    } // init
	
	private void editColor() {
		Color color = (Color) m_input.get();
		
		final ColorPicker colorPicker = new ColorPicker();
        colorPicker.setValue(color);
        
        final Text text = new Text("Choose colour for " + m_input.getName());
        text.setFont(Font.font ("Verdana", 20));
        text.setFill(colorPicker.getValue());
        
        Dialog dlg = new Dialog<>();
        dlg.getDialogPane().getChildren().add(colorPicker);
        dlg.getDialogPane().getButtonTypes().add(ButtonType.CANCEL);
        dlg.getDialogPane().getButtonTypes().add(ButtonType.OK);
        
        Optional<ButtonType> result = dlg.showAndWait();
      
        if (result.toString().toLowerCase().contains("ok")) {
        	color = colorPicker.getValue();
        	m_input.setValue(color, m_beastObject);
        }
    	button.setStyle("-fx-background-color: " + color + ";");
	}
}

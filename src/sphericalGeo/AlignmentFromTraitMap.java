package sphericalGeo;

import java.util.ArrayList;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.Input.Validate;
import beast.base.evolution.alignment.Alignment;


@Description("Creates and alignment of continuous data, such as a geographic location, from a trait set")
public class AlignmentFromTraitMap extends Alignment {
	
	public Input<TreeTraitMap> traitInput = new Input<>("traitMap", "trait map to be interpreted as single site alignment");

	TreeTraitMap traitMap;

	public AlignmentFromTraitMap() {
		dataTypeInput.setRule(Validate.OPTIONAL);
		taxonSetInput.setRule(Validate.OPTIONAL);
		stateCountInput.setRule(Validate.OPTIONAL);
		stripInvariantSitesInput.setRule(Validate.OPTIONAL);
		sequenceInput.setRule(Validate.OPTIONAL);
	}

	@Override
    public void initAndValidate() {
		
		if (taxonSetInput.get() != null && taxonSetInput.get().getTaxonCount() == 0) {
			taxonSetInput.setValue(null, this);
		}
    	traitMap = traitInput.get();
    	patternIndex = new int[0];
        counts = new ArrayList<>();
    	if (traitMap == null) { // assume we are in beauti
    		return;
    	}
    	m_dataType = userDataTypeInput.get();
        if (!(m_dataType instanceof LocationDataType)) {
        	throw new IllegalArgumentException("Data type must be a LocationDataType, not " + m_dataType.getClass().getName());
        }

        taxaNames = new ArrayList<>();
        for (String name : traitMap.treeInput.get().getTaxonset().asStringList()) {
        	taxaNames.add(name);
        }
        
        if (traitMap.value.get() == null || traitMap.value.get().matches("^\\s*$")) {
        	// prevent initialisation when in beauti
        	patternIndex = new int[1];
            return;
        }
        

        stateCounts = new ArrayList<>();
        for (@SuppressWarnings("unused") String s : taxaNames) {
        	stateCounts.add(m_dataType.getStateCount());
        }

    }
	
	public TreeTraitMap getTraitMap() {
		return traitMap;
	}
	
	@Override
	public int getSiteCount() {
		return 1;
	}


}

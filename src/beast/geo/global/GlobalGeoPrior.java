package beast.geo.global;


import java.lang.reflect.InvocationTargetException;
import java.lang.reflect.Method;
import java.util.List;
import java.util.Random;

import beast.core.Description;
import beast.core.Distribution;
import beast.core.Input;
import beast.core.Input.Validate;
import beast.core.State;
import beast.evolution.likelihood.GenericTreeLikelihood;
import beast.evolution.tree.Node;
import beast.evolution.tree.TreeInterface;

@Description("Prior on global language tree for the global language tree project")
public class GlobalGeoPrior extends Distribution {
	public Input<GenericTreeLikelihood> locationsInput = new Input<>("likelihood", "approx. trait likelihood containing locations", Validate.REQUIRED);
	public Input<Boolean> hasNorthPoleRestrictionInput = new Input<>("hasNorthPoleRestriction", "flag to indicate internal nodes are not allowed over 70 degrees latitude" , true); 
			
	GenericTreeLikelihood locations;
	TreeInterface tree;
	boolean hasNorthPoleRestriction;
	Method getPositions;
	
	
	@Override
	public void initAndValidate() {
		locations = locationsInput.get();
		try {
			getPositions = locations.getClass().getMethod("getPositions");
		} catch (NoSuchMethodException | SecurityException e) {
			throw new IllegalArgumentException(e);
		}
		
		tree = locations.treeInput.get();
		hasNorthPoleRestriction = hasNorthPoleRestrictionInput.get();
	}
	
	final static double SINGLE_CONSTRAINT_CROSSING_PENALTY = -1000000.0;
	
	@Override
	public double calculateLogP() {
		logP = 0;
		
		double[][] positions;
		try {
			positions = (double [][]) getPositions.invoke(locations);
		} catch (IllegalAccessException | IllegalArgumentException | InvocationTargetException e) {
			throw new RuntimeException(e);
		}//.getPositions();
		
		// count nr of trans-atlantic branches at -30 degrees longitude
		for (Node node : tree.getNodesAsArray()) {
			if (!node.isRoot()) {
				double long0 = positions[node.getNr()][1];
				while (long0 > 180.0) {
					long0 -= 360.0;
				}
				double long1 = positions[node.getParent().getNr()][1];
				while (long1 > 180.0) {
					long1 -= 360.0;
				}
				logP += getPenalty(long0, 0, long1, 0);
			}
			if (hasNorthPoleRestriction && !node.isLeaf()) {
				// discourage internal nodes above 70 latitude
				double lat0 = positions[node.getNr()][0];
				if (lat0 > 70) {
					logP -= (lat0 - 70) * 1000;
				}
			}
		}
		return logP;
	}
	
    public double getCurrentLogP() {
        try {
			logP = calculateLogP();
		} catch (Exception e) {
			e.printStackTrace();
		}
        return logP;
    }


	public static double getPenalty(double long0, double lat0, double long1, double lat1) {
		if (long0 > -150 &&  long0 < -30 && long1 > -30 && long1 < 90) {
			return SINGLE_CONSTRAINT_CROSSING_PENALTY;
		}
		if (long1 > -150 &&  long1 < -30 && long0 > -30 && long0 < 90) {
			return SINGLE_CONSTRAINT_CROSSING_PENALTY;
		}
		return 0.0;
	}

	@Override
	public List<String> getArguments() {return null;}

	@Override
	public List<String> getConditions() {return null;}

	@Override
	public void sample(State state, Random random) {}

}

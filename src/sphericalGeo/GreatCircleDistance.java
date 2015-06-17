package sphericalGeo;


import beast.core.BEASTObject;
import beast.core.Description;
import beast.core.Input;
import beast.core.Input.Validate;
import beast.evolution.alignment.distance.Distance;
import beast.evolution.tree.TreeInterface;
//import beast.evolution.tree.TreeTraitMap;

@Description("Distance between points on a globe") 
public class GreatCircleDistance extends BEASTObject implements Distance {
	public Input<TreeTraitMap> traitInput = new Input<TreeTraitMap>("trait", "trait specifying latitude/longitude locations", Validate.REQUIRED);
			
	TreeTraitMap trait;
	TreeInterface tree;

	final public static double EARTHRADIUS = 6371; // mean radius, according to http://en.wikipedia.org/wiki/Earth_radius
	
	@Override
	public void initAndValidate() throws Exception {
		trait = traitInput.get();
		tree = trait.treeInput.get();
	}

	@Override
	public double pairwiseDistance(int taxon1, int taxon2) {
		double [] loc1 = trait.getTrait(tree, tree.getNode(taxon1));
		double [] loc2 = trait.getTrait(tree, tree.getNode(taxon2));
		return pairwiseDistance(loc1, loc2);
	}
	
	public static double pairwiseDistance(double [] start, double [] end) {
		if (start[0] == end[0] && start[1] == end[1]) {
			return 0.0;
		}
		
		double latitude1 = start[0];
		double longitude1 = start[1];
		double theta1 = (latitude1)*Math.PI/180.0;
		if (longitude1 < 0) longitude1 += 360;
		double phi1 = longitude1 * Math.PI/180;

		double latitude2 = end[0];
		double longitude2 = end[1];
		double theta2 = (latitude2)*Math.PI/180.0;
		if (longitude2 < 0) longitude2 += 360;
		double phi2 = longitude2 * Math.PI/180;
		
		double Deltalambda = phi2 - phi1;
		
		double angle = Math.acos(Math.sin(theta1)*Math.sin(theta2)+Math.cos(theta1) * Math.cos(theta2) * Math.cos(Deltalambda)); 

		//double inverseVariance = 10;
        //double logP =  0.5 * Math.log(angle * Math.sin(angle)) - 0.5 * angle*angle * inverseVariance;
        // double logP = Math.log(Math.sqrt(angle * Math.sin(angle)) * inverseVariance) - 0.5 * angle*angle * inverseVariance;
        //double logP = 0.5 * Math.log(angle * Math.sin(angle)) + 0.5 * Math.log(inverseVariance) - 0.5 * angle*angle * inverseVariance;
        //return logP;

		return angle * EARTHRADIUS;
	}

}

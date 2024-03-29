package sphericalGeo;


import beast.base.core.BEASTObject;
import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.Input.Validate;
import beast.base.evolution.distance.Distance;
import beast.base.evolution.tree.TreeInterface;
//import beast.base.evolution.tree.TreeTraitMap;

@Description("Distance between points on a globe") 
public class GreatCircleDistance extends BEASTObject implements Distance {
	public Input<TreeTraitMap> traitInput = new Input<>("trait", "trait specifying latitude/longitude locations", Validate.REQUIRED);
			
	TreeTraitMap trait;
	TreeInterface tree;

	final public static double EARTHRADIUS = 6371; // mean radius, according to http://en.wikipedia.org/wiki/Earth_radius
	
	@Override
	public void initAndValidate() {
		trait = traitInput.get();
		tree = trait.treeInput.get();
	}

	// returns distance in km
	@Override
	public double pairwiseDistance(int taxon1, int taxon2) {
		double [] loc1 = trait.getTrait(tree, tree.getNode(taxon1));
		double [] loc2 = trait.getTrait(tree, tree.getNode(taxon2));
		return pairwiseDistance(loc1, loc2);
	}
	
	// returns distance in km
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

	public static void main(String[] args) {
		System.out.println("vancouver - berlin " + pairwiseDistance(new double[]{11.8, -15.5}, new double[]{10.76262, 106.660172}));
		System.out.println("galapagos - nairobi " + pairwiseDistance(new double[]{-0.777259,-91.142578}, new double[]{-1.280423, 36.816311}));
		
		
		System.out.println("vancouver - berlin " + pairwiseDistance(new double[]{49.246292, -123.116226}, new double[]{52.516266,13.377775}));
	}
}

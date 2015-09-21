package sphericalGeo;

import java.util.HashMap;
import java.util.Map;
import java.util.Set;

import beast.evolution.operators.DistanceProvider;
import org.apache.commons.math3.util.FastMath;

import beast.core.BEASTObject;
import beast.core.Description;
import beast.core.Input;
import beast.core.Input.Validate;
import beast.evolution.alignment.Alignment;
import beast.evolution.operators.AttachOperator;
import beast.evolution.tree.Node;
import beast.evolution.tree.TreeInterface;

@Description("Provide distance between geographical locations -- to be used by AttachOperator")
public class CircleDistanceProvider extends BEASTObject implements DistanceProvider {
	enum Method {
		DISTANCE("distance"),
		SQRT("sqrt"),
		ARC("arc");

		Method(final String name) {
			this.ename = name;
		}

		public String toString() {
			return ename;
		}

		private final String ename;
	}
    public Input<TreeInterface> treeInput = new Input<TreeInterface>("tree", "phylogenetic beast.tree with sequence data in the leafs", Validate.REQUIRED);
    public Input<Alignment> dataInput = new Input<Alignment>("data", "sequence data for the beast.tree", Validate.REQUIRED);

	public Input<Transformer> transformerInput = new Input<Transformer>("transformer","landscape transformer to capture some inheterogenuity in the diffusion process");
	public Input<Method> distMethod = new Input<Method>("method", "for calculating distance between clade positions (for operator weights). sqrt takes " +
	         "square root of distance (default distance)",  Method.DISTANCE, Method.values());

	private TreeInterface tree;
	private double [][] position;
    private Method distanceMethod;

	@Override
	public void initAndValidate() throws Exception {
		distanceMethod = distMethod.get();
		tree = treeInput.get();
		
		// initialise leaf positions
		position = new double[tree.getNodeCount()][2];
		AlignmentFromTraitMap data = (AlignmentFromTraitMap) dataInput.get();
		TreeTraitMap traitMap = data.getTraitMap(); 
		Node [] nodes = tree.getNodesAsArray();
		for (int i = 0; i < tree.getLeafNodeCount(); i++) {
			position[i] = traitMap.getTrait(tree, nodes[i]);
			if (transformerInput.get() != null) {
				position[i] = transformerInput.get().project(position[i][0], position[i][1]);
			}
		}

	}
	
	class LocationData implements DistanceProvider.Data {
		final static int MAX_ITER = 0;
		final static double MIN_EPSILON = 0.001;

		double[] position;
		int weight;

		public LocationData(double[] pos) {
			position = pos;
			weight = 1;
		}
		public LocationData() {
			position = new double[3];
	        weight = 0;
		}
	}
	
	@Override
	public Map<String, DistanceProvider.Data> init(Set<String> taxa) {
		final HashMap<String, DistanceProvider.Data> m = new HashMap<>();
		int count = 0;
		for (int i = 0; i < tree.getLeafNodeCount(); i++) {
			Node n = tree.getNode(i);
			final String taxon = n.getID();
			if( taxa.contains(taxon) ) {
				final double[] xyz = SphericalDiffusionModel.spherical2Cartesian(position[i][0], position[i][1]);
				m.put(taxon, new LocationData(xyz));
				count += 1;
			}
		}
		if( count != taxa.size() ) {
			return null;
		}
		return m;
	}

//	@Override
//	public AttachOperatorNew.DistanceProvider.Data combine(AttachOperatorNew.DistanceProvider.Data x1, AttachOperatorNew.DistanceProvider.Data x2) {
//		LocationData d1 = (LocationData) x1;
//		LocationData d2 = (LocationData) x2;
//		LocationData r = new LocationData();
//
//		final double w = r.weight = d1.weight + d2.weight;
//        assert d1.weight >= 0 &&  d2.weight >= 0 && w > 0;
//		for(int i = 0; i < 3; ++i) {
//			r.position[i] = (d1.position[i] * d1.weight + d2.position[i] * d2.weight) / w;
//		}
//		return r;
//	}

    @Override
    public Data empty() {
        return new LocationData();
    }

    @Override
    public void clear(Data d) {
        ((LocationData)d).weight = 0;
    }

    @Override
    public void update(Data info, Data with) {
        LocationData d1 = (LocationData) info;
        LocationData d2 = (LocationData) with;
        assert d1.weight >= 0 &&  d2.weight > 0;

        if( d1.weight == 0 ) {
            System.arraycopy(d2.position, 0, d1.position, 0, 3);
            d1.weight = d2.weight;
        } else {
            final int w = d1.weight + d2.weight;
            for (int i = 0; i < 3; ++i) {
                d1.position[i] = (d1.position[i] * d1.weight + d2.position[i] * d2.weight) / w;
            }
            d1.weight = w;
        }
        assert d1.weight > 0;
    }

    @Override
    public double dist(DistanceProvider.Data info1, DistanceProvider.Data info2) {
        LocationData d1 = (LocationData) info1;
        LocationData d2 = (LocationData) info2;
        double s = 0;
        for(int k = 0; k < 3; ++k) {
            double x = (d1.position[k] - d2.position[k]);
            s += x*x;
        }
        s = (s == 0) ? 1e-8 : s;
        switch (distanceMethod) {
            case DISTANCE: break;
            case SQRT: s = Math.sqrt(s); break;
            case ARC: s =  FastMath.asin(FastMath.sqrt(s) / 2); break;
        }
        return s;
    }

}

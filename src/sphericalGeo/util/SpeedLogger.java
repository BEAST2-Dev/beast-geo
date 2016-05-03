package sphericalGeo.util;

import java.io.PrintStream;

import beast.core.BEASTObject;
import beast.core.Loggable;
import beast.core.Param;
import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;
import sphericalGeo.GreatCircleDistance;
import sphericalGeo.LocationProvider;

public class SpeedLogger extends BEASTObject implements Loggable {
	LocationProvider locations;
	Tree tree;
	
	public SpeedLogger() {}
	public SpeedLogger(@Param(name="locations", description="position provider used to calculate speed") LocationProvider locations,
			@Param(name="tree",description="tree providing branch lengths to determine total time") Tree tree) {
		this.locations = locations;
		this.tree = tree;
	}

	public LocationProvider getLocations() {
		return locations;
	}
	public void setLocations(LocationProvider locations) {
		this.locations = locations;
	}
	public Tree getTree() {
		return tree;
	}
	public void setTree(Tree tree) {
		this.tree = tree;
	}
	
	@Override
	public void initAndValidate() {
	}

	@Override
	public void init(PrintStream out) {
		out.append((getID() != null ? getID() : "speed") + "\t");
	}

	@Override
	public void log(int sample, PrintStream out) {
		double sumOfTime = 0;
		double sumOfDistance = 0;
		for (Node node : tree.getNodesAsArray()) {
			if (!node.isRoot()) {
				sumOfTime += node.getLength();
				double [] start = locations.getPosition(node.getParent().getNr());
				double [] end = locations.getPosition(node.getNr());
				sumOfDistance += GreatCircleDistance.pairwiseDistance(start, end);
			}
		}
		double speed = sumOfDistance / sumOfTime;
		out.append(speed + "\t");
	}

	@Override
	public void close(PrintStream out) {
	}

} // SpeedLogger

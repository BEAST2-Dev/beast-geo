package sphericalGeo.util.treeset;

import java.io.IOException;

import beast.base.evolution.tree.Tree;

abstract public class TreeSet {
	public abstract boolean hasNext();
	public abstract Tree next() throws IOException;
	public abstract void reset() throws IOException;
	public abstract int size();
} 

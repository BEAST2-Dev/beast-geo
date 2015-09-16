package sphericalGeo.scapetoad;

import java.util.Vector;

public interface StatusTracker {
	
	public void updateRunningStatus(final int progress, final String label1, final String label2);

	public void setComputationError (String title, String message, String stackTrace);

	default public void goToFinishedPanel() {}
	
	default public Vector getSimultaneousLayers () { return new Vector();}
	
	default public Vector getConstrainedDeformationLayers() { return new Vector();}

}

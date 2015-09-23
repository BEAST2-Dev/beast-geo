package sphericalGeo.scapetoad;


import java.awt.Color;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.text.DecimalFormat;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import javax.xml.parsers.DocumentBuilderFactory;

import org.w3c.dom.Node;
import org.w3c.dom.NodeList;

import com.vividsolutions.jts.geom.Coordinate;
import com.vividsolutions.jts.geom.Envelope;
import com.vividsolutions.jts.geom.GeometryFactory;
import com.vividsolutions.jts.geom.LinearRing;
import com.vividsolutions.jts.geom.MultiPolygon;
import com.vividsolutions.jts.geom.Polygon;
import com.vividsolutions.jump.feature.AttributeType;
import com.vividsolutions.jump.feature.BasicFeature;
import com.vividsolutions.jump.feature.Feature;
import com.vividsolutions.jump.feature.FeatureCollection;
import com.vividsolutions.jump.feature.FeatureDataset;
import com.vividsolutions.jump.feature.FeatureSchema;
import com.vividsolutions.jump.io.IllegalParametersException;
import com.vividsolutions.jump.workbench.model.Layer;
import com.vividsolutions.jump.workbench.model.LayerManager;

import sphericalGeo.Transformer;
import sphericalGeo.scapetoad.core.Cartogram;
import sphericalGeo.scapetoad.core.CartogramGrid;
import sphericalGeo.scapetoad.core.CartogramLayer;
import sphericalGeo.scapetoad.core.IOManager;
import sphericalGeo.scapetoad.gui.AppContext;
import beast.core.BEASTObject;
import beast.core.Input;
import beast.core.Input.Validate;
import beast.core.util.Log;

public class ScapeToadTransfomer extends BEASTObject implements StatusTracker, Transformer {
	final static String DEFAULT_ATT_IDENTIFIER = "NAME";
	final static String ATTRIBUTE = "weight";

	public Input<File> shapeFileInput = new Input<File>("shapeFile", "file containing regions in .shp or .kml format", Validate.REQUIRED);

	public Input<String> attIdentifierInput = new Input<String>("attIdentifier", "name of attribute in shape file that identifies a region", DEFAULT_ATT_IDENTIFIER);
	public Input<String> weightInput = new Input<String>("value", "comma separated list of 'id=value' pairs where values representing relative size for each regions with ID from the shape file");
	public Input<String> weightIdentifiedInput = new Input<String>("weightIdentifier", "name of attribute in shape file that identifies a weight", Validate.XOR, weightInput);
	
	public Input<File> cartogramFileInput = new Input<File>("cartogram", "File containing cartogram. If specified and file does not exist, the cartogram will be calculated and stored there.");
	

	public Input<Boolean> isDensityInput = new Input<>("isDensity", "if true, weights are interpreted as densities, otherwise they are interpreted as mass", true);
	public Input<Integer> deformationAmountInput = new Input<>("deformationAmount", "Defines the amount of deformation. This is an integer value between" +
			" 0 and 100. The default value is 50", 50);
	public Input<Boolean> isAdvancedInput = new Input<>("advanced", "if true the advanced parameters should be taken in account, otherwise these are estimated.", false);
		

	public Input<Integer> gridSizeXInput = new Input<>("gridSizeX", "A first grid is applied to the main transformation layer. "
			+ "This rectangular grid is defined by the number of columns. "
			+ "Higher numbers produce denser grids and thus a better cartogram quality. "
			+ "However, a denser grid also implies a longer treatment. "
			+ "Default 200, ignored if advanced=false", 200);
	public Input<Integer> gridSizeYInput = new Input<>("gridSizeY", "As gridSizeX but for number of rows, default 200, ignored if advanced=false", 200);

	public Input<Integer> diffusionGridSizeInput = new Input<>("diffusionGridSize", "A second grid is applied to the main transformation layer. "
			+ "This square grid is defined by the number of rows. "
			+ "Denser grids imply again a better cartogram quality but longer computation times."
			+ "must be a power of 2, default 128, ignored if advanced=false", 128);

	public Input<Integer> diffusionIterationsInput = new Input<>("diffusionIterations", "The second grid is transformed with the Gastner/Newman diffusion "
			+ "algorithm, which can be run several times to obtain a higher transformation quality. "
			+ "Higher numbers of iterations also imply longer treatment times."
			+ "Default 3, ignored if advanced=false", 3);
	
	public Input<Integer> gridLayerSizeInput = new Input<>("gridLayerSize", "The visualisation grid layer size. "
			+ "This is the grid which is produced for visual effect only in the files /tmp/scapetoad.svg and /tmp/deformationLayer.svg. "
			+ "Default 100.", 100);
	
	
	
	CartogramGrid cartogramGrid;
	Map<String,Double> map;
	
	@Override
	public void initAndValidate() throws Exception {
		// input sanity checks
		if (deformationAmountInput.get() < 0 || deformationAmountInput.get() > 100) {
			throw new RuntimeException("deformationAmount must be between 0 and 100");
		}
		if (gridSizeXInput.get() < 0) {
			throw new RuntimeException("gridSizeX must be between positive");
		}
		if (gridSizeYInput.get() < 0) {
			throw new RuntimeException("gridSizeY must be between positive");
		}
		int i = diffusionGridSizeInput.get();
		do {
			if (i % 2 != 0) {
				throw new RuntimeException("diffusionGridSize must be a power of 2");
			}
			i /= 2;
		} while (i > 1);
		if (gridLayerSizeInput.get() < 0) {
			throw new RuntimeException("gridLayerSize must be between positive");
		}
		
		
        // parse region values and store in map
		if (weightInput.get() != null) {
			map = new HashMap<String, Double>();
			String ts = weightInput.get();
			String [] traits = ts.split(",");
	        for (String trait : traits) {
	            trait = trait.replaceAll("\\s+", " ");
	            String[] sStrs = trait.split("=");
	            if (sStrs.length != 2) {
	                throw new Exception("could not parse trait: " + trait);
	            }
	            String ID = sStrs[0].trim();
	            String value = sStrs[1].trim();
	            double dvalue = Double.parseDouble(value);
		        map.put(ID, dvalue);
	        }
		}

		
		
		// Create a new layer manager.
		AppContext.layerManager = new LayerManager();
		AppContext.layerManager.addCategory("Original layers");

		Layer lyr = createLayer(shapeFileInput.get().getPath());
		
		File file = null;
		if (cartogramFileInput.get() != null) {
			file = cartogramFileInput.get();
		}
		if (file != null && file.exists()) {
			FileInputStream fis = new FileInputStream(file);
			ObjectInputStream in = new ObjectInputStream(fis);
			cartogramGrid = (CartogramGrid) in.readObject();
			in.close();
		} else {
			Cartogram cartogram;
			cartogram = calcCartogram(lyr);
			
			updateRunningStatus(975, "Writing file /tmp/deformationLayer.svg", "");	
			Layer deformationLAyer = AppContext.layerManager.getLayer("Deformation grid");
			Layer[] lyrs = new Layer[1];
			lyrs[0] = deformationLAyer;
			IOManager.writeSvg(lyrs, "/tmp/deformationLayer.svg" );
			
			updateRunningStatus(990, "Writing file /tmp/toadmap.svg", "");	
			Layer[] lyrs2 = new Layer[2];
			lyrs2[0] = deformationLAyer;
			lyrs2[1] = lyr;
			IOManager.writeSvg(lyrs2, "/tmp/toadmap.svg" );

			updateRunningStatus(1000, "Done!", "");	
			
			if (file != null) {
				FileOutputStream fos = new FileOutputStream(file);
				ObjectOutputStream out = new ObjectOutputStream(fos);
				out.writeObject(cartogram.getGrid());
				out.close();
			}
			cartogramGrid = cartogram.getGrid();
		}
	}

	@Override
	public double [] project(double lat0, double long0) {
		double [] p = cartogramGrid.projectPoint(long0, lat0);
		return new double[]{p[1], p[0]};
	}
	
	@Override
	public double [] projectInverse(double lat0, double long0) {
		double [][] mNodeX = cartogramGrid.getXCoordinates();
		double [][] mNodeY = cartogramGrid.getYCoordinates();
		Coordinate [] points = new Coordinate[5];
		for (int i = 0; i < 5; i++) {
			points[i] = new Coordinate();
		}
		for (int i = 0; i < mNodeX.length - 1; i++) {
			for (int j = 0; j < mNodeX[0].length - 1; j++) {
				points[0].x = mNodeX[i][j];
				points[0].y = mNodeY[i][j];
				points[1].x = mNodeX[i+1][j];
				points[1].y = mNodeY[i+1][j];
				points[2].x = mNodeX[i+1][j+1];
				points[2].y = mNodeY[i+1][j+1];
				points[3].x = mNodeX[i][j+1];
				points[3].y = mNodeY[i][j+1];
				points[4].x = mNodeX[i][j];
				points[4].y = mNodeY[i][j];
				if (insideQuadrangle(long0, lat0, points)) {
					double [] portions = sphericalGeo.scapetoad.core.Geometry.reverseprojection(long0, lat0, 
							points[0].x, points[0].y, 
							points[1].x, points[1].y, 
							points[2].x, points[2].y, 
							points[3].x, points[3].y
							);
					Envelope envelope = cartogramGrid.envelope(); 
					double lat0org = envelope.getMinY() + envelope.getHeight() * (portions[1] + j) / mNodeX[0].length;    
					double long0org = envelope.getMinX() + envelope.getWidth() * (portions[0] + i) / mNodeX.length;
					return new double[]{lat0org, long0org};
				}
			}
		}
		return null; 
	}
	
    /**
     * Return true if the given point is contained inside the boundary.
     * See: http://www.ecse.rpi.edu/Homepages/wrf/Research/Short_Notes/pnpoly.html
     * @param test The point to check
     * @return true if the point is inside the boundary, false otherwise
     *
     */
    public boolean insideQuadrangle(double testx, double testy, Coordinate [] points) {
      int i;
      int j;
      boolean result = false;
      for (i = 0, j = points.length - 1; i < points.length; j = i++) {
        if ((points[i].y > testy) != (points[j].y > testy)) {
        	if (testx < (points[j].x - points[i].x) * (testy - points[i].y) / (points[j].y-points[i].y) + points[i].x) {
        		result = !result;
        	}
         }
      }
      return result;
    }

	private Layer createLayer(String path) throws Exception {
		Log.warning.println("Processing file " + shapeFileInput.get().getPath());
		if (path.toLowerCase().endsWith(".shp")) {
			Layer lyr = IOManager.readShapefile(shapeFileInput.get().getPath());
			return lyr;
		} else if (path.toLowerCase().endsWith(".kml")) {
			// Define a layer fill color.
			Color lc = Color.GREEN;

			// Read the Shape file.
			String ln = IOManager.fileNameFromPath(path);;
			FeatureCollection fc = readKML(path);
			
			// If there is no category "Original layers", we add one.
			if (AppContext.layerManager.getCategory("Original layers") == null)
				AppContext.layerManager.addCategory("Original layers");
			
			// Add the layer to the "Original layers" category.
			Layer lyr = new Layer(ln, lc, fc, AppContext.layerManager);
			lyr = AppContext.layerManager.addLayer("Original layers", ln, fc);
			return lyr;
		} else {
			throw new RuntimeException("Unrecoginised file type for shapeFile -- should be .shp or .kml");
		}
	}
	
    /**
     * Main method to read a shapefile.  Most of the work is done in the org.geotools.* package.
     *
     *@param dp 'InputFile' or 'DefaultValue' to specify output .shp file.
     *
     */
    private FeatureCollection readKML(String kmlfile)
        throws IllegalParametersException, Exception {
    	Map<String,List<List<Double>>> map = parseKML(kmlfile);

    	//List<com.vividsolutions.jts.geom.Geometry> geometryList = new ArrayList<>();
        GeometryFactory factory = new GeometryFactory();
        
        FeatureSchema fs = new FeatureSchema();
        fs.addAttribute(attIdentifierInput.get(), AttributeType.STRING);
        fs.addAttribute("density", AttributeType.DOUBLE);
        fs.addAttribute("GEOMETRY", AttributeType.GEOMETRY);

        FeatureCollection featureCollection = new FeatureDataset(fs);

        for (String name : map.keySet()) {
        	List<List<Double>> polygons = map.get(name);
        	List<Polygon> geoPolygons = new ArrayList<>();
        	for (List<Double> polygon : polygons) {
        		Coordinate [] coordinates = new Coordinate[polygon.size()/2];
        		for (int i = 0; i < coordinates.length; i++) {
        			coordinates[i] = new Coordinate(polygon.get(i*2), polygon.get(i*2+1));
        		}
        		LinearRing ring = factory.createLinearRing(coordinates);
        		Polygon geoPolygon = factory.createPolygon(ring, null);
        		geoPolygons.add(geoPolygon);
        	}
        	MultiPolygon multiPolygon = factory.createMultiPolygon(geoPolygons.toArray(new Polygon[]{}));
    		//geometryList.add(multiPolygon);

    		Feature feature = new BasicFeature(fs);
    		feature.setAttribute(0, name);
    		feature.setAttribute(1, map.get(name));
            feature.setGeometry(multiPolygon);
            featureCollection.add(feature);
    	}
    	
        return featureCollection;
    }

	private Map<String,List<List<Double>>> parseKML(String path) throws Exception {
		DocumentBuilderFactory factory = DocumentBuilderFactory.newInstance();
		factory.setValidating(false);
		org.w3c.dom.Document doc = factory.newDocumentBuilder().parse(path);
		doc.normalize();

		Map<String, List<List<Double>>> coordinates = new HashMap();

		// grab 'coordinates' elements out of the KML file
		NodeList oPlacemarks = doc.getElementsByTagName("Placemark");
		for (int iNode = 0; iNode < oPlacemarks.getLength(); iNode++) {
			Node oPlaceMark = oPlacemarks.item(iNode);
			String id = getChild(oPlaceMark, "name").getTextContent();
			List<List<Double>> polygons = new ArrayList<List<Double>>();
			getPolygons(oPlaceMark, polygons);
			coordinates.put(id, polygons);
		}
		return coordinates;
	}

	private void getPolygons(Node node, List<List<Double>> polygons) {
		NodeList children = node.getChildNodes();
		for (int i = 0; i < children.getLength(); i++) {
			Node child = children.item(i);
			if (child.getNodeType() == Node.ELEMENT_NODE) {
				if (child.getNodeName().equals("coordinates")) {
					String sCoordinates = child.getTextContent();
					List<Double> polygon = new ArrayList<>();
					String[] sStrs = sCoordinates.split("\\s+");
					for (String sStr : sStrs) {
						if (sStr.contains(",")) {
							String[] sCoords = sStr.split(",");
							polygon.add(Double.parseDouble(sCoords[0].trim()));
							polygon.add(Double.parseDouble(sCoords[1].trim()));
						}
					}
					polygons.add(polygon);
				} else {
					getPolygons(child, polygons);
				}
			}
		}
	}

	private Node getChild(Node node, String name) {
		NodeList children = node.getChildNodes();
		for (int i = 0; i < children.getLength(); i++) {
			Node child = children.item(i);
			if (child.getNodeType() == Node.ELEMENT_NODE && child.getNodeName().equals(name)) {
				return child;
			}
		}
		return null;
	}



	Cartogram calcCartogram(Layer lyr) throws Exception {
		// Get the name of the selected layer.
		String selectedLayer = lyr.getName();
		
		// Get the name of the selected attribute.
		CartogramLayer.addAttribute(lyr, ATTRIBUTE, AttributeType.DOUBLE);
		
		
        // set up attribute
        String attIdentifier = attIdentifierInput.get();
        String weightIdentifier = weightIdentifiedInput.get();
		for (Object o : lyr.getFeatureCollectionWrapper().getFeatures()) {
			Feature feature = (Feature) o;
			String id = null;
			try {
				id = (String) feature.getAttribute(attIdentifier);
			} catch (IllegalArgumentException e) {
				String options = "";
				for (int i = 0; i < feature.getSchema().getAttributeCount(); i++) {
					options += feature.getSchema().getAttributeName(i) + " ";
				}
				throw new RuntimeException("The attribute '" + attIdentifier + "' does not exist in the shape file.\n"+
						"Either add them to the regions, or use another attIdentifier.\n"+
						"Available features: " + options);
			}
			id = id.trim();
			if (map != null) {
				if (map.containsKey(id)) {
					double weight = map.get(id);
					feature.setAttribute(ATTRIBUTE, weight);
					map.remove(id);
				} else {
					Log.warning.println("WARNING: No value in map for " + id);
					feature.setAttribute(ATTRIBUTE, 0.0);
				}
			} else {
				double weight = 0;
				try {
					weight = Double.parseDouble(feature.getAttribute(weightIdentifier).toString());
				} catch (IllegalArgumentException e) {
					String options = "";
					for (int i = 0; i < feature.getSchema().getAttributeCount(); i++) {
						options += feature.getSchema().getAttributeName(i) + " ";
					}
					throw new RuntimeException("The weight attribute '" + weightIdentifier + "' does not exist in the shape file.\n"+
							"Either add them to the regions, or use another attIdentifier.\n"+
							"Available features: " + options);
				}
				feature.setAttribute(ATTRIBUTE, weight);
			}
		}
		
		if (map != null) {
			for (String key : map.keySet()) {
				Log.warning.println("WARNING: No region in map for " + key);
			}
		}
		
		String selectedAttribute = ATTRIBUTE; //mCartogramWizard.getCartogramAttributeName();
		
		// Get the attribute type (population or density value).
		boolean isDensityValue = isDensityInput.get();
		
		
		//mCartogramWizard.setMissingValue(""/*mCartogramWizard.getPanelTwo().getMissingValue()*/);
		
		
		// Create a new cartogram instance and set the parameters.
		AppContext.setStatusTracker(this);
		Cartogram cg = new Cartogram(this);
		cg.setLayerManager(AppContext.layerManager);
		cg.setMasterLayer(selectedLayer);
		cg.setMasterAttribute(selectedAttribute);
		cg.setMasterAttributeIsDensityValue(isDensityValue);
		cg.setMissingValue(""/*mCartogramWizard.getMissingValue()*/);
		cg.setSlaveLayers(null /*mCartogramWizard.getSimultaneousLayers()*/);
		cg.setConstrainedDeformationLayers(null/*mCartogramWizard.getConstrainedDeformationLayers()*/);
			
		
		cg.setAmountOfDeformation(deformationAmountInput.get()/*mCartogramWizard.getAmountOfDeformation()*/);
		
		cg.setAdvancedOptionsEnabled(isAdvancedInput.get()/*mCartogramWizard.getAdvancedOptionsEnabled()*/);

		cg.setGridSize(gridSizeXInput.get(), gridSizeYInput.get());
		
		cg.setDiffusionGridSize(diffusionGridSizeInput.get()/*mCartogramWizard.getDiffusionGridSize()*/);

		cg.setDiffusionIterations(diffusionIterationsInput.get()/*mCartogramWizard.getDiffusionIterations()*/);
		
		// Set the parameters for the deformation grid layer.
		cg.setCreateGridLayer(true/*mCartogramWizard.getCreateGridLayer()*/);
		cg.setGridLayerSize(gridLayerSizeInput.get()/*mCartogramWizard.getDeformationGridSize()*/);
		
		// Set the parameters for the legend layer.
		// We have to estimate the legend values.
		if (isDensityValue)
			cg.setCreateLegendLayer(false);
		else
		{
			cg.setCreateLegendLayer(true);
		}
		
		//CartogramWizard.setCartogram(cg);
		
		// Start the cartogram computation.
		cg.construct();
		cg.finished();
		
		return cg;
	}

	
	@Override
	public void updateRunningStatus(int progress, String label1, String label2) {
		NumberFormat formatter = new DecimalFormat("#0.0"); 
		Log.warning.println(formatter.format(progress/10.0) + "%: " + (label2.length() > 0 ? label2 : label1));
	}

	@Override
	public void setComputationError(String title, String message, String stackTrace) {
		Log.err.println(title);
		Log.err.println(message);
		Log.err.println(stackTrace);
	}

	
	
	public static void main(String[] args) throws Exception {
//		ScapeToadTransfomer st = new ScapeToadTransfomer();
//		st.initByName("shapeFile", "data/x.shp", "attIdentifier", "Name", 
//				"value", "Netherlands=17,Belgium=14,Luxembourg=1,France=50,Germany=80"
//				//,"cartogram", "data/cartogram.ser"
//				);

		ScapeToadTransfomer st2 = new ScapeToadTransfomer();
//		st2.initByName("shapeFile", "data/x.kml", "attIdentifier", "Name", 
//				"value", "Netherlands=17,Belgium=14,Luxembourg=1,France=50,Germany=80"
//				//,"cartogram", "data/cartogram2.ser"
//				);
		st2.initByName("shapeFile","/Users/remco/ie/geo/TM_WORLD_BORDERS_SIMPL-0.3.shp" ,
	            	"attIdentifier", "NAME", 
	            	"weightIdentifier", "POP2005" ,
	            	"cartogram", "/Users/remco/ie/geo/cartogram-population.ser");
		double [] p = st2.project(37.9683334850582, 23.5258329215863);
		double [] p0 = st2.projectInverse(p[0], p[1]);
	}
	
}

package sphericalGeo;

import beast.core.Description;
import beast.evolution.datatype.DataType;

@Description("Datatype for representing geographic locations in 2 dimension, latitude x longitude")
public class LocationDataType extends DataType.Base {
	
	@Override
	public void initAndValidate() throws Exception {
		// nothing to do
	}
	
	public int getDimension() {
		return 2;
	}
	
	/** by default assume comma separated list of values **/
	public double[] string2values(String sValue) throws Exception {
		String [] strs = sValue.trim().split(",");
		if (strs.length != 2) {
			throw new Exception ("Expected 2 comma separated numbers, not " + strs.length);
		}
		double [] values = new double [strs.length];
		for (int i = 0; i < strs.length; i++) {
			try {
			values[i] = Double.parseDouble(strs[i].trim());
			} catch (NumberFormatException e) {
				throw new Exception("String is not a comma separated list of numbers");
			}
		}
		return values;
	}
	
	
	@Override
	public String getTypeDescription() {
		// TODO Auto-generated method stub
		return null;
	}


}

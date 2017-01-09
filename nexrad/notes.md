For the nexrad analysis here is some psedo code:

1. For each hourly time step for each day
	1. Estimate kriged surface from hourly rain gauge data	
		2. For each watershed, get average of surface
		3. For each watershed, get the average of the nexrad raster
		4. Compare the two to get the difference between the rainfall estimated from the rain gauges and the rainfall estimated from the radar


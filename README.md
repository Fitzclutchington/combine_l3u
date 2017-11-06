#Hourly L3U Concatenation

This program concatenates a day of granule per x minute data into granule per
hour data.  The output file contains the following data:

* Sea Surface Temperature
* Satellite Zenith Angle
* Retrieval Time
* Day Mask

The SST and Day mask data are taken dirrectly from the input L3U files, while 
the SZA and retrieval times are unquantized and remapped from corresponding L2P
files.

To run:
`python combine_l3u.py input_l3u_folder input_l2p_folder output_concatenated_l3u_folder`


dir = "Z:/BIOP-HRM/guiet/Deconvolved/STED/";

image_name ="test-STED03_5d1eff185885c_hrm" ;

run("Bio-Formats Importer", "open="+dir+image_name+".ics color_mode=Default rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT");
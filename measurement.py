class Measurements: 
    def __init__(self):
        self.measures_dic = {
            "ID" : ["count", "Morphology"],
            "x center [micron]": ["x-center coordinate [µm]", "Spatial Parameters"],
            "y center [micron]": ["y-center coordinate [µm]", "Spatial Parameters"],
            "z center [micron]": ["z-center coordinate [µm]", "Spatial Parameters"],
            "Volume [voxel]": ["Cilium volume [voxel]", "Morphology"],
            "Volume [micron^3]": [r"Cilium volume [$\mu m^3$]", "Morphology"],
            "# surface voxels": ["Cilium surface [voxels]", "Morphology"],
            "Surface [micron^2]": [r"Cilium surface voxels [$\mu m^2$]", "Morphology"],
            "Shape complexity index": ["Shape complexity [a.u.]", "Morphology"],
            "Sphere radius [micron]": ["Sphere radius [µm]", "Morphology"],
            "Maximum span [micron]": ["Maximum span [µm]", "Morphology"],
            "A: Colocalized volume [micron^3] (if channel in input image was background-removed)": [r"Colocalized volume [$\mu m^3$]", "Channel A/B:"],
            "A: Colocalized volume [% total volume] (if channel in input image was background-removed)": ["Colocalized volume [%]", "Channel A/B:"],
            "B: Colocalized volume [micron^3] (if channel in input image was background-removed)": [r"Colocalized volume [$\mu m^3$]", "Channel A/B:"],
            "B: Colocalized volume [% total volume] (if channel in input image was background-removed)": ["Colocalized volume [%]", "Channel A/B:"],
            "A: Colocalized compared to BG volume [micron^3]": [r"Colocalized compared to BG volume [$\mu m^3$]", "Channel A/B:"],
            "A: Colocalized compared to BG volume [% total volume]": ["Colocalized compared to BG volume [%]", "Channel A/B:"],
            "B: Colocalized compared to BG volume [micron^3]": [r"Colocalized compared to BG volume [$\mu m^3$]", "Channel A/B:"],
            "B: Colocalized compared to BG volume [% total volume]": ["Colocalized compared to BG volume [%]", "Channel A/B:"],
            "minimum intensity (in reconstruction channel)": ["Minimum intensity [a.u.]", "Cilia Markers"],
            "maximum intensity (in reconstruction channel)": ["Maximum intensity [a.u.]", "Cilia Markers"],
            "average intensity of the 10% of voxels with highest intensity (in reconstruction channel)": ["Average intensity [a.u.]", "Cilia Markers"],
            "average intensity (in reconstruction channel)": ["Average intensity [a.u.]", "Cilia Markers"],
            "SD of intensity (in reconstruction channel)": ["Standard deviation of intensity [a.u.]", "Cilia Markers"],
            "minimum A intensity": ["Minimum intensity [a.u.]", "Channel A/B:"],
            "maximum A intensity": ["Maximum intensity [a.u.]", "Channel A/B:"],
            "average A intensity of the 10% of voxels with highest A intensity": ["Average intensity [a.u.]", "Channel A/B:"],
            "average A intensity": ["Average intensity [a.u.]", "Channel A/B:"],
            "SD of A intensity": ["SD of intensity [a.u.]", "Channel A/B:"],
            "minimum B intensity": ["Minimum intensity [a.u.]", "Channel A/B:"],
            "maximum B intensity": ["Maximum intensity [a.u.]", "Channel A/B:"],
            "average B intensity of the 10% of voxels with highest B intensity": ["Average intensity [a.u.]", "Channel A/B:"],
            "average B intensity": ["Average intensity [a.u.]", "Channel A/B:"],
            "SD of B intensity": ["Standard deviation of intensity [a.u.]", "Channel A/B:"],
            "# of found skeletons (quality parameter)" : ["Count", "Morphology"],
            "# branches (quality parameter)" : ["Count", "Morphology"],
            "tree length [micron] (quality parameter)" : ["Ciliary length [µm]", "Morphology"],
            "cilia length [micron] (largest shortest path of largest skeleton)" : ["Ciliary length [µm]", "Morphology"],
            "orientation vector x [micron] (vector from first to last skeleton point)" : ["Orientation vector x [µm]", "Spatial Parameters"],
            "orientation vector y [micron] (vector from first to last skeleton point)" : ["Orientation vector y [µm]", "Spatial Parameters"],
            "orientation vector z [micron] (vector from first to last skeleton point)" : ["Orientation vector z [µm]", "Spatial Parameters"],
            "cilia bending index (arc length of cilium / eucledian distance of first and last skeleton point)" : ["Cilia Bending Index", "Morphology"],
            "Intensity threshold A" : ["Intensity threshold [a.u.]", "Channel A/B:"],
            "Intensity threshold B" : ["Intensity threshold [a.u.]", "Channel A/B:"],
            "Intensity threshold Basal Stain" : ["Intensity threshold [a.u.]", "Channel A/B:"],
            "Integrated A intensity" : ["Integrated intensity [a.u.]", "Channel A/B:"],
            "Average A intensity on center line" : ["Average intensity [a.u.]", "Channel A/B:"],
            "Integrated B intensity" : ["Integrated intensity [a.u.]", "Channel A/B:"],
            "Average B intensity on center line" : ["Average intensity [a.u.]", "Channel A/B:"],
            "A: Colocalized on centerline compared to BG volume [micron]" : [r"Colocalized volume [$\mu m^3$]" , "Channel A/B:"],
            "A: Colocalized on centerline compared to BG volume [% total length]" : ["Colocalized volume [%]", "Channel A/B:"],
            "B: Colocalized on centerline compared to BG volume [micron]" : [r"Colocalized volume [$\mu m^3$]", "Channel A/B:"],
            "B: Colocalized on centerline compared to BG volume [% total length]" : ["Colocalized volume [%]", "Channel A/B:"]
        }   

    def get_measurement_dic(self):
        dic_temp = {key: values[0] for key, values in self.measures_dic.items()}
        return dic_temp
    
    def get_measurement_category(self):
        dic_temp = {key: values[1] for key, values in self.measures_dic.items()}
        return dic_temp 
   
    def get_measurement_list(self):
        measurement_list = list(self.measures_dic.keys())
        measurement_list += ["Group", "Replicate"]
        return list(measurement_list)
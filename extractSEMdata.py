import os

file_list = os.listdir("/Users/nbense/NASA/BRAILLE/SEM_Updated_All/hdr_files/")
metadata_list = []
with open("/Users/nbense/NASA/BRAILLE/SEM_Updated_All/SEM_Metadata.csv", "w") as outputfile:
    outputfile.write("Filename"+","+
        "Date"+","+
        "Description"+","+
        "Device"+","+
        "FullUserName"+","+
        "ImageStripSize"+","+
        "Magnification"+","+
        "Note"+","+
        "PixelSizeX"+","+
        "PixelSizeY"+","+
        "ProjectName"+","+
        "SerialNumber"+","+
        "Sign"+","+
        "TagRevision"+","+
        "Time"+","+
        "UserName"+","+
        "Version"+","+
        "3DBeamTiltX"+","+
        "3DBeamTiltY"+","+
        "BeamIntensityIndex"+","+
        "ChamberPressure"+","+
        "Detector"+","+
        "DwellTime"+","+
        "EmissionCurrent"+","+
        "Gun"+","+
        "GunShiftX"+","+
        "GunShiftY"+","+
        "GunTiltX"+","+
        "GunTiltY"+","+
        "HV"+","+
        "IMLCenteringX"+","+
        "IMLCenteringY"+","+
        "ImageShiftX"+","+
        "ImageShiftY"+","+
        "LUTGamma"+","+
        "LUTMaximum"+","+
        "LUTMinimum"+","+
        "MixingMode"+","+
        "OBJCenteringX"+","+
        "OBJCenteringY"+","+
        "PrimaryDetectorGain"+","+
        "PrimaryDetectorOffset"+","+
        "ScanMode"+","+
        "ScanRotation"+","+
        "ScanSpeed"+","+
        "SpecimenCurrent"+","+
        "SpotSize"+","+
        "StageRotation"+","+
        "StageTilt"+","+
        "StageX"+","+
        "StageY"+","+
        "StageZ"+","+
        "StigmatorX"+","+
        "StigmatorY"+","+
        "SystemPressure"+","+
        "TiltCorrection"+","+
        "WD"+","+
        "AccFrames"+"\n")
    for file in file_list:
        with open("/Users/nbense/NASA/BRAILLE/SEM_Updated_All/hdr_files/"+ file, "r") as hdrfile:
            file_metadata = ["NA"] * 58
            file_metadata[0] = file.replace("-tif.hdr", ".tif")
            for line in hdrfile:
                if line.startswith(("Date=", '\t')):
                    if len(line.split('=')) > 1:
                        file_metadata[1] = line.split('=')[1]
                    else:
                        file_metadata[1] = "NA"
                elif line.startswith(("Description=", '\t')):
                    if len(line.split('=')) > 1:
                        file_metadata[2] = line.split('=')[1]
                    else:
                        file_metadata[2] = "NA"
                elif line.startswith(("Device=", '\t')):
                    if len(line.split('=')) > 1:
                        file_metadata[3] = line.split('=')[1]
                    else:
                        file_metadata[3] = "NA"
                elif line.startswith(("FullUserName=", '\t')):
                    if len(line.split('=')) > 1:
                        file_metadata[4] = line.split('=')[1]
                    else:
                        file_metadata[4] = "NA"
                elif line.startswith(("ImageStripSize=", '\t')):
                    if len(line.split('=')) > 1:
                        file_metadata[5] = float(line.split('=')[1])
                    else:
                        file_metadata[5] = "NA"
                elif line.startswith(("Magnification=", '\t')):
                    if len(line.split('=')) > 1:
                        file_metadata[6] = float(line.split('=')[1])
                    else:
                        file_metadata[6] = "NA"
                elif line.startswith(("Note=", '\t')):
                    if len(line.split('=')) > 1:
                        file_metadata[7] = line.split('=')[1]
                    else:
                        file_metadata[7] = "NA"
                elif line.startswith(("PixelSizeX=", '\t')):
                    if len(line.split('=')) > 1:
                        file_metadata[8] = float(line.split('=')[1])
                    else:
                        file_metadata[8] = "NA"
                elif line.startswith(("PixelSizeY=", '\t')):
                    if len(line.split('=')) > 1:
                        file_metadata[9] = float(line.split('=')[1])
                    else:
                        file_metadata[9] = "NA"
                elif line.startswith(("ProjectName=", '\t')):
                    if len(line.split('=')) > 1:
                        file_metadata[10] = line.split('=')[1]
                    else:
                        file_metadata[10] = "NA"
                elif line.startswith(("SerialNumber=", '\t')):
                    if len(line.split('=')) > 1:
                        file_metadata[11] = line.split('=')[1]
                    else:
                        file_metadata[11] = "NA"
                elif line.startswith(("Sign=", '\t')):
                    if len(line.split('=')) > 1:
                        file_metadata[12] = line.split('=')[1]
                    else:
                        file_metadata[12] = "NA"
                elif line.startswith(("TagRevision=", '\t')):
                    if len(line.split('=')) > 1:
                        file_metadata[13] = line.split('=')[1]
                    else:
                        file_metadata[13] = "NA"
                elif line.startswith(("Time=", '\t')):
                    if len(line.split('=')) > 1:
                        file_metadata[14] = line.split('=')[1]
                    else:
                        file_metadata[14] = "NA"
                elif line.startswith(("UserName=", '\t')):
                    if len(line.split('=')) > 1:
                        file_metadata[15] = line.split('=')[1]
                    else:
                        file_metadata[15] = "NA"
                elif line.startswith(("Version=", '\t')):
                    if len(line.split('=')) > 1:
                        file_metadata[16] = line.split('=')[1]
                    else:
                        file_metadata[16] = "NA"
                elif line.startswith(("DBeamTiltX=", '\t')):
                    if len(line.split('=')) > 1:
                        file_metadata[17] = float(line.split('=')[1])
                    else:
                        file_metadata[17] = "NA"
                elif line.startswith(("DBeamTiltY=", '\t')):
                    if len(line.split('=')) > 1:
                        file_metadata[18] = float(line.split('=')[1])
                    else:
                        file_metadata[18] = "NA"
                elif line.startswith(("BeamIntensityIndex=", '\t')):
                    if len(line.split('=')) > 1:
                        file_metadata[19] = float(line.split('=')[1])
                    else:
                        file_metadata[19] = "NA"
                elif line.startswith(("ChamberPressure=", '\t')):
                    if len(line.split('=')) > 1:
                        file_metadata[20] = float(line.split('=')[1])
                    else:
                        file_metadata[20] = "NA"
                elif line.startswith(("Detector=", '\t')):
                    if len(line.split('=')) > 1:
                        file_metadata[21] = line.split('=')[1]
                    else:
                        file_metadata[21] = "NA"
                elif line.startswith(("DwellTime=", '\t')):
                    if len(line.split('=')) > 1:
                        file_metadata[22] = float(line.split('=')[1])
                    else:
                        file_metadata[22] = "NA"
                elif line.startswith(("EmissionCurrent=", '\t')):
                    if len(line.split('=')) > 1:
                        file_metadata[23] = float(line.split('=')[1])
                    else:
                        file_metadata[23] = "NA"
                elif line.startswith(("Gun=", '\t')):
                    if len(line.split('=')) > 1:
                        file_metadata[24] = line.split('=')[1]
                    else:
                        file_metadata[24] = "NA"
                elif line.startswith(("GunShiftX=", '\t')):
                    if len(line.split('=')) > 1:
                        file_metadata[25] = float(line.split('=')[1])
                    else:
                        file_metadata[25] = "NA"
                elif line.startswith(("GunShiftY=", '\t')):
                    if len(line.split('=')) > 1:
                        file_metadata[26] = float(line.split('=')[1])
                    else:
                        file_metadata[26] = "NA"
                elif line.startswith(("GunTiltX", '\t')):
                    if len(line.split('=')) > 1:
                        file_metadata[27] = float(line.split('=')[1])
                    else:
                        file_metadata[27] = "NA"
                elif line.startswith(("GunTiltY", '\t')):
                    if len(line.split('=')) > 1:
                        file_metadata[28] = float(line.split('=')[1])
                    else:
                        file_metadata[28] = "NA"
                elif line.startswith(("HV=", '\t')):
                    if len(line.split('=')) > 1:
                        file_metadata[29] = float(line.split('=')[1])
                    else:
                        file_metadata[29] = "NA"
                elif line.startswith(("IMLCenteringX=", '\t')):
                    if len(line.split('=')) > 1:
                        file_metadata[30] = float(line.split('=')[1])
                    else:
                        file_metadata[30] = "NA"
                elif line.startswith(("IMLCenteringY=", '\t')):
                    if len(line.split('=')) > 1:
                        file_metadata[31] = float(line.split('=')[1])
                    else:
                        file_metadata[31] = "NA"
                elif line.startswith(("ImageShiftX=", '\t')):
                    if len(line.split('=')) > 1:
                        file_metadata[32] = float(line.split('=')[1])
                    else:
                        file_metadata[32] = "NA"
                elif line.startswith(("ImageShiftY=", '\t')):
                    if len(line.split('=')) > 1:
                        file_metadata[33] = float(line.split('=')[1])
                    else:
                        file_metadata[33] = "NA"
                elif line.startswith(("LUTGamma=", '\t')):
                    if len(line.split('=')) > 1:
                        file_metadata[34] = float(line.split('=')[1])
                    else:
                        file_metadata[34] = "NA"
                elif line.startswith(("LUTMaximum=", '\t')):
                    if len(line.split('=')) > 1:
                        file_metadata[35] = float(line.split('=')[1])
                    else:
                        file_metadata[35] = "NA"
                elif line.startswith(("LUTMinimum=", '\t')):
                    if len(line.split('=')) > 1:
                        file_metadata[36] = float(line.split('=')[1])
                    else:
                        file_metadata[36] = "NA"
                elif line.startswith(("MixingMode=", '\t')):
                    if len(line.split('=')) > 1:
                        file_metadata[37] = float(line.split('=')[1])
                    else:
                        file_metadata[37] = "NA"
                elif line.startswith(("OBJCenteringX=", '\t')):
                    if len(line.split('=')) > 1:
                        file_metadata[38] = float(line.split('=')[1])
                    else:
                        file_metadata[38] = "NA"
                elif line.startswith(("OBJCenteringY=", '\t')):
                    if len(line.split('=')) > 1:
                        file_metadata[39] = float(line.split('=')[1])
                    else:
                        file_metadata[39] = "NA"
                elif line.startswith(("PrimaryDetectorGain=", '\t')):
                    if len(line.split('=')) > 1:
                        file_metadata[40] = float(line.split('=')[1])
                    else:
                        file_metadata[40] = "NA"
                elif line.startswith(("PrimaryDetectorOffset=", '\t')):
                    if len(line.split('=')) > 1:
                        file_metadata[41] = float(line.split('=')[1])
                    else:
                        file_metadata[41] = "NA"
                elif line.startswith(("ScanMode=", '\t')):
                    if len(line.split('=')) > 1:
                        file_metadata[42] = line.split('=')[1]
                    else:
                        file_metadata[42] = "NA"
                elif line.startswith(("ScanRotation=", '\t')):
                    if len(line.split('=')) > 1:
                        file_metadata[43] = float(line.split('=')[1])
                    else:
                        file_metadata[43] = "NA"
                elif line.startswith(("ScanSpeed=", '\t')):
                    if len(line.split('=')) > 1:
                        file_metadata[44] = float(line.split('=')[1])
                    else:
                        file_metadata[44] = "NA"
                elif line.startswith(("SpecimenCurrent=", '\t')):
                    if len(line.split('=')) > 1:
                        file_metadata[45] = float(line.split('=')[1])
                    else:
                        file_metadata[45] = "NA"
                elif line.startswith(("SpotSize=", '\t')):
                    if len(line.split('=')) > 1:
                        file_metadata[46] = float(line.split('=')[1])
                    else:
                        file_metadata[46] = "NA"
                elif line.startswith(("StageRotation=", '\t')):
                    if len(line.split('=')) > 1:
                        file_metadata[47] = float(line.split('=')[1])
                    else:
                        file_metadata[47] = "NA"
                elif line.startswith(("StageTilt=", '\t')):
                    if len(line.split('=')) > 1:
                        file_metadata[48] = float(line.split('=')[1])
                    else:
                        file_metadata[48] = "NA"
                elif line.startswith(("StageX=", '\t')):
                    if len(line.split('=')) > 1:
                        file_metadata[49] = float(line.split('=')[1])
                    else:
                        file_metadata[49] = "NA"
                elif line.startswith(("StageY=", '\t')):
                    if len(line.split('=')) > 1:
                        file_metadata[50] = float(line.split('=')[1])
                    else:
                        file_metadata[50] = "NA"
                elif line.startswith(("StageZ=", '\t')):
                    if len(line.split('=')) > 1:
                        file_metadata[51] = float(line.split('=')[1])
                    else:
                        file_metadata[51] = "NA"
                elif line.startswith(("StigmatorX=", '\t')):
                    if len(line.split('=')) > 1:
                        file_metadata[52] = float(line.split('=')[1])
                    else:
                        file_metadata[52] = "NA"
                elif line.startswith(("StigmatorY=", '\t')):
                    if len(line.split('=')) > 1:
                        file_metadata[53] = float(line.split('=')[1])
                    else:
                        file_metadata[53] = "NA"
                elif line.startswith(("SystemPressure=", '\t')):
                    if len(line.split('=')) > 1:
                        file_metadata[54] = float(line.split('=')[1])
                    else:
                        file_metadata[54] = "NA"
                elif line.startswith(("TiltCorrection=", '\t')):
                    if len(line.split('=')) > 1:
                        file_metadata[55] = float(line.split('=')[1])
                    else:
                        file_metadata[55] = "NA"
                elif line.startswith(("WD=", '\t')):
                    if len(line.split('=')) > 1:
                        file_metadata[56] = float(line.split('=')[1])
                    else:
                        file_metadata[56] = "NA"
                elif line.startswith(("AccFrames=", '\t')):
                    if len(line.split('=')) > 1:
                        file_metadata[57] = float(line.split('=')[1])
                    else:
                        file_metadata[57] = "NA"
            file_metadata[0:5] = [x.strip() for x in file_metadata[0:5]]
            file_metadata[7] = file_metadata[7].strip()
            file_metadata[10:17] = [x.strip() for x in file_metadata[10:17]]
            file_metadata[21] = file_metadata[21].strip()
            file_metadata[24] = file_metadata[24].strip()
            file_metadata[42] = file_metadata[42].strip()
        outputfile.write(",".join(str(x) for x in file_metadata)+"\n")


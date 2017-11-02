def StrCount(str):
    """
    Help on function StrCount in module CalMorph:
        
    StrCount(str)
    Gives the number of streams corresponding to each Strahler order(Strahler, 1957).
    
    Parameters (Params)
    ----------
    str : Stream network file. It is created using the 'Stream to Feature' tool in the ArcGIS Spatial analyst -> Hydrology tool box. It is an Esri shape file. The extension of the file name is '.shp'. It has a number of fields like 'FID', 'Shape', 'ARCID', 'GRID_CODE' etc.
    
    A review of the str can be found at: http://desktop.arcgis.com/en/arcmap/10.3/tools/spatial-analyst-toolbox/an-overview-of-the-hydrology-tools.htm
    
    Returns
    --------
    out : a python directory object. The keys of the directory are the Strahler orders, whereas item corresponding to each key is their individual stream frequency.
    
    See also
    --------
    StrLength, Length_Maxordr, AreaEst, MorphRatio
    
    Notes
    -----
    1. This function is compatible with ArcGIS 10.1 or higher version.
    2. Name of the input parameters can be prefixed with the address of the directory or it can be stand-alone. In case of fully prefixed filename, it should be noted that in python '\' is a reserved character. So it is advised to use '\\' or '/' in the place of '\', while writing the input parameters.
    
    Examples
    ---------
    >>> from CalMorph import StrCount
    >>> str = 'filename' # The name of the file can be stand-alone or can be prefixed by the address of the directory
    >>> n = StrCount(str) # n is a dictionary
    >>> n1 = n[1] # n1 = stream frequency of order 1. '1' is one of the keys of the dictionary n
    """   
    
    # Import the necessary libraries
    import arcpy
    import numpy as np
    # -----------------------------------------------------------------------------
    
    # Checking if the input parameter type is consistent
    desc = arcpy.Describe(str)
    type_param = desc.dataType
    
    if type_param != 'FeatureLayer' and type_param != 'ShapeFile':
        print('Please check the type of your parameter param1')
        raise TypeError('param1 -> Expecting a FeatureLayer or Shapefile')
        # print('Please check the type of your parameter')
    
    number_of_features = int(arcpy.GetCount_management(str).getOutput(0)) # Get the feature count of the shape file
    
    fields = arcpy.ListFields(str) # List of fields in the attribute table of a shape file. For instance it can be (in order) ->  'FID', 'Shape', 'ARCID', 'GRID_CODE', 'FROM_NODE', 'TO_NODE' and other user defined fields
    
    number_of_fields = len(fields) - 3 # The function does not require the first 3 fields listed in the comments for the previous line
    
    name_of_fields = [] # Stores the name of the relevant fields
    for i in range(number_of_fields):
        name_of_fields.append(fields[i + 3].name)
    
    # creating a 2d numpy array 'ar1' containg the values of the relevant fields for all the features
    ar1 = np.zeros((number_of_features, number_of_fields))
    j = 0
    cursor = arcpy.SearchCursor(str)
    for row in cursor:
        for k in range(number_of_fields):
            ar1[j, k] = row.getValue(fields[k + 3].name)
        j = j + 1
        
    st_order_column = ar1[:, 0] # Accessing the first column of the array, which stores the GRID_CODE field (Strahler order of the features)
    st_order = np.unique(st_order_column)
    
    str_count = {} # Stores the total number of streams for each order
    
    for j1 in st_order:
        str_count[j1] = 0 # Set the dictionary key values
    
    #--------------------------------------------------------------------------------
    # Algorithm part for calculating the number of streams for each strahler order
    for j2 in st_order:
        ar2 = ar1[ar1[:, 0] == j2]
        if j2 == 1:
            str_count[j2] = ar2.shape[0]
        else:
            for j3 in range(ar2.shape[0]):
                indx = ar1[:, 1] == ar2[j3, 2]
                if len(ar1[indx]) == 0:
                    str_count[j2] += 1
                elif ar1[indx, 0] > j2:
                    str_count[j2] += 1
    
    #-------------------------------------------------------------------------------
    return(str_count)

def StrLength(str):
    """
    Help on function StrLength in module CalMorph:
        
    StrLength(str)
    Gives the average stream length corresponding to each Strahler order(Strahler, 1957).
    
    Parameters (Params)
    ----------
    str : Stream network file. It is created using the 'Stream to Feature' tool in the ArcGIS Spatial analyst -> Hydrology tool box. It is an Esri shape file. The extension of the file name is '.shp'. It has a number of fields like 'FID', 'Shape', 'ARCID', 'GRID_CODE' etc, along with 'length', which has to be calculated using the calculate geometry method inside ArcGIS. The shape file has to be in projected coordinate system (UTM), to facilitate the calculation of length of each feature in the shape file.
    
    A review of the str can be found at: http://desktop.arcgis.com/en/arcmap/10.3/tools/spatial-analyst-toolbox/an-overview-of-the-hydrology-tools.htm
    
    Returns
    --------
    out : a python directory object. The keys of the directory are the Strahler orders, whereas item corresponding to each key is their individual average stream length, which can be calculated for each Strahler order as:
        
        average stream length = total stream length / stream frequency
    
    See also
    --------
    StrCount, Length_Maxordr, AreaEst, MorphRatio
    
    Notes
    -----
    1. This function is compatible with ArcGIS 10.1 or higher version.
    2. Name of the input parameters can be prefixed with the address of the directory or it can be stand-alone. In case of fully prefixed filename, it should be noted that in python '\' is a reserved character. So it is advised to use '\\' or '/' in the place of '\', while writing the input parameters.
    
    Examples
    ---------
    >>> from CalMorph import StrCount
    >>> str = 'filename' # The name of the file can be stand-alone or can be prefixed by the address of the directory
    >>> l = StrLength(str) # l is a dictionary
    >>> l1 = l[1] # l1 = average stream length of Strahler order 1. '1' is one of the keys of the dictionary l
    """
    
    
    # Import the necessary libraries
    import arcpy
    import CalMorph
    import numpy as np
    # ----------------------------------------------------------------------------
    
    desc = arcpy.Describe(str)
    type_param = desc.dataType
    SR = desc.SpatialReference.type
    
    # Checking if the input parameter type is consistent
    if type_param != 'FeatureLayer' and type_param != 'ShapeFile':
        print('Please check the type of your parameter param1')
        raise TypeError('param1 -> Expecting a FeatureLayer')
    
    # Checking the projection system
    if SR != 'Projected':
        print('Please check the projection system of param1')
        raise TypeError('param1 was expected to be in projected coordinate system')
    # ----------------------------------------------------------------------------
    
    fields = arcpy.ListFields(str)
    number_of_features = int(arcpy.GetCount_management(str).getOutput(0))
    number_of_fields = len(fields) - 3
        
    ar1 = np.zeros((number_of_features, number_of_fields))
    j = 0
    cursor = arcpy.SearchCursor(str)
    for row in cursor:
        for k in range(number_of_fields):
            ar1[j, k] = row.getValue(fields[k + 3].name)
        j = j + 1
    
    # -----------------------------------------------------------------------------
    # Main algorithm part for average length calculation
    str_count = CalMorph.StrCount(str)
    str_length = {} # Stores the average length for each stream order
    for j2 in str_count.keys():
        ar2 = ar1[ar1[:, 0] == j2]
        total_length = sum(ar2[:, 3])
        str_length[j2] = float(total_length) / str_count[j2]
    
    # ----------------------------------------------------------------------------
    return str_length

def Length_Maxordr(str):
    """
    Help on function Length_Maxordr in module CalMorph:
        
    Length_Maxordr(str)
    Gives the length of the stream of maximum  Strahler order(Strahler, 1957).
    
    Parameters (Params)
    ----------
    str : Stream network file. It is created using the 'Stream to Feature' tool in the ArcGIS Spatial analyst -> Hydrology tool box. It is an Esri shape file. The extension of the file name is '.shp'. It has a number of fields like 'FID', 'Shape', 'ARCID', 'GRID_CODE' etc, along with 'length', which has to be calculated using the calculate geometry method inside ArcGIS. The shape file has to be in projected coordinate system (UTM), to facilitate the calculation of length of each feature in the shape file.
    
    A review of the str can be found at: http://desktop.arcgis.com/en/arcmap/10.3/tools/spatial-analyst-toolbox/an-overview-of-the-hydrology-tools.htm
    
    Returns
    --------
    out : a numpy.float64 object.
        
    
    See also
    --------
    StrCount, StrLength, AreaEst, MorphRatio
    
    Notes
    -----
    1. This function is compatible with ArcGIS 10.1 or higher version.
    2. Name of the input parameters can be prefixed with the address of the directory or it can be stand-alone. In case of fully prefixed filename, it should be noted that in python '\' is a reserved character. So it is advised to use '\\' or '/' in the place of '\', while writing the input parameters.
    
    Examples
    ---------
    >>> from CalMorph import Length_Maxordr
    >>> str = 'filename' # The name of the file can be stand-alone or can be prefixed by the address of the directory
    >>> lmax = Length_Maxordr(str)
    """
    
    # Importing modules
    import arcpy
    import numpy as np
    
    # ---------------------------------------------------------------------------
    desc = arcpy.Describe(str)
    type_param = desc.dataType
    SR = desc.SpatialReference.type
    
    # Checking if the input parameter type is consistent
    if type_param != 'FeatureLayer' and type_param != 'ShapeFile':
        print('Please check the type of your parameter param1')
        raise TypeError('param1 -> Expecting a FeatureLayer')
    
    # Checking the projection system
    if SR != 'Projected':
        print('Please check the projection system of param1')
        raise TypeError('param1 was expected to be in projected coordinate system')
    # ----------------------------------------------------------------------------
    
    
    fields = arcpy.ListFields(str)
    number_of_features = int(arcpy.GetCount_management(str).getOutput(0))
    number_of_fields = len(fields) - 3
                          
    ar1 = np.zeros((number_of_features, number_of_fields))
    j = 0
    cursor = arcpy.SearchCursor(str)
    for row in cursor:
        for k in range(number_of_fields):
            ar1[j, k] = row.getValue(fields[k + 3].name)
        j = j + 1
    
    st_order = np.unique(ar1[:, 0])
    max_ord = max(st_order)
    ar2 = ar1[ar1[:, 0] == max_ord]
    Lw = sum(ar2[:, 3])
    
    # -----------------------------------------------------------------------------
    return Lw

def AreaEst(FlowDir, FlowAcc, StrOrdr, str):
    """
    Help on function AreaEst in module CalMorph:
        
    AreaEst(FlowDir, FlowAcc, StrOrdr, str)
    Gives the average discharge area corresponding to each Strahler order(Strahler, 1957).
    
    Parameters(Params)
    ----------
    FlowDir : Flow direction raster file. It is generated using the 'Flow direction' tool in the ArcGIS Spatial analyst -> Hydrology tool box.
    
    FlowAcc : Flow Accumulation raster file. It is generated using the 'Flow accumulation' tool in the ArcGIS Spatial analyst -> Hydrology tool box.
    
    StrOrdr : Stream order raster file. It is generated using the 'Stream order' tool in the ArcGIS Spatial analyst -> Hydrology tool box.
    
    str : Stream network file. It is created using the 'Stream to Feature' tool in the ArcGIS Spatial analyst -> Hydrology tool box. It is an Esri shape file. The extension of the file name is '.shp'. It has a number of fields like 'FID', 'Shape', 'ARCID', 'GRID_CODE' etc.
    
    A review of the paremeters can be found at: http://desktop.arcgis.com/en/arcmap/10.3/tools/spatial-analyst-toolbox/an-overview-of-the-hydrology-tools.htm
    
    Returns
    --------
    out : a python directory object. The keys of the directory are the Strahler orders, whereas item corresponding to each key is their individual average discharge area, which can be calculated for each Strahler order as:
        
        average discharge area = total discharge area / stream frequency
    
    See also
    --------
    StrCount, StrLength, Length_Maxordr, MorphRatio
    
    Notes
    -----
    1. This function is compatible with ArcGIS 10.1 or higher version.
    2. Name of the input parameters can be prefixed with the address of the directory or it can be stand-alone. In case of fully prefixed filename, it should be noted that in python '\' is a reserved character. So it is advised to use '\\' or '/' in the place of '\', while writing the input parameters.
    
    Examples
    ---------
    >>> from CalMorph import AreaEst
    >>> A = AreaEst(FlowDir, FlowAcc, StrOrdr, str) # A is a dictionary
    >>> A1 = A[1] # A1 = average discharge area of Strahler order 1. '1' is one of the keys of the dictionary A
    """
    
    # Importing the modules
    import arcpy
    import numpy as np
    import CalMorph
    
    # ---------------------------------------------------------------------------
    desc_param1 = arcpy.Describe(FlowDir)
    desc_param2 = arcpy.Describe(FlowAcc)
    desc_param3 = arcpy.Describe(StrOrdr)
    desc_param4 = arcpy.Describe(str)
    
    SR_param4 = desc_param4.SpatialReference.type
    
    type_param1 = desc_param1.dataType
    if type_param1 != 'RasterLayer' and type_param1 != 'RasterDataset':
        print('Please check the type of your parameter param1')
        raise TypeError('param1 -> Expecting a Raster Layer')
    
    type_param2 = desc_param2.dataType
    if type_param2 != 'RasterLayer' and type_param2 != 'RasterDataset':
        print('Please check the type of your parameter param2')
        raise TypeError('param2 -> Expecting a Raster Layer')
        
    type_param3 = desc_param3.dataType
    if type_param3 != 'RasterLayer' and type_param3 != 'RasterDataset':
        print('Please check the type of your parameter param3')
        raise TypeError('param3 -> Expecting a Raster Layer')
        
    type_param4 = desc_param4.dataType
    if type_param4 != 'FeatureLayer' and type_param4 != 'ShapeFile':
        print('Please check the type of your parameter param4')
        raise TypeError('param4 -> Expecting a Feature Layer')
        
    if SR_param4 != 'Projected':
        print('Please check the projection system of param1')
        raise TypeError('param4 was expected to be in projected coordinate system')

    # ------------------------------------------------------------------------------
    
    str_count = CalMorph.StrCount(str)
    
    data_FlowDir = arcpy.RasterToNumPyArray(FlowDir)
    data_FlowDir[(data_FlowDir < 1) | (data_FlowDir > 128)] = 0
    
    data_FlowAcc = arcpy.RasterToNumPyArray(FlowAcc)
    
    data_StrOrder = arcpy.RasterToNumPyArray(StrOrdr)
    data_StrOrder[(data_StrOrder < 1) | (data_StrOrder > 50)] = 0
    
    m = data_StrOrder.shape[0]
    n = data_StrOrder.shape[1]
    
    st = data_StrOrder.flatten()
    fdir = data_FlowDir.flatten()
    facc = data_FlowAcc.flatten()
    
    max_order = st.max()
    dic_area = {i:0 for i in range(1,max_order + 1)}
    max_area = facc.max()
    dic_area[max(dic_area)] = max_area
    
    fdir_new = np.zeros(len(fdir))
    dic_fdir = {1: 1, 2: n+1, 4: n, 8: n-1, 16: -1, 32: -n-1, 64: -n, 128: -n+1, 0: 0}
    for i in range(len(fdir)):
        fdir_new[i] = i + dic_fdir[fdir[i]]
        
    # --------------------------------------------------------------------------
    # Main algorithm for calculation of average area
    for k in range(n, n*(m-1)):
        if (k % n != 0) & ((k + 1) % n != 0) & (st[k] > 0):
            list1 = [k+1,k-1,k+n,k-n,k+n-1,k+n+1,k-n+1,k-n-1]
            for j in list1:
                if (int(fdir_new[k]) == j) & (st[j] > st[k]):
                    dic_area[st[k]] = dic_area[st[k]] + facc[k]
    avg_area = {}
    for key in dic_area.keys():
        avg_area[key] = dic_area[key] / str_count[key]
    
    # --------------------------------------------------------------------------
    return avg_area

def MorphRatio(FlowDir, FlowAcc, StrOrdr, str):
    """
    Help on function MorphRatio in module CalMorph:
        
    MorphRatio(FlowDir, FlowAcc, StrOrdr, str)
    Gives the Horton's ratios (Bifurcation ratio Rb, Length ratio Rl, Area ratio Ra) for the drainage basin(Horton, 1945; Schumm, 1956; Strahler, 1957).
    
    Parameters
    ----------
    FlowDir : Flow direction raster file. It is generated using the 'Flow direction' tool in the ArcGIS Spatial analyst -> Hydrology tool box.
    
    FlowAcc : Flow Accumulation raster file. It is generated using the 'Flow accumulation' tool in the ArcGIS Spatial analyst -> Hydrology tool box.
    
    StrOrdr : Stream order raster file. It is generated using the 'Stream order' tool in the ArcGIS Spatial analyst -> Hydrology tool box.
    
    str : Stream network file. It is created using the 'Stream to Feature' tool in the ArcGIS Spatial analyst -> Hydrology tool box. It is an Esri shape file. The extension of the file name is '.shp'. It has a number of fields like 'FID', 'Shape', 'ARCID', 'GRID_CODE' etc, along with 'length', which has to be calculated using the calculate geometry method inside ArcGIS. The shape file has to be in projected coordinate system (UTM), to facilitate the calculation of length of each feature in the shape file.
    
    A review of the paremeters can be found at: http://desktop.arcgis.com/en/arcmap/10.3/tools/spatial-analyst-toolbox/an-overview-of-the-hydrology-tools.htm
    
    Returns
    --------
    out : a python directory object. The keys of the directory are 'Rb', 'Rl', 'Ra' and their corresponding items are the ratios.
    
    See also
    --------
    StrCount, StrLength, Length_Maxordr, AreaEst
    
    Notes
    -----
    1. This function is compatible with ArcGIS 10.1 or higher version.
    2. Name of the input parameters can be prefixed with the address of the directory or it can be stand-alone. In case of fully prefixed filename, it should be noted that in python '\' is a reserved character. So it is advised to use '\\' or '/' in the place of '\', while writing the input parameters.
    
    Examples
    ---------
    >>> from CalMorph import MorphRatio
    >>> R = MorphRatio(FlowDir, FlowAcc, StrOrdr, str) # R is a dictionary
    >>> Rb = R['Rb'] # Rb = Bifurcation ratio of the drainage basin. 'Rb' is one of the key of the dirctory R. Similarly Rl and Ra can be retrieved from the dictionary
    """
    # Importin the modules
    import arcpy
    import CalMorph
    import numpy as np
    
    # -----------------------------------------------------------------------------
    desc_param1 = arcpy.Describe(FlowDir)
    desc_param2 = arcpy.Describe(FlowAcc)
    desc_param3 = arcpy.Describe(StrOrdr)
    desc_param4 = arcpy.Describe(str)
    
    SR_param4 = desc_param4.SpatialReference.type
    
    type_param1 = desc_param1.dataType
    if type_param1 != 'RasterLayer' and type_param1 != 'RasterDataset':
        print('Please check the type of your parameter param1')
        raise TypeError('param1 -> Expecting a Raster Layer')
    
    type_param2 = desc_param2.dataType
    if type_param2 != 'RasterLayer' and type_param2 != 'RasterDataset':
        print('Please check the type of your parameter param2')
        raise TypeError('param2 -> Expecting a Raster Layer')
        
    type_param3 = desc_param3.dataType
    if type_param3 != 'RasterLayer' and type_param3 != 'RasterDataset':
        print('Please check the type of your parameter param3')
        raise TypeError('param3 -> Expecting a Raster Layer')
        
    type_param4 = desc_param4.dataType
    if type_param4 != 'FeatureLayer' and type_param4 != 'ShapeFile':
        print('Please check the type of your parameter param4')
        raise TypeError('param4 -> Expecting a Feature Layer')
        
    if SR_param4 != 'Projected':
        print('Please check the projection system of param1')
        raise TypeError('param4 was expected to be in projected coordinate system')
        
    # ---------------------------------------------------------------------------
    
    # Main algorithm part to calculate the ratios
    number_of_str = CalMorph.StrCount(str)
    avg_str_length = CalMorph.StrLength(str)
    avg_area = CalMorph.AreaEst(FlowDir, FlowAcc, StrOrdr, str)
    
    ar_orders = np.zeros(len(avg_area.keys()))
    ar_number_of_str = np.zeros(len(ar_orders))
    ar_avgLength = np.zeros(len(ar_orders))
    ar_avgArea = np.zeros(len(ar_orders))
    for i in range(len(avg_area.keys())):
        ar_orders[i] = avg_area.keys()[i]
        ar_number_of_str[i] = number_of_str[ar_orders[i]]
        ar_avgLength[i] = avg_str_length[ar_orders[i]]
        ar_avgArea[i] = avg_area[ar_orders[i]]
        
    A = np.vstack([ar_orders, np.ones(len(ar_orders))]).T
    m1, c1 = np.linalg.lstsq(A, np.log10(ar_number_of_str))[0]
    m2, c2 = np.linalg.lstsq(A, np.log10(ar_avgLength))[0]
    m3, c3 = np.linalg.lstsq(A, np.log10(ar_avgArea))[0]
    
    dic_morphRatio = {}
    dic_morphRatio['Rb'] = 10 ** abs(m1)
    dic_morphRatio['Rl'] = 10 ** m2
    dic_morphRatio['Ra'] = 10 ** m3
    
    # -----------------------------------------------------------------------------
    return dic_morphRatio

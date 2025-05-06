def from_170k_to_59k(img, data, return_concat_hemis=False, return_59k_mask=False):
    """
    Transform 170k data into 59k data by retaining only the left and right cortex and medial wall vertices.

    Parameters
    ----------
    img : the cifti image of your 170k 
    data : your 170k data
    return_concat_hemis : if False, will return two arrays for the two hemispheres 
                          (2 arrays with dimensions (time x 59k vertices))
                          
                          if True, will return the concatenation of hemi-L and hemi-R 
                          (1 array with dimensions (time x 118 vertices))
                          
    return_59k_mask : if True, will return a mask where True corresponds to cortex vertices and
                      False to medial wall vertices.
    

    Returns
    -------
    result : dict
        A dictionary containing the following keys:
        - 'concatenated_data': numpy array, stacked data for the requested hemisphere(s). 2-dimensional array (time x vertices).
        - 'mask_59k': optional numpy array, mask where True corresponds to cortex vertices and False to medial wall vertices for 59k data.
    """

    import numpy as np 

    brain_model = img.header.get_axis(1)
    surf_name_L = 'CIFTI_STRUCTURE_CORTEX_LEFT'
    surf_name_R = 'CIFTI_STRUCTURE_CORTEX_RIGHT'

    result = {}

    for structure_name, data_indices, model in brain_model.iter_structures(): 
        if structure_name == surf_name_L: 
            data_L = data.T[data_indices]
            vtx_indices_L = model.vertex 

            # include inter hemi vertex
            surf_data_L = np.zeros((vtx_indices_L.max() + 1,) + data_L.shape[1:], dtype=data_L.dtype)
            surf_data_L[vtx_indices_L] = data_L    
            surf_data_L = surf_data_L.T 

            # Know where are inter hemi vertex
            mask_L = np.any(surf_data_L != 0, axis=0)

        elif structure_name == surf_name_R: 
            data_R = data.T[data_indices]
            vtx_indices_R = model.vertex 

            # include inter hemi vertex
            surf_data_R = np.zeros((vtx_indices_R.max() + 1,) + data_R.shape[1:], dtype=data_R.dtype)
            surf_data_R[vtx_indices_R] = data_R
            surf_data_R = surf_data_R.T 

            # Know where are inter hemi vertex
            mask_R = np.any(surf_data_R != 0, axis=0)

    brain_mask_59k = np.concatenate((mask_L, mask_R))
    
    if return_concat_hemis:
        result['data_concat'] = np.concatenate((surf_data_L, surf_data_R), axis=1)
    else:
        result['data_L'] = surf_data_L
        result['data_R'] = surf_data_R

    if return_59k_mask:
        result['mask_59k'] = brain_mask_59k


    return result


def from_59k_to_170k(data_59k, brain_mask_59k):
    """
    Transform 59k data into 170k data by filling non-59k vertices with numpy.nan.

    Parameters
    ----------
    
    data_59k : The 59k data you want to transform into 170k.
    brain_mask_59k : 59k brain mask output from from_170k_to_59k.
    
    Returns
    -------
    The transformed data in 170k format with non-59k vertices filled with numpy.nan.
    """
    import numpy as np
    n_vertex_170k = 170494
    # mask 59k data to optain only cortex vertex (and not medial wall vertices ) 
    data_54k = data_59k[:,brain_mask_59k]

    # create an 170k full nan array
    nan_object = np.full((data_54k.shape[0], n_vertex_170k - data_54k.shape[1]), np.nan)
    data_170k_final = np.concatenate((data_54k, nan_object), axis=1)
    
    return data_170k_final


def from_91k_to_32k(img, data, return_concat_hemis=False, return_32k_mask=False):
    """
    Transform 91k data into 32k data by retaining only the left and right cortex and medial wall vertices.

    Parameters
    ----------
    img : the cifti image of your 91k 
    data : your 91k data
    return_concat_hemis : if False, will return two arrays for the two hemispheres 
                          (2 arrays with dimensions (time x 32k vertices))
                          
                          if True, will return the concatenation of hemi-L and hemi-R 
                          (1 array with dimensions (time x 65k vertices))
                          
    return_32k_mask : if True, will return a mask where True corresponds to cortex vertices and
                      False to medial wall vertices.
    

    Returns
    -------
    result : dict
        A dictionary containing the following keys:
        - 'concatenated_data': numpy array, stacked data for the requested hemisphere(s). 2-dimensional array (time x vertices).
        - 'mask_32k': optional numpy array, mask where True corresponds to cortex vertices and False to medial wall vertices for 32k data.
    """

    import numpy as np 

    brain_model = img.header.get_axis(1)
    surf_name_L = 'CIFTI_STRUCTURE_CORTEX_LEFT'
    surf_name_R = 'CIFTI_STRUCTURE_CORTEX_RIGHT'

    result = {}

    for structure_name, data_indices, model in brain_model.iter_structures(): 
        if structure_name == surf_name_L: 
            data_L = data.T[data_indices]
            vtx_indices_L = model.vertex 

            # include inter hemi vertex
            surf_data_L = np.zeros((vtx_indices_L.max() + 1,) + data_L.shape[1:], dtype=data_L.dtype)
            surf_data_L[vtx_indices_L] = data_L    
            surf_data_L = surf_data_L.T 

            # Know where are inter hemi vertex
            mask_L = np.any(surf_data_L != 0, axis=0)

        elif structure_name == surf_name_R: 
            data_R = data.T[data_indices]
            vtx_indices_R = model.vertex 

            # include inter hemi vertex
            surf_data_R = np.zeros((vtx_indices_R.max() + 1,) + data_R.shape[1:], dtype=data_R.dtype)
            surf_data_R[vtx_indices_R] = data_R
            surf_data_R = surf_data_R.T 

            # Know where are inter hemi vertex
            mask_R = np.any(surf_data_R != 0, axis=0)

    brain_mask_32k = np.concatenate((mask_L, mask_R))
    
    if return_concat_hemis:
        result['data_concat'] = np.concatenate((surf_data_L, surf_data_R), axis=1)
    else:
        result['data_L'] = surf_data_L
        result['data_R'] = surf_data_R

    if return_32k_mask:
        result['mask_32k'] = brain_mask_32k

    return result


def from_32k_to_91k(data_32k, brain_mask_32k):
    """
    Transform 32k data into 91k data by filling non-32k vertices with numpy.nan.

    Parameters
    ----------
    
    data_32k : The 32k data you want to transform into 91k.
    brain_mask_32k : 32k brain mask output from from_91k_to_32k.
    
    Returns
    -------
    The transformed data in 91k format with non-32k vertices filled with numpy.nan.
    """
    import numpy as np
    n_vertex_91k = 91282
    # mask 59k data to optain only cortex vertex (and not medial wall vertices ) 
    data_30k = data_32k[:,brain_mask_32k]

    # create an 170k full nan array
    nan_object = np.full((data_30k.shape[0], n_vertex_91k - data_30k.shape[1]), np.nan)
    data_91k_final = np.concatenate((data_30k, nan_object), axis=1)
    
    return data_91k_final
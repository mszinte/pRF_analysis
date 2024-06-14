def load_gifti_image(gii_fn): 
    """
    load a gifti image

    Parameters
    ----------
    gii_fn : gifti single hemisphere filename
    
    Returns
    -------
    data_hemi : numpy stacked data
                2 dim (time x vertices)    
    """
    # import 
    import nibabel as nb
    import numpy as np
    img_hemi = nb.load(gii_fn)
    data_hemi = [x.data for x in img_hemi.darrays]
    data_hemi = np.vstack(data_hemi)

    return img_hemi, data_hemi

def make_gifti_image(source_img, data_to_write) : 
    """
    make a gifti image with data

    Parameters
    ----------
    source_img : image from with new data derives
    data_to_write : data you want to write on your image 
                    numpy darray 2 dim (time x vertices)
    
    Returns
    -------
    The gifti image of your data with the same strucure the source image 
    
    """
    # import 
    import nibabel as nb
    
    # gete source img informations
    header = source_img.header
    meta = source_img.meta
    file_map = source_img.file_map
    labeltable = source_img.labeltable
    
    # initialise the final image 
    final_img = nb.gifti.GiftiImage(header=header, meta=meta,
                                    file_map=file_map,
                                    labeltable=labeltable)
    
    # fill final image
    for i in range(data_to_write.shape[0]):
        data = data_to_write[i,:]
        darray = nb.gifti.GiftiDataArray(data, datatype='NIFTI_TYPE_FLOAT32',
                                         intent=source_img.darrays[i].intent, 
                                         meta=source_img.darrays[i].meta, 
                                         coordsys=source_img.darrays[i].coordsys)
        final_img.add_gifti_data_array(darray)
    
    return final_img
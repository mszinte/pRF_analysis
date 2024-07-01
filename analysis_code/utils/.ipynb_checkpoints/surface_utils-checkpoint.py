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

def make_gifti_image(data, source_img) : 
    """
    make a gifti image with data

    Parameters
    ----------
    data_to_write : data you want to write on your image 
                    numpy darray 2 dim (time x vertices)
    source_img : image from with new data derives

    
    Returns
    -------
    The gifti image of your data with the same strucure of the source image 
    
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
    for i in range(data.shape[0]):
        time_point = data[i,:]
        darray = nb.gifti.GiftiDataArray(time_point, datatype='NIFTI_TYPE_FLOAT32',
                                         intent=source_img.darrays[i].intent, 
                                         meta=source_img.darrays[i].meta, 
                                         coordsys=source_img.darrays[i].coordsys)
        final_img.add_gifti_data_array(darray)
    
    return final_img

def make_cifti_image(data, source_img):
    """
    make cifti image with data using a template CIFTI file.
    
    Parameters:
    -----------
    - data: numpy.ndarray
        The data array to be saved. In shape (trs, nr_voxels)
    - source_img: str
        The path to the template CIFTI file for header information.
    
    Returns:
    --------
    - Cifti2Image Object
    
    Example usage:
    --------------
    cor_final_new = np.random.random((170494, 1))
    template_file = 'sub-02_ses-01_task-pRF_run-01_space-fsLR_den-170k_bold_dct.dtseries.nii'
    output_file = 'test_correlation.dtseries.nii'
    save_cifti(cor_final_new, template_file, output_file)
    """
    import nibabel as nb
    # Load the template CIFTI file to get header information

    # Create a new header with the necessary axes
    ax_0 = nb.cifti2.SeriesAxis(start=0, step=1, size=data.shape[0])
    ax_1 = source_img.header.get_axis(1)
    new_header = nb.cifti2.Cifti2Header.from_axes((ax_0, ax_1))
    
    extra = source_img.extra
    file_map = source_img.file_map
    dtype = source_img.get_data_dtype()

    # Create a CIFTI image with the new data and header
    img = nb.Cifti2Image(dataobj=data, 
                         header=new_header,
                         extra=extra, 
                         file_map=file_map,
                         dtype=dtype)
    
    return img

def load_surface(fn): 
    """
    load a surface image inndependently if it's CIFTI or GIFTI

    Parameters
    ----------
    fn : surface filename

    Returns
    -------
    img : your image 
    data : a np.array of the data from your imahe. 2 dim (time x vertices)   
    """
    
    # import 
    import nibabel as nb
    
    if fn.endswith('.gii'):
        img, data = load_gifti_image(fn)

    elif fn.endswith('.nii'):
        img = nb.load(fn)
        data = img.get_fdata()
            
    else:
         raise ValueError("The type of fn is neither Cifti2Image nor GiftiImage")

    return img, data

def make_surface_image(data, source_img):
    """
    write a surface image inndependently if it's CIFTI or GIFTI

    Parameters
    ----------
    data_to_write : data you want to write on your image 
                    numpy darray 2 dim (time x vertices)
    source_img : image from with new data derives

    
    Returns
    -------
    The surface image of your data with the same strucure of the source image 
    
    """
    
    import nibabel as nb

    if type(source_img) == nb.cifti2.cifti2.Cifti2Image:
        img = make_cifti_image(data=data, source_img=source_img)
        

    elif type(source_img) == nb.gifti.gifti.GiftiImage:
        img = make_gifti_image(data=data, source_img=source_img)
        
    else:
         raise ValueError("The type of source_img is neither Cifti2Image nor GiftiImage")
         
    return img
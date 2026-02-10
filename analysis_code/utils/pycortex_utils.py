def get_rois(subject, surf_format, atlas_name, mask=True, rois=None, hemis=None):
    """
    Load ROI masks stored as .npz files from a pycortex subject database.

    The ROI files are expected to be dictionaries where:
        - keys are ROI names (str)
        - values are boolean numpy arrays
          (True = vertex belongs to the ROI, False otherwise)

    Depending on the arguments, the function can return:
        - full-brain ROI masks
        - left or right hemisphere ROI masks
        - both hemispheres simultaneously

    Parameters
    ----------
    subject : str
        Pycortex subject name.
        Must be:
            - 'sub-hcp59k' for surf_format in ['170k', '59k']
            - 'sub-hcp32k' for surf_format in ['91k', '32k']

    surf_format : str
        Surface format of the ROI masks.
        Allowed values are:
            ['170k', '59k', '91k', '32k', 'fsnative']

    atlas_name : str
        Name of the atlas used to define the ROIs.
        Allowed values are:
            ['mmp', 'group-mmp', 'drawn']

    mask : bool, optional (default=True)
        If True, return boolean masks.
        If False, return vertex indices (np.where(mask)[0]).
        Indices are relative to the requested space
        (full brain or hemisphere).

    rois : list of str or None, optional (default=None)
        List of ROI names to load.
        If None, all ROIs found in the file are returned.

    hemis : None, str, or list of str, optional (default=None)
        Hemisphere selection:
            - None : load full-brain ROIs
            - 'hemi-L' : load left hemisphere ROIs
            - 'hemi-R' : load right hemisphere ROIs
            - ['hemi-L', 'hemi-R'] : load both hemispheres

    Returns
    -------
    dict
        ROI masks or vertex indices, depending on input arguments.

    Notes
    -----
    Expected file naming conventions:

    - Full brain:
        {subject}_{surf_format}_rois-{atlas_name}.npz

    - Single hemisphere:
        {subject}_{surf_format}_hemi-L_rois-{atlas_name}.npz
        {subject}_{surf_format}_hemi-R_rois-{atlas_name}.npz

    All files must be located in:
        <pycortex_db>/<subject>/rois/
    """
    import cortex
    import numpy as np
    import os

    # --------------------------------------------------
    # Input checks
    # --------------------------------------------------
    if surf_format not in ['170k', '59k', '91k', '32k', 'fsnative']:
        raise ValueError("Invalid value for 'surf_format'.")
    if atlas_name not in ['mmp', 'group-mmp', 'drawn']:
        raise ValueError("Invalid value for 'atlas_name'.")
    if surf_format in ['170k', '59k'] and subject != 'sub-hcp59k':
        raise ValueError("For this format subject should be sub-hcp59k")
    if surf_format in ['91k', '32k'] and subject != 'sub-hcp32k':
        raise ValueError("For this format subject should be sub-hcp32k")

    # --------------------------------------------------
    # Paths
    # --------------------------------------------------
    db_dir = cortex.database.default_filestore
    rois_dir = os.path.join(db_dir, subject, 'rois')

    # --------------------------------------------------
    # Full brain
    # --------------------------------------------------
    if hemis is None:
        rois_dict_brain = dict(
            np.load(
                '{}/{}_{}_rois-{}.npz'.format(
                    rois_dir, subject, surf_format, atlas_name
                )
            )
        )

        if rois is not None:
            rois_dict_brain = {
                roi_name: roi_mask
                for roi_name, roi_mask in rois_dict_brain.items()
                if roi_name in rois
            }

        if not mask:
            rois_dict_brain = {
                roi_name: np.where(roi_mask)[0]
                for roi_name, roi_mask in rois_dict_brain.items()
            }

        return rois_dict_brain

    # --------------------------------------------------
    # Hemispheres
    # --------------------------------------------------
    if isinstance(hemis, str):
        hemis = [hemis]

    if not isinstance(hemis, (list, tuple)):
        raise ValueError(
            "Invalid value for 'hemis'. It should be 'hemi-L', 'hemi-R' or ['hemi-L', 'hemi-R']."
        )

    for hemi in hemis:
        if hemi not in ['hemi-L', 'hemi-R']:
            raise ValueError(
                "Invalid value for 'hemis'. It should be 'hemi-L', 'hemi-R' or ['hemi-L', 'hemi-R']."
            )

    rois_dict_by_hemi = {}

    for hemi in hemis:
        rois_dict_hemi = dict(
            np.load(
                '{}/{}_{}_{}_rois-{}.npz'.format(
                    rois_dir, subject, surf_format, hemi, atlas_name
                )
            )
        )

        if rois is not None:
            rois_dict_hemi = {
                roi_name: roi_mask
                for roi_name, roi_mask in rois_dict_hemi.items()
                if roi_name in rois
            }

        if not mask:
            rois_dict_hemi = {
                roi_name: np.where(roi_mask)[0]
                for roi_name, roi_mask in rois_dict_hemi.items()
            }

        rois_dict_by_hemi[hemi] = rois_dict_hemi

    if len(rois_dict_by_hemi) == 1:
        return rois_dict_by_hemi[hemis[0]]

    return rois_dict_by_hemi         



# def data_from_rois(fn, subject, rois):
#     """
#     Load a surface, and returne vertex only data from the specified ROIs
#     ----------
#     fn : surface filename
#     subject : subject 
#     rois : list of rois you want extract
#     filter_rois : bool, optional
#         If True (default), filters out NaN vertices from the ROI data.
#         If False, returns all vertices in the ROIs without NaN filtering.
    
#     Returns
#     -------
#     img : the image load from fn   
#     data_roi : numpy rois data 
#               2 dim (time x vertices from all the rois)  
              
#     roi_idx : indices of the rois vertices 
    
    
#     data_hemi : numpy stacked data
#                 2 dim (time x vertices)    
#     """
#     import cortex
#     import numpy as np
#     from surface_utils import load_surface

#     # Import data
#     img, data = load_surface(fn=fn)
#     len_data = data.shape[1]
    
#     # Get regions of interest (ROIs) mask
#     if fn.endswith('.gii'):
#         roi_verts = cortex.get_roi_verts(subject=subject, roi=rois, mask=True)
#     elif fn.endswith('.nii'):
#         surf_size = '170k' if len_data > 60000 else '59k'
#         roi_verts = load_rois_atlas(atlas_name='mmp', 
#                                     surf_size=surf_size, 
#                                     return_hemis=False,
#                                     rois=rois, 
#                                     mask=True)

#     # Create a brain mask
#     # na_vertices = np.where(np.isnan(data).any(axis=0))[0]
#     brain_mask = np.any(list(roi_verts.values()), axis=0)
    
#     # Filter out NaN vertices
#     na_vertices = np.isnan(data).any(axis=0)
    
#     # create a hemi mask  
#     if 'hemi-L' in fn:
#         hemi_mask = brain_mask[:len_data]
#         for i, na_vertex in enumerate(na_vertices):
#             hemi_mask[i] = not na_vertex and hemi_mask[i]
        
#     elif 'hemi-R' in fn: 
#         hemi_mask = brain_mask[-len_data:]
#         for i, na_vertex in enumerate(na_vertices):
#             hemi_mask[i] = not na_vertex and hemi_mask[i]
#     else: 
#         hemi_mask = brain_mask
#         for i, na_vertex in enumerate(na_vertices):
#             hemi_mask[i] = not na_vertex and hemi_mask[i]
    
#     # Get indices of regions of interest (ROIs)
#     roi_idx = np.where(hemi_mask)[0]
    
#     # Extract data corresponding to regions of interest (ROIs)
#     data_roi = data[:, hemi_mask]

        
#     return img, data, data_roi, roi_idx


def load_surface_pycortex(L_fn=None, R_fn=None, brain_fn=None, return_img=False, 
                          return_hemi_len=False, return_59k_mask=False, return_source_data=False):
    """
    Load a surface image independently if it's CIFTI or GIFTI, and return 
    concatenated data from the left and right cortex if data are GIFTI and 
    a decomposition from 170k vertices to 59k verticices if data are CIFTI.

    Parameters
    ----------
    L_fn : gifti left hemisphere filename
    R_fn : gifti right hemisphere filename
    brain_fn : brain data in cifti format
    return_img : whether to include img in the return
    return_hemi_len : whether to include hemisphere lengths in the return
    return_59k_mask : whether to include a mask corresponding to cortex vertices 
                      (True) or medial wall vertices (False) for 59k data
    return_source_data : whether to include the source data in the return (both for GIFTI and CIFTI)
    
    Returns
    -------
    result : dict
        A dictionary containing the following keys:
        - 'data_concat': numpy array, stacked data of the two hemispheres. 2-dimensional array (time x vertices).
        - 'img_L': optional numpy array, surface image data for the left hemisphere.
        - 'img_R': optional numpy array, surface image data for the right hemisphere.
        - 'len_L': optional int, length of the left hemisphere data.
        - 'len_R': optional int, length of the right hemisphere data.
        - 'mask_59k': optional numpy array, mask where True corresponds to cortex vertices and False to medial wall vertices for 59k data.
        - 'source_data_L': optional numpy array, source data for the left hemisphere (only available if return_source_data is True).
        - 'source_data_R': optional numpy array, source data for the right hemisphere (only available if return_source_data is True).
        - 'source_data': optional numpy array, source data for the entire brain (only available if return_source_data is True).
    """
    
    import numpy as np
    from surface_utils import load_surface
    from cifti_utils import from_170k_to_59k
    
    result = {}

    if L_fn and R_fn: 
        img_L, data_L = load_surface(L_fn)
        len_L = data_L.shape[1]
        img_R, data_R = load_surface(R_fn)
        len_R = data_R.shape[1]
        data_concat = np.concatenate((data_L, data_R), axis=1)
        result['data_concat'] = data_concat
        if return_img:
            result['img_L'] = img_L
            result['img_R'] = img_R
        if return_hemi_len:
            result['len_L'] = len_L
            result['len_R'] = len_R
        if return_source_data:
            result['source_data_L'] = data_L
            result['source_data_R'] = data_R

    elif brain_fn:
        img, data = load_surface(brain_fn)
        result.update(from_170k_to_59k(img=img, 
                                        data=data, 
                                        return_concat_hemis=True, 
                                        return_59k_mask=return_59k_mask))

        if return_img:
            result['img'] = img
        if return_source_data:
            result['source_data'] = data

    return result

def make_image_pycortex(data, 
                        maps_names=None,
                        img_L=None, 
                        img_R=None, 
                        lh_vert_num=None, 
                        rh_vert_num=None, 
                        img=None, 
                        brain_mask_59k=None):
    """
    Make a Cifti or Gifti image with data imported by PyCortex. This means that Gifti data 
    will be split by hemisphere, and Cifti data will be transformed back into 170k size.

    Parameters:
    - data: numpy array, your data.
    - maps_names: list of strings, optional, names for the mapped data.
    - img_L: Gifti Surface, left hemisphere surface object.
    - img_R: Gifti Surface, right hemisphere surface object.
    - lh_vert_num: int, number of vertices in the left hemisphere.
    - rh_vert_num: int, number of vertices in the right hemisphere.
    - img: Cifti Surface, source volume for mapping onto the surface.
    - brain_mask_59k: numpy array, optional, brain mask for 59k vertices (output of the from_170k_to_59k function).

    Returns:
    If mapping onto separate hemispheres (img_L and img_R provided):
    - new_img_L: Gifti img, new surface representing data on the left hemisphere.
    - new_img_R: Gifti img, new surface representing data on the right hemisphere.

    If mapping onto a single hemisphere (img provided):
    - new_img: Cifti img, new surface representing data on 170k size.
    """
    from cifti_utils import from_59k_to_170k
    from surface_utils import make_surface_image 
  
    

    if img_L and img_R: 
        data_L = data[:,:lh_vert_num]
        data_R = data[:,-rh_vert_num:]

        new_img_L = make_surface_image(data_L, img_L, maps_names=maps_names)
        new_img_R = make_surface_image(data_R, img_R, maps_names=maps_names)
        return new_img_L, new_img_R
        
    elif img:
        data_170k = from_59k_to_170k(data_59k=data, 
                                     brain_mask_59k=brain_mask_59k)
                
        new_img = make_surface_image(data=data_170k, 
                                     source_img=img, 
                                     maps_names=maps_names)
        return new_img

def set_pycortex_config_file(cortex_folder):

    # Import necessary modules
    import os
    import sys
    import cortex
    from pathlib import Path

    # Get pycortex config file location
    pycortex_config_file  = cortex.options.usercfg

    # Define the filestore and colormaps path
    filestore_line = 'filestore={}/db/\n'.format(cortex_folder)
    colormaps_line = 'colormaps={}/colormaps/\n'.format(cortex_folder)
    
    # Check if path correct
    with open(pycortex_config_file, 'r') as fileIn:
        for line in fileIn:
            if 'filestore' in line:
                if line==filestore_line: correct_filestore = True
                else: correct_filestore = False
            elif 'colormaps' in line:
                if line==colormaps_line: correct_colormaps = True
                else: correct_colormaps = False
                    
    # Change config file
    if correct_filestore==False or correct_colormaps==False:

        # Create name of new config file that will be written
        new_pycortex_config_file = pycortex_config_file[:-4] + '_new.cfg'
    
        # Create the new config file
        Path(new_pycortex_config_file).touch()
    
        # Write back the two lines
        with open(pycortex_config_file, 'r') as fileIn:
            with open(new_pycortex_config_file, 'w') as fileOut:
                for line in fileIn:
                    if 'filestore' in line:
                        fileOut.write(filestore_line)
                    elif 'colormaps' in line:
                        fileOut.write(colormaps_line)
                    else:
                        fileOut.write(line)
                        
        # Renames the original config file
        os.rename(new_pycortex_config_file, pycortex_config_file)
        sys.exit('Pycortex config file changed: please restart your code')

    return None

def draw_cortex(subject, data, vmin, vmax, description, cortex_type='VolumeRGB', cmap='Viridis',\
                cbar = 'discrete', cmap_dict=None, cmap_steps=255, xfmname=None, \
                alpha=None, depth=1, thick=1, height=1024, sampler='nearest',\
                with_curvature=True, with_labels=False, with_colorbar=False,\
                with_borders=False, curv_brightness=0.95, curv_contrast=0.05, add_roi=False,\
                roi_name='empty', col_offset=0, zoom_roi=None, zoom_hem=None, zoom_margin=0.0, cbar_label='', \
                overlay_fn='None'):
    """
    Plot brain data onto a previously saved flatmap.
    
    Parameters
    ----------
    subject             : subject id (e.g. 'sub-001')
    xfmname             : xfm transform
    data                : the data you would like to plot on a flatmap
    cmap                : colormap that shoudl be used for plotting
    cmap_dict           : colormap dict of label and color for personalized colormap
    vmins               : minimal values of 1D 2D colormap [0] = 1D, [1] = 2D
    vmaxs               : minimal values of 1D/2D colormap [0] = 1D, [1] = 2D
    description         : plot title
    cortex_type         : cortex function to create the volume (VolumeRGB, Volume2D, VertexRGB)
    cbar                : color bar layout
    cbar_label          : colorbar label
    cmap_steps          : number of colormap bins
    alpha               : alpha map
    depth               : Value between 0 and 1 for how deep to sample the surface for the flatmap (0 = gray/white matter boundary, 1 = pial surface)
    thick               : Number of layers through the cortical sheet to sample. Only applies for pixelwise = True
    height              : Height of the image to render. Automatically scales the width for the aspect of the subject's flatmap
    sampler             : Name of sampling function used to sample underlying volume data. Options include 'trilinear', 'nearest', 'lanczos'
    with_curvature      : Display the rois, labels, colorbar, annotated flatmap borders, or cross-hatch dropout?
    with_labels         : Display labels?
    with_colorbar       : Display pycortex colorbar?
    with_borders        : Display borders?
    curv_brightness     : Mean brightness of background. 0 = black, 1 = white, intermediate values are corresponding grayscale values.
    curv_contrast       : Contrast of curvature. 1 = maximal contrast (black/white), 0 = no contrast (solid color for curvature equal to curvature_brightness).
    add_roi             : add roi -image- to overlay.svg
    roi_name            : roi name
    col_offset          : colormap offset between 0 and 1
    zoom_roi            : name of the roi on which to zoom on
    zoom_hem            : hemifield fo the roi zoom
    zoom_margin         : margin in mm around the zoom
    overlay_fn          : file name of the overlay file (e.g. 'overlay_rois-drawn.svg')
    
    Returns
    -------
    braindata - pycortex volumr or vertex file
    """
    
    import cortex
    import numpy as np
    import matplotlib.pyplot as plt
    import matplotlib.colors as colors
    from matplotlib import cm
    import matplotlib as mpl
    import ipdb
    
    deb = ipdb.set_trace
    
    # define colormap
    try: base = plt.cm.get_cmap(cmap)
    except: base = cortex.utils.get_cmap(cmap)

    

    if overlay_fn == 'None': 
        overlay_file = None
    else:
        # define overlay path
        pycortex_config_file  = cortex.options.usercfg
        with open(pycortex_config_file, 'r') as fileIn:
            for line in fileIn:
                if 'filestore' in line:
                    db_path=line[10:-2]
        overlay_file = f"{db_path}/{subject}/{overlay_fn}"


    if '_alpha' in cmap: base.colors = base.colors[1,:,:]
    val = np.linspace(0, 1, cmap_steps, endpoint=False)
    
    colmap = colors.LinearSegmentedColormap.from_list('my_colmap', base(val), N=cmap_steps)

    
    if cortex_type=='VolumeRGB':
        # convert data to RGB
        vrange = float(vmax) - float(vmin)
        norm_data = ((data-float(vmin))/vrange)*cmap_steps
        mat = colmap(norm_data.astype(int))*255.0
        alpha = alpha*255.0

        # define volume RGB
        braindata = cortex.VolumeRGB(channel1 = mat[...,0].T.astype(np.uint8),
                                     channel2 = mat[...,1].T.astype(np.uint8),
                                     channel3 = mat[...,2].T.astype(np.uint8),
                                     alpha = alpha.T.astype(np.uint8),
                                     subject = subject,
                                     xfmname = xfmname)
    elif cortex_type=='Volume2D':
        braindata = cortex.Volume2D(dim1 = data.T,
                                 dim2 = alpha.T,
                                 subject = subject,
                                 xfmname = xfmname,
                                 description = description,
                                 cmap = cmap,
                                 vmin = vmin[0],
                                 vmax = vmax[0],
                                 vmin2 = vmin[1],
                                 vmax2 = vmax[1])
    elif cortex_type=='VertexRGB':
        
        # convert data to RGB
        vrange = float(vmax) - float(vmin)
        norm_data = ((data-float(vmin))/vrange)*cmap_steps
        mat = colmap(norm_data.astype(int))*255.0
        alpha = alpha*255.0
        
        # define Vertex RGB
        braindata = cortex.VertexRGB( red = mat[...,0].astype(np.uint8),
                                      green = mat[...,1].astype(np.uint8),
                                      blue = mat[...,2].astype(np.uint8),
                                      subject = subject,
                                      alpha = alpha.astype(np.uint8))
        braindata = braindata.blend_curvature(alpha)
        
    elif cortex_type=='Vertex':
        
        # define Vertex 
        braindata = cortex.Vertex(data = data,
                                 subject = subject,
                                 description = description,
                                 cmap = cmap,
                                 vmin = vmin,
                                 vmax = vmax)

    braindata_fig = cortex.quickshow(braindata = braindata,
                                     depth = depth,
                                     thick = thick,
                                     height = height,
                                     sampler = sampler,
                                     with_curvature = with_curvature,
                                     nanmean = True,
                                     overlay_file=overlay_file,
                                     with_labels = with_labels,
                                     with_colorbar = with_colorbar,
                                     with_borders = with_borders,
                                     curvature_brightness = curv_brightness,
                                     curvature_contrast = curv_contrast)
    if cbar == 'polar':
        try: base = plt.cm.get_cmap(cmap)
        except: base = cortex.utils.get_cmap(cmap)
        val = np.arange(1,cmap_steps+1)/cmap_steps - (1/(cmap_steps*2))
        val = np.fmod(val+col_offset,1)
        cbar_axis = braindata_fig.add_axes([0.5, 0.07, 0.8, 0.2], projection='polar')
        norm = colors.Normalize(0, 2*np.pi)
        t = np.linspace(0,2*np.pi,200,endpoint=True)
        r = [0,1]
        rg, tg = np.meshgrid(r,t)
        im = cbar_axis.pcolormesh(t, r, tg.T,norm=norm, cmap=colmap)
        cbar_axis.set_yticklabels([])
        cbar_axis.set_xticklabels([])
        cbar_axis.set_theta_zero_location("W")
        cbar_axis.spines['polar'].set_visible(False)

    elif cbar == 'ecc':
        colorbar_location = [0.5, 0.07, 0.8, 0.2]
        n = 200
        cbar_axis = braindata_fig.add_axes(colorbar_location, projection='polar')
        t = np.linspace(0,2*np.pi, n)
        r = np.linspace(0,1, n)
        rg, tg = np.meshgrid(r,t)
        c = tg
        im = cbar_axis.pcolormesh(t, r, c, norm = mpl.colors.Normalize(0, 2*np.pi), cmap=colmap)
        cbar_axis.tick_params(pad=1,labelsize=15)
        cbar_axis.spines['polar'].set_visible(False)
        box = cbar_axis.get_position()
        cbar_axis.set_yticklabels([])
        cbar_axis.set_xticklabels([])
        axl = braindata_fig.add_axes([0.97*box.xmin,0.5*(box.ymin+box.ymax), box.width/600,box.height*0.5])
        axl.spines['top'].set_visible(False)
        axl.spines['right'].set_visible(False)
        axl.spines['bottom'].set_visible(False)
        axl.yaxis.set_ticks_position('right')
        axl.xaxis.set_ticks_position('none')
        axl.set_xticklabels([])
        axl.set_yticklabels(np.linspace(vmin, vmax, 3),size = 'x-large')
        axl.set_ylabel('$dva$\t\t', rotation=0, size='x-large')
        axl.yaxis.set_label_coords(box.xmax+30,0.4)
        axl.patch.set_alpha(0.5)

    elif cbar == 'discrete':
        colorbar_location= [0.8, 0.05, 0.1, 0.05]
        cmaplist = [colmap(i) for i in range(colmap.N)]
        bounds = np.linspace(vmin, vmax, cmap_steps + 1)  
        bounds_label = np.linspace(vmin, vmax, 3)
        norm = mpl.colors.BoundaryNorm(bounds, colmap.N)
        cbar_axis = braindata_fig.add_axes(colorbar_location)
        cb = mpl.colorbar.ColorbarBase(cbar_axis, cmap=colmap, norm=norm, ticks=bounds_label, boundaries=bounds,orientation='horizontal')
        cb.set_label(cbar_label,size='x-large')

    elif cbar == '2D':
        cbar_axis = braindata_fig.add_axes([0.8, 0.05, 0.15, 0.15])
        base = cortex.utils.get_cmap(cmap)
        cbar_axis.imshow(np.dstack((base.colors[...,0], base.colors[...,1], base.colors[...,2],base.colors[...,3])))
        cbar_axis.set_xticks(np.linspace(0,255,3))
        cbar_axis.set_yticks(np.linspace(0,255,3))
        cbar_axis.set_xticklabels(np.linspace(vmin[0],vmax[0],3))
        cbar_axis.set_yticklabels(np.linspace(vmax[1],vmin[1],3))
        cbar_axis.set_xlabel(cbar_label[0], size='x-large')
        cbar_axis.set_ylabel(cbar_label[1], size='x-large')
        
    elif cbar == 'discrete_personalized':
        colorbar_location = [0.05, 0.02, 0.04, 0.3]
        cbar_axis = braindata_fig.add_axes(colorbar_location)
        norm = mpl.colors.BoundaryNorm(np.linspace(0, len(cmap_dict), len(cmap_dict)+1),
                                       len(cmap_dict))
        
        cb = mpl.colorbar.ColorbarBase(cbar_axis,
                                       cmap=base.reversed(),
                                       norm=norm, 
                                       ticks=(np.arange(0, len(cmap_dict), 1) + 0.5),
                                       orientation='vertical')
        cb.set_ticklabels(list(reversed(cmap_dict.keys())))
        cb.ax.tick_params(size=0, labelsize=15) 
        
    elif cbar == 'glm':
        
        val = np.linspace(0, 1, cmap_steps + 1, endpoint=False)

        # Exclure les valeurs proches du blanc
        val = val[val > 0.25]
        
        colmapglm = colors.LinearSegmentedColormap.from_list('my_colmap', base(val), N=len(val))
        colorbar_location = [0.85, 0.02, 0.04, 0.2]
        bounds_label = ['Both','Saccade','Pursuit']  
        bounds = np.linspace(vmin, vmax, colmap.N) 
        ticks_positions = [0.5, 1.5, 2.5]  
        norm = mpl.colors.BoundaryNorm(bounds, colmap.N)
        cbar_axis = braindata_fig.add_axes(colorbar_location)
        cb = mpl.colorbar.ColorbarBase(cbar_axis, cmap=colmapglm.reversed(), norm=norm, ticks=ticks_positions, orientation='vertical')
        cb.set_ticklabels(bounds_label)
        cb.ax.tick_params(size=0,labelsize=20) 
    elif cbar == 'stats':
        # colmap = colors.LinearSegmentedColormap.from_list('my_colmap', base(val), N=cmap_steps)
        val = np.linspace(0, 1, cmap_steps + 1, endpoint=False)
        val = val[val > 0.13]
        
        colmapglm = colors.LinearSegmentedColormap.from_list('my_colmap', base(val), N=len(val))
        colorbar_location = [0.05, 0.02, 0.04, 0.2]
        bounds_label = ['pursuit', 'saccade', 'pursuit_and_saccade', 'vision', 'vision_and_pursuit', 'vision_and_saccade', 'vision_and_saccade_and_pursuite']  
        bounds = np.linspace(vmin, vmax, colmap.N) 
        ticks_positions = [6.5, 5.5, 4.5, 3.5, 2.5, 1.5, 0.5] 
        norm = mpl.colors.BoundaryNorm(bounds, colmap.N)
        cbar_axis = braindata_fig.add_axes(colorbar_location)
        cb = mpl.colorbar.ColorbarBase(cbar_axis, cmap=colmapglm.reversed(), norm=norm, ticks=ticks_positions, orientation='vertical')
        cb.set_ticklabels(bounds_label)
        cb.ax.tick_params(size=0,labelsize=20) 
    
    # add to overlay
    if add_roi == True:
        cortex.utils.add_roi(   data = braindata,
                                name = roi_name,
                                open_inkscape = False,
                                add_path = False,
                                depth = depth,
                                thick = thick,
                                sampler = sampler,
                                with_curvature = with_curvature,
                                with_colorbar = with_colorbar,
                                with_borders = with_borders,
                                curvature_brightness = curv_brightness,
                                curvature_contrast = curv_contrast)

    return braindata

def create_colormap(cortex_dir, colormap_name, colormap_dict, recreate=False):
    """
    Add a 1 dimensional colormap in pycortex dataset
    
    Parameters
    ----------
    cortex_colormaps_dir:          directory of the pycortex dataset colormaps folder
    colormap_name:                 name of the colormap
    colormap_dict:                 dict containing keys of the color refence and tuple of the color list
    
    Returns
    -------
    colormap PNG in pycortex dataset colormaps folder
    """
    from PIL import Image
    import os
    import ipdb
    deb=ipdb.set_trace

    colormap_fn = '{}/colormaps/{}.png'.format(cortex_dir, colormap_name)
    
    # Create image of the colormap
    if (os.path.isfile(colormap_fn) == False) or recreate: 
        image = Image.new("RGB", (len(colormap_dict), 1))
        i = 0
        for color in colormap_dict.values():
            image.putpixel((i, 0), color)
            i +=1
        print('Saving new colormap: {}'.format(colormap_fn))
        image.save(colormap_fn)
        

    return None

# Check and setup pycortex directory structure
def setup_pycortex_dirs(cortex_dir):
    import os
    import urllib.request
    import json

    """
    - Create cortex/colormaps and cortex/db directories if missing
    - Download pycortex colormaps if not present
    - Download HCP pycortex subjects (sub-hcp32k, sub-hcp59k) from GitHub if missing
    """

    colormaps_dir = os.path.join(cortex_dir, "colormaps")
    db_dir = os.path.join(cortex_dir, "db")

    # Create directories if they don't exist
    os.makedirs(colormaps_dir, exist_ok=True)
    os.makedirs(db_dir, exist_ok=True)

    # ------------------------------------------------------------------
    # 1) Pycortex colormaps
    # ------------------------------------------------------------------
    if not os.listdir(colormaps_dir):
        print("Downloading colormaps from pycortex GitHub...")
        api_url = "https://api.github.com/repos/gallantlab/pycortex/contents/filestore/colormaps"

        try:
            with urllib.request.urlopen(api_url) as response:
                files = json.loads(response.read())

            # Download each colormap file
            for file_info in files:
                if file_info["type"] == "file":
                    dst = os.path.join(colormaps_dir, file_info["name"])
                    print(f"  Downloading {file_info['name']}...")
                    urllib.request.urlretrieve(file_info["download_url"], dst)

            print("Colormaps downloaded successfully.")
        except Exception as e:
            print(f"Warning: Could not download colormaps: {e}")
    else:
        print("Colormaps directory already contains files.")

    # ------------------------------------------------------------------
    # 2) HCP pycortex subjects (sub-hcp32k / sub-hcp59k)
    # ------------------------------------------------------------------
    def download_github_dir(api_url, local_dir):
        """Recursively download a GitHub directory using the GitHub API."""
        os.makedirs(local_dir, exist_ok=True)

        with urllib.request.urlopen(api_url) as response:
            items = json.loads(response.read())

        for item in items:
            local_path = os.path.join(local_dir, item["name"])

            if item["type"] == "file":
                print(f"  Downloading {local_path}")
                urllib.request.urlretrieve(item["download_url"], local_path)

            elif item["type"] == "dir":
                download_github_dir(item["url"], local_path)

    github_base = (
        "https://api.github.com/repos/ulascombes/"
        "HCP_pycortex_subjects/contents/data/cortex/db"
    )

    for subject in ["sub-hcp1.6mm", "sub-hcp2.0mm"]:
        local_subject_dir = os.path.join(db_dir, subject)

        if not os.path.exists(local_subject_dir):
            print(f"Downloading {subject} from GitHub...")
            api_url = f"{github_base}/{subject}"

            try:
                download_github_dir(api_url, local_subject_dir)
                print(f"{subject} downloaded successfully.")
            except Exception as e:
                print(f"Warning: Could not download {subject}: {e}")
        else:
            print(f"{subject} already exists.")

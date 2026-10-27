"""
-----------------------------------------------------------------------------------------
compute_css_pcm.py
-----------------------------------------------------------------------------------------
Goal of the script:
Compute population cortical magnification and add to derivatives
Note:
CM is computed using the geodesic distances (mm) of vertices located within a radius on
the cortical surface set by pcm_surf_type (see vertex_pcm_rad) and restricted by the ROI boundaries
Optimization (same computation, same outputs):
- surfaces and ROI vertices are loaded once per format / ROI method
- geodesic patches (get_geodesic_patch) depend only on the surface and the vertex:
  they are computed once per vertex and reused across averaging methods and tasks
- ROI exclusion is vectorized (boolean mask instead of python list comprehension)
[CORRECTION] fichier entierement reorganise pour la vitesse (boucles) : comparer avec l'original
par diff ; les seules corrections de fond sont marquees #7 et #24
-----------------------------------------------------------------------------------------
Input(s):
sys.argv[1]: main project directory
sys.argv[2]: project name (correspond to directory)
sys.argv[3]: subject name (e.g. sub-01)
sys.argv[4]: analysis name (e.g. prf)
sys.argv[5]: group (e.g. 327)
-----------------------------------------------------------------------------------------
Output(s):
New brain volume in derivative nifti file
-----------------------------------------------------------------------------------------
To run:
1. cd to function
>> cd ~/projects/pRF_analysis/analysis_code/postproc/prf/postfit
2. run python command
>> python compute_pcm.py [main directory] [project name] [subject] [analysis name] [group]
-----------------------------------------------------------------------------------------
Exemple:
cd ~/projects/pRF_analysis/analysis_code/postproc/prf/postfit
python compute_css_pcm.py /scratch/mszinte/data RetinoMaps sub-01 prf 327
python compute_css_pcm.py /scratch/mszinte/data RetinoMaps template_avg prf 327
-----------------------------------------------------------------------------------------
Written by Martin Szinte (martin.szinte@gmail.com)
and Uriel Lascombes (uriel.lascombes@laposte.net)
-----------------------------------------------------------------------------------------
"""
# Stop warnings
import warnings
warnings.filterwarnings("ignore")

# Debug
import ipdb
deb = ipdb.set_trace

# General imports
import os
# [CORRECTION #74] 1 thread BLAS par processus : les patchs geodesiques sont calcules en parallele (un processus par coeur)
# AVANT : (ligne ajoutee)
for _var in ['OMP_NUM_THREADS', 'MKL_NUM_THREADS', 'OPENBLAS_NUM_THREADS']: os.environ.setdefault(_var, '1')
import sys
# [CORRECTION #74] calcul parallele des patchs geodesiques
# AVANT : (ligne ajoutee)
import multiprocessing
import cortex
import datetime
import numpy as np
import nibabel as nb

# Personal iports
sys.path.append("{}/../../../utils".format(os.getcwd()))
from settings_utils import load_settings
from surface_utils import make_surface_image
from maths_utils import  median_subject_template, weighted_nan_median
from pycortex_utils import set_pycortex_config_file, load_surface_pycortex, get_rois, make_image_pycortex

# [CORRECTION #74] fonction executee par chaque coeur : patch geodesique d'un vertex (code deplace depuis la boucle des ROI)
# AVANT : (fonction ajoutee ; le calcul etait fait vertex par vertex dans la boucle des ROI, sur un seul coeur)
def geo_patch_worker(key):
    hemi, surf_idx = key
    surf = geo_surfs[hemi]
    try:
        geo_patch = surf.get_geodesic_patch(radius=vert_dist_th, vertex=surf_idx)
        patch_idx = np.where(geo_patch['vertex_mask'])[0].astype(np.int32)
        patch_dist = np.asarray(geo_patch['geodesic_distance'])
    # [CORRECTION #24] en cas d'erreur, patch vide propre (l'original ecrivait dans geo_patch, non defini ou reste du vertex precedent)
    # AVANT : except Exception as e:
    # AVANT : print("Vertex #{}: error: {} within {} mm".format(vert_idx, e, vert_dist_th))
    # AVANT : geo_patch['vertex_mask'] = np.zeros(surf.pts.shape[0]).astype(bool)
    # AVANT : geo_patch['geodesic_distance'] = []
    except Exception as e:
        print("Vertex #{} ({}): error: {} within {} mm".format(surf_idx, hemi, e, vert_dist_th))
        patch_idx = np.array([], dtype=np.int32)
        patch_dist = np.array([])
    return patch_idx, patch_dist


# Inputs
main_dir = sys.argv[1]
project_dir = sys.argv[2]
subject = sys.argv[3]
analysis_name = sys.argv[4]
group = sys.argv[5]

# Load settings
base_dir = os.path.abspath(os.path.join(os.getcwd(), "../../../../"))
general_settings_path = os.path.join(base_dir, project_dir, "settings.yml")
analysis_settings_path = os.path.join(base_dir, project_dir, f"{analysis_name}-analysis.yml")
settings = load_settings([general_settings_path, analysis_settings_path])
analysis_info = settings[0]

vert_dist_th = analysis_info['vertex_pcm_rad']
formats = analysis_info['formats']
maps_names_css_stats = analysis_info['maps_names_css_stats']
maps_names_pcm = analysis_info['maps_names_pcm']
subjects = analysis_info['subjects']
task_names = analysis_info['analysis_task_names']
preproc_prep = analysis_info['preproc_prep']
filtering = analysis_info['filtering']
normalization = analysis_info['normalization']
avg_methods = analysis_info['avg_methods']
rois_methods = analysis_info['rois_methods']
pycortex_subject_template = analysis_info['pycortex_subject_template']
averaging_templates = analysis_info['averaging_templates']
formats_conversion = analysis_info['formats_conversion']
output_folder = analysis_info["output_folder"]
dm_name = analysis_info["dm_name"]

# [CORRECTION #7] surface des distances geodesiques lue dans prf-analysis.yml (pcm_surf_type) ; valeur retenue : 'fiducial'
# [CORRECTION #7] (surface 3D entre substance blanche et piale) : la surface aplatie deforme les distances
# AVANT : cortex.db.get_surf(pycortex_subject, "flat") ecrit en dur plus bas
# Surface used for geodesic distances
surf_type = analysis_info['pcm_surf_type']

# Set pycortex db and colormaps
cortex_dir = "{}/{}/derivatives/pp_data/cortex".format(main_dir, project_dir)
set_pycortex_config_file(cortex_dir)

# Derivatives and stats idx
# compute duration
start_time = datetime.datetime.now()

# template_avg exception
if subject != 'template_avg':

    for format_ in formats:

        # define directories
        prf_dir = f"{main_dir}/{project_dir}/derivatives/pp_data/{subject}/{format_}/{output_folder}"
        prf_deriv_dir = "{}/prf_derivatives".format(prf_dir)

        # Format settings
        if format_ == 'fsnative':
            pycortex_subject = subject
            surf_format = format_
        elif format_ == '170k':
            pycortex_subject = pycortex_subject_template
            surf_format = formats_conversion[format_]

        # Get surfaces for each hemisphere (loaded once per format)
        surfs = [cortex.polyutils.Surface(*d) for d in cortex.db.get_surf(pycortex_subject, surf_type)]
        surf_lh, surf_rh = surfs[0], surfs[1]

        # [CORRECTION #74] surfaces visibles par les processus paralleles ; nombre de coeurs = --cpus-per-task du sbatch
        # AVANT : (lignes ajoutees)
        geo_surfs = {'lh': surf_lh, 'rh': surf_rh}
        nb_procs = int(os.environ.get('SLURM_CPUS_PER_TASK', 1))

        # [CORRECTION #74] cache {(hemi, surf_idx): (patch idx, distances geodesiques)} : le patch ne depend ni des donnees
        # ni de la ROI (le masque de ROI est applique a la lecture), donc il est reutilise pour toutes les ROI,
        # methodes de ROI, methodes de moyennage et taches du format
        # AVANT : geo_cache = {} (remis a zero a chaque methode de ROI, cle (roi, hemi, surf_idx), voir plus bas)
        geo_cache = {}

        # Get the vertices number per hemisphere
        lh_vert_num, rh_vert_num = surf_lh.pts.shape[0], surf_rh.pts.shape[0]
        vert_num = lh_vert_num + rh_vert_num

        # Loop over rois-methods
        rois_methods_format = rois_methods[format_]
        for rois_method_format in rois_methods_format:

            # get rois
            if rois_method_format == 'rois-drawn':
                rois = analysis_info[rois_method_format]
            elif rois_method_format == 'rois-group-mmp':
                rois = list(analysis_info[rois_method_format].keys())

            # Get a dict with the surface vertices contained in each ROI (once per ROI method)
            print('Get a dict with the surface vertices contained in each ROI...')
            roi_verts_dict = get_rois(pycortex_subject,
                                      surf_format=surf_format,
                                      rois_type=rois_method_format,
                                      mask=False,
                                      rois=rois,
                                      hemis=None)

            # ROI vertex indices and ROI boolean masks per hemisphere (once per ROI method)
            roi_hemi_info = {}
            for roi in rois:
                roi_vert_lh_idx = roi_verts_dict[roi][roi_verts_dict[roi] < lh_vert_num]
                roi_vert_rh_idx = roi_verts_dict[roi][roi_verts_dict[roi] >= lh_vert_num]
                roi_surf_lh_idx = roi_vert_lh_idx
                roi_surf_rh_idx = roi_vert_rh_idx - lh_vert_num

                roi_mask_lh = np.zeros(lh_vert_num, dtype=bool)
                roi_mask_lh[roi_surf_lh_idx] = True
                roi_mask_rh = np.zeros(rh_vert_num, dtype=bool)
                roi_mask_rh[roi_surf_rh_idx] = True

                roi_hemi_info[roi] = {'lh': (roi_vert_lh_idx, roi_surf_lh_idx, roi_mask_lh),
                                      'rh': (roi_vert_rh_idx, roi_surf_rh_idx, roi_mask_rh)}

            # [CORRECTION #74] cache deplace au niveau du format (voir plus haut)
            # AVANT : # Geodesic patches cache: {(roi, hemi, surf_idx): (neighbor idx in roi, geodesic distances)}
            # AVANT : # geodesic patches do not depend on the data, they are computed once per vertex
            # AVANT : # and reused for all averaging methods and tasks
            # AVANT : geo_cache = {}

            for avg_method in avg_methods:
                if 'loo' in avg_method:
                    maps_names_css = analysis_info['maps_names_css_loo']
                else:
                    maps_names_css = analysis_info['maps_names_css']

                # Maps settings
                for idx, col_name in enumerate(maps_names_css + maps_names_css_stats):
                    exec("{}_idx = idx".format(col_name))
                if 'loo' in avg_method: rsq_idx2use = prf_loo_rsq_idx
                else: rsq_idx2use = prf_rsq_idx

                for task_name in task_names:

                    # [CORRECTION #73] temps ecoule affiche pour suivre la duree de chaque etape
                    # AVANT : print(f'{avg_method} - {format_} - {rois_method_format} - {task_name}')
                    print(f'{avg_method} - {format_} - {rois_method_format} - {task_name} (elapsed {datetime.datetime.now() - start_time})')
                    if format_ == 'fsnative':
                        # Derivatives
                        deriv_fn_L = f'{prf_deriv_dir}/{subject}_task-{task_name}_hemi-L_{preproc_prep}_{filtering}_{normalization}_{avg_method}_{analysis_name}-css{dm_name}_deriv.func.gii'
                        deriv_fn_R = f'{prf_deriv_dir}/{subject}_task-{task_name}_hemi-R_{preproc_prep}_{filtering}_{normalization}_{avg_method}_{analysis_name}-css{dm_name}_deriv.func.gii'
                        # [CORRECTION #73] methode de moyennage sans fichiers (ex : concat non fitte) ignoree au lieu de faire planter le job
                        # AVANT : (ligne ajoutee)
                        if not (os.path.isfile(deriv_fn_L) and os.path.isfile(deriv_fn_R)):
                            print(f'[SKIP] missing file(s): {deriv_fn_L}')
                            continue
                        print(f'Loading:\t{deriv_fn_L}\n\t\t{deriv_fn_R}')
                        results = load_surface_pycortex(L_fn=deriv_fn_L,
                                                        R_fn=deriv_fn_R,
                                                        return_img=True)
                        deriv_mat = results['data_concat']
                        img_L = results['img_L']
                        img_R = results['img_R']

                        # Stats
                        stats_fn_L = f'{prf_deriv_dir}/{subject}_task-{task_name}_hemi-L_{preproc_prep}_{filtering}_{normalization}_{avg_method}_{analysis_name}-css{dm_name}_stats.func.gii'
                        stats_fn_R = f'{prf_deriv_dir}/{subject}_task-{task_name}_hemi-R_{preproc_prep}_{filtering}_{normalization}_{avg_method}_{analysis_name}-css{dm_name}_stats.func.gii'
                        print(f'Loading:\t{stats_fn_L}\n\t\t{stats_fn_R}')
                        stats_results = load_surface_pycortex(L_fn=stats_fn_L,
                                                              R_fn=stats_fn_R)
                        stats_mat = stats_results['data_concat']

                    elif format_ == '170k':
                        # Derivatives
                        deriv_fn = f'{prf_deriv_dir}/{subject}_task-{task_name}_{preproc_prep}_{filtering}_{normalization}_{avg_method}_{analysis_name}-css{dm_name}_deriv.dtseries.nii'
                        # [CORRECTION #73] methode de moyennage sans fichiers (ex : concat non fitte) ignoree au lieu de faire planter le job
                        # AVANT : (ligne ajoutee)
                        if not os.path.isfile(deriv_fn):
                            print(f'[SKIP] missing file: {deriv_fn}')
                            continue
                        print(f'Loading:\t{deriv_fn}')
                        results = load_surface_pycortex(brain_fn=deriv_fn,
                                                        return_img=True,
                                                        return_59k_mask=True,
                                                        return_source_data=True)
                        deriv_mat = results['data_concat']
                        mask_59k = results['mask_59k']
                        deriv_mat_170k = results['source_data']
                        img = results['img']

                        # Stats
                        stats_fn = f'{prf_deriv_dir}/{subject}_task-{task_name}_{preproc_prep}_{filtering}_{normalization}_{avg_method}_{analysis_name}-css{dm_name}_stats.dtseries.nii'
                        print(f'Loading:\t{stats_fn}')
                        stats_results = load_surface_pycortex(brain_fn=stats_fn)
                        stats_mat = stats_results['data_concat']

                    # Combine deriv and stats mat
                    deriv_mat = np.concatenate((deriv_mat, stats_mat))

                    # Threshold data
                    deriv_mat_th = deriv_mat
                    amp_down = deriv_mat_th[amplitude_idx,...] > 0
                    rsq_down = deriv_mat_th[rsq_idx2use,...] >= analysis_info['rsqr_th']
                    size_th_down = deriv_mat_th[prf_size_idx,...] >= analysis_info['size_th'][0]
                    size_th_up = deriv_mat_th[prf_size_idx,...] <= analysis_info['size_th'][1]
                    ecc_th_down = deriv_mat_th[prf_ecc_idx,...] >= analysis_info['ecc_th'][0]
                    ecc_th_up = deriv_mat_th[prf_ecc_idx,...] <= analysis_info['ecc_th'][1]
                    n_th_down = deriv_mat_th[prf_n_idx,...] >= analysis_info['n_th'][0]
                    n_th_up = deriv_mat_th[prf_n_idx,...] <= analysis_info['n_th'][1]
                    if analysis_info['stats_th'] == 0.05: stats_th_down = deriv_mat_th[corr_pvalue_5pt_idx,...] <= 0.05
                    elif analysis_info['stats_th'] == 0.01: stats_th_down = deriv_mat_th[corr_pvalue_1pt_idx,...] <= 0.01
                    # [CORRECTION #37] else ajoute : toute autre valeur de stats_th laissait la variable indefinie (NameError plus loin)
                    # AVANT : (ligne(s) ajoutee(s))
                    else: raise ValueError(f"stats_th must be 0.05 or 0.01 (got {analysis_info['stats_th']})")
                    all_th = np.array((amp_down, rsq_down, size_th_down, size_th_up, ecc_th_down,
                                       ecc_th_up,n_th_down, n_th_up,stats_th_down))
                    deriv_mat[rsq_idx2use, np.logical_and.reduce(all_th)==False]=0

                    # Derivatives settings
                    vert_rsq_data = deriv_mat[rsq_idx2use, ...]
                    vert_x_data = deriv_mat[prf_x_idx, ...]
                    vert_y_data = deriv_mat[prf_y_idx, ...]

                    # Split per hemisphere (once per task, not per ROI)
                    hemi_data = {'lh': (surf_lh, vert_rsq_data[:lh_vert_num], vert_x_data[:lh_vert_num], vert_y_data[:lh_vert_num]),
                                 'rh': (surf_rh, vert_rsq_data[lh_vert_num:], vert_x_data[lh_vert_num:], vert_y_data[lh_vert_num:])}

                    # [CORRECTION #74] patchs geodesiques manquants calcules en parallele sur nb_procs coeurs,
                    # pour tous les vertex utilises (rsq > 0) de toutes les ROI, avant la boucle des ROI
                    # AVANT : (lignes ajoutees ; calcul sequentiel vertex par vertex dans la boucle des ROI)
                    todo = []
                    for roi in rois:
                        for hemi in ['lh', 'rh']:
                            roi_surf_idx = roi_hemi_info[roi][hemi][1]
                            todo += [(hemi, int(i)) for i in roi_surf_idx[hemi_data[hemi][1][roi_surf_idx] > 0]]
                    todo = [key for key in dict.fromkeys(todo) if key not in geo_cache]
                    print(f'Geodesic patches to compute: {len(todo)} on {nb_procs} core(s) (elapsed {datetime.datetime.now() - start_time})')
                    if len(todo) > 0:
                        with multiprocessing.get_context('fork').Pool(nb_procs) as pool:
                            for key, patch in zip(todo, pool.imap(geo_patch_worker, todo, chunksize=20)):
                                geo_cache[key] = patch
                    print(f'Geodesic patches done (elapsed {datetime.datetime.now() - start_time})')

                    # Create empty results
                    vert_cm = np.zeros((4,vert_num))*np.nan
                    for roi in rois:
                        for hemi in ['lh', 'rh']:
                            surf, vert_rsq, vert_x, vert_y = hemi_data[hemi]
                            roi_vert_idx, roi_surf_idx, roi_mask = roi_hemi_info[roi][hemi]

                            print('ROI -> {} / Hemisphere -> {}'.format(roi, hemi))
                            for vert_idx, surf_idx in zip(roi_vert_idx, roi_surf_idx):
                                if vert_rsq[surf_idx] > 0:

                                    # [CORRECTION #74] patch lu dans le cache (calcule en parallele plus haut) puis restreint a la ROI
                                    # AVANT : # Get geodesic distances (mm), computed once per vertex
                                    # AVANT : cache_key = (roi, hemi, surf_idx)
                                    # AVANT : if cache_key not in geo_cache:
                                    # AVANT : try:
                                    # AVANT : geo_patch = surf.get_geodesic_patch(radius=vert_dist_th,
                                    # AVANT : vertex=surf_idx)
                                    # AVANT : patch_idx = np.where(geo_patch['vertex_mask'])[0]
                                    # AVANT : patch_dist = np.asarray(geo_patch['geodesic_distance'])
                                    # AVANT : except Exception as e:
                                    # AVANT : print("Vertex #{}: error: {} within {} mm".format(vert_idx, e, vert_dist_th))
                                    # AVANT : patch_idx = np.array([], dtype=int)
                                    # AVANT : patch_dist = np.array([])
                                    # AVANT : # Exclude vextex out of roi
                                    # AVANT : in_roi = roi_mask[patch_idx]
                                    # AVANT : geo_cache[cache_key] = (patch_idx[in_roi], patch_dist[in_roi])
                                    # AVANT : neigh_idx, neigh_dist = geo_cache[cache_key]
                                    patch_idx, patch_dist = geo_cache[(hemi, int(surf_idx))]
                                    # Exclude vextex out of roi
                                    in_roi = roi_mask[patch_idx]
                                    neigh_idx, neigh_dist = patch_idx[in_roi], patch_dist[in_roi]

                                    # Compute the n neighbor of each vertex
                                    n_neighbor = neigh_idx.shape[0]
                                    vert_cm[0,vert_idx] = n_neighbor

                                    if n_neighbor > 1:

                                        # Get prf parameters of vertices (center and surround) in geodesic distance threshold
                                        srd = neigh_idx != surf_idx
                                        srd_idx, srd_dist = neigh_idx[srd], neigh_dist[srd]
                                        vert_ctr_x, vert_ctr_y = vert_x[surf_idx], vert_y[surf_idx]
                                        vert_srd_x = vert_x[srd_idx]
                                        vert_srd_y = vert_y[srd_idx]
                                        vert_srd_rsq = vert_rsq[srd_idx]

                                        # median
                                        # Compute median geodesic distance
                                        vert_geo_dist_median = weighted_nan_median(srd_dist, vert_srd_rsq)

                                        # Compute prf center suround distance (deg)
                                        vert_prf_dist_median_array = np.sqrt((vert_ctr_x - vert_srd_x)**2 + (vert_ctr_y - vert_srd_y)**2)

                                        # Compute median of prf center suround distance
                                        vert_prf_dist_median = weighted_nan_median(vert_prf_dist_median_array, vert_srd_rsq)

                                        # Compute cortical magnification in mm/deg (surface distance / pRF positon distance)
                                        vert_cm[1, vert_idx] = vert_geo_dist_median/vert_prf_dist_median

                                        # export median geodesic and prf distance
                                        vert_cm[2, vert_idx] = vert_geo_dist_median
                                        vert_cm[3, vert_idx] = vert_prf_dist_median

                    deriv_mat_new = vert_cm

                    # Save data
                    if format_ == 'fsnative':
                        # Save as image
                        new_img_L, new_img_R  = make_image_pycortex(data=deriv_mat_new,
                                                                    maps_names=maps_names_pcm,
                                                                    img_L=img_L,
                                                                    img_R=img_R,
                                                                    lh_vert_num=lh_vert_num,
                                                                    rh_vert_num=rh_vert_num,
                                                                    img=None,
                                                                    brain_mask_59k=None)

                        deriv_fn_L = f'{prf_deriv_dir}/{subject}_task-{task_name}_hemi-L_{preproc_prep}_{filtering}_{normalization}_{avg_method}_{rois_method_format}_{analysis_name}-css{dm_name}_pcm.func.gii'
                        deriv_fn_R = f'{prf_deriv_dir}/{subject}_task-{task_name}_hemi-R_{preproc_prep}_{filtering}_{normalization}_{avg_method}_{rois_method_format}_{analysis_name}-css{dm_name}_pcm.func.gii'
                        # [CORRECTION #73] le chemin est deja complet : l'ancien print affichait le dossier deux fois (fichier bien sauve, affichage seulement)
                        # AVANT : print('Saving {}/{}'.format(prf_deriv_dir, deriv_fn_L))
                        print('Saving {}'.format(deriv_fn_L))
                        nb.save(new_img_L, deriv_fn_L)
                        # [CORRECTION #73] le chemin est deja complet : l'ancien print affichait le dossier deux fois (fichier bien sauve, affichage seulement)
                        # AVANT : print('Saving {}/{}'.format(prf_deriv_dir, deriv_fn_R))
                        print('Saving {}'.format(deriv_fn_R))
                        nb.save(new_img_R, deriv_fn_R)

                    elif format_ == '170k':
                        # Save as image
                        new_img = make_image_pycortex(data=deriv_mat_new,
                                              maps_names=None,
                                              img_L=None,
                                              img_R=None,
                                              lh_vert_num=None,
                                              rh_vert_num=None,
                                              img=img,
                                              brain_mask_59k=mask_59k)

                        deriv_fn = f'{prf_deriv_dir}/{subject}_task-{task_name}_{preproc_prep}_{filtering}_{normalization}_{avg_method}_{rois_method_format}_{analysis_name}-css{dm_name}_pcm.dtseries.nii'
                        # [CORRECTION #73] le chemin est deja complet : l'ancien print affichait le dossier deux fois (fichier bien sauve, affichage seulement)
                        # AVANT : print('Saving {}/{}'.format(prf_deriv_dir, deriv_fn))
                        print('Saving {}'.format(deriv_fn))
                        nb.save(new_img, deriv_fn)

            # [CORRECTION #74] cache conserve pour la methode de ROI suivante (libere a la fin du format)
            # AVANT : # Free memory before next ROI method
            # AVANT : del geo_cache

        # [CORRECTION #74] liberation memoire avant le format suivant
        # AVANT : (ligne ajoutee)
        del geo_cache

# template_avg median
elif subject == 'template_avg':
    for averaging_template_name, averaging_template_format in averaging_templates.items():
        print('{}, Median corr across subject...'. format(averaging_template_name))

        for task_name in task_names:
            for avg_method in avg_methods:

                rois_methods_format = rois_methods['{}'.format(averaging_template_format)]
                for rois_method_format in rois_methods_format:

                    # find all the subject prf stats
                    prf_pcm_fns = []
                    for subject in subjects:
                        prf_deriv_dir = f"{main_dir}/{project_dir}/derivatives/pp_data/{subject}/{averaging_template_format}/{output_folder}/prf_derivatives"
                        prf_pcm_fns += [f"{prf_deriv_dir}/{subject}_task-{task_name}_{preproc_prep}_{filtering}_{normalization}_{avg_method}_{rois_method_format}_{analysis_name}-css{dm_name}_pcm.dtseries.nii"]

                    # [CORRECTION #73] methode de moyennage absente chez les sujets (ex : concat) ignoree au lieu de faire planter le job
                    # AVANT : (ligne ajoutee)
                    if not all(os.path.isfile(fn) for fn in prf_pcm_fns):
                        print(f'[SKIP] missing subject pcm file(s) for {avg_method} - {rois_method_format}')
                        continue

                     # Computing  across subject
                    img, data_pcm_median = median_subject_template(fns=prf_pcm_fns)

                    # Export results
                    template_pcm_dir = f"{main_dir}/{project_dir}/derivatives/pp_data/{averaging_template_name}/{averaging_template_format}/{output_folder}/prf_derivatives"
                    os.makedirs(template_pcm_dir, exist_ok=True)

                    template_pcm_fn = f"{template_pcm_dir}/{averaging_template_name}_task-{task_name}_{preproc_prep}_{filtering}_{normalization}_{avg_method}_{rois_method_format}_{analysis_name}-css{dm_name}_pcm.dtseries.nii"
                    print("saving: {}".format(template_pcm_fn))
                    template_pcm_img = make_surface_image(data=data_pcm_median,
                                                          source_img=img,
                                                          maps_names=maps_names_pcm)
                    nb.save(template_pcm_img, template_pcm_fn)

# Print duration
end_time = datetime.datetime.now()
print("\nStart time:\t{start_time}\nEnd time:\t{end_time}\nDuration:\t{dur}".format(
    start_time=start_time,
    end_time=end_time,
    dur=end_time - start_time))

# Define permission cmd
print('Changing files permissions in {}/{}'.format(main_dir, project_dir))
os.system("chmod -Rf 771 {}/{}".format(main_dir, project_dir))
os.system("chgrp -Rf {} {}/{}".format(group, main_dir, project_dir))

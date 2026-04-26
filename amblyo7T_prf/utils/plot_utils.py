# General imports
import numpy as np
import pandas as pd
import os, sys
import ipdb 
deb = ipdb.set_trace

# Figure imports
import plotly.graph_objects as go
from plotly.subplots import make_subplots

# Personal import
sys.path.append("{}/../../../../analysis_code/utils".format(os.getcwd()))
from maths_utils import weighted_regression


def weighted_deming_regression(x, y, weights=None):
    """
    Weighted Deming (orthogonal) regression.
    Minimizes perpendicular distances, symmetric when x and y are swapped.

    Parameters
    ----------
    x, y : array-like
        Data arrays
    weights : array-like, optional
        Weights for each point (e.g. mean R²)

    Returns
    -------
    slope, intercept : float
    """
    x = np.array(x)
    y = np.array(y)

    if weights is None:
        weights = np.ones(len(x))
    weights = np.array(weights)
    weights = weights / weights.sum()

    # Weighted means
    x_mean = np.sum(weights * x)
    y_mean = np.sum(weights * y)

    # Weighted variances and covariance
    sxx = np.sum(weights * (x - x_mean) ** 2)
    syy = np.sum(weights * (y - y_mean) ** 2)
    sxy = np.sum(weights * (x - x_mean) * (y - y_mean))

    # Deming slope (assumes equal error variance in x and y, lambda=1)
    slope = (syy - sxx + np.sqrt((syy - sxx) ** 2 + 4 * sxy ** 2)) / (2 * sxy)
    intercept = y_mean - slope * x_mean

    return slope, intercept


def plotly_template(template_specs):
    """
    Define the template for plotly
    Parameters
    ----------
    template_specs : dict
        dictionary containing specific figure settings
    Returns
    -------
    fig_template : plotly.graph_objs.layout._template.Template
        Template for plotly figure
    """
    fig_template = go.layout.Template()

    # Violin plots
    fig_template.data.violin = [go.Violin(
                                    box_visible=False,
                                    points=False,
                                    line_color="rgba(0, 0, 0, 1)",
                                    line_width=template_specs['rois_plot_width'],
                                    width=0.8,
                                    hoveron='violins',
                                    meanline_visible=False,
                                    showlegend=False)]

    # Barpolar
    fig_template.data.barpolar = [go.Barpolar(
                                    marker_line_color="rgba(0,0,0,1)",
                                    marker_line_width=template_specs['rois_plot_width'],
                                    showlegend=False)]

    # Pie plots
    fig_template.data.pie = [go.Pie(textposition=["inside", "none"],
                                    marker_line_width=0,
                                    rotation=0,
                                    direction="clockwise",
                                    hole=0.4,
                                    sort=False)]

    # Layout
    fig_template.layout = (go.Layout(
                                    font_family=template_specs['font'],
                                    font_size=template_specs['axes_font_size'],
                                    plot_bgcolor=template_specs['bg_col'],
                                    xaxis_visible=True,
                                    xaxis_linewidth=template_specs['axes_width'],
                                    xaxis_color=template_specs['axes_color'],
                                    xaxis_showgrid=False,
                                    xaxis_ticks="outside",
                                    xaxis_ticklen=8,
                                    xaxis_tickwidth=template_specs['axes_width'],
                                    xaxis_title_font_family=template_specs['font'],
                                    xaxis_title_font_size=template_specs['title_font_size'],
                                    xaxis_tickfont_family=template_specs['font'],
                                    xaxis_tickfont_size=template_specs['axes_font_size'],
                                    xaxis_zeroline=False,
                                    yaxis_visible=True,
                                    yaxis_linewidth=template_specs['axes_width'],
                                    yaxis_color=template_specs['axes_color'],
                                    yaxis_showgrid=False,
                                    yaxis_ticks="outside",
                                    yaxis_ticklen=8,
                                    yaxis_tickwidth=template_specs['axes_width'],
                                    yaxis_tickfont_family=template_specs['font'],
                                    yaxis_tickfont_size=template_specs['axes_font_size'],
                                    yaxis_title_font_family=template_specs['font'],
                                    yaxis_title_font_size=template_specs['title_font_size'],
                                    yaxis_zeroline=False))

    fig_template.layout.annotationdefaults = go.layout.Annotation(
                                    font_color=template_specs['axes_color'],
                                    font_family=template_specs['font'],
                                    font_size=template_specs['title_font_size'])

    return fig_template


def corr_plot(tsv_dir, subject, fn_spec_combined, figure_info, rsq2use):
    """
    Make correlation plot across eyes (FE/LE x-axis vs AE/RE y-axis).
    One figure with rows = ROIs, columns = corr_params.
    KDE contour from pre-computed TSV, regression line from regression TSV.

    Parameters
    ----------
    tsv_dir : str
        Directory containing the per-parameter correlation TSVs and merged TSV
    subject : str
        Subject name (e.g. sub-17, group-patient)
    fn_spec_combined : str
        Filename spec string (without parameter suffix)
    figure_info : dict
        Dictionary with figure settings from analysis_info
    rsq2use : str
        R² column to use for weighting (e.g. 'prf_rsq' or 'prf_loo_rsq')

    Returns
    -------
    fig : plotly figure object
    """

    # Axis label settings per parameter
    corr_plot_settings = {
        'prf_rsq':    {'axes': 'pRF R<sup>2</sup>',   'tick_step': 0.2},
        'prf_x':      {'axes': 'pRF x coord. (dva)',   'tick_step': 2.5},
        'prf_y':      {'axes': 'pRF y coord. (dva)',   'tick_step': 2.5},
        'prf_size':   {'axes': 'pRF size (dva)',        'tick_step': 1},
        'prf_ecc':    {'axes': 'pRF ecc. (dva)',        'tick_step': 2},
        'pcm_median': {'axes': 'pRF CM (mm/dva)',       'tick_step': 5}
    }

    # Template
    template_specs = dict(axes_color="rgba(0, 0, 0, 1)",
                          axes_width=2,
                          axes_font_size=15,
                          bg_col="rgba(255, 255, 255, 1)",
                          font='Arial',
                          title_font_size=15,
                          rois_plot_width=1.5)
    fig_template = plotly_template(template_specs)

    # Figure layout settings
    roi_colors = figure_info['roi_colors']
    fig_margin = figure_info['rois_fig_margin']
    rois = figure_info['rois']
    corr_params = figure_info['corr_params']
    corr_param_ranges = figure_info['corr_param_ranges']
    rows, cols = len(rois), len(corr_params)

    rois_hor_spacing = figure_info['rois_hor_spacing']
    rois_ver_spacing = figure_info['rois_ver_spacing']
    rois_plot_height = figure_info['rois_plot_height']
    rois_plot_width = figure_info['rois_plot_width']
    subject_group = figure_info['subject_group']

    fig_height = rois_plot_height * rows + fig_margin[1] + fig_margin[3] + (rois_ver_spacing * (rows - 1))
    fig_width = rois_plot_width * cols + fig_margin[0] + fig_margin[2] + (rois_hor_spacing * (cols - 1))
    hor_spacing = rois_hor_spacing / (fig_width - fig_margin[0] - fig_margin[2])
    ver_spacing = rois_ver_spacing / (fig_height - fig_margin[1] - fig_margin[3])

    fig = make_subplots(rows=rows, cols=cols, print_grid=False,
                        horizontal_spacing=hor_spacing,
                        vertical_spacing=ver_spacing)

    # Axis labels: x = FE/LE (reference), y = AE/RE (outcome)
    if subject_group == 'patient':
        x_label_prefix = 'FE'
        y_label_prefix = 'AE'
    else:
        x_label_prefix = 'LE'
        y_label_prefix = 'RE'

    for l, corr_param in enumerate(corr_params):

        param_range = corr_param_ranges[corr_param]
        axis_settings = corr_plot_settings[corr_param]

        # Load KDE TSV
        kde_tsv_fn = "{}/{}_{}_{}-corr.tsv".format(tsv_dir, subject, fn_spec_combined, corr_param)
        df_kde = pd.read_table(kde_tsv_fn, sep="\t")

        # Build grid coordinates
        x_grid = np.sort(df_kde['x_grid'].unique())
        y_grid = np.sort(df_kde['y_grid'].unique())

        # Load regression TSV
        reg_tsv_fn = "{}/{}_{}_{}-regression.tsv".format(tsv_dir, subject, fn_spec_combined, corr_param)
        if os.path.isfile(reg_tsv_fn):
            df_reg = pd.read_table(reg_tsv_fn, sep='\t')
        else:
            df_reg = None

        for j, roi in enumerate(rois):

            roi_color = roi_colors[roi]

            # Parse roi_color to rgba for colorscale
            rgb = roi_color.replace('rgb(', '').replace(')', '').split(',')
            r, g, b = int(rgb[0]), int(rgb[1]), int(rgb[2])

            df_roi = df_kde.loc[df_kde.roi == roi].copy()
            df_roi = df_roi.sort_values(['y_grid', 'x_grid'])
            density = df_roi['density'].values.reshape(len(y_grid), len(x_grid))

            # Normalize density to [0, 1] for colorscale
            d_max = np.nanmax(density)
            if d_max > 0:
                density_norm = density / d_max
            else:
                density_norm = density

            # Identity line
            line_xy = np.linspace(param_range[0], param_range[1], 50)
            fig.add_trace(go.Scatter(x=line_xy, y=line_xy, mode='lines',
                                     line=dict(color='rgba(0,0,0,0.5)', width=2, dash='dash'),
                                     showlegend=False),
                          row=j + 1, col=l + 1)

            # KDE contour
            fig.add_trace(go.Contour(
                x=x_grid,
                y=y_grid,
                z=density_norm,
                colorscale=[[0, f'rgba({r},{g},{b},0)'],
                            [1, f'rgba({r},{g},{b},1)']],
                showscale=False,
                contours=dict(coloring='fill',
                              showlines=False,
                              start=0.05,
                              end=1.0,
                              size=0.05),
                zmin=0, zmax=1),
                row=j + 1, col=l + 1)

            # Regression line from regression TSV
            if df_reg is not None:
                df_reg_roi = df_reg.loc[df_reg.roi == roi]
                if len(df_reg_roi) > 0:
                    slope = df_reg_roi['slope'].values[0]
                    intercept = df_reg_roi['intercept'].values[0]
                    if not np.isnan(slope):
                        reg_x = np.linspace(param_range[0], param_range[1], 50)
                        reg_y = slope * reg_x + intercept
                        fig.add_trace(go.Scatter(x=reg_x, y=reg_y, mode='lines',
                                                 line=dict(color=roi_color, width=3),
                                                 showlegend=False),
                                      row=j + 1, col=l + 1)

            # ROI label annotation
            r = param_range
            annotation = go.layout.Annotation(
                x=r[1] - 0.05 * (r[1] - r[0]),
                y=r[0] + 0.05 * (r[1] - r[0]),
                text=roi, xanchor='right', yanchor='bottom',
                showarrow=False, font_color=roi_color,
                font_family=template_specs['font'],
                font_size=template_specs['axes_font_size'])
            fig.add_annotation(annotation, row=j + 1, col=l + 1)

            # Axis labels: x only on last row, y on all rows
            x_title = f"{x_label_prefix} – {axis_settings['axes']}" if j == len(rois) - 1 else ''
            fig.update_xaxes(title_text=x_title,
                             range=param_range,
                             tickmode='linear',
                             tick0=param_range[0],
                             dtick=axis_settings['tick_step'],
                             showline=True,
                             row=j + 1, col=l + 1)
            fig.update_yaxes(title_text=f"{y_label_prefix} – {axis_settings['axes']}",
                             range=param_range,
                             tickmode='linear',
                             tick0=param_range[0],
                             dtick=axis_settings['tick_step'],
                             showline=True,
                             row=j + 1, col=l + 1)

    fig.update_layout(height=fig_height,
                      width=fig_width,
                      showlegend=False,
                      template=fig_template,
                      margin_l=fig_margin[0],
                      margin_t=fig_margin[1],
                      margin_r=fig_margin[2],
                      margin_b=fig_margin[3])

    return fig

def eyes_ecc_size_pcm_plot(df, figure_info, rsq2use):
    """
    Make scatter plots of eccentricity vs pRF size (row 1), eccentricity vs pCM (row 2),
    eccentricity vs mean R² (row 3), and eccentricity vs n_vertex (row 4),
    separately for AE-RE and FE-LE eye conditions. One column per ROI.
    FE-LE: filled square markers + solid regression line.
    AE-RE: white-filled square markers with ROI color contour + dashed regression line.
    Legend per subplot: ROI + eye condition label (e.g. V1 FE, V1 AE).
 
    Parameters
    ----------
    df : dataframe with columns roi, eye_condition, num_bins, prf_ecc_bins,
         prf_size_bins_median, prf_size_bins_ci_upper_bound, prf_size_bins_ci_lower_bound,
         prf_pcm_bins_median, prf_pcm_bins_ci_upper_bound, prf_pcm_bins_ci_lower_bound,
         {rsq2use}_bins_median, {rsq2use}_bins_ci_upper_bound, {rsq2use}_bins_ci_lower_bound,
         n_vert_bins (+ ci bounds for group)
    figure_info : dict with figure settings
    rsq2use : str, R² column to use for weighting
 
    Returns
    -------
    fig : plotly figure
    """
    from maths_utils import weighted_regression
 
    # General figure settings
    template_specs = dict(axes_color="rgba(0, 0, 0, 1)",
                          axes_width=2,
                          axes_font_size=15,
                          bg_col="rgba(255, 255, 255, 1)",
                          font='Arial',
                          title_font_size=15,
                          rois_plot_width=1.5)
    fig_template = plotly_template(template_specs)
 
    rois = figure_info['rois']
    roi_colors = figure_info['roi_colors']
    fig_margin = figure_info['rois_fig_margin']
    rois_hor_spacing = figure_info['rois_hor_spacing']
    rois_ver_spacing = figure_info['rois_ver_spacing']
    rois_plot_height = figure_info['rois_plot_height']
    rois_plot_width = figure_info['rois_plot_width']
    ecc_size_axis = figure_info['ecc_size_axis']
    ecc_size_max = figure_info['ecc_size_max']
    ecc_pcm_max = figure_info['ecc_pcm_max']
    ecc_rsq_range = figure_info['ecc_rsq_range']
    ecc_n_vert_range = figure_info['ecc_n_vert_range']
    ecc_range = figure_info['ecc_range']
    size_range = figure_info['size_range']
    pcm_range = figure_info['pmc_range']
    subject_group = figure_info['subject_group']
 
    rows, cols = 4, len(rois)
    fig_height = rois_plot_height * rows + fig_margin[1] + fig_margin[3] + (rois_ver_spacing * (rows - 1))
    fig_width = rois_plot_width * cols + fig_margin[0] + fig_margin[2] + (rois_hor_spacing * (cols - 1))
    hor_spacing = rois_hor_spacing / (fig_width - fig_margin[0] - fig_margin[2])
    ver_spacing = rois_ver_spacing / (fig_height - fig_margin[1] - fig_margin[3])
 
    fig = make_subplots(rows=rows, cols=cols, print_grid=False,
                        horizontal_spacing=hor_spacing,
                        vertical_spacing=ver_spacing)
 
    # Eye condition labels
    if subject_group == 'patient':
        eye_labels = {'AE-RE': 'AE', 'FE-LE': 'FE'}
    else:
        eye_labels = {'AE-RE': 'RE', 'FE-LE': 'LE'}
 
    eye_conditions = ['FE-LE', 'AE-RE']
 
    y_titles = {1: 'pRF size (dva)', 2: 'pRF CM (mm/dva)',
                3: 'pRF R<sup>2</sup>', 4: 'Vertices per bin'}
    y_ranges = {1: size_range, 2: pcm_range,
                3: ecc_rsq_range, 4: ecc_n_vert_range}
    x_range = ecc_range
 
    for j, roi in enumerate(rois):
        roi_color = roi_colors[roi]
        roi_color_opac = f"rgba{roi_color[3:-1]}, 0.15)"
 
        for eye_condition in eye_conditions:
            df_roi = df.loc[(df.roi == roi) & (df.eye_condition == eye_condition)]
            eye_label = eye_labels[eye_condition]
            legend_label = f'{roi} {eye_label}'
            legendgroup = f'{roi}_{eye_condition}'
 
            # FE-LE: filled marker, solid line
            # AE-RE: white-filled marker with ROI contour, dashed line
            if eye_condition == 'FE-LE':
                marker_style = dict(color=roi_color, symbol='square', size=8,
                                    line=dict(color=roi_color, width=3))
                line_dash = 'solid'
            else:
                marker_style = dict(color='white', symbol='square', size=8,
                                    line=dict(color=roi_color, width=3))
                line_dash = 'dash'
 
            ecc_median = np.array(df_roi.prf_ecc_bins)
            size_median = np.array(df_roi.prf_size_bins_median)
            r2_median = np.array(df_roi[f'{rsq2use}_bins_median'])
            size_upper = np.array(df_roi.prf_size_bins_ci_upper_bound)
            size_lower = np.array(df_roi.prf_size_bins_ci_lower_bound)
 
            # --- ROW 1: pRF size vs eccentricity ---
            valid = ~np.isnan(ecc_median) & ~np.isnan(size_median)
            if valid.sum() > 1:
                slope, intercept = weighted_regression(
                    ecc_median[valid], size_median[valid], r2_median[valid], model='linear')
                line_x = np.linspace(ecc_median[valid][0], ecc_median[valid][-1], 50)
                line = slope * line_x + intercept
 
                fig.add_trace(go.Scatter(x=line_x, y=line, mode='lines',
                                         name=legend_label, legendgroup=legendgroup,
                                         showlegend=False,
                                         line=dict(color=roi_color, width=3, dash=line_dash)),
                              row=1, col=j + 1)
 
                valid_up = ~np.isnan(size_upper)
                valid_lo = ~np.isnan(size_lower)
                slope_up, intercept_up = weighted_regression(
                    ecc_median[valid_up], size_upper[valid_up], r2_median[valid_up], model='linear') if valid_up.sum() > 1 else (slope, intercept)
                slope_lo, intercept_lo = weighted_regression(
                    ecc_median[valid_lo], size_lower[valid_lo], r2_median[valid_lo], model='linear') if valid_lo.sum() > 1 else (slope, intercept)
 
                fig.add_trace(go.Scatter(
                    x=np.concatenate([line_x, line_x[::-1]]),
                    y=np.concatenate([slope_up * line_x + intercept_up,
                                      (slope_lo * line_x + intercept_lo)[::-1]]),
                    mode='lines', fill='toself', fillcolor=roi_color_opac,
                    line=dict(color=roi_color_opac, width=0),
                    legendgroup=legendgroup, showlegend=False),
                    row=1, col=j + 1)
 
            # Markers only — no connecting line
            fig.add_trace(go.Scatter(
                x=ecc_median, y=size_median,
                mode='markers',
                error_y=dict(type='data', array=size_upper - size_median,
                             arrayminus=size_median - size_lower,
                             visible=True, thickness=3, width=0, color=roi_color),
                marker=marker_style,
                legendgroup=legendgroup, showlegend=False),
                row=1, col=j + 1)
 
            # --- ROW 2: pCM vs eccentricity ---
            pcm_median = np.array(df_roi.prf_pcm_bins_median)
            pcm_upper = np.array(df_roi.prf_pcm_bins_ci_upper_bound)
            pcm_lower = np.array(df_roi.prf_pcm_bins_ci_lower_bound)
 
            valid_pcm = ~np.isnan(ecc_median) & ~np.isnan(pcm_median)
            try:
                if valid_pcm.sum() > 1:
                    slope_pcm, intercept_pcm = weighted_regression(
                        ecc_median[valid_pcm], pcm_median[valid_pcm], r2_median[valid_pcm], model='pcm')
                    line_x_pcm = np.linspace(ecc_median[valid_pcm][0], ecc_median[valid_pcm][-1], 50)
                    line_pcm = 1 / (slope_pcm * line_x_pcm + intercept_pcm)
 
                    fig.add_trace(go.Scatter(x=line_x_pcm, y=line_pcm, mode='lines',
                                             legendgroup=legendgroup, showlegend=False,
                                             line=dict(color=roi_color, width=3, dash=line_dash)),
                                  row=2, col=j + 1)
 
                    valid_up_pcm = ~np.isnan(pcm_upper)
                    valid_lo_pcm = ~np.isnan(pcm_lower)
                    slope_up_pcm, intercept_up_pcm = weighted_regression(
                        ecc_median[valid_up_pcm], pcm_upper[valid_up_pcm], r2_median[valid_up_pcm], model='pcm') if valid_up_pcm.sum() > 1 else (slope_pcm, intercept_pcm)
                    slope_lo_pcm, intercept_lo_pcm = weighted_regression(
                        ecc_median[valid_lo_pcm], pcm_lower[valid_lo_pcm], r2_median[valid_lo_pcm], model='pcm') if valid_lo_pcm.sum() > 1 else (slope_pcm, intercept_pcm)
 
                    fig.add_trace(go.Scatter(
                        x=np.concatenate([line_x_pcm, line_x_pcm[::-1]]),
                        y=np.concatenate([1 / (slope_up_pcm * line_x_pcm + intercept_up_pcm),
                                          (1 / (slope_lo_pcm * line_x_pcm + intercept_lo_pcm))[::-1]]),
                        mode='lines', fill='toself', fillcolor=roi_color_opac,
                        line=dict(color=roi_color_opac, width=0),
                        legendgroup=legendgroup, showlegend=False),
                        row=2, col=j + 1)
 
            except RuntimeError:
                print(f"Fit failed for ROI: {roi}, eye: {eye_condition}")
 
            # Markers only — no connecting line
            fig.add_trace(go.Scatter(
                x=ecc_median, y=pcm_median,
                mode='markers',
                error_y=dict(type='data', array=pcm_upper - pcm_median,
                             arrayminus=pcm_median - pcm_lower,
                             visible=True, thickness=3, width=0, color=roi_color),
                marker=marker_style,
                legendgroup=legendgroup, showlegend=False),
                row=2, col=j + 1)
 
            # --- ROW 3: mean R² per bin --- connected dots, line width=3
            r2_median_row = np.array(df_roi[f'{rsq2use}_bins_median'])
            r2_upper = np.array(df_roi[f'{rsq2use}_bins_ci_upper_bound'])
            r2_lower = np.array(df_roi[f'{rsq2use}_bins_ci_lower_bound'])
 
            fig.add_trace(go.Scatter(
                x=ecc_median, y=r2_median_row,
                mode='markers+lines',
                error_y=dict(type='data', array=r2_upper - r2_median_row,
                             arrayminus=r2_median_row - r2_lower,
                             visible=True, thickness=3, width=0, color=roi_color),
                marker=marker_style,
                line=dict(color=roi_color, width=3, dash=line_dash),
                legendgroup=legendgroup, showlegend=False),
                row=3, col=j + 1)
 
            # --- ROW 4: n_vertex per bin --- connected dots, line width=3
            n_vert = np.array(df_roi['n_vert_bins'])
            if 'n_vert_bins_ci_upper_bound' in df_roi.columns:
                n_vert_upper = np.array(df_roi['n_vert_bins_ci_upper_bound'])
                n_vert_lower = np.array(df_roi['n_vert_bins_ci_lower_bound'])
            elif 'n_vert_bins_ci_upper' in df_roi.columns:
                n_vert_upper = np.array(df_roi['n_vert_bins_ci_upper'])
                n_vert_lower = np.array(df_roi['n_vert_bins_ci_lower'])
            else:
                n_vert_upper = n_vert
                n_vert_lower = n_vert
 
            fig.add_trace(go.Scatter(
                x=ecc_median, y=n_vert,
                mode='markers+lines',
                error_y=dict(type='data', array=n_vert_upper - n_vert,
                             arrayminus=n_vert - n_vert_lower,
                             visible=True, thickness=3, width=0, color=roi_color),
                marker=marker_style,
                line=dict(color=roi_color, width=3, dash=line_dash),
                legendgroup=legendgroup, showlegend=False),
                row=4, col=j + 1)
 
        # Legend per subplot — line with single dot in middle, left of text
        legend_x_text = x_range[1]
        legend_line_len = 0.15 * (x_range[1] - x_range[0])
        legend_gap = 0.45 * (x_range[1] - x_range[0])
        legend_x_right = legend_x_text - legend_gap
        legend_x_left = legend_x_right - legend_line_len
        legend_x_mid = (legend_x_left + legend_x_right) / 2
        legend_y_start = y_ranges[1][1] - 0.05 * (y_ranges[1][1] - y_ranges[1][0])
        legend_y_step = 0.12 * (y_ranges[1][1] - y_ranges[1][0])
 
        for k, eye_condition in enumerate(eye_conditions):
            eye_label = eye_labels[eye_condition]
            legend_y = legend_y_start - k * legend_y_step
            legend_label = f'{eye_label} - {roi}'
 
            if eye_condition == 'FE-LE':
                leg_marker = dict(color=roi_color, symbol='square', size=8,
                                  line=dict(color=roi_color, width=3))
                leg_line_dash = 'solid'
            else:
                leg_marker = dict(color='white', symbol='square', size=8,
                                  line=dict(color=roi_color, width=3))
                leg_line_dash = 'dash'
 
            # Line with single dot in the middle
            fig.add_trace(go.Scatter(
                x=[legend_x_left, legend_x_mid, legend_x_right],
                y=[legend_y, legend_y, legend_y],
                mode='lines+markers',
                marker={**leg_marker, **{'size': [0, 8, 0]}},
                line=dict(color=roi_color, width=3, dash=leg_line_dash),
                showlegend=False),
                row=1, col=j + 1)
 
            fig.add_annotation(
                x=legend_x_text,
                y=legend_y,
                text=legend_label,
                xanchor='right', yanchor='middle',
                showarrow=False,
                font_color=roi_color,
                font_family=template_specs['font'],
                font_size=template_specs['axes_font_size'],
                row=1, col=j + 1)
 
        # Axes
        for row_idx in range(1, rows + 1):
            x_title = 'pRF eccentricity (dva)' if row_idx == rows else ''
            fig.update_xaxes(title_text=x_title,
                             range=x_range,
                             showline=True, row=row_idx, col=j + 1)
            y_title = y_titles[row_idx] if j == 0 else ''
            fig.update_yaxes(title_text=y_title,
                             range=y_ranges[row_idx],
                             showline=True, nticks=6,
                             row=row_idx, col=j + 1)
            # For row 4 (vertices), format ticks in k
            if row_idx == 4:
                tick_vals = list(range(ecc_n_vert_range[0], ecc_n_vert_range[1] + 1,
                                       ecc_n_vert_range[1] // 4))
                tick_text = [f'{v//1000}k' if v >= 1000 else str(v) for v in tick_vals]
                fig.update_yaxes(tickvals=tick_vals, ticktext=tick_text,
                                 row=row_idx, col=j + 1)
 
    fig.update_layout(height=fig_height,
                      width=fig_width,
                      showlegend=False,
                      template=fig_template,
                      margin_l=fig_margin[0],
                      margin_t=fig_margin[1],
                      margin_r=fig_margin[2],
                      margin_b=fig_margin[3])
 
    return fig
    
def eyes_active_vert_plot(df, figure_info, format_):
    """
    Make grouped bar plots of significant vertices per ROI for AE-RE and FE-LE eye conditions.
    Two subplots: FDR 0.05 and FDR 0.01.
    AE-RE: solid bars. FE-LE: x-pattern bars.
    Legend: two gray bar entries inside each subplot top-right.

    Parameters
    ----------
    df : dataframe with columns roi, eye_condition, n_vert_tot,
         n_vert_corr_pvalue_5pt, n_vert_corr_pvalue_1pt, ratio_5pt, ratio_1pt
    figure_info : dict with figure settings
    format_ : str, data format (e.g. 'fsnative', '170k')

    Returns
    -------
    fig : grouped bar plot
    """

    # General figure settings
    template_specs = dict(axes_color="rgba(0, 0, 0, 1)",
                          axes_width=2,
                          axes_font_size=15,
                          bg_col="rgba(255, 255, 255, 1)",
                          font='Arial',
                          title_font_size=15,
                          rois_plot_width=1.5)
    fig_template = plotly_template(template_specs)

    rois = figure_info['rois']
    roi_colors = figure_info['roi_colors']
    fig_margin = figure_info['rois_fig_margin']
    rois_hor_spacing = figure_info['rois_hor_spacing']
    rois_ver_spacing = figure_info['rois_ver_spacing']
    bar_width = figure_info['rois_bar_width']
    rois_plot_height = figure_info['rois_plot_height']
    subject_group = figure_info['subject_group']

    rows, cols = 1, 2
    fig_height = rois_plot_height * rows + fig_margin[1] + fig_margin[3] + (rois_ver_spacing * (rows - 1))
    fig_width = bar_width * cols * len(rois) + fig_margin[0] + fig_margin[2] + (rois_hor_spacing * (cols - 1))
    hor_spacing = rois_hor_spacing / (fig_width - fig_margin[0] - fig_margin[2])
    ver_spacing = rois_ver_spacing / (fig_height - fig_margin[1] - fig_margin[3])

    # Eye condition labels depend on subject group
    if subject_group == 'patient':
        eye_labels = {'AE-RE': 'AE', 'FE-LE': 'FE'}
    else:
        eye_labels = {'AE-RE': 'RE', 'FE-LE': 'LE'}

    # Pattern per eye condition
    eye_patterns = {'FE-LE': '', 'AE-RE': '/'}
    eye_conditions = ['FE-LE', 'AE-RE']

    fig = make_subplots(rows=rows, cols=cols,
                        subplot_titles=['FDR threshold = 0.05', 'FDR threshold = 0.01'],
                        vertical_spacing=ver_spacing,
                        horizontal_spacing=hor_spacing)

    if format_ == 'fsnative':
        range_val = figure_info['active_vert_fsnative_range']
    elif format_ == '170k':
        range_val = figure_info['active_vert_170k_range']

    for col_idx, (fdr_label, ratio_col, sig_col) in enumerate([
            ('FDR 0.05', 'ratio_5pt', 'n_vert_corr_pvalue_5pt'),
            ('FDR 0.01', 'ratio_1pt', 'n_vert_corr_pvalue_1pt')], start=1):

        legend_ref = 'legend' if col_idx == 1 else 'legend2'

        for eye_condition in eye_conditions:
            df_eye = df.loc[df.eye_condition == eye_condition]
            eye_label = eye_labels[eye_condition]
            pattern_shape = eye_patterns[eye_condition]
            bar_colors = [roi_colors[roi] for roi in df_eye['roi']]

            # Total vertices — transparent background
            fig.add_trace(go.Bar(
                x=df_eye['roi'],
                y=df_eye['n_vert_tot'],
                name=eye_label,
                text=(df_eye[ratio_col] * 100).astype(int).astype(str) + '%',
                textposition='outside',
                textangle=-90,
                textfont=dict(size=15),
                offsetgroup=eye_condition,
                marker=dict(color=bar_colors,
                            opacity=0.15),
                showlegend=False,
                legend=legend_ref),
                row=1, col=col_idx)

            # Significant vertices — solid/patterned bar
            fig.add_trace(go.Bar(
                x=df_eye['roi'],
                y=df_eye[sig_col],
                name=eye_label,
                offsetgroup=eye_condition,
                marker=dict(color=bar_colors,
                            pattern=dict(shape=pattern_shape,
                                         fgcolor='rgba(0,0,0,1)',
                                         bgcolor='rgba(0,0,0,0)',
                                         size=8)),
                showlegend=False,
                legend=legend_ref),
                row=1, col=col_idx)

        # Dummy bar traces for legend — one per subplot using legend/legend2
        for eye_condition in eye_conditions:
            eye_label = eye_labels[eye_condition]
            pattern_shape = eye_patterns[eye_condition]
            fig.add_trace(go.Bar(
                x=[None], y=[None],
                name=eye_label,
                marker=dict(color='rgba(128,128,128,1)',
                            pattern=dict(shape=pattern_shape,
                                         fgcolor='rgba(0,0,0,1)',
                                         bgcolor='rgba(0,0,0,0.3)',
                                         size=8)),
                showlegend=True,
                legend=legend_ref),
                row=1, col=col_idx)

    fig.update_xaxes(showline=True, ticklen=0, linecolor='rgba(255,255,255,0)')
    fig.update_yaxes(range=range_val, showline=True, nticks=10,
                     title_text='Number of vertex')

    legend_style = dict(orientation='h',
                        font_family=template_specs['font'],
                        font_size=template_specs['axes_font_size'],
                        x=0.44, y=0.95,
                        xanchor='right',
                        yanchor='top',
                        bgcolor='rgba(0,0,0,0)',
                        borderwidth=0,
                        traceorder='normal',
                        itemwidth=30)

    legend2_style = {**legend_style, 'x': 1}

    fig.update_layout(barmode='group',
                      bargap=0.3,
                      bargroupgap=0.05,
                      height=fig_height,
                      width=fig_width,
                      template=fig_template,
                      legend=legend_style,
                      legend2=legend2_style,
                      margin_l=fig_margin[0],
                      margin_t=fig_margin[1],
                      margin_r=fig_margin[2],
                      margin_b=fig_margin[3])

    return fig
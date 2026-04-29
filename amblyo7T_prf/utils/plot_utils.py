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

def eyes_ecc_size_pcm_plot(df, figure_info, rsq2use, df_ctrl=None):
    """
    Make scatter plots of eccentricity vs pRF size (row 1), eccentricity vs pCM (row 2),
    eccentricity vs mean R² (row 3), and eccentricity vs n_vertex (row 4),
    separately for AE-RE and FE-LE eye conditions. One column per ROI.
    FE-LE: filled square markers + solid regression line.
    AE-RE: white-filled square markers with ROI color contour + dashed regression line.

    For control subjects (individual and group-control): a CTRL condition is also plotted
    from df itself (eye_condition == 'CTRL'). Filled gray square, solid line, gray shade
    goes from light (first ROI) to dark (last ROI).

    For group-patient: CTRL is plotted from df_ctrl (loaded from group-control TSV).
    Same gray style. Patient individual subjects: no CTRL plotted.

    Legend per subplot: ROI + eye condition label (e.g. V1 FE, V1 AE, V1 CTRL).

    Parameters
    ----------
    df : dataframe
    figure_info : dict with figure settings
    rsq2use : str, R² column to use for weighting
    df_ctrl : dataframe or None

    Returns
    -------
    fig : plotly figure
    """
    from maths_utils import weighted_regression

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

    if subject_group == 'patient':
        eye_labels = {'AE-RE': 'AE', 'FE-LE': 'FE'}
    else:
        eye_labels = {'AE-RE': 'RE', 'FE-LE': 'LE'}

    eye_conditions = ['FE-LE', 'AE-RE']

    plot_ctrl = False
    df_ctrl_source = None
    if subject_group == 'control':
        if 'CTRL' in df['eye_condition'].values:
            plot_ctrl = True
            df_ctrl_source = df
    elif subject_group == 'patient' and df_ctrl is not None:
        plot_ctrl = True
        df_ctrl_source = df_ctrl

    y_titles = {1: 'pRF size (dva)', 2: 'pRF CM (mm/dva)',
                3: 'pRF R<sup>2</sup>', 4: 'Vertices per bin'}
    y_ranges = {1: size_range, 2: pcm_range, 3: ecc_rsq_range}
    x_range = ecc_range

    ctrl_grays = ['rgb(210,210,210)', 'rgb(180,180,180)', 'rgb(150,150,150)',
                  'rgb(120,120,120)', 'rgb(90,90,90)', 'rgb(60,60,60)',
                  'rgb(40,40,40)']

    for j, roi in enumerate(rois):
        roi_color = roi_colors[roi]
        roi_color_opac = f"rgba{roi_color[3:-1]}, 0.15)"

        ctrl_color = ctrl_grays[j % len(ctrl_grays)]
        ctrl_color_opac = ctrl_color.replace('rgb(', 'rgba(').replace(')', ', 0.15)')

        # --- CTRL traces ---
        if plot_ctrl:
            legendgroup_ctrl = f'{roi}_CTRL'
            df_roi_ctrl = df_ctrl_source.loc[
                (df_ctrl_source.roi == roi) & (df_ctrl_source.eye_condition == 'CTRL')]

            ctrl_marker_style = dict(color=ctrl_color, symbol='square', size=8,
                                     line=dict(color=ctrl_color, width=3))

            ecc_ctrl = np.array(df_roi_ctrl.prf_ecc_bins)
            r2_ctrl = np.array(df_roi_ctrl[f'{rsq2use}_bins_median'])

            # ROW 1: pRF size
            size_ctrl = np.array(df_roi_ctrl.prf_size_bins_median)
            size_ctrl_upper = np.array(df_roi_ctrl.prf_size_bins_ci_upper_bound)
            size_ctrl_lower = np.array(df_roi_ctrl.prf_size_bins_ci_lower_bound)

            valid_ctrl = ~np.isnan(ecc_ctrl) & ~np.isnan(size_ctrl)
            if valid_ctrl.sum() > 1:
                slope_c, intercept_c = weighted_regression(
                    ecc_ctrl[valid_ctrl], size_ctrl[valid_ctrl], r2_ctrl[valid_ctrl], model='linear')
                line_x_c = np.linspace(ecc_ctrl[valid_ctrl][0], ecc_ctrl[valid_ctrl][-1], 50)
                line_c = slope_c * line_x_c + intercept_c

                fig.add_trace(go.Scatter(x=line_x_c, y=line_c, mode='lines',
                                         legendgroup=legendgroup_ctrl, showlegend=False,
                                         line=dict(color=ctrl_color, width=3, dash='solid')),
                              row=1, col=j + 1)

                valid_up_c = ~np.isnan(size_ctrl_upper)
                valid_lo_c = ~np.isnan(size_ctrl_lower)
                slope_up_c, intercept_up_c = weighted_regression(
                    ecc_ctrl[valid_up_c], size_ctrl_upper[valid_up_c], r2_ctrl[valid_up_c], model='linear') if valid_up_c.sum() > 1 else (slope_c, intercept_c)
                slope_lo_c, intercept_lo_c = weighted_regression(
                    ecc_ctrl[valid_lo_c], size_ctrl_lower[valid_lo_c], r2_ctrl[valid_lo_c], model='linear') if valid_lo_c.sum() > 1 else (slope_c, intercept_c)

                fig.add_trace(go.Scatter(
                    x=np.concatenate([line_x_c, line_x_c[::-1]]),
                    y=np.concatenate([slope_up_c * line_x_c + intercept_up_c,
                                      (slope_lo_c * line_x_c + intercept_lo_c)[::-1]]),
                    mode='lines', fill='toself', fillcolor=ctrl_color_opac,
                    line=dict(color=ctrl_color_opac, width=0),
                    legendgroup=legendgroup_ctrl, showlegend=False),
                    row=1, col=j + 1)

            fig.add_trace(go.Scatter(
                x=ecc_ctrl, y=size_ctrl,
                mode='markers',
                error_y=dict(type='data', array=size_ctrl_upper - size_ctrl,
                             arrayminus=size_ctrl - size_ctrl_lower,
                             visible=True, thickness=3, width=0, color=ctrl_color),
                marker=ctrl_marker_style,
                legendgroup=legendgroup_ctrl, showlegend=False),
                row=1, col=j + 1)

            # ROW 2: pCM
            pcm_ctrl = np.array(df_roi_ctrl.prf_pcm_bins_median)
            pcm_ctrl_upper = np.array(df_roi_ctrl.prf_pcm_bins_ci_upper_bound)
            pcm_ctrl_lower = np.array(df_roi_ctrl.prf_pcm_bins_ci_lower_bound)

            valid_pcm_ctrl = ~np.isnan(ecc_ctrl) & ~np.isnan(pcm_ctrl)
            try:
                if valid_pcm_ctrl.sum() > 1:
                    slope_pcm_c, intercept_pcm_c = weighted_regression(
                        ecc_ctrl[valid_pcm_ctrl], pcm_ctrl[valid_pcm_ctrl],
                        r2_ctrl[valid_pcm_ctrl], model='pcm')
                    line_x_pcm_c = np.linspace(ecc_ctrl[valid_pcm_ctrl][0], ecc_ctrl[valid_pcm_ctrl][-1], 50)
                    line_pcm_c = 1 / (slope_pcm_c * line_x_pcm_c + intercept_pcm_c)

                    fig.add_trace(go.Scatter(x=line_x_pcm_c, y=line_pcm_c, mode='lines',
                                             legendgroup=legendgroup_ctrl, showlegend=False,
                                             line=dict(color=ctrl_color, width=3, dash='solid')),
                                  row=2, col=j + 1)

                    valid_up_pcm_c = ~np.isnan(pcm_ctrl_upper)
                    valid_lo_pcm_c = ~np.isnan(pcm_ctrl_lower)
                    slope_up_pcm_c, intercept_up_pcm_c = weighted_regression(
                        ecc_ctrl[valid_up_pcm_c], pcm_ctrl_upper[valid_up_pcm_c],
                        r2_ctrl[valid_up_pcm_c], model='pcm') if valid_up_pcm_c.sum() > 1 else (slope_pcm_c, intercept_pcm_c)
                    slope_lo_pcm_c, intercept_lo_pcm_c = weighted_regression(
                        ecc_ctrl[valid_lo_pcm_c], pcm_ctrl_lower[valid_lo_pcm_c],
                        r2_ctrl[valid_lo_pcm_c], model='pcm') if valid_lo_pcm_c.sum() > 1 else (slope_pcm_c, intercept_pcm_c)

                    fig.add_trace(go.Scatter(
                        x=np.concatenate([line_x_pcm_c, line_x_pcm_c[::-1]]),
                        y=np.concatenate([1 / (slope_up_pcm_c * line_x_pcm_c + intercept_up_pcm_c),
                                          (1 / (slope_lo_pcm_c * line_x_pcm_c + intercept_lo_pcm_c))[::-1]]),
                        mode='lines', fill='toself', fillcolor=ctrl_color_opac,
                        line=dict(color=ctrl_color_opac, width=0),
                        legendgroup=legendgroup_ctrl, showlegend=False),
                        row=2, col=j + 1)

            except RuntimeError:
                print(f"Fit failed for ROI: {roi}, CTRL")

            fig.add_trace(go.Scatter(
                x=ecc_ctrl, y=pcm_ctrl,
                mode='markers',
                error_y=dict(type='data', array=pcm_ctrl_upper - pcm_ctrl,
                             arrayminus=pcm_ctrl - pcm_ctrl_lower,
                             visible=True, thickness=3, width=0, color=ctrl_color),
                marker=ctrl_marker_style,
                legendgroup=legendgroup_ctrl, showlegend=False),
                row=2, col=j + 1)

            # ROW 3: R²
            r2_ctrl_row = np.array(df_roi_ctrl[f'{rsq2use}_bins_median'])
            r2_ctrl_upper = np.array(df_roi_ctrl[f'{rsq2use}_bins_ci_upper_bound'])
            r2_ctrl_lower = np.array(df_roi_ctrl[f'{rsq2use}_bins_ci_lower_bound'])

            fig.add_trace(go.Scatter(
                x=ecc_ctrl, y=r2_ctrl_row,
                mode='markers+lines',
                error_y=dict(type='data', array=r2_ctrl_upper - r2_ctrl_row,
                             arrayminus=r2_ctrl_row - r2_ctrl_lower,
                             visible=True, thickness=3, width=0, color=ctrl_color),
                marker=ctrl_marker_style,
                line=dict(color=ctrl_color, width=3, dash='solid'),
                legendgroup=legendgroup_ctrl, showlegend=False),
                row=3, col=j + 1)

            # ROW 4: n_vertex
            n_vert_ctrl = np.array(df_roi_ctrl['n_vert_bins'])
            if 'n_vert_bins_ci_upper_bound' in df_roi_ctrl.columns:
                n_vert_ctrl_upper = np.array(df_roi_ctrl['n_vert_bins_ci_upper_bound'])
                n_vert_ctrl_lower = np.array(df_roi_ctrl['n_vert_bins_ci_lower_bound'])
            elif 'n_vert_bins_ci_upper' in df_roi_ctrl.columns:
                n_vert_ctrl_upper = np.array(df_roi_ctrl['n_vert_bins_ci_upper'])
                n_vert_ctrl_lower = np.array(df_roi_ctrl['n_vert_bins_ci_lower'])
            else:
                n_vert_ctrl_upper = n_vert_ctrl
                n_vert_ctrl_lower = n_vert_ctrl

            fig.add_trace(go.Scatter(
                x=ecc_ctrl, y=n_vert_ctrl,
                mode='markers+lines',
                error_y=dict(type='data', array=n_vert_ctrl_upper - n_vert_ctrl,
                             arrayminus=n_vert_ctrl - n_vert_ctrl_lower,
                             visible=True, thickness=3, width=0, color=ctrl_color),
                marker=ctrl_marker_style,
                line=dict(color=ctrl_color, width=3, dash='solid'),
                legendgroup=legendgroup_ctrl, showlegend=False),
                row=4, col=j + 1)

        for eye_condition in eye_conditions:
            df_roi = df.loc[(df.roi == roi) & (df.eye_condition == eye_condition)]
            eye_label = eye_labels[eye_condition]
            legend_label = f'{roi} {eye_label}'
            legendgroup = f'{roi}_{eye_condition}'

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

            fig.add_trace(go.Scatter(
                x=ecc_median, y=pcm_median,
                mode='markers',
                error_y=dict(type='data', array=pcm_upper - pcm_median,
                             arrayminus=pcm_median - pcm_lower,
                             visible=True, thickness=3, width=0, color=roi_color),
                marker=marker_style,
                legendgroup=legendgroup, showlegend=False),
                row=2, col=j + 1)

            # --- ROW 3: mean R² per bin ---
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

            # --- ROW 4: n_vertex per bin ---
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

        # --- Legend per subplot (in row 1, size subplot) ---
        legend_conditions = list(eye_conditions)
        if plot_ctrl:
            legend_conditions.append('CTRL')

        legend_x_text = x_range[1]
        legend_line_len = 0.15 * (x_range[1] - x_range[0])
        legend_gap = 0.55 * (x_range[1] - x_range[0])
        legend_x_right = legend_x_text - legend_gap
        legend_x_left = legend_x_right - legend_line_len
        legend_x_mid = (legend_x_left + legend_x_right) / 2
        legend_y_start = size_range[1] - 0.05 * (size_range[1] - size_range[0])
        legend_y_step = 0.12 * (size_range[1] - size_range[0])

        for k, cond in enumerate(legend_conditions):
            legend_y = legend_y_start - k * legend_y_step

            if cond == 'CTRL':
                eye_label = 'CTRL'
                leg_marker = dict(color=ctrl_color, symbol='square', size=8,
                                  line=dict(color=ctrl_color, width=3))
                leg_line_color = ctrl_color
                leg_line_dash = 'solid'
            elif cond == 'FE-LE':
                eye_label = eye_labels['FE-LE']
                leg_marker = dict(color=roi_color, symbol='square', size=8,
                                  line=dict(color=roi_color, width=3))
                leg_line_color = roi_color
                leg_line_dash = 'solid'
            else:
                eye_label = eye_labels['AE-RE']
                leg_marker = dict(color='white', symbol='square', size=8,
                                  line=dict(color=roi_color, width=3))
                leg_line_color = roi_color
                leg_line_dash = 'dash'

            legend_label = f'{eye_label} - {roi}'

            fig.add_trace(go.Scatter(
                x=[legend_x_left, legend_x_right],
                y=[legend_y, legend_y],
                mode='lines',
                line=dict(color=leg_line_color, width=3, dash=leg_line_dash),
                showlegend=False),
                row=1, col=j + 1)

            fig.add_trace(go.Scatter(
                x=[legend_x_mid],
                y=[legend_y],
                mode='markers',
                marker=leg_marker,
                showlegend=False),
                row=1, col=j + 1)

            fig.add_annotation(
                x=legend_x_text,
                y=legend_y,
                text=legend_label,
                xanchor='right', yanchor='middle',
                showarrow=False,
                font_color=leg_line_color,
                font_family=template_specs['font'],
                font_size=template_specs['axes_font_size'],
                row=1, col=j + 1)

        # --- Axes rows 1-3: fixed ranges ---
        for row_idx in range(1, 4):
            fig.update_xaxes(title_text='',
                             range=x_range,
                             showline=True, row=row_idx, col=j + 1)
            y_title = y_titles[row_idx] if j == 0 else ''
            fig.update_yaxes(title_text=y_title,
                             range=y_ranges[row_idx],
                             showline=True, nticks=6,
                             row=row_idx, col=j + 1)

        # --- Row 4: x-axis with title, y-axis fully auto (Plotly decides range and ticks) ---
        fig.update_xaxes(title_text='pRF eccentricity (dva)',
                         range=x_range,
                         showline=True, row=4, col=j + 1)

        y_title_row4 = y_titles[4] if j == 0 else ''
        fig.update_yaxes(title_text=y_title_row4,
                         showline=True,
                         row=4, col=j + 1)

    fig.update_layout(height=fig_height,
                      width=fig_width,
                      showlegend=False,
                      template=fig_template,
                      margin_l=fig_margin[0],
                      margin_t=fig_margin[1],
                      margin_r=fig_margin[2],
                      margin_b=fig_margin[3])

    return fig

def ecc_comp_plot(df, df_stats, figure_info, eye_conditions, show_stats=True):
    """
    Plot pRF parameters per ROI grouped by eccentricity category (all / foveal / peripheral).
    5 rows (params) × 3 columns (ecc categories).
    x-axis: categorical ROI labels, axis line and ticks hidden.

    Marker style:
      FE / LE : filled square, ROI color
      AE / RE : white square with ROI color border
      CTRL    : filled gray square

    Legend: top-right of first subplot (row 1 col 1).
             7 ROI-colored filled square markers side by side, no lines.

    Error bars: vertex-level 2.5/97.5 CI (individual) or subject-level CI (group).
                n_vert: no error bar for individual, Plotly auto-scales y-axis.

    Significance lines (group only): horizontal brackets per ROI at ~90% of y-axis max.

    Parameters
    ----------
    df           : TSV dataframe
    df_stats     : stats TSV (None for individual subjects)
    figure_info  : dict with figure settings
    eye_conditions : list of eye_condition values to plot
    show_stats   : bool

    Returns
    -------
    fig : plotly figure
    """
    # ---------------------------------------------------------------------------------
    # Settings
    # ---------------------------------------------------------------------------------
    template_specs = dict(axes_color="rgba(0, 0, 0, 1)",
                          axes_width=2,
                          axes_font_size=15,
                          bg_col="rgba(255, 255, 255, 1)",
                          font='Arial',
                          title_font_size=15,
                          rois_plot_width=1.5)
    fig_template = plotly_template(template_specs)

    rois          = figure_info['rois']
    roi_colors    = figure_info['roi_colors']
    subject_group = figure_info['subject_group']
    fig_width     = figure_info['ecc_comp_fig_width']
    fig_height    = figure_info['ecc_comp_fig_height']
    marker_size   = figure_info['ecc_comp_marker_size']
    roi_spacing   = figure_info['ecc_comp_roi_x_spacing']
    marker_offsets = figure_info['ecc_comp_marker_offset']

    # Enforce fixed plotting order: FE-LE → AE-RE → CTRL
    fixed_order    = ['FE-LE', 'AE-RE', 'CTRL']
    eye_conditions = [e for e in fixed_order if e in eye_conditions]

    # Select offset based on number of conditions
    marker_offset = marker_offsets[0] if len(eye_conditions) <= 2 else marker_offsets[1]

    # Eye condition labels
    if subject_group == 'patient':
        eye_labels = {'AE-RE': 'AE', 'FE-LE': 'FE', 'CTRL': 'CTRL'}
    else:
        eye_labels = {'AE-RE': 'RE', 'FE-LE': 'LE', 'CTRL': 'CTRL'}

    ctrl_grays = ['rgb(210,210,210)', 'rgb(180,180,180)', 'rgb(150,150,150)',
                  'rgb(120,120,120)', 'rgb(90,90,90)', 'rgb(60,60,60)',
                  'rgb(40,40,40)']

    def get_marker(eye_cond, roi_color, roi_idx=0):
        if eye_cond == 'FE-LE':
            return dict(color=roi_color, symbol='square', size=marker_size,
                        line=dict(color=roi_color, width=3))
        elif eye_cond == 'AE-RE':
            return dict(color='white', symbol='square', size=marker_size,
                        line=dict(color=roi_color, width=3))
        else:
            ctrl_color = ctrl_grays[roi_idx % len(ctrl_grays)]
            return dict(color=ctrl_color, symbol='square', size=marker_size,
                        line=dict(color=ctrl_color, width=3))

    def get_err_color(eye_cond, roi_color, roi_idx=0):
        if eye_cond == 'CTRL':
            return ctrl_grays[roi_idx % len(ctrl_grays)]
        return roi_color

    # ---------------------------------------------------------------------------------
    # Parameters and y-axis settings
    # Per param × per ecc_category range from figure_info
    # n_vert: no range (Plotly auto-scales)
    # ---------------------------------------------------------------------------------
    ecc_categories = ['all', 'foveal', 'peripheral']
    col_titles     = ['All pRF', 'Foveal pRF (ecc ≤ 2.5 dva)', 'Peripheral pRF (ecc > 2.5 dva)']

    params_settings = [
        {'param': 'prf_rsq',
         'y_title': 'pRF R²',
         'y_range': {'all':        figure_info['ecc_comp_rsq_all_range'],
                     'foveal':     figure_info['ecc_comp_rsq_foveal_range'],
                     'peripheral': figure_info['ecc_comp_rsq_peripheral_range']}},
        {'param': 'prf_size',
         'y_title': 'pRF size (dva)',
         'y_range': {'all':        figure_info['ecc_comp_size_all_range'],
                     'foveal':     figure_info['ecc_comp_size_foveal_range'],
                     'peripheral': figure_info['ecc_comp_size_peripheral_range']}},
        {'param': 'prf_ecc',
         'y_title': 'pRF ecc (dva)',
         'y_range': {'all':        figure_info['ecc_comp_ecc_all_range'],
                     'foveal':     figure_info['ecc_comp_ecc_foveal_range'],
                     'peripheral': figure_info['ecc_comp_ecc_peripheral_range']}},
        {'param': 'pcm_median',
         'y_title': 'pRF CM (mm/dva)',
         'y_range': {'all':        figure_info['ecc_comp_pcm_all_range'],
                     'foveal':     figure_info['ecc_comp_pcm_foveal_range'],
                     'peripheral': figure_info['ecc_comp_pcm_peripheral_range']}},
        {'param': 'n_vert',
         'y_title': 'N vertices',
         'y_range': None},  # Plotly auto-scales
    ]

    rows, cols = len(params_settings), len(ecc_categories)

    fig = make_subplots(rows=rows, cols=cols,
                        print_grid=False,
                        vertical_spacing=0.07,
                        horizontal_spacing=0.06)

    # Column titles — centered on subplot domains, not bold
    # Build after fig is created so domains are available
    col_domain_centers = []
    for c_idx in range(1, cols + 1):
        axis_key = f'xaxis{c_idx}' if c_idx > 1 else 'xaxis'
        domain   = fig.layout[axis_key].domain
        if domain is None or len(domain) == 0:
            # fallback to equal spacing
            col_domain_centers.append((c_idx - 0.5) / cols)
        else:
            col_domain_centers.append((domain[0] + domain[1]) / 2)

    for c_idx, col_title in enumerate(col_titles):
        fig.add_annotation(
            x=col_domain_centers[c_idx],
            y=1.02,
            xref='paper', yref='paper',
            text=col_title,
            showarrow=False,
            font=dict(size=template_specs['title_font_size'],
                      family=template_specs['font']),
            xanchor='center', yanchor='bottom')

    # ---------------------------------------------------------------------------------
    # Fake continuous x-axis
    # ---------------------------------------------------------------------------------
    n_rois      = len(rois)
    roi_centers = {roi: i * roi_spacing for i, roi in enumerate(rois)}
    # cond_offsets: FE-LE leftmost, AE-RE middle, CTRL rightmost
    n_conds = len(eye_conditions)
    if n_conds == 1:
        cond_offsets = {eye_conditions[0]: 0.0}
    elif n_conds == 2:
        cond_offsets = {eye_conditions[0]: -marker_offset,
                        eye_conditions[1]:  marker_offset}
    else:
        cond_offsets = {eye_conditions[0]: -marker_offset,
                        eye_conditions[1]:  0.0,
                        eye_conditions[2]:  marker_offset}

    x_min     = -roi_spacing * 0.5
    x_max     = (n_rois - 1) * roi_spacing + roi_spacing * 0.5
    tick_vals = [roi_centers[roi] for roi in rois]

    # Significance comparisons ordered top→bottom, matching marker positions
    if len(eye_conditions) == 2:
        sig_comparisons = [('FE-LE', 'AE-RE')]
    else:
        sig_comparisons = [
            ('FE-LE', 'CTRL'),
            ('AE-RE', 'CTRL'),
            ('FE-LE', 'AE-RE'),
        ]

    # ---------------------------------------------------------------------------------
    # Plot
    # ---------------------------------------------------------------------------------
    for row_idx, ps in enumerate(params_settings, 1):
        param   = ps['param']
        y_title = ps['y_title']
        col_med = 'n_vert_median' if param == 'n_vert' else f'{param}_median'
        col_lo  = 'n_vert_ci_lo'  if param == 'n_vert' else f'{param}_ci_lo'
        col_hi  = 'n_vert_ci_hi'  if param == 'n_vert' else f'{param}_ci_hi'

        for col_idx, ecc_cat in enumerate(ecc_categories, 1):
            y_range  = ps['y_range'][ecc_cat] if ps['y_range'] is not None else None
            y_max    = y_range[1] if y_range is not None else None
            sig_y_top  = y_max * 0.92  if y_max is not None else None
            sig_y_step = y_max * 0.06  if y_max is not None else None

            for j, roi in enumerate(rois):
                roi_color = roi_colors[roi]
                x_center  = roi_centers[roi]

                for eye_cond in eye_conditions:
                    x_pos = x_center + cond_offsets[eye_cond]

                    df_row = df.loc[
                        (df.roi == roi) &
                        (df.ecc_category == ecc_cat) &
                        (df.eye_condition == eye_cond)]

                    if len(df_row) == 0:
                        continue

                    y_val = df_row[col_med].values[0]
                    y_lo  = df_row[col_lo].values[0] if col_lo in df_row.columns else np.nan
                    y_hi  = df_row[col_hi].values[0] if col_hi in df_row.columns else np.nan
                    has_ci = not (np.isnan(y_lo) or np.isnan(y_hi))

                    marker  = get_marker(eye_cond, roi_color, roi_idx=j)
                    err_col = get_err_color(eye_cond, roi_color, roi_idx=j)

                    fig.add_trace(go.Scatter(
                        x=[x_pos], y=[y_val],
                        mode='markers',
                        marker=marker,
                        error_y=dict(
                            type='data',
                            array=[y_hi - y_val] if has_ci else [0],
                            arrayminus=[y_val - y_lo] if has_ci else [0],
                            visible=has_ci,
                            thickness=2, width=0,
                            color=err_col),
                        showlegend=False),
                        row=row_idx, col=col_idx)

                # Significance lines
                if show_stats and df_stats is not None and sig_y_top is not None:
                    for sig_lvl, (cond_a, cond_b) in enumerate(sig_comparisons):
                        if cond_a not in eye_conditions or cond_b not in eye_conditions:
                            continue

                        stat_row = df_stats.loc[
                            (df_stats.roi == roi) &
                            (df_stats.ecc_category == ecc_cat) &
                            (df_stats.param == param) &
                            (df_stats.cond_A == cond_a) &
                            (df_stats.cond_B == cond_b)]

                        if len(stat_row) == 0:
                            continue

                        stars  = stat_row['stars'].values[0]
                        x_a    = x_center + cond_offsets[cond_a]
                        x_b    = x_center + cond_offsets[cond_b]
                        sig_y  = sig_y_top - sig_lvl * sig_y_step
                        line_color = 'rgba(0,0,0,0.7)'

                        fig.add_trace(go.Scatter(
                            x=[x_a, x_a, x_b, x_b],
                            y=[sig_y - sig_y_step * 0.15,
                               sig_y, sig_y,
                               sig_y - sig_y_step * 0.15],
                            mode='lines',
                            line=dict(color=line_color, width=1.5),
                            showlegend=False),
                            row=row_idx, col=col_idx)

                        fig.add_annotation(
                            x=(x_a + x_b) / 2,
                            y=sig_y + sig_y_step * 0.05,
                            text=stars,
                            showarrow=False,
                            font=dict(size=13 if stars != 'ns' else 11,
                                      family=template_specs['font'],
                                      color=line_color),
                            xanchor='center', yanchor='bottom',
                            row=row_idx, col=col_idx)

            # Axes
            x_title      = 'ROI' if row_idx == rows else ''
            y_title_show = y_title if col_idx == 1 else ''

            # x-axis: hide line and ticks, keep ROI labels
            fig.update_xaxes(
                range=[x_min, x_max],
                tickvals=tick_vals,
                ticktext=rois,
                tickangle=0,
                ticklen=0,
                linecolor='rgba(255,255,255,0)',
                title_text=x_title,
                row=row_idx, col=col_idx)

            # y-axis
            if y_range is not None:
                fig.update_yaxes(range=y_range, showline=True, nticks=6,
                                 title_text=y_title_show,
                                 row=row_idx, col=col_idx)
            else:
                # n_vert: let Plotly auto-scale
                fig.update_yaxes(showline=True,
                                 title_text=y_title_show,
                                 row=row_idx, col=col_idx)

    # ---------------------------------------------------------------------------------
    # Legend — top-right of first subplot (row 1, col 1)
    # Row 1: 7 filled squares (ROI colors)    → "FE / LE"
    # Row 2: 7 white squares (ROI color border) → "AE / RE"
    # Row 3: 7 gray squares (light→dark)       → "CTRL" (only if in eye_conditions)
    # ---------------------------------------------------------------------------------
    first_y_range  = params_settings[0]['y_range']['all']
    y_top          = first_y_range[1] - 0.04 * (first_y_range[1] - first_y_range[0])
    legend_row_gap = 0.10 * (first_y_range[1] - first_y_range[0])
    legend_x_step  = marker_offset * 2.2
    # Anchor the legend block to the right of the subplot
    legend_x_end   = x_max - 1.5 * roi_spacing
    legend_x_start = legend_x_end - (len(rois) - 1) * legend_x_step

    legend_label_x = legend_x_end + 0.3 * roi_spacing

    # Determine which rows to draw
    legend_rows = []
    if 'FE-LE' in eye_conditions:
        legend_rows.append(('FE-LE', eye_labels['FE-LE']))
    if 'AE-RE' in eye_conditions:
        legend_rows.append(('AE-RE', eye_labels['AE-RE']))
    if 'CTRL' in eye_conditions:
        legend_rows.append(('CTRL', 'CTRL'))

    for col_leg in range(1, cols + 1):
        for row_k, (eye_cond, row_label) in enumerate(legend_rows):
            legend_y = y_top - row_k * legend_row_gap

            for j, roi in enumerate(rois):
                roi_color  = roi_colors[roi]
                legend_x   = legend_x_start + j * legend_x_step
                leg_marker = get_marker(eye_cond, roi_color, roi_idx=j)

                fig.add_trace(go.Scatter(
                    x=[legend_x], y=[legend_y],
                    mode='markers',
                    marker=leg_marker,
                    showlegend=False),
                    row=1, col=col_leg)

            # Row label at the right end
            fig.add_annotation(
                x=legend_label_x,
                y=legend_y,
                text=row_label,
                xanchor='left', yanchor='middle',
                showarrow=False,
                font=dict(color='rgba(0,0,0,1)',
                          family=template_specs['font'],
                          size=template_specs['axes_font_size']),
                row=1, col=col_leg)

    fig.update_layout(
        height=fig_height,
        width=fig_width,
        template=fig_template,
        showlegend=False,
        margin_l=100, margin_t=120, margin_r=50, margin_b=80)

    return fig
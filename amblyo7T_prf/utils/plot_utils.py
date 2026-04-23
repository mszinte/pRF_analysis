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
    Make correlation plot across eyes (AE/RE vs FE/LE).
    One figure with rows = ROIs, columns = corr_params.

    Parameters
    ----------
    tsv_dir : str
        Directory containing the per-parameter correlation TSVs
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

    # Axis settings per parameter
    corr_plot_settings = {
        'prf_rsq':    {'axes': 'pRF R<sup>2</sup>',
                       'range': [0, 0.8], 'tick_step': 0.2, 'n_ticks': 5},
        'prf_x':      {'axes': 'pRF x coord. (dva)',
                       'range': [-5, 5], 'tick_step': 2.5, 'n_ticks': 5},
        'prf_y':      {'axes': 'pRF y coord. (dva)',
                       'range': [-5, 5], 'tick_step': 2.5, 'n_ticks': 5},
        'prf_size':   {'axes': 'pRF size (dva)',
                       'range': [0, 4], 'tick_step': 1, 'n_ticks': 5},
        'prf_ecc':    {'axes': 'pRF ecc. (dva)',
                       'range': [0, 10], 'tick_step': 2, 'n_ticks': 6},
        'pcm_median': {'axes': 'pRF CM (mm/dva)',
                       'range': [0, 30], 'tick_step': 5, 'n_ticks': 7}
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
    rows, cols = len(rois), len(corr_params)

    rois_hor_spacing = figure_info['rois_hor_spacing']
    rois_ver_spacing = figure_info['rois_ver_spacing']
    rois_plot_height = figure_info['rois_plot_height']
    rois_plot_width = figure_info['rois_plot_width']
    subject_group = figure_info['subject_group']
    corr_bin_eye = figure_info['corr_bin_eye']
    corr_other_eye = 'FE-LE' if corr_bin_eye == 'AE-RE' else 'AE-RE'

    fig_height = rois_plot_height * rows + fig_margin[1] + fig_margin[3] + (rois_ver_spacing * (rows - 1))
    fig_width = rois_plot_width * cols + fig_margin[0] + fig_margin[2] + (rois_hor_spacing * (cols - 1))
    hor_spacing = rois_hor_spacing / (fig_width - fig_margin[0] - fig_margin[2])
    ver_spacing = rois_ver_spacing / (fig_height - fig_margin[1] - fig_margin[3])

    fig = make_subplots(rows=rows, cols=cols, print_grid=False,
                        horizontal_spacing=hor_spacing,
                        vertical_spacing=ver_spacing)

    # Axis labels depend on subject group and binning eye
    if subject_group == 'patient':
        eye_labels = {'AE-RE': 'AE', 'FE-LE': 'FE'}
    else:
        eye_labels = {'AE-RE': 'RE', 'FE-LE': 'LE'}
    x_label_prefix = eye_labels[corr_bin_eye]
    y_label_prefix = eye_labels[corr_other_eye]

    for l, corr_param in enumerate(corr_params):

        # Load per-parameter TSV
        tsv_fn = "{}/{}_{}_{}-corr.tsv".format(tsv_dir, subject, fn_spec_combined, corr_param)
        df = pd.read_table(tsv_fn, sep="\t")

        for j, roi in enumerate(rois):

            roi_color = roi_colors[roi]
            axis_settings = corr_plot_settings[corr_param]

            df_roi = df.loc[(df.roi == roi)]

            # X-axis eye (corr_bin_eye)
            param_x_median = np.array(df_roi[f'{corr_param}_{corr_bin_eye}_median'])
            param_x_upper_bound = np.array(df_roi[f'{corr_param}_{corr_bin_eye}_ci_upper_bound'])
            param_x_lower_bound = np.array(df_roi[f'{corr_param}_{corr_bin_eye}_ci_lower_bound'])

            # Y-axis eye (corr_other_eye)
            param_y_median = np.array(df_roi[f'{corr_param}_{corr_other_eye}_median'])
            param_y_upper_bound = np.array(df_roi[f'{corr_param}_{corr_other_eye}_ci_upper_bound'])
            param_y_lower_bound = np.array(df_roi[f'{corr_param}_{corr_other_eye}_ci_lower_bound'])

            # Binned R² median
            r2_median = np.array(df_roi[f'{rsq2use}_median'])

            # Reference line range
            line_x = np.linspace(axis_settings['range'][0], axis_settings['range'][1], 50)

            # Diagonal reference (identity line)
            fig.add_trace(go.Scatter(x=line_x, y=line_x, mode='lines',
                                     line=dict(color='rgba(0, 0, 0, 1)', width=2, dash='dash'),
                                     showlegend=False),
                          row=j + 1, col=l + 1)

            # Weighted linear regression
            slope, intercept = weighted_regression(param_x_median, param_y_median, r2_median, model='linear')
            line = slope * line_x + intercept
            fig.add_trace(go.Scatter(x=line_x, y=line, mode='lines',
                                     line=dict(color=roi_color, width=3),
                                     showlegend=False),
                          row=j + 1, col=l + 1)

            # Marker size scaled by R² (normalized per ROI)
            min_size, max_size = 5, 15
            r2_min = np.nanmin(r2_median)
            r2_max = np.nanmax(r2_median)
            r2_normalized = (r2_median - r2_min) / (r2_max - r2_min + 1e-8)
            marker_size = min_size + (r2_normalized * (max_size - min_size))

            # Scatter with error bars
            fig.add_trace(go.Scatter(
                x=param_x_median,
                y=param_y_median,
                mode='markers',
                error_x=dict(type='data',
                             array=param_x_upper_bound - param_x_median,
                             arrayminus=param_x_median - param_x_lower_bound,
                             visible=True, thickness=3, width=0, color=roi_color),
                error_y=dict(type='data',
                             array=param_y_upper_bound - param_y_median,
                             arrayminus=param_y_median - param_y_lower_bound,
                             visible=True, thickness=3, width=0, color=roi_color),
                marker=dict(color=roi_color, symbol='square', size=marker_size,
                            opacity=1.0, line=dict(color=roi_color, width=3)),
                showlegend=False),
                row=j + 1, col=l + 1)

            # ROI label annotation (top-right corner)
            r = axis_settings['range']
            annotation = go.layout.Annotation(
                x=r[1] - 0.05 * (r[1] - r[0]),
                y=r[1] - 0.9 * (r[1] - r[0]),
                text=roi, xanchor='right',
                showarrow=False, font_color=roi_color,
                font_family=template_specs['font'],
                font_size=template_specs['axes_font_size'])
            fig.add_annotation(annotation, row=j + 1, col=l + 1)

            # Axis labels: x only on last row, y on all rows
            x_title = f"{x_label_prefix} – {axis_settings['axes']}" if j == len(rois) - 1 else ''
            fig.update_xaxes(title_text=x_title,
                             range=axis_settings['range'],
                             tickmode='linear',
                             tick0=axis_settings['range'][0],
                             dtick=axis_settings['tick_step'],
                             showline=True,
                             row=j + 1, col=l + 1)
            fig.update_yaxes(title_text=f"{y_label_prefix} – {axis_settings['axes']}",
                             range=axis_settings['range'],
                             tickmode='linear',
                             tick0=axis_settings['range'][0],
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

            # Total vertices — transparent background with pattern
            
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
                            opacity=0.15,
                            ),
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
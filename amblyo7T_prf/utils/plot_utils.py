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
                                    yaxis_visible=True,
                                    yaxis_linewidth=template_specs['axes_width'],
                                    yaxis_color=template_specs['axes_color'],
                                    yaxis_showgrid=False,
                                    yaxis_ticks="outside"))

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
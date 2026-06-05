# General imports
import numpy as np

# Figure imports
import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots

# Debug
import ipdb
deb = ipdb.set_trace

def plotly_template(template_specs):
    """
    Define the template for plotly
    Parameters
    ----------
    template_specs : dict
        dictionary contain specific figure settings
    
    Returns
    -------
    fig_template : plotly.graph_objs.layout._template.Template
        Template for plotly figure
    """
    fig_template=go.layout.Template()

    # Violin plots
    fig_template.data.violin = [go.Violin(
                                    box_visible=False,
                                    points=False,
                                    # opacity=1,
                                    line_color= "rgba(0, 0, 0, 1)",
                                    line_width=template_specs['rois_plot_width'],
                                    width=0.8,
                                    #marker_symbol='x',
                                    #marker_opacity=1,
                                    hoveron='violins',
                                    meanline_visible=False,
                                    # meanline_color="rgba(0, 0, 0, 1)",
                                    # meanline_width=template_specs['rois_plot_width'],
                                    showlegend=False,
                                    )]

    # Barpolar
    fig_template.data.barpolar = [go.Barpolar(
                                    marker_line_color="rgba(0,0,0,1)",
                                    marker_line_width=template_specs['rois_plot_width'], 
                                    showlegend=False, 
                                    )]
    # Pie plots
    fig_template.data.pie = [go.Pie(textposition=["inside","none"],
                                    # marker_line_color=['rgba(0,0,0,1)','rgba(255,255,255,0)'],
                                    marker_line_width=0,#[template_specs['rois_plot_width'],0],
                                    rotation=0,
                                    direction="clockwise",
                                    hole=0.4,
                                    sort=False,
                                    )]

    # Layout
    fig_template.layout = (go.Layout(# general
                                    font_family=template_specs['font'],
                                    font_size=template_specs['axes_font_size'],
                                    plot_bgcolor=template_specs['bg_col'],

                                    # # x axis
                                    xaxis_visible=True,
                                    xaxis_linewidth=template_specs['axes_width'],
                                    xaxis_color= template_specs['axes_color'],
                                    xaxis_showgrid=False,
                                    xaxis_ticks="outside",
                                    xaxis_ticklen=8,
                                    xaxis_tickwidth = template_specs['axes_width'],
                                    xaxis_title_font_family=template_specs['font'],
                                    xaxis_title_font_size=template_specs['title_font_size'],
                                    xaxis_tickfont_family=template_specs['font'],
                                    xaxis_tickfont_size=template_specs['axes_font_size'],
                                    xaxis_zeroline=False,
                                    xaxis_zerolinecolor=template_specs['axes_color'],
                                    xaxis_zerolinewidth=template_specs['axes_width'],
                                    # xaxis_range=[0,1],
                                    xaxis_hoverformat = '.1f',
                                    
                                    # y axis
                                    yaxis_visible=True,
                                    yaxis_linewidth=template_specs['axes_width'],
                                    yaxis_color= template_specs['axes_color'],
                                    yaxis_showgrid=False,
                                    yaxis_ticks="outside",
                                    yaxis_ticklen=8,
                                    yaxis_tickwidth = template_specs['axes_width'],
                                    yaxis_tickfont_family=template_specs['font'],
                                    yaxis_tickfont_size=template_specs['axes_font_size'],
                                    yaxis_title_font_family=template_specs['font'],
                                    yaxis_title_font_size=template_specs['title_font_size'],
                                    yaxis_zeroline=False,
                                    yaxis_zerolinecolor=template_specs['axes_color'],
                                    yaxis_zerolinewidth=template_specs['axes_width'],
                                    yaxis_hoverformat = '.1f',

                                    # bar polar
                                    polar_radialaxis_visible = False,
                                    polar_radialaxis_showticklabels=False,
                                    polar_radialaxis_ticks='',
                                    polar_angularaxis_visible = False,
                                    polar_angularaxis_showticklabels = False,
                                    polar_angularaxis_ticks = ''
                                    ))

    # Annotations
    fig_template.layout.annotationdefaults = go.layout.Annotation(
                                    font_color=template_specs['axes_color'],
                                    font_family=template_specs['font'],
                                    font_size=template_specs['title_font_size'])

    return fig_template


def nCSF_roi_active_vert_plot(df, figure_info, format):
    """
    Make bar plots of each roi number of vertex and the corresponding significative activer vertex for nCSF  
    
    Parameters
    ----------
    df : dataframe for corresponding plot
    figure_info : dict with figure settings
    format : format of data to define axis size
    
    Returns
    -------
    fig : bar plot
    """
    
    # General figure settings
    template_specs = dict(axes_color="rgba(0, 0, 0, 1)",
                          axes_width=2,
                          axes_font_size=15,
                          bg_col="rgba(255, 255, 255, 1)",
                          font='Arial',
                          title_font_size=15,
                          rois_plot_width=1.5)
    
    # General figure settings
    fig_template = plotly_template(template_specs)    
    rows = 1    
    cols = 2
    rois = figure_info['rois']
    roi_colors = figure_info['roi_colors']
    fig_margin = figure_info['rois_fig_margin']
    rois_hor_spacing = figure_info['rois_hor_spacing']
    rois_ver_spacing = figure_info['rois_ver_spacing']
    bar_width = figure_info['rois_bar_width']
    rois_plot_height = figure_info['rois_plot_height']
    
    fig_height = rois_plot_height * rows + fig_margin[1] + fig_margin[3] + (rois_ver_spacing * (rows-1))
    fig_width = bar_width * cols * len(rois) + fig_margin[0] + fig_margin[2] + (rois_hor_spacing * (cols-1))
    hor_spacing = rois_hor_spacing / (fig_width - fig_margin[0] - fig_margin[2])
    ver_spacing = rois_ver_spacing / (fig_height - fig_margin[1] - fig_margin[3])

    # colors 
    roi_colors = list(roi_colors.values())

    # General settings
    fig = make_subplots(rows=rows, 
                        cols=cols, 
                        subplot_titles=['FDR threshold = 0.05', 'FDR threshold = 0.01'],
                        vertical_spacing=ver_spacing, 
                        horizontal_spacing=hor_spacing
                       )
    
    # FDR 0.05 
    # All vertices
    fig.add_trace(go.Bar(x=df['roi'], 
                         y=df['n_vert_tot'], 
                         text=(df['ratio_5pt']*100).astype(int).astype(str) + '%',
                         textposition='outside',
                         textangle=-60,
                         showlegend=False, 
                         marker=dict(color=roi_colors, opacity=0.2)),
                 row=1, col=1)
 
    # Significant vertices
    fig.add_trace(go.Bar(x=df['roi'], 
                         y=df['n_vert_corr_pvalue_5pt'], 
                         showlegend=False, 
                         marker=dict(color=roi_colors)),
                 row=1, col=1)
    
    
    # FDR 0.01
    # All vertices
    fig.add_trace(go.Bar(x=df['roi'], 
                         y=df['n_vert_tot'], 
                         text=(df['ratio_1pt']*100).astype(int).astype(str) + '%',
                         textposition='outside',
                         textangle=-60,
                         showlegend=False, 
                         marker=dict(color=roi_colors, opacity=0.2)),
                 row=1, col=2)
 
    # Significant vertices
    fig.add_trace(go.Bar(x=df['roi'], 
                         y=df['n_vert_corr_pvalue_1pt'], 
                         showlegend=False, 
                         marker=dict(color=roi_colors)),
                 row=1, col=2)

    # Define parameters
    fig.update_xaxes(showline=True, 
                     ticklen=0, 
                     linecolor=('rgba(255,255,255,0)'))      

    if format == 'fsnative':
        range_val = figure_info['active_vert_fsnative_range']
    elif format == '170k':
        range_val = figure_info['active_vert_170k_range']
    
    fig.update_yaxes(range=range_val, 
                     showline=True, 
                     nticks=10, 
                     title_text='Number of vertex',secondary_y=False)
    
    fig.update_layout(barmode='overlay',
                      height=fig_height, 
                      width=fig_width, 
                      template=fig_template,
                      margin_l=fig_margin[0], 
                      margin_t=fig_margin[1], 
                      margin_r=fig_margin[2], 
                      margin_b=fig_margin[3]
                     )

    # Return outputs
    return fig

def ncsf_violins_plot(df, figure_info, rsq2use):
    """
    Make violins plots for pRF loo_r2, size, n and pcm

    Parameters
    ----------
    df : dataframe
    figure_info : dict with figure settings
    rsq2use : rsquare value to use
    
    Returns
    -------
    fig : violins plot
    """
    
    # General figure settings
    template_specs = dict(axes_color="rgba(0, 0, 0, 1)",
                          axes_width=2,
                          axes_font_size=15,
                          bg_col="rgba(255, 255, 255, 1)",
                          font='Arial',
                          title_font_size=15,
                          rois_plot_width=1.5)
    
    # General figure settings
    fig_template = plotly_template(template_specs)
    rows = 3
    cols = 3
    rois = figure_info['rois']
    roi_colors = figure_info['roi_colors']
    fig_margin = figure_info['rois_fig_margin']
    rois_hor_spacing = figure_info['rois_hor_spacing']
    rois_ver_spacing = figure_info['rois_ver_spacing']
    bar_width = figure_info['rois_bar_width']
    rois_plot_height = figure_info['rois_plot_height']
    
    fig_height = rois_plot_height * rows + fig_margin[1] + fig_margin[3] + (rois_ver_spacing * (rows-1))
    fig_width = bar_width * cols * len(rois) + fig_margin[0] + fig_margin[2] + (rois_hor_spacing * (cols-1))
    hor_spacing = rois_hor_spacing / (fig_width - fig_margin[0] - fig_margin[2])
    ver_spacing = rois_ver_spacing / (fig_height - fig_margin[1] - fig_margin[3])
    
    fig = make_subplots(rows=rows, 
                        cols=cols, 
                        print_grid=False, 
                        vertical_spacing=ver_spacing,
                        horizontal_spacing=hor_spacing)

    for j, roi in enumerate(rois):
        
        df_roi = df.loc[(df.roi == roi)]
        
        # nCSF r2 or loor2    
        fig.add_trace(go.Violin(x=df_roi.roi[df_roi.roi==roi], 
                                y=df_roi[rsq2use], 
                                name=roi, 
                                opacity=1,
                                showlegend=False, 
                                points=False, 
                                spanmode='manual', 
                                span=figure_info['violin_ncsf_rsq_range'],
                                scalemode='width', 
                                fillcolor=roi_colors[roi],
                                line_color=roi_colors[roi]), 
                      row=1, col=1)
                
        # SFp
        fig.add_trace(go.Violin(x=df_roi.roi[df_roi.roi==roi], 
                                y=df_roi.SFp, 
                                name=roi, 
                                opacity=1,
                                showlegend=False, 
                                points=False, 
                                spanmode='manual', 
                                span=figure_info['violin_SFp_range'],
                                scalemode='width', 
                                fillcolor=roi_colors[roi],
                                line_color=roi_colors[roi]), 
                      row=1, col=2)
        
        # CSp
        fig.add_trace(go.Violin(x=df_roi.roi[df_roi.roi==roi], 
                                y=df_roi.CSp, 
                                name=roi, 
                                opacity=1,
                                showlegend=False, 
                                points=False, 
                                spanmode='manual', 
                                span=figure_info['violin_CSp_range'],
                                scalemode='width', 
                                fillcolor=roi_colors[roi],
                                line_color=roi_colors[roi]), 
                      row=1, col=3)
        
        # width_l
        fig.add_trace(go.Violin(x=df_roi.roi[df_roi.roi==roi], 
                                y=df_roi.width_l, 
                                name=roi, 
                                opacity=1,
                                showlegend=False, 
                                points=False, 
                                spanmode='manual', 
                                span=figure_info['violin_width_l_range'],
                                scalemode='width', 
                                fillcolor=roi_colors[roi],
                                line_color=roi_colors[roi]), 
                      row=2, col=1)
        
        # width_r
        fig.add_trace(go.Violin(x=df_roi.roi[df_roi.roi==roi], 
                                y=df_roi.width_r, 
                                name=roi, 
                                opacity=1,
                                showlegend=False, 
                                points=False, 
                                spanmode='manual', 
                                span=figure_info['violin_width_r_range'],
                                scalemode='width', 
                                fillcolor=roi_colors[roi],
                                line_color=roi_colors[roi]), 
                      row=2, col=2)
            
        # crf_exp
        fig.add_trace(go.Violin(x=df_roi.roi[df_roi.roi==roi], 
                        y=df_roi.crf_exp, 
                        name=roi, 
                        opacity=1,
                        showlegend=False, 
                        points=False, 
                        spanmode='manual', 
                        span=figure_info['violin_crf_exp_range'],
                        scalemode='width', 
                        fillcolor=roi_colors[roi],
                        line_color=roi_colors[roi]), 
              row=2, col=3)
        
        # auc
        fig.add_trace(go.Violin(x=df_roi.roi[df_roi.roi==roi], 
                                y=df_roi.auc, 
                                name=roi, 
                                opacity=1,
                                showlegend=False, 
                                points=False, 
                                spanmode='manual', 
                                span=figure_info['violin_auc_range'],
                                scalemode='width', 
                                fillcolor=roi_colors[roi],
                                line_color=roi_colors[roi]), 
                      row=3, col=1)
        

        
        # normalize_auc
        fig.add_trace(go.Violin(x=df_roi.roi[df_roi.roi==roi], 
                                y=df_roi.normalize_auc, 
                                name=roi, 
                                opacity=1,
                                showlegend=False, 
                                points=False, 
                                spanmode='manual', 
                                span=figure_info['violin_normalize_auc_range'],
                                scalemode='width', 
                                fillcolor=roi_colors[roi],
                                line_color=roi_colors[roi]), 
                      row=3, col=2)
        
        # SFmax
        fig.add_trace(go.Violin(x=df_roi.roi[df_roi.roi==roi], 
                        y=df_roi.SFmax, 
                        name=roi, 
                        opacity=1,
                        showlegend=False, 
                        points=False, 
                        spanmode='manual', 
                        span=figure_info['violin_SFmax_range'],
                        scalemode='width', 
                        fillcolor=roi_colors[roi],
                        line_color=roi_colors[roi]), 
              row=3, col=3)
        


    fig.update_yaxes(showline=True, range=figure_info['violin_ncsf_rsq_range'], 
                     title_text=rsq2use, row=1, col=1)    
    fig.update_yaxes(showline=True, range=figure_info['violin_SFp_range'], 
                     title_text="SFp", row=1, col=2)
    fig.update_yaxes(showline=True, range=figure_info['violin_CSp_range'], 
                     title_text="CSp", row=1, col=3)

    fig.update_yaxes(showline=True, range=figure_info['violin_width_l_range'], 
                     title_text="width_l", row=2, col=1)
    fig.update_yaxes(showline=True, range=figure_info['violin_width_r_range'], 
                     title_text="width_r", row=2, col=2)
    fig.update_yaxes(showline=True, range=figure_info['violin_crf_exp_range'], 
                     title_text="crf_exp", row=2, col=3)
    
    fig.update_yaxes(showline=True, range=figure_info['violin_auc_range'], 
                     title_text="auc", row=3, col=1)
    fig.update_yaxes(showline=True, range=figure_info['violin_normalize_auc_range'], 
                     title_text="normalize_auc", row=3, col=2)
    fig.update_yaxes(showline=True, range=figure_info['violin_SFmax_range'], 
                     title_text="SFmax", row=3, col=3)
        
        
    fig.update_xaxes(showline=True, ticklen=0, linecolor='rgba(255,255,255,0)')
        
    fig.update_layout(height=fig_height, 
                      width=fig_width, 
                      showlegend=False,
                      legend=dict(orientation="h", 
                                  font_family=template_specs['font'],
                                  font_size=template_specs['axes_font_size'],
                                  y=1.1, 
                                  yanchor='top', 
                                  xanchor='left', 
                                  traceorder='normal', 
                                  itemwidth=30), 
                      template=fig_template,
                      margin_l=fig_margin[0], 
                      margin_t=fig_margin[1], 
                      margin_r=fig_margin[2], 
                      margin_b=fig_margin[3]
                     )

    return fig

def ncsf_params_median_plot(df, figure_info, rsq2use):

    template_specs = dict(axes_color="rgba(0, 0, 0, 1)",
                          axes_width=2,
                          axes_font_size=15,
                          bg_col="rgba(255, 255, 255, 1)",
                          font='Arial',
                          title_font_size=15,
                          rois_plot_width=1.5)

    fig_template = plotly_template(template_specs)

    rows, cols = 3, 3
    rois = figure_info['rois']
    roi_colors = figure_info['roi_colors']
    fig_margin = figure_info['rois_fig_margin']
    rois_hor_spacing = figure_info['rois_hor_spacing']
    rois_ver_spacing = figure_info['rois_ver_spacing']
    bar_width = figure_info['rois_bar_width']
    rois_plot_height = figure_info['rois_plot_height']

    fig_height = rois_plot_height * rows + fig_margin[1] + fig_margin[3] + (rois_ver_spacing * (rows - 1))
    fig_width = bar_width * cols * len(rois) + fig_margin[0] + fig_margin[2] + (rois_hor_spacing * (cols - 1))
    hor_spacing = rois_hor_spacing / (fig_height - fig_margin[1] - fig_margin[3])
    ver_spacing = rois_ver_spacing / (fig_height - fig_margin[1] - fig_margin[3])

    fig = make_subplots(rows=rows,
                        cols=cols,
                        print_grid=False,
                        vertical_spacing=ver_spacing,
                        horizontal_spacing=hor_spacing)

    for j, roi in enumerate(rois):

        df_roi = df.loc[df.roi == roi]
        # rsq2use
        fig.add_trace(go.Scatter(
            x=[roi], y=tuple(df_roi[f'{rsq2use}_weighted_median']),
            mode='markers',
            error_y=dict(type='data',
                         array=[df_roi[f'{rsq2use}_ci_up'] - df_roi[f'{rsq2use}_weighted_median']],
                         arrayminus=[df_roi[f'{rsq2use}_weighted_median'] - df_roi[f'{rsq2use}_ci_down']],
                         visible=True, thickness=3, width=0, color=roi_colors[roi]),
            marker=dict(symbol="square", color=roi_colors[roi], size=12,
                        line=dict(color=roi_colors[roi], width=3)),
            name=roi, showlegend=False),
            row=1, col=1)

        # SFp
        fig.add_trace(go.Scatter(
            x=[roi], y=tuple(df_roi['SFp_weighted_median']),
            mode='markers',
            error_y=dict(type='data',
                         array=[df_roi['SFp_ci_up'] - df_roi['SFp_weighted_median']],
                         arrayminus=[df_roi['SFp_weighted_median'] - df_roi['SFp_ci_down']],
                         visible=True, thickness=3, width=0, color=roi_colors[roi]),
            marker=dict(symbol="square", color=roi_colors[roi], size=12,
                        line=dict(color=roi_colors[roi], width=3)),
            name=roi, showlegend=False),
            row=1, col=2)

        # CSp
        fig.add_trace(go.Scatter(
            x=[roi], y=tuple(df_roi['CSp_median']),
            mode='markers',
            error_y=dict(type='data',
                         array=[df_roi['CSp_ci_up'] - df_roi['CSp_median']],
                         arrayminus=[df_roi['CSp_median'] - df_roi['CSp_ci_down']],
                         visible=True, thickness=3, width=0, color=roi_colors[roi]),
            marker=dict(symbol="square", color=roi_colors[roi], size=12,
                        line=dict(color=roi_colors[roi], width=3)),
            name=roi, showlegend=False),
            row=1, col=3)

        # width_l
        fig.add_trace(go.Scatter(
            x=[roi], y=tuple(df_roi['width_l_weighted_median']),
            mode='markers',
            error_y=dict(type='data',
                         array=[df_roi['width_l_ci_up'] - df_roi['width_l_weighted_median']],
                         arrayminus=[df_roi['width_l_weighted_median'] - df_roi['width_l_ci_down']],
                         visible=True, thickness=3, width=0, color=roi_colors[roi]),
            marker=dict(symbol="square", color=roi_colors[roi], size=12,
                        line=dict(color=roi_colors[roi], width=3)),
            name=roi, showlegend=False),
            row=2, col=1)

        # width_r
        fig.add_trace(go.Scatter(
            x=[roi], y=tuple(df_roi['width_r_weighted_median']),
            mode='markers',
            error_y=dict(type='data',
                         array=[df_roi['width_r_ci_up'] - df_roi['width_r_weighted_median']],
                         arrayminus=[df_roi['width_r_weighted_median'] - df_roi['width_r_ci_down']],
                         visible=True, thickness=3, width=0, color=roi_colors[roi]),
            marker=dict(symbol="square", color=roi_colors[roi], size=12,
                        line=dict(color=roi_colors[roi], width=3)),
            name=roi, showlegend=False),
            row=2, col=2)

        # crf_exp
        fig.add_trace(go.Scatter(
            x=[roi], y=tuple(df_roi['crf_exp_weighted_median']),
            mode='markers',
            error_y=dict(type='data',
                         array=[df_roi['crf_exp_ci_up'] - df_roi['crf_exp_weighted_median']],
                         arrayminus=[df_roi['crf_exp_weighted_median'] - df_roi['crf_exp_ci_down']],
                         visible=True, thickness=3, width=0, color=roi_colors[roi]),
            marker=dict(symbol="square", color=roi_colors[roi], size=12,
                        line=dict(color=roi_colors[roi], width=3)),
            name=roi, showlegend=False),
            row=2, col=3)

        # auc
        fig.add_trace(go.Scatter(
            x=[roi], y=tuple(df_roi['auc_weighted_median']),
            mode='markers',
            error_y=dict(type='data',
                         array=[df_roi['auc_ci_up'] - df_roi['auc_weighted_median']],
                         arrayminus=[df_roi['auc_weighted_median'] - df_roi['auc_ci_down']],
                         visible=True, thickness=3, width=0, color=roi_colors[roi]),
            marker=dict(symbol="square", color=roi_colors[roi], size=12,
                        line=dict(color=roi_colors[roi], width=3)),
            name=roi, showlegend=False),
            row=3, col=1)

        # normalize_auc
        fig.add_trace(go.Scatter(
            x=[roi], y=tuple(df_roi['normalize_auc_weighted_median']),
            mode='markers',
            error_y=dict(type='data',
                         array=[df_roi['normalize_auc_ci_up'] - df_roi['normalize_auc_weighted_median']],
                         arrayminus=[df_roi['normalize_auc_weighted_median'] - df_roi['normalize_auc_ci_down']],
                         visible=True, thickness=3, width=0, color=roi_colors[roi]),
            marker=dict(symbol="square", color=roi_colors[roi], size=12,
                        line=dict(color=roi_colors[roi], width=3)),
            name=roi, showlegend=False),
            row=3, col=2)

        # SFmax
        fig.add_trace(go.Scatter(
            x=[roi], y=tuple(df_roi['SFmax_weighted_median']),
            mode='markers',
            error_y=dict(type='data',
                         array=[df_roi['SFmax_ci_up'] - df_roi['SFmax_weighted_median']],
                         arrayminus=[df_roi['SFmax_weighted_median'] - df_roi['SFmax_ci_down']],
                         visible=True, thickness=3, width=0, color=roi_colors[roi]),
            marker=dict(symbol="square", color=roi_colors[roi], size=12,
                        line=dict(color=roi_colors[roi], width=3)),
            name=roi, showlegend=False),
            row=3, col=3)

    fig.update_yaxes(showline=True, range=figure_info['violin_ncsf_rsq_range'], 
                     title_text=rsq2use, row=1, col=1)    
    fig.update_yaxes(showline=True, range=figure_info['violin_SFp_range'], 
                     title_text="SFp", row=1, col=2)
    fig.update_yaxes(showline=True, range=figure_info['violin_CSp_range'], 
                     title_text="CSp", row=1, col=3)

    fig.update_yaxes(showline=True, range=figure_info['violin_width_l_range'], 
                     title_text="width_l", row=2, col=1)
    fig.update_yaxes(showline=True, range=figure_info['violin_width_r_range'], 
                     title_text="width_r", row=2, col=2)
    fig.update_yaxes(showline=True, range=figure_info['violin_crf_exp_range'], 
                     title_text="crf_exp", row=2, col=3)
    
    fig.update_yaxes(showline=True, range=figure_info['violin_auc_range'], 
                     title_text="auc", row=3, col=1)
    fig.update_yaxes(showline=True, range=figure_info['violin_normalize_auc_range'], 
                     title_text="normalize_auc", row=3, col=2)
    fig.update_yaxes(showline=True, range=figure_info['violin_SFmax_range'], 
                     title_text="SFmax", row=3, col=3)

    fig.update_xaxes(showline=True, ticklen=0, linecolor='rgba(255,255,255,0)')

    fig.update_layout(height=fig_height,
                      width=fig_width,
                      showlegend=False,
                      legend=dict(orientation="h",
                                  font_family=template_specs['font'],
                                  font_size=template_specs['axes_font_size'],
                                  y=1.1,
                                  yanchor='top',
                                  xanchor='left',
                                  traceorder='normal',
                                  itemwidth=30),
                      template=fig_template,
                      margin_l=fig_margin[0],
                      margin_t=fig_margin[1],
                      margin_r=fig_margin[2],
                      margin_b=fig_margin[3])

    return fig

def ecc_SFp_plot(df, figure_info, rsq2use):
    """
    Make scatter plot for linear relationship between eccentricity and nCSF SFp

    Parameters
    ----------
    df : A data dataframe
    figure_info : dict with figure settings
    rsq2use : rsquare to use
    
    Returns
    -------
    fig : eccentricy as a function of size plot
    """

    # General figure settings
    template_specs = dict(axes_color="rgba(0, 0, 0, 1)",
                          axes_width=2,
                          axes_font_size=15,
                          bg_col="rgba(255, 255, 255, 1)",
                          font='Arial',
                          title_font_size=15,
                          rois_plot_width=1.5)
    
    # General figure settings
    fig_template = plotly_template(template_specs)
    roi_colors = figure_info['roi_colors']
    fig_margin = figure_info['rois_fig_margin']
    rois_groups = figure_info['rois_groups_plot']
    rows, cols = 1, len(rois_groups)
    rois_hor_spacing = figure_info['rois_hor_spacing']
    rois_ver_spacing = figure_info['rois_ver_spacing']
    rois_plot_height = figure_info['rois_plot_height']
    rois_plot_width = figure_info['rois_plot_width']    
    ecc_axis = figure_info['ecc_SFp_x_axis']
    SFp_axis = figure_info['ecc_SFp_y_axis']
    max_SFp = figure_info['ecc_SFp_max'][1]
    
    fig_height = rois_plot_height * rows + fig_margin[1] + fig_margin[3] + (rois_ver_spacing * (rows-1))
    fig_width = rois_plot_width * cols + fig_margin[0] + fig_margin[2] + (rois_hor_spacing * (cols-1))
    hor_spacing = rois_hor_spacing / (fig_width - fig_margin[0] - fig_margin[2])
    ver_spacing = rois_ver_spacing / (fig_height - fig_margin[1] - fig_margin[3])

    # General settings
    fig = make_subplots(rows=rows, cols=cols, print_grid=False, 
                       horizontal_spacing=hor_spacing,
                       vertical_spacing=ver_spacing)
    
    for l, line_label in enumerate(rois_groups):
        for j, roi in enumerate(line_label):
            
            # Parametring colors
            roi_color = roi_colors[roi]
            
            # Get data
            df_roi = df.loc[(df.roi == roi)]
            ecc_median = np.array(df_roi.prf_ecc_bins)
            SFp_median = np.array(df_roi.SFp_bins_median)
            SFp_upper_bound = np.array(df_roi.SFp_bins_ci_upper_bound)
            SFp_lower_bound = np.array(df_roi.SFp_bins_ci_lower_bound)
            
            # Markers
            fig.add_trace(go.Scatter(x=ecc_median, 
                                     y=SFp_median, 
                                     mode='markers', 
                                     error_y=dict(type='data', 
                                                  array=SFp_upper_bound - SFp_median, 
                                                  arrayminus=SFp_median - SFp_lower_bound,
                                                  visible=True, 
                                                  thickness=3, 
                                                  width=0, 
                                                  color=roi_color),
                                      marker=dict(color=roi_color,
                                                  symbol='square',
                                                  size=8, 
                                                  line=dict(color=roi_color, 
                                                            width=3)), 
                                      showlegend=False), 
                          row=1, col=l + 1)
            
            # Add legend
            annotation = go.layout.Annotation(x=1, y=max_SFp-(j*0.1*max_SFp), text=roi, xanchor='left',
                                              showarrow=False, font_color=roi_color, 
                                              font_family=template_specs['font'],
                                              font_size=template_specs['axes_font_size'],
                                             )
            fig.add_annotation(annotation, row=1, col=l+1)

        # Set axis
        fig.update_xaxes(title_text='pRF eccentricity (dva)', range=[ecc_axis[0], ecc_axis[1]], showline=True)
        fig.update_yaxes(title_text='nCSF SFp (cycles/dva)', range=[SFp_axis[0], SFp_axis[1]], showline=True)
        fig.update_layout(height=fig_height, 
                          width=fig_width, 
                          showlegend=False, 
                          template=fig_template,
                          margin_l=fig_margin[0], 
                          margin_t=fig_margin[1], 
                          margin_r=fig_margin[2], 
                          margin_b=fig_margin[3]
                         )
        
    return fig

def ecc_auc_plot(df, figure_info, rsq2use):
    """
    Make scatter plot for linear relationship between eccentricity and nCSF auc

    Parameters
    ----------
    df : A data dataframe
    figure_info : dict with figure settings
    rsq2use : rsquare to use
    
    Returns
    -------
    fig : eccentricy as a function of size plot
    """
    
    # General figure settings
    template_specs = dict(axes_color="rgba(0, 0, 0, 1)",
                          axes_width=2,
                          axes_font_size=15,
                          bg_col="rgba(255, 255, 255, 1)",
                          font='Arial',
                          title_font_size=15,
                          rois_plot_width=1.5)
    
    # General figure settings
    fig_template = plotly_template(template_specs)
    max_auc = figure_info['ecc_auc_max'][1]
    roi_colors = figure_info['roi_colors']
    fig_margin = figure_info['rois_fig_margin']
    rois_groups = figure_info['rois_groups_plot']
    rows, cols = 1, len(rois_groups)
    rois_hor_spacing = figure_info['rois_hor_spacing']
    rois_ver_spacing = figure_info['rois_ver_spacing']
    rois_plot_height = figure_info['rois_plot_height']
    rois_plot_width = figure_info['rois_plot_width']    
    ecc_axis = figure_info['ecc_SFp_x_axis']
    auc_axis = figure_info['ecc_auc_y_axis']
    
    fig_height = rois_plot_height * rows + fig_margin[1] + fig_margin[3] + (rois_ver_spacing * (rows-1))
    fig_width = rois_plot_width * cols + fig_margin[0] + fig_margin[2] + (rois_hor_spacing * (cols-1))
    hor_spacing = rois_hor_spacing / (fig_width - fig_margin[0] - fig_margin[2])
    ver_spacing = rois_ver_spacing / (fig_height - fig_margin[1] - fig_margin[3])

    # General settings
    
    fig = make_subplots(rows=rows, cols=cols, print_grid=False, 
                       horizontal_spacing=hor_spacing,
                       vertical_spacing=ver_spacing)
    
    for l, line_label in enumerate(rois_groups):
        for j, roi in enumerate(line_label):
            
            # Parametring colors
            roi_color = roi_colors[roi]
            
            # Get data
            df_roi = df.loc[(df.roi == roi)]
            ecc_median = np.array(df_roi.prf_ecc_bins)
            auc_median = np.array(df_roi.auc_bins_median)
            auc_upper_bound = np.array(df_roi.auc_bins_ci_upper_bound)
            auc_lower_bound = np.array(df_roi.auc_bins_ci_lower_bound)
            
            # Markers
            fig.add_trace(go.Scatter(x=ecc_median, 
                                     y=auc_median, 
                                     mode='markers', 
                                     error_y=dict(type='data', 
                                                  array=auc_upper_bound - auc_median, 
                                                  arrayminus=auc_median - auc_lower_bound,
                                                  visible=True, 
                                                  thickness=3, 
                                                  width=0, 
                                                  color=roi_color),
                                      marker=dict(color=roi_color,
                                                  symbol='square',
                                                  size=8, 
                                                  line=dict(color=roi_color, 
                                                            width=3)), 
                                      showlegend=False), 
                          row=1, col=l + 1)
            
            # Add legend
            annotation = go.layout.Annotation(x=1, y=max_auc - (j*0.1*max_auc), text=roi, xanchor='left',
                                              showarrow=False, font_color=roi_color, 
                                              font_family=template_specs['font'],
                                              font_size=template_specs['axes_font_size'],
                                             )
            fig.add_annotation(annotation, row=1, col=l+1)

        # Set axis
        fig.update_xaxes(title_text='pRF eccentricity (dva)', range=[ecc_axis[0], ecc_axis[1]], showline=True)
        fig.update_yaxes(title_text='nCSF auc', range=[auc_axis[0], auc_axis[1]], showline=True)
        fig.update_layout(height=fig_height, 
                          width=fig_width, 
                          showlegend=False, 
                          template=fig_template,
                          margin_l=fig_margin[0], 
                          margin_t=fig_margin[1], 
                          margin_r=fig_margin[2], 
                          margin_b=fig_margin[3]
                         )
        
    return fig

def size_SFp_plot(df, figure_info, rsq2use):
    """
    Make scatter plot for linear relationship between size and SFp

    Parameters
    ----------
    df : A data dataframe
    figure_info : dict with figure settings
    rsq2use : rsquare to use
    
    Returns
    -------
    fig : eccentricy as a function of size plot
    """
    

    # General figure settings
    template_specs = dict(axes_color="rgba(0, 0, 0, 1)",
                          axes_width=2,
                          axes_font_size=15,
                          bg_col="rgba(255, 255, 255, 1)",
                          font='Arial',
                          title_font_size=15,
                          rois_plot_width=1.5)
    
    # General figure settings
    fig_template = plotly_template(template_specs)
    roi_colors = figure_info['roi_colors']
    fig_margin = figure_info['rois_fig_margin']
    rois_groups = figure_info['rois_groups_plot']
    rows, cols = 1, len(rois_groups)
    rois_hor_spacing = figure_info['rois_hor_spacing']
    rois_ver_spacing = figure_info['rois_ver_spacing']
    rois_plot_height = figure_info['rois_plot_height']
    rois_plot_width = figure_info['rois_plot_width']    
    max_SFp = figure_info['ecc_SFp_max'][1]
    size_axis = figure_info['size_SFp_x_axis']
    SFp_axis = figure_info['ecc_SFp_y_axis']
    
    
    
    fig_height = rois_plot_height * rows + fig_margin[1] + fig_margin[3] + (rois_ver_spacing * (rows-1))
    fig_width = rois_plot_width * cols + fig_margin[0] + fig_margin[2] + (rois_hor_spacing * (cols-1))
    hor_spacing = rois_hor_spacing / (fig_width - fig_margin[0] - fig_margin[2])
    ver_spacing = rois_ver_spacing / (fig_height - fig_margin[1] - fig_margin[3])

    # General settings
    
    fig = make_subplots(rows=rows, cols=cols, print_grid=False, 
                       horizontal_spacing=hor_spacing,
                       vertical_spacing=ver_spacing)
    
    for l, line_label in enumerate(rois_groups):
        for j, roi in enumerate(line_label):
            
            # Parametring colors
            roi_color = roi_colors[roi]
            
            # Get data
            df_roi = df.loc[(df.roi == roi)]
            size_median = np.array(df_roi.prf_size_bins)
            SFp_median = np.array(df_roi.sfp_bins_median)
            SFp_upper_bound = np.array(df_roi.sfp_bins_ci_upper_bound)
            SFp_lower_bound = np.array(df_roi.sfp_bins_ci_lower_bound)
            # Markers
            fig.add_trace(go.Scatter(x=size_median, 
                                     y=SFp_median, 
                                     mode='markers', 
                                     error_y=dict(type='data', 
                                                  array=SFp_upper_bound - SFp_median, 
                                                  arrayminus=SFp_median - SFp_lower_bound,
                                                  visible=True, 
                                                  thickness=3, 
                                                  width=0, 
                                                  color=roi_color),
                                      marker=dict(color=roi_color,
                                                  symbol='square',
                                                  size=8, 
                                                  line=dict(color=roi_color, 
                                                            width=3)), 
                                      showlegend=False), 
                          row=1, col=l + 1)
            
            # Add legend
            annotation = go.layout.Annotation(x=1, y=max_SFp-(j*0.1*max_SFp), text=roi, xanchor='left',
                                              showarrow=False, font_color=roi_color, 
                                              font_family=template_specs['font'],
                                              font_size=template_specs['axes_font_size'],
                                             )
            fig.add_annotation(annotation, row=1, col=l+1)

        # Set axis
        fig.update_xaxes(title_text='pRF size (dva)', range=[size_axis[0], size_axis[1]], showline=True)
        fig.update_yaxes(title_text='nCSF SFp (cyles/dva)', range=[SFp_axis[0], SFp_axis[1]], showline=True)
        fig.update_layout(height=fig_height, 
                          width=fig_width, 
                          showlegend=False, 
                          template=fig_template,
                          margin_l=fig_margin[0], 
                          margin_t=fig_margin[1], 
                          margin_r=fig_margin[2], 
                          margin_b=fig_margin[3]
                         )
        
    return fig

def size_auc_plot(df, figure_info, rsq2use):
    """
    Make scatter plot for linear relationship between size and auc

    Parameters
    ----------
    df : A data dataframe
    figure_info : dict with figure settings
    rsq2use : rsquare to use
    
    Returns
    -------
    fig : eccentricy as a function of size plot
    """
    
    # General figure settings
    template_specs = dict(axes_color="rgba(0, 0, 0, 1)",
                          axes_width=2,
                          axes_font_size=15,
                          bg_col="rgba(255, 255, 255, 1)",
                          font='Arial',
                          title_font_size=15,
                          rois_plot_width=1.5)
    
    # General figure settings
    fig_template = plotly_template(template_specs)
    roi_colors = figure_info['roi_colors']
    fig_margin = figure_info['rois_fig_margin']
    rois_groups = figure_info['rois_groups_plot']
    rows, cols = 1, len(rois_groups)
    rois_hor_spacing = figure_info['rois_hor_spacing']
    rois_ver_spacing = figure_info['rois_ver_spacing']
    rois_plot_height = figure_info['rois_plot_height']
    rois_plot_width = figure_info['rois_plot_width']    
    max_auc = figure_info['ecc_auc_max'][1]
    size_axis = figure_info['size_SFp_x_axis']
    auc_axis = figure_info['ecc_auc_y_axis']
    
    fig_height = rois_plot_height * rows + fig_margin[1] + fig_margin[3] + (rois_ver_spacing * (rows-1))
    fig_width = rois_plot_width * cols + fig_margin[0] + fig_margin[2] + (rois_hor_spacing * (cols-1))
    hor_spacing = rois_hor_spacing / (fig_width - fig_margin[0] - fig_margin[2])
    ver_spacing = rois_ver_spacing / (fig_height - fig_margin[1] - fig_margin[3])

    # General settings
    
    fig = make_subplots(rows=rows, cols=cols, print_grid=False, 
                       horizontal_spacing=hor_spacing,
                       vertical_spacing=ver_spacing)
    
    for l, line_label in enumerate(rois_groups):
        for j, roi in enumerate(line_label):
            
            # Parametring colors
            roi_color = roi_colors[roi]
            
            # Get data
            df_roi = df.loc[(df.roi == roi)]
            size_median = np.array(df_roi.prf_size_bins)
            auc_median = np.array(df_roi.auc_bins_median)
            auc_upper_bound = np.array(df_roi.auc_bins_ci_upper_bound)
            auc_lower_bound = np.array(df_roi.auc_bins_ci_lower_bound)
           
            # Markers
            fig.add_trace(go.Scatter(x=size_median, 
                                     y=auc_median, 
                                     mode='markers', 
                                     error_y=dict(type='data', 
                                                  array=auc_upper_bound - auc_median, 
                                                  arrayminus=auc_median - auc_lower_bound,
                                                  visible=True, 
                                                  thickness=3, 
                                                  width=0, 
                                                  color=roi_color),
                                      marker=dict(color=roi_color,
                                                  symbol='square',
                                                  size=8, 
                                                  line=dict(color=roi_color, 
                                                            width=3)), 
                                      showlegend=False), 
                          row=1, col=l + 1)
            
            # Add legend
            annotation = go.layout.Annotation(x=1, y=max_auc-(j*0.1*max_auc), text=roi, xanchor='left',
                                              showarrow=False, font_color=roi_color, 
                                              font_family=template_specs['font'],
                                              font_size=template_specs['axes_font_size'],
                                             )
            fig.add_annotation(annotation, row=1, col=l+1)

        # Set axis
        fig.update_xaxes(title_text='pRF size (dva)', range=[size_axis[0], size_axis[1]], showline=True)
        fig.update_yaxes(title_text='nCSF auc', range=[auc_axis[0], auc_axis[1]], showline=True)
        fig.update_layout(height=fig_height, 
                          width=fig_width, 
                          showlegend=False, 
                          template=fig_template,
                          margin_l=fig_margin[0], 
                          margin_t=fig_margin[1], 
                          margin_r=fig_margin[2], 
                          margin_b=fig_margin[3]
                         )
        
    return fig

def remove_second_page(pdf_path):
    """
    Remove the second page from a PDF file if it exists.
    Parameters
    ----------
    pdf_path : str
        Path to the PDF file.
    """
    import pikepdf
    with pikepdf.open(pdf_path, allow_overwriting_input=True) as pdf:
        if len(pdf.pages) > 1:
            del pdf.pages[1]  # Remove the second page
        pdf.save(pdf_path)
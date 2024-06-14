# Figure imports
import plotly.graph_objects as go
import plotly.express as px
from plotly.subplots import make_subplots
import numpy as np
import pandas as pd
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
    import plotly.graph_objects as go
    fig_template=go.layout.Template()

    # Violin plots
    fig_template.data.violin = [go.Violin(
                                    box_visible=False,
                                    points=False,
                                    # opacity=1,
                                    line_color= "rgba(0, 0, 0, 1)",
                                    line_width=template_specs['plot_width'],
                                    width=0.8,
                                    #marker_symbol='x',
                                    #marker_opacity=1,
                                    hoveron='violins',
                                    meanline_visible=False,
                                    # meanline_color="rgba(0, 0, 0, 1)",
                                    # meanline_width=template_specs['plot_width'],
                                    showlegend=False,
                                    )]

    # Barpolar
    fig_template.data.barpolar = [go.Barpolar(
                                    marker_line_color="rgba(0,0,0,1)",
                                    marker_line_width=template_specs['plot_width'], 
                                    showlegend=False, 
                                    )]
    # Pie plots
    fig_template.data.pie = [go.Pie(textposition=["inside","none"],
                                    # marker_line_color=['rgba(0,0,0,1)','rgba(255,255,255,0)'],
                                    marker_line_width=0,#[template_specs['plot_width'],0],
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


def prf_roi_area(df_roi_area, fig_width, fig_height, roi_colors):
    """
    Make bar plots of each roi area and the corresponding significative area of pRF  
    
    Parameters
    ----------
    df_roi_area : dataframe for corresponding plot
    fig_width : figure width in pixels
    fig_height : figure height in pixels
    roi_colors : list of rgb colors for plotly
    
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
                          plot_width=1.5)
    
    # General figure settings
    fig_template = plotly_template(template_specs)
    
    # General settings
    fig = make_subplots(rows=1, 
                        cols=2, 
                        subplot_titles=['FDR threshold = 0.05', 'FDR threshold = 0.01'],
                       )
    
    # FDR 0.01 
    # All vertices
    fig.add_trace(go.Bar(x=df_roi_area.roi, 
                         y=df_roi_area.vert_area, 
                         text=(df_roi_area['ratio_corr_pvalue_5pt']*100).astype(int).astype(str) + '%',
                         textposition='outside',
                         textangle=-60,
                         showlegend=False, 
                         marker=dict(color=roi_colors, opacity=0.2)),
                 row=1, col=1)
 
    # Significant vertices
    fig.add_trace(go.Bar(x=df_roi_area.roi, 
                         y=df_roi_area.vert_area_corr_pvalue_5pt, 
                         showlegend=False, 
                         marker=dict(color=roi_colors)),
                 row=1, col=1)
    
    
    # FDR 0.01 
    # All vertices
    fig.add_trace(go.Bar(x=df_roi_area.roi, 
                         y=df_roi_area.vert_area, 
                         text=(df_roi_area['ratio_corr_pvalue_1pt']*100).astype(int).astype(str) + '%',
                         textposition='outside',
                         textangle=-60,
                         showlegend=False, 
                         marker=dict(color=roi_colors, opacity=0.1)),
                 row=1, col=2)
    
    # Significant vertices
    fig.add_trace(go.Bar(x=df_roi_area.roi, 
                         y=df_roi_area.vert_area_corr_pvalue_1pt,
                         showlegend=False, 
                         marker=dict(color=roi_colors)),
                 row=1, col=2)

    # Define parameters
    fig.update_xaxes(showline=True, 
                     ticklen=0, 
                     linecolor=('rgba(255,255,255,0)'))      
    
    fig.update_yaxes(range=[0,100], 
                     showline=True, 
                     nticks=10, 
                     title_text='Surface area (cm<sup>2</sup>)',secondary_y=False)
    
    fig.update_layout(barmode='overlay',
                      height=fig_height, 
                      width=fig_width, 
                      template=fig_template,
                      margin_l=100, 
                      margin_r=50, 
                      margin_t=50, 
                      margin_b=50,
                     )

    # Return outputs
    return fig

def prf_violins_plot(df_violins, fig_width, fig_height, rois, roi_colors):
    """
    Make violins plots for pRF loo_r2, size, n and pcm

    Parameters
    ----------
    df_violins : dataframe
    fig_width : figure width in pixels
    fig_height : figure height in pixels
    rois : list of rois
    roi_colors : list of rgb colors for plotly
    
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
                          plot_width=1.5)
    
    # General figure settings
    fig_template = plotly_template(template_specs)

    rows, cols = 2,2
    fig = make_subplots(rows=rows, 
                        cols=cols, 
                        print_grid=False, 
                        vertical_spacing=0.08, 
                        horizontal_spacing=0.05)


    for j, roi in enumerate(rois):
        
        df = df_violins.loc[(df_violins.roi == roi)]
        
        # pRF loo r2
        fig.add_trace(go.Violin(x=df.roi[df.roi==roi], 
                                y=df.prf_loo_r2, 
                                name=roi, 
                                opacity=1,
                                showlegend=False, 
                                legendgroup='loo', 
                                points=False, 
                                spanmode='manual', 
                                span=[0, 1],
                                scalemode='width', 
                                fillcolor=roi_colors[j],
                                line_color=roi_colors[j]), 
                      row=1, col=1)
                
        # pRF size
        fig.add_trace(go.Violin(x=df.roi[df.roi==roi], 
                                y=df.prf_size, 
                                name=roi, 
                                opacity=1,
                                showlegend=False, 
                                legendgroup='size', 
                                points=False, 
                                spanmode='manual', 
                                span=[0, 20],
                                scalemode='width', 
                                fillcolor=roi_colors[j],
                                line_color=roi_colors[j]), 
                      row=1, col=2)
        
        # # pRF n
        # fig.add_trace(go.Violin(x=df.roi[df.roi==roi], 
        #                         y=df.prf_n, 
        #                         name=roi, 
        #                         opacity=1,
        #                         showlegend=False, 
        #                         legendgroup='n', 
        #                         points=False, 
        #                         scalemode='width', 
        #                         fillcolor=roi_colors[j],
        #                         line_color=roi_colors[j]), 
        #               row=2, col=1)

        # pRF ecc
        fig.add_trace(go.Violin(x=df.roi[df.roi==roi], 
                                y=df.prf_ecc, 
                                name=roi, 
                                opacity=1,
                                showlegend=False, 
                                legendgroup='n', 
                                points=False, 
                                spanmode='manual', 
                                span=[0, 20],
                                scalemode='width', 
                                fillcolor=roi_colors[j],
                                line_color=roi_colors[j]), 
                      row=2, col=1)
        
        # pcm
        fig.add_trace(go.Violin(x=df.roi[df.roi==roi], 
                                y=df.pcm_median, 
                                name=roi, 
                                opacity=1,
                                showlegend=False, 
                                legendgroup='pcm', 
                                points=False, 
                                spanmode='manual', 
                                span=[0, 50],
                                scalemode='width', 
                                fillcolor=roi_colors[j],
                                line_color=roi_colors[j]), 
                      row=2, col=2)
        
        # Set axis titles only for the left-most column and bottom-most row
        fig.update_yaxes(showline=True, 
                         range=[0, 1],
                         nticks=10, 
                         title_text='pRF LOO R<sup>2</sup>',
                         row=1, col=1)
        
        fig.update_yaxes(showline=True, 
                         range=[0, 20], 
                         nticks=5, 
                         title_text='pRF size (dva)', 
                         row=1, col=2)
        
        # fig.update_yaxes(showline=True, 
        #                  range=[0, 2], 
        #                  nticks=5, 
        #                  title_text='pRF n', 
        #                  row=2, col=1)

        fig.update_yaxes(showline=True, 
                         range=[0, 20], 
                         nticks=5, 
                         title_text='pRF eccentricity (dva)', 
                         row=2, col=1)
        
        fig.update_yaxes(showline=True, 
                         range=[0, 50],
                         nticks=10, 
                         title_text='pRF pCM (mm/dva)', 
                         row=2, col=2)
        
        fig.update_xaxes(showline=True, 
                         ticklen=0, 
                         linecolor=('rgba(255,255,255,0)'))
        
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
                      margin_l=100, 
                      margin_r=50, 
                      margin_t=100, 
                      margin_b=100)

    return fig

def prf_params_avg_plot(df_params_avg, fig_width, fig_height, rois, roi_colors):
    """
    Make parameters average plots for pRF loo_r2, size, n and pcm

    Parameters
    ----------
    df_params_avg : dataframe
    fig_width : figure width in pixels
    fig_height : figure height in pixels
    rois : list of rois
    roi_colors : list of rgb colors for plotly
    
    Returns
    -------
    fig : parameters average plot
    """
    
    # General figure settings
    template_specs = dict(axes_color="rgba(0, 0, 0, 1)",
                          axes_width=2,
                          axes_font_size=15,
                          bg_col="rgba(255, 255, 255, 1)",
                          font='Arial',
                          title_font_size=15,
                          plot_width=1.5)
    
    # General figure settings
    fig_template = plotly_template(template_specs)

    rows, cols = 2,2
    fig = make_subplots(rows=rows, 
                        cols=cols, 
                        print_grid=False, 
                        vertical_spacing=0.08, 
                        horizontal_spacing=0.05)
    
    for j, roi in enumerate(rois):
        
        df = df_params_avg.loc[(df_params_avg.roi == roi)]

        weighted_median = df.prf_loo_r2_weighted_median
        ci_up = df.prf_loo_r2_ci_up
        ci_down = df.prf_loo_r2_ci_down
        

        fig.add_trace(go.Scatter(x=[roi],
                                 y=tuple(weighted_median),
                                 mode='markers', 
                                 name=roi,
                                 error_y=dict(type='data', 
                                              array=[ci_up-weighted_median], 
                                              arrayminus=[weighted_median-ci_down],
                                              visible=True, 
                                              thickness=3,
                                              width=0, 
                                              color=roi_colors[j]),
                                 marker=dict(symbol="square",
                                             color=roi_colors[j],
                                             size=12, 
                                             line=dict(color=roi_colors[j], 
                                                       width=3)),
                                 legendgroup='loo',
                                 showlegend=False), 
                          row=1, col=1)
        
        # pRF size
        weighted_median = df.prf_size_weighted_median
        ci_up = df.prf_size_ci_up
        ci_down = df.prf_size_ci_down
        
        fig.add_trace(go.Scatter(x=[roi],
                                 y=tuple(weighted_median),
                                 mode='markers', 
                                 name=roi,
                                 error_y=dict(type='data', 
                                              array=[ci_up-weighted_median], 
                                              arrayminus=[weighted_median-ci_down],
                                              visible=True, 
                                              thickness=3,
                                              width=0, 
                                              color=roi_colors[j]),
                                 marker=dict(symbol="square",
                                             color=roi_colors[j],
                                             size=12, 
                                             line=dict(color=roi_colors[j], 
                                                       width=3)),
                                 legendgroup='size',
                                 showlegend=False), 
                          row=1, col=2)
                
        # # pRF n
        # weighted_median = df.prf_n_weighted_median
        # ci_up = df.prf_n_ci_up
        # ci_down = df.prf_n_ci_down
        
        # fig.add_trace(go.Scatter(x=[roi],
        #                          y=tuple(weighted_median),
        #                          mode='markers', 
        #                          name=roi,
        #                          error_y=dict(type='data', 
        #                                       array=[ci_up-weighted_median], 
        #                                       arrayminus=[weighted_median-ci_down],
        #                                       visible=True, 
        #                                       thickness=3,
        #                                       width=0, 
        #                                       color=roi_colors[j]),
        #                          marker=dict(symbol="square",
        #                                      color=roi_colors[j],
        #                                      size=12, 
        #                                      line=dict(color=roi_colors[j], 
        #                                                width=3)),
        #                          legendgroup='n',
        #                          showlegend=False), 
        #                   row=2, col=1)
        
        # pRF ecc
        weighted_median = df.prf_ecc_weighted_median
        ci_up = df.prf_ecc_ci_up
        ci_down = df.prf_ecc_ci_down
        
        fig.add_trace(go.Scatter(x=[roi],
                                 y=tuple(weighted_median),
                                 mode='markers', 
                                 name=roi,
                                 error_y=dict(type='data', 
                                              array=[ci_up-weighted_median], 
                                              arrayminus=[weighted_median-ci_down],
                                              visible=True, 
                                              thickness=3,
                                              width=0, 
                                              color=roi_colors[j]),
                                 marker=dict(symbol="square",
                                             color=roi_colors[j],
                                             size=12, 
                                             line=dict(color=roi_colors[j], 
                                                       width=3)),
                                 legendgroup='ecc',
                                 showlegend=False), 
                          row=2, col=1)
        
        # pcm
        weighted_median = df.pcm_median_weighted_median
        ci_up = df.pcm_median_ci_up
        ci_down = df.pcm_median_ci_down
        
        fig.add_trace(go.Scatter(x=[roi],
                                 y=tuple(weighted_median),
                                 mode='markers', 
                                 name=roi,
                                 error_y=dict(type='data', 
                                              array=[ci_up-weighted_median], 
                                              arrayminus=[weighted_median-ci_down],
                                              visible=True, 
                                              thickness=3,
                                              width=0, 
                                              color=roi_colors[j]),
                                 marker=dict(symbol="square",
                                             color=roi_colors[j],
                                             size=12, 
                                             line=dict(color=roi_colors[j], 
                                                       width=3)),
                                 legendgroup='pcm',
                                 showlegend=False), 
                          row=2, col=2)

        
        # Set axis titles only for the left-most column and bottom-most row
        fig.update_yaxes(showline=True, 
                         range=[0, 1],
                         nticks=10, 
                         title_text='pRF LOO R<sup>2</sup>',
                         row=1, col=1)
        
        fig.update_yaxes(showline=True, 
                         range=[0, 20], 
                         nticks=5, 
                         title_text='pRF size (dva)', 
                         row=1, col=2)
        
        # fig.update_yaxes(showline=True, 
        #                  range=[0, 2], 
        #                  nticks=5, 
        #                  title_text='pRF n', 
        #                  row=2, col=1)

        fig.update_yaxes(showline=True, 
                         range=[0, 20], 
                         nticks=5, 
                         title_text='pRF eccentricity (dva)', 
                         row=2, col=1)
        
        fig.update_yaxes(showline=True, 
                         range=[0, 20], 
                         nticks=10, 
                         title_text='pRF pCM (mm/dva)', 
                         row=2, col=2)
        
        fig.update_xaxes(showline=True, 
                         ticklen=0, 
                         linecolor=('rgba(255,255,255,0)'))
        
    fig.update_layout(height=fig_height, 
                      width=fig_width, 
                      legend=dict(orientation="h", 
                                  font_family=template_specs['font'],
                                  font_size=template_specs['axes_font_size'],
                                  y=1.1, 
                                  yanchor='top', 
                                  xanchor='left', 
                                  traceorder='normal', 
                                  itemwidth=30), 
                      template=fig_template,
                      margin_l=100, 
                      margin_r=50, 
                      margin_t=100, 
                      margin_b=100)

    return fig

def prf_ecc_size_plot(df_ecc_size, fig_width, fig_height, rois, roi_colors, plot_groups, max_ecc):
    """
    Make scatter plot for linear relationship between eccentricity and size

    Parameters
    ----------
    df_ecc_size : A data dataframe
    fig_width : figure width in pixels
    fig_height : figure height in pixels
    rois : list of rois
    roi_colors : list of rgb colors for plotly
    plot_groups : groups of roi to plot together
    max_ecc : maximum eccentricity 
    
    Returns
    -------
    fig : eccentricy as a function of size plot
    """
    
    from maths_utils import weighted_regression

    # General figure settings
    template_specs = dict(axes_color="rgba(0, 0, 0, 1)",
                          axes_width=2,
                          axes_font_size=15,
                          bg_col="rgba(255, 255, 255, 1)",
                          font='Arial',
                          title_font_size=15,
                          plot_width=1.5)
    
    # General figure settings
    fig_template = plotly_template(template_specs)

    # General settings
    rows, cols = 1, len(plot_groups)
    fig = make_subplots(rows=rows, cols=cols, print_grid=False)
    
    for l, line_label in enumerate(plot_groups):
        for j, roi in enumerate(line_label):
            
            # Parametring colors
            roi_color = roi_colors[j + l * 3]
            roi_color_opac = f"rgba{roi_color[3:-1]}, 0.15)"
            
            # Get data
            df = df_ecc_size.loc[(df_ecc_size.roi == roi)]
            ecc_median = np.array(df.prf_ecc_bins)
            size_median = np.array(df.prf_size_bins_median)
            r2_median = np.array(df.prf_loo_r2_bins_median)
            size_upper_bound = np.array(df.prf_size_bins_ci_upper_bound)
            size_lower_bound = np.array(df.prf_size_bins_ci_lower_bound)
            
            # Linear regression
            slope, intercept = weighted_regression(ecc_median, size_median, r2_median, model='linear')
            slope_upper, intercept_upper = weighted_regression(ecc_median[np.where(~np.isnan(size_upper_bound))], 
                                                               size_upper_bound[~np.isnan(size_upper_bound)], 
                                                               r2_median[np.where(~np.isnan(size_upper_bound))], 
                                                               model='linear')
            slope_lower, intercept_lower = weighted_regression(ecc_median[np.where(~np.isnan(size_lower_bound))], 
                                                               size_lower_bound[~np.isnan(size_lower_bound)], 
                                                               r2_median[np.where(~np.isnan(size_lower_bound))], 
                                                               model='linear')

            line_x = np.linspace(ecc_median[0], ecc_median[-1], 50)
            line = slope * line_x + intercept
            line_upper = slope_upper * line_x + intercept_upper
            line_lower = slope_lower * line_x + intercept_lower

            fig.add_trace(go.Scatter(x=line_x, y=line, mode='lines', name=roi, legendgroup=roi, 
                                      line=dict(color=roi_color, width=3), showlegend=False), 
                          row=1, col=l+1)

            # Error area
            fig.add_trace(go.Scatter(x=np.concatenate([line_x, line_x[::-1]]), 
                                      y=np.concatenate([list(line_upper), list(line_lower[::-1])]), 
                                      mode='lines', fill='toself', fillcolor=roi_color_opac, 
                                      line=dict(color=roi_color_opac, width=0), showlegend=False), 
                          row=1, col=l+1)

            # Markers
            fig.add_trace(go.Scatter(x=ecc_median, 
                                     y=size_median, mode='markers', 
                                     error_y=dict(type='data', 
                                                  array=size_upper_bound - size_median, 
                                                  arrayminus=size_median - size_lower_bound,
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
            annotation = go.layout.Annotation(x=1, y=max_ecc-j*1.5, text=roi, xanchor='left',
                                              showarrow=False, font_color=roi_color, 
                                              font_family=template_specs['font'],
                                              font_size=template_specs['axes_font_size'],
                                             )
            fig.add_annotation(annotation, row=1, col=l+1)

        # Set axis titles only for the left-most column and bottom-most row
        fig.update_yaxes(title_text='pRF size (dva)', row=1, col=1)
        fig.update_xaxes(title_text='pRF eccentricity (dva)', range=[0, max_ecc], showline=True, row=1, col=l+1)
        fig.update_yaxes(range=[0, max_ecc], showline=True)
        fig.update_layout(height=fig_height, width=fig_width, showlegend=False, template=fig_template,
                         margin_l=100, margin_r=50, margin_t=50, margin_b=100)
        
    return fig

def prf_ecc_pcm_plot(df_ecc_pcm, fig_width, fig_height, rois, roi_colors, plot_groups, max_ecc):
    """
    Make scatter plot for relationship between eccentricity and pCM

    Parameters
    ----------
    df_ecc_pcm : dataframe for the plot
    fig_width : figure width in pixels
    fig_height : figure height in pixels
    rois : list of rois
    roi_colors : list of rgb colors for plotly
    plot_groups : groups of roi to plot together
    max_ecc : maximum eccentricity
    
    Returns
    -------
    fig : eccentricy as a function of pcm plot
    """

    from maths_utils import weighted_regression
    
    # General figure settings
    template_specs = dict(axes_color="rgba(0, 0, 0, 1)",
                          axes_width=2,
                          axes_font_size=15,
                          bg_col="rgba(255, 255, 255, 1)",
                          font='Arial',
                          title_font_size=15,
                          plot_width=1.5)
    
    # General figure settings
    fig_template = plotly_template(template_specs)

    # General settings
    rows, cols = 1, len(plot_groups)
    fig = make_subplots(rows=rows, cols=cols, print_grid=False)
    
    for l, line_label in enumerate(plot_groups):
        for j, roi in enumerate(line_label):

            # Parametring colors
            roi_color = roi_colors[j + l * 3]
            roi_color_opac = f"rgba{roi_color[3:-1]}, 0.15)"
            
            # Get data
            df = df_ecc_pcm.loc[(df_ecc_pcm.roi == roi)]
            ecc_median = np.array(df.prf_ecc_bins)
            pcm_median = np.array(df.prf_pcm_bins_median)
            r2_median = np.array(df.prf_loo_r2_bins_median)
            pcm_upper_bound = np.array(df.prf_pcm_bins_ci_upper_bound)
            pcm_lower_bound = np.array(df.prf_pcm_bins_ci_lower_bound)
            
            # Linear regression
            slope, intercept = weighted_regression(ecc_median, pcm_median, r2_median, model='pcm')
            
            slope_upper, intercept_upper = weighted_regression(ecc_median[~np.isnan(pcm_upper_bound)], 
                                                               pcm_upper_bound[~np.isnan(pcm_upper_bound)], 
                                                               r2_median[~np.isnan(pcm_upper_bound)], 
                                                               model='pcm')
            
            slope_lower, intercept_lower = weighted_regression(ecc_median[~np.isnan(pcm_lower_bound)], 
                                                               pcm_lower_bound[~np.isnan(pcm_lower_bound)], 
                                                               r2_median[~np.isnan(pcm_lower_bound)], 
                                                               model='pcm')

            line_x = np.linspace(ecc_median[0], ecc_median[-1], 50)
            line = 1 / (slope * line_x + intercept)
            line_upper = 1 / (slope_upper * line_x + intercept_upper)
            line_lower = 1 / (slope_lower * line_x + intercept_lower)

            fig.add_trace(go.Scatter(x=line_x, 
                                     y=line, 
                                     mode='lines', 
                                     name=roi, 
                                     legendgroup=roi, 
                                     line=dict(color=roi_color, width=3), 
                                     showlegend=False), 
                          row=1, col=l+1)

            # Error area
            fig.add_trace(go.Scatter(x=np.concatenate([line_x, line_x[::-1]]),
                                      y=np.concatenate([list(line_upper), list(line_lower[::-1])]), 
                                      mode='lines', fill='toself', fillcolor=roi_color_opac, 
                                      line=dict(color=roi_color_opac, width=0), showlegend=False), 
                          row=1, col=l+1)

            # Markers
            fig.add_trace(go.Scatter(x=ecc_median, 
                                     y=pcm_median, 
                                     mode='markers', 
                                     error_y=dict(type='data', 
                                                  array=pcm_upper_bound - pcm_median, 
                                                  arrayminus=pcm_median - pcm_lower_bound,
                                                  visible=True, 
                                                  thickness=3, 
                                                  width=0, 
                                                  color=roi_color),
                                     marker=dict(color=roi_color, 
                                                 symbol='square',
                                                 size=8, line=dict(color=roi_color,
                                                                   width=3)), 
                                     showlegend=False), 
                          row=1, col=l + 1)
            
            # Add legend
            annotation = go.layout.Annotation(x=12, y=(max_ecc+15)-j*3, text=roi, xanchor='left',
                                              showarrow=False, font_color=roi_color, 
                                              font_family=template_specs['font'],
                                              font_size=template_specs['axes_font_size'],
                                             )
            fig.add_annotation(annotation, row=1, col=l+1)

        # Set axis titles only for the left-most column and bottom-most row
        fig.update_yaxes(title_text='pRF cortical magn. (mm/dva)', row=1, col=1)
        fig.update_xaxes(title_text='pRF eccentricity (dva)', range=[0, max_ecc], showline=True, row=1, col=l+1)
        fig.update_yaxes(range=[0, max_ecc + 15], showline=True)
        fig.update_layout(height=fig_height, width=fig_width, showlegend=False, template=fig_template,
                         margin_l=100, margin_r=50, margin_t=50, margin_b=100)
        
    return fig

def prf_polar_angle_plot(df_polar_angle, fig_width, fig_height, rois, roi_colors, num_polar_angle_bins) :    
    """
    Make polar angle distribution plots
    
    Parameters
    ----------
    df_polar_angle : polar angle dataframe
    fig_width : figure width in pixels
    fig_height : figure height in pixels
    rois : list of rois
    roi_colors : list of rgb colors for plotly
    num_bins : bins for the polar angle 
     
    Returns
    -------
    figs : a list of three figures
    hemispheres : a list of corresponding hemispheres
    """
    # General figure settings
    template_specs = dict(axes_color="rgba(0, 0, 0, 1)",
                          axes_width=2,
                          axes_font_size=15,
                          bg_col="rgba(255, 255, 255, 1)",
                          font='Arial',
                          title_font_size=15,
                          plot_width=1.5)
    
    # General figure settings
    fig_template = plotly_template(template_specs)

    # General settings
    rows, cols = 1, len(rois)
    specs = [[{'type': 'polar'}] * cols]
    
    figs = []
    hemispheres = []
    hemis = ['hemi-L', 'hemi-R', 'hemi-LR']
    for i, hemi in enumerate(hemis):
        fig = make_subplots(rows=rows, cols=cols, print_grid=False, specs=specs)
            
        for j, roi in enumerate(rois):
            if j == 0: showlegend = True
            else: showlegend = False
    
            # Parts of polar angles and number of voxels in each part
            df = df_polar_angle.loc[(df_polar_angle.roi==roi) & (df_polar_angle.hemi==hemi)]
            
            # barpolar
            fig.add_trace(go.Barpolar(r=df.loo_r2_sum, 
                                      theta=df.theta_slices, 
                                      marker_color=roi_colors[j], 
                                      width=360/(num_polar_angle_bins),
                                      marker_line_color='white', 
                                      marker_line_width=3, 
                                      opacity=1,
                                      showlegend=True,
                                      name=roi, 
                                     ), 
                          row=1, col=j+1)
    
        # Define parameters
        fig.update_polars(angularaxis=dict(visible=False), 
                          radialaxis=dict(visible=False))
        
        fig.update_layout(height=fig_height, 
                          width=fig_width, 
                          legend=dict(orientation="h", 
                                  font_family=template_specs['font'],
                                  font_size=template_specs['axes_font_size'],
                                  y=1.1, 
                                  yanchor='top', 
                                  xanchor='left', 
                                  traceorder='normal', 
                                  itemwidth=30), 
                          template=fig_template,
                          margin_l=50, 
                          margin_r=50, 
                          margin_t=50, 
                          margin_b=50)
                          
        figs.append(fig)
        hemispheres.append(hemi)
        
    return figs, hemispheres

def prf_contralaterality_plot(df_contralaterality, fig_height, fig_width, rois, roi_colors):
    """
    Make contralaterality pie plot
    
    Parameters
    ----------
    df_contralaterality : dataframe
    fig_width : figure width in pixels
    fig_height : figure height in pixels
    rois : list of rois
    roi_colors : list of rgb colors for plotly
     
    Returns
    -------
    fig : contralaterality figure
    """

    # General figure settings
    template_specs = dict(axes_color="rgba(0, 0, 0, 1)",
                          axes_width=2,
                          axes_font_size=15,
                          bg_col="rgba(255, 255, 255, 1)",
                          font='Arial',
                          title_font_size=15,
                          plot_width=1.5)
    
    # General figure settings
    fig_template = plotly_template(template_specs)

    # General settings
    rows, cols = 1, len(rois)
    specs = [[{'type': 'pie'}] * cols]    
    fig = make_subplots(rows=rows, cols=cols, print_grid=False, specs=specs)
    
    for j, roi in enumerate(rois):

        df = df_contralaterality.loc[df_contralaterality.roi==roi]
        percentage_total = np.array(df.contralaterality_prct)
        percentage_rest = 1 - percentage_total
        percentage_total = percentage_total.tolist()
        percentage_rest = percentage_rest.tolist()
        values = [percentage_total[0], percentage_rest[0]]


        fig.add_trace(go.Pie(values=values,
                             marker=dict(colors=[roi_colors[j], 'white'])
                            ),
                      row=1, col=j+1)

    # Define parameters
    fig.update_layout(height=fig_height, 
                      width=fig_width, 
                      showlegend=False,
                      template=fig_template,
                      margin_l=50, 
                      margin_r=50, 
                      margin_t=50, 
                      margin_b=50)
    
    return fig 

def prf_distribution_plot(df_distribution, fig_height, fig_width, rois, roi_colors, screen_side):
    """
    Make prf distribution contour plot
    
    Parameters
    ----------
    df_distribution : dataframe
    fig_width : figure width in pixels
    fig_height : figure height in pixels
    rois : list of rois
    roi_colors : list of rgb colors for plotly
    screen_side: mesh screen side (square) im dva (e.g. 20 dva from -10 to 10 dva)
     
    Returns
    -------
    fig : distribution figure
    """
    
    # Template settings
    template_specs = dict(axes_color="rgba(0, 0, 0, 1)",
                          axes_width=2,
                          axes_font_size=15,
                          bg_col="rgba(255, 255, 255, 1)",
                          font='Arial',
                          title_font_size=15,
                          plot_width=1.5)
    fig_template = plotly_template(template_specs)
    
    # General figure settings
    rows, cols = 1, len(rois)
    line_width = 1
    contour_width = 0.5
    
    figs = []
    hemispheres = []
    hemis = ['hemi-L', 'hemi-R', 'hemi-LR']
    for i, hemi in enumerate(hemis):  
        fig = make_subplots(rows=rows ,cols=cols)
        for j, roi in enumerate(rois) :
            # Make df roi
            df_roi = df_distribution.loc[(df_distribution.roi == roi) & (df_distribution.hemi == hemi)]

            # make the two dimensional mesh for z dimension
            gauss_z_tot = df_roi.drop(columns=['roi', 'x', 'y', 'hemi']).values
            
            # Contour plot
            fig.add_trace(go.Contour(x=df_roi.x, 
                                     y=df_roi.y, 
                                     z=gauss_z_tot, 
                                     colorscale=[[0, 'white'],[0.25, 'white'], [1, roi_colors[j]]],  
                                     showscale=False,  
                                     line=dict(color='black', width=contour_width),  
                                     # ncontours=6, 
                                     contours=dict(coloring='fill', 
                                                   start=0, 
                                                   end=1, 
                                                   size=0.25, 
                                                   showlines=True)
                                     ),row=1, col=j+1)
            
            # x line
            fig.add_trace(go.Scatter(x=[0,0],
                                     y=[-20,20],
                                     mode='lines',
                                     line=dict(dash='2px',color='rgba(0, 0, 0, 0.6)', width=line_width)
                                    ),row=1, col=j+1)
            # y line
            fig.add_trace(go.Scatter(x=[-20,20], 
                                     y=[0,0], 
                                     mode='lines', 
                                     line=dict(dash='2px',color='rgba(0, 0, 0, 0.6)', width=line_width)),row=1, col=j+1)
            
            # # square
            # fig.add_shape(type="rect", 
            #               x0=-10, 
            #               y0=-10, 
            #               x1=10, 
            #               y1=10, 
            #               line=dict(dash='2px',color='black', width=line_width),row=1, col=j+1)
            
        fig.update_xaxes(range=[-20,20], color= ('rgba(255,255,255,0)'))
        fig.update_yaxes(range=[-20,20], color= ('rgba(255,255,255,0)'))
        
        # Define parameters
        fig.update_layout(height=fig_height, 
                          width=fig_width, 
                          showlegend=False,
                          template=fig_template,
                          margin_l=10, 
                          margin_r=10, 
                          margin_t=100, 
                          margin_b=100)
        figs.append(fig)
        hemispheres.append(hemi)
    return figs, hemispheres

def prf_barycentre_plot(df_barycentre, fig_height, fig_width, rois, roi_colors, screen_side):
    """
    Make prf barycentre plot
    
    Parameters
    ----------
    df_barycentre : dataframe
    fig_width : figure width in pixels
    fig_height : figure height in pixels
    rois : list of rois
    roi_colors : list of rgb colors for plotly
    screen_side: mesh screen side (square) im dva (e.g. 20 dva from -10 to 10 dva)
     
    Returns
    -------
    fig : barycentre figure
    """
    
    # Template settings
    template_specs = dict(axes_color="rgba(0, 0, 0, 1)",
                          axes_width=2,
                          axes_font_size=15,
                          bg_col="rgba(255, 255, 255, 1)",
                          font='Arial',
                          title_font_size=15,
                          plot_width=1.5)
    fig_template = plotly_template(template_specs)
    
    # General figure settings
    line_width = 1
    fig = go.Figure()
    hemis = ['hemi-L', 'hemi-R']
    for i, hemi in enumerate(hemis): 
        if hemi=='hemi-L': symbol, showlegend = 'square' , True
        elif hemi=='hemi-R': symbol, showlegend = 'circle' , False
        for j, roi in enumerate(rois) :
            # Make df roi
            df_roi = df_barycentre.loc[(df_barycentre.roi == roi) & (df_barycentre.hemi == hemi)]    
    
            # barycentre position
            fig.add_trace(go.Scatter(x=df_roi.barycentre_x, 
                                     y=df_roi.barycentre_y, 
                                     mode='markers', 
                                     name = roi,
                                     marker=dict(symbol=symbol, 
                                                 color=roi_colors[j], 
                                                 size=12),
                                     error_x=dict(type='data', 
                                                  array=[df_roi.upper_ci_x - df_roi.barycentre_x], 
                                                  arrayminus=[df_roi.barycentre_x - df_roi.lower_ci_x],
                                                  visible=True, 
                                                  thickness=3, 
                                                  width=0, 
                                                  color=roi_colors[j]),
                                     error_y=dict(type='data', 
                                                  array=[df_roi.upper_ci_y - df_roi.barycentre_y], 
                                                  arrayminus=[df_roi.barycentre_y - df_roi.lower_ci_y],
                                                  visible=True, 
                                                  thickness=3, 
                                                  width=0, 
                                                  color=roi_colors[j]),
                                 showlegend=showlegend))
        # Center lignes
        fig.add_trace(go.Scatter(x=[0,0], 
                                 y=[-screen_side, screen_side], 
                                 mode='lines', 
                                 showlegend=False, 
                                 line=dict(dash='2px',color='grey', width=line_width)))
        
        fig.add_trace(go.Scatter(x=[-screen_side,screen_side], 
                                 y=[0,0], 
                                 mode='lines', 
                                 showlegend=False,
                                 line=dict(dash='2px',color='grey', width=line_width)))
        
        # Add squares 
        for position in [2,4,6,8,10]:
            fig.add_shape(type="rect", 
                          x0=-position, 
                          y0=-position, 
                          x1=position, 
                          y1=position, 
                          line=dict(dash='2px',color='grey', width=line_width))
        # Add annotations 
        fig.add_trace(go.Scatter(x=[0, 0, 0, 0, 0], 
                                 y=[2.2, 4.2, 6.2, 8.2, 10.2], 
                                 showlegend=False, 
                                 text=["2 dva", 
                                       "4 dva", 
                                       "6 dva",
                                       "8 dva", 
                                       "10 dva"], 
                                 mode="text", 
                                 textfont=dict(size=10)))
    
    fig.update_yaxes(range=[-12,12],color= ('rgba(255,255,255,0)'))
    fig.update_xaxes(range=[-12,12],color= ('rgba(255,255,255,0)'))

    # Define parameters
    fig.update_layout(height=fig_height, 
                      width=fig_width, 
                      showlegend=True,
                      template=fig_template,
                      margin_l=570, 
                      margin_r=570, 
                      margin_t=50, 
                      margin_b=50)
        
    return fig

def categories_proportions_roi_plot(data, subject, fig_height, fig_width):
    data = data.copy()
    filtered_data = data[data['stats_final'] != 'non_responding']
    
    # Sort categories
    categories_order = ['vision', 'vision_and_pursuit_and_saccade', 'pursuit_and_saccade', 'vision_and_saccade', 'vision_and_pursuit', 'saccade', 'pursuit']
    filtered_data['stats_final'] = pd.Categorical(filtered_data['stats_final'], categories=categories_order, ordered=True)
    filtered_data = filtered_data.sort_values(['rois', 'stats_final'])

    
    #  Defines colors settings 
    roi_colors = px.colors.sequential.Sunset[:4] + px.colors.sequential.Rainbow[:]
    stats_categories_colors = list(reversed(px.colors.qualitative.D3))[2:]
    
    # To write make the percent visible only for choose categories 
    percent_color =  {'pursuit': 'rgba(255,255,255,0)', 
                      'saccade': 'rgba(255,255,255,0)', 
                      'pursuit_and_saccade': 'rgba(0, 0, 0, 1)', 
                      'vision': 'rgba(0, 0, 0, 1)', 
                      'vision_and_pursuit': 'rgba(255,255,255,0)', 
                      'vision_and_saccade': 'rgba(255,255,255,0)', 
                      'vision_and_pursuit_and_saccade': 'rgba(0, 0, 0, 1)'}
    
    categorie_color_map = {'pursuit': 'rgba(255,255,255,0)', 
                           'saccade': 'rgba(255,255,255,0)', 
                           'pursuit_and_saccade': stats_categories_colors[3], 
                           'vision': stats_categories_colors[4], 
                           'vision_and_pursuit': 'rgba(255,255,255,0)', 
                           'vision_and_saccade': 'rgba(255,255,255,0)', 
                           'vision_and_pursuit_and_saccade': stats_categories_colors[7]}
    
    rois = pd.unique(data.rois)
    #  Make the subplot
    # fig_height, fig_width = 300, 1920
    rows = 2 
    cols =len(rois)
    specs = [[{'type': 'domain'}] * cols,  [{'type': 'xy'}] * cols]
    
    
    fig = make_subplots(rows=rows, cols=cols, print_grid=False, specs=specs, row_heights=[1,0.05])
    
    
    for i, roi in enumerate(rois):
        df_rois = filtered_data.loc[filtered_data.rois == roi]
        #  Colors for categories 
        categorie_colors = [categorie_color_map[label] for label in df_rois.stats_final]
        #  Colors for the percentages 
        percentage_colors = [percent_color[label] for label in df_rois.stats_final]
        
        
        fig.add_trace(go.Pie(labels=df_rois.stats_final, 
                             values=df_rois.vertex_surf, 
                             showlegend=False, 
                             sort=False,
                             textinfo='percent',
                             textposition='inside',
                             direction='clockwise',
                             name= roi,    
                             marker=dict(colors=categorie_colors),
                             insidetextfont=dict(color=percentage_colors),
                             hole=0.3), 
                      row=1, col=i+1)
        
        fig.add_annotation(text=roi, 
                           yshift =10,
                           showarrow=False, 
                           font=dict(size=13,color=roi_colors[i]), 
                           row=2, col=i+1)
    
    
        
    fig.update_layout(height=fig_height, 
                      width=fig_width,
                      template='simple_white')  
    
    fig.update_xaxes(showline=True, 
                     ticklen=0, 
                     linecolor=('rgba(255,255,255,0)'), 
                     tickfont=dict(color='rgba(255,255,255,0)'))
    
    fig.update_yaxes(showline=True, 
                     ticklen=0, 
                     linecolor=('rgba(255,255,255,0)'), 
                     tickfont=dict(color='rgba(255,255,255,0)'))
    
    return fig

def surface_rois_categories_plot(data, subject, fig_height, fig_width):   
    data = data.copy()
    #  Defines colors settings 
    roi_colors = px.colors.sequential.Sunset[:4] + px.colors.sequential.Rainbow[:]
    stats_categories_colors = list(reversed(px.colors.qualitative.D3))[2:]
    
    
    
    
    categorie_color_map = {'non_responding': stats_categories_colors[0], 
                           'pursuit': stats_categories_colors[1], 
                           'saccade': stats_categories_colors[2], 
                           'pursuit_and_saccade': stats_categories_colors[3], 
                           'vision': stats_categories_colors[4], 
                           'vision_and_pursuit': stats_categories_colors[5], 
                           'vision_and_saccade': stats_categories_colors[6], 
                           'vision_and_pursuit_and_saccade': stats_categories_colors[7]}
    
    
    #  grpup df 
    group_df_rois = data.groupby(['rois'], sort=False)['vertex_surf'].sum().reset_index()
    group_rois_categories = data.groupby(['rois', 'stats_final'], sort=False)['vertex_surf'].sum().reset_index()
    
    #  Make subplot 
    # fig_height, fig_width = 1080, 1920
    fig = make_subplots(specs=[[{"secondary_y": True}]])

    
    
    fig.add_trace(go.Bar(x=group_df_rois.rois, 
                         y=group_df_rois.vertex_surf, 
                         showlegend=False, 
                         marker=dict(color=roi_colors, opacity=0.1)), 
                  secondary_y=False)
    
    #Choose categories to plot
    stats_categories = ['vision',  'vision_and_pursuit_and_saccade','pursuit_and_saccade']
    for categorie in stats_categories:
        df = group_rois_categories.loc[group_rois_categories.stats_final == categorie]
    
        fig.add_trace(go.Bar(x=df.rois, 
                             y=df.vertex_surf, 
                             name=categorie,  
                             legendgroup=categorie, 
                             marker_color=categorie_color_map[categorie]), 
                      secondary_y=True) 
    
    fig.update_layout(yaxis2=dict(overlaying='y',
                                  side='right',
                                  range=[0, 5000],  
                                  showticklabels=False, 
                                  ticklen=0, 
                                  linecolor=('rgba(255,255,255,0)')))
    
    
    
    fig.update_xaxes(showline=True, 
                     ticklen=0, 
                     linecolor=('rgba(255,255,255,0)'), 
                     tickfont=dict(size=12))      
    
    fig.update_yaxes(range=[0,5000], 
                     nticks=5, 
                     title_text='Surface in mm<sup>2</sup>',secondary_y=False)
    
    fig.update_layout(height=fig_height, 
                      width=fig_width,
                      barmode='stack',
                      showlegend=True, 
                      template='simple_white')  
    
    
    return fig 


def surface_rois_all_categories_plot(data, subject, fig_height, fig_width):  
    data = data.copy()
    #  Defines colors settings 
    stats_categories_colors = list(reversed(px.colors.qualitative.D3))[2:]
    
    
    
    
    categorie_color_map = {'non_responding': stats_categories_colors[0], 
                           'pursuit': stats_categories_colors[1], 
                           'saccade': stats_categories_colors[2], 
                           'pursuit_and_saccade': stats_categories_colors[3], 
                           'vision': stats_categories_colors[4], 
                           'vision_and_pursuit': stats_categories_colors[5], 
                           'vision_and_saccade': stats_categories_colors[6], 
                           'vision_and_pursuit_and_saccade': stats_categories_colors[7]}
    
    #  grpup df 
    group_df = data.groupby(['rois', 'stats_final'], sort=False)['vertex_surf'].sum().reset_index()
    
    #  Figure settings
    # fig_height, fig_width = 1080, 1920
    fig = go.Figure()
    
    
    stats_categories= ['non_responding', 'vision','vision_and_pursuit_and_saccade', 'pursuit_and_saccade', 'pursuit', 'saccade', 'vision_and_pursuit', 'vision_and_saccade']
    for categorie in stats_categories:
        df = group_df.loc[group_df.stats_final == categorie]
    
        fig.add_trace(go.Bar(x=df.rois, 
                             y=df.vertex_surf, 
                             name=categorie,  
                             legendgroup=categorie, 
                             marker_color=categorie_color_map[categorie])) 
    
    fig.update_xaxes(showline=True, 
                     ticklen=0, 
                     linecolor=('rgba(255,255,255,0)'), 
                     tickfont=dict(size=12))      
    
    fig.update_yaxes(range=[0,5000], 
                     nticks=5, 
                     title_text='Surface in mm<sup>2</sup>')
    
    fig.update_layout(height=fig_height, 
                      width=fig_width,
                      barmode='stack',
                      showlegend=True, 
                      template='simple_white')  
    
    return fig 
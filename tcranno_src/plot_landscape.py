from comut import comut
import palettable
import pandas as pd
import matplotlib.pyplot as plt
from numpy import *
from pandas import Series
import matplotlib.gridspec as gridspec
import argparse

class CoMut_tcranno(comut.CoMut): # modify the plot_comut function to hide the xticklabels
    def plot_comut(self, fig=None, spec=None, x_padding=0, y_padding=0,tri_padding=0, heights=None, hspace=0.2, 
                   subplot_hspace=None,widths=None, wspace=0.2, structure=None, figsize=(10,6)):
        if structure is None:
            structure = [[plot] for plot in self._plots]
        if heights is None:
            heights = {}
        num_subplots = len(structure)

        # get height structure based on input heights
        heights = self._get_height_spec(structure, heights)

        # calculate height of plots for gridspeccing. Heights are
        # reversed to match CoMut plotting (bottom to top)
        plot_heights = [sum(height) for height in heights][::-1]

        # create figure if none given
        if fig is None:
            fig = plt.figure(figsize=figsize)

        # make default widths and determine location of CoMut. If widths give,
        # just calculate location of CoMut
        if widths is None:
            widths, comut_idx = self._get_default_widths_and_comut_loc()
        else:
            _, comut_idx = self._get_default_widths_and_comut_loc()

        # number of cols is equal to size of widths
        num_cols = len(widths)

        # create gridspec if none given
        if spec is None:
            spec = gridspec.GridSpec(ncols=num_cols, nrows=num_subplots, figure=fig,
                                     height_ratios=plot_heights, width_ratios=widths,
                                     hspace=hspace, wspace=wspace)

        # otherwise, create gridspec in spec
        else:
            spec = gridspec.GridSpecFromSubplotSpec(ncols=num_cols, nrows=num_subplots,
                                                    height_ratios=plot_heights, width_ratios=widths,
                                                    hspace=hspace, wspace=wspace, subplot_spec=spec)

        # plot each plot in subplots
        for i, (plot, height) in enumerate(zip(structure, heights)):
            # subplots share an x axis with first plot
            if i == 0:
                sharex = None
                first_plot = plot[0]
            else:
                sharex = self.axes[first_plot]

            # if only one plot in subplot, just add subplot and plot
            if len(plot) == 1:
                plot_name = plot[0]
                ax = fig.add_subplot(spec[num_subplots - i - 1, comut_idx], sharex=sharex)
                ax = self._plot_data_on_axis(ax=ax, plot_name=plot_name, x_padding=x_padding, y_padding=y_padding, tri_padding=tri_padding)

                # extract all sideplots on this axis
                side_plots = self._side_plots[plot_name]

                # identify the locations of each sideplot and plot from inward -> outward
                left_idx, right_idx = 1, 1
                for side_name, side_plot in side_plots.items():
                    position = side_plot['position']
                    if position == 'left':
                        sideplot_idx = comut_idx - left_idx
                        left_idx += 1
                    elif position == 'right':
                        sideplot_idx = comut_idx + right_idx
                        right_idx += 1
                    side_ax = fig.add_subplot(spec[num_subplots - i - 1, sideplot_idx])
                    side_ax = self._plot_side_bar_data(side_ax, side_name, y_padding=y_padding,
                                                       **side_plot)
                    side_ax.set_ylim(ax.get_ylim())
            else:
                num_plots = len(plot)
                height = height[::-1]
                subplot_spec = gridspec.GridSpecFromSubplotSpec(ncols=1, nrows=num_plots,hspace=subplot_hspace, 
                        subplot_spec=spec[num_subplots - i - 1, comut_idx],height_ratios=height)
                for j, plot_name in enumerate(plot):
                    ax = fig.add_subplot(subplot_spec[num_plots - j - 1, 0], sharex=sharex)
                    ax = self._plot_data_on_axis(ax=ax, plot_name=plot_name, x_padding=x_padding, y_padding=y_padding, 
                                                 tri_padding=tri_padding)
                    if self._side_plots[plot_name]:
                        raise ValueError('Side bar plot for {} cannot be created. '
                                         'Plots within a subplot cannot have a side plot.'.format(plot_name))
        # add x axis labels to the bottom-most axis, make it visible
        self.axes[first_plot].set_xticks([])
        self.figure = fig
        return self

def plot_tcr2ept(tcr2tcr_df, items, outprefix, extension):
    items0 = items.copy()
    items = []
    for each in items0:
        ept_names = (each.split(' [')[0]).split(',')
        item = ept_names[0]
        items.append(item+' ['+each.split(' [')[1])
    #print(items)
    cm = tcr2tcr_df[tcr2tcr_df.Rank=='Complete_Match']
    pm = tcr2tcr_df[tcr2tcr_df['Index'].str.endswith('-1')]
    s = cm['Matched_Epitope(s) [Matched_Antigen->Matched_Organism]'].str.split('; ').apply(Series, 1).stack()
    s.index = s.index.droplevel(-1)
    s.name = 'category'
    del cm['Matched_Epitope(s) [Matched_Antigen->Matched_Organism]']
    cm = cm.join(s)
    cm_df=pd.DataFrame(list(zip(cm.Index.tolist(),cm.category.tolist(),['CM']*len(cm),cm.Frequency.tolist(),[0]*len(cm))),columns=['sample','category','value','CM','PM'])
    cm_df=cm_df[cm_df.category!='']
    cm_df.drop_duplicates(inplace=True)
    cm_df.sort_values(by=['CM','sample'], ascending=False, inplace=True)
    s = pm['Matched_Epitope(s) [Matched_Antigen->Matched_Organism]'].str.split('; ').apply(Series, 1).stack()
    s.index = s.index.droplevel(-1)
    s.name = 'category'
    del pm['Matched_Epitope(s) [Matched_Antigen->Matched_Organism]']
    pm = pm.join(s)
    pm_df=pd.DataFrame(list(zip(pm.Index.tolist(),pm.category.tolist(),['PM']*len(pm),[0]*len(pm),pm.Frequency.tolist())),columns=['sample','category','value','CM','PM'])
    pm_df=pm_df[pm_df.category!='']
    pm_df.drop_duplicates(inplace=True)
    pm_df.sort_values(by=['PM','sample'], ascending=False, inplace=True)
    df = pd.concat([cm_df,pm_df],ignore_index=True)
    df=df[df['category'].isin(items)]
    plot_landscape(df, items, 'tcr2ept', outprefix, extension)
    
def plot_tcr2ag(tcr2tcr_df, items, outprefix, extension):
    cm = tcr2tcr_df[tcr2tcr_df.Rank=='Complete_Match']
    pm = tcr2tcr_df[tcr2tcr_df['Index'].str.endswith('-1')]
    s = cm['Matched_Epitope(s) [Matched_Antigen->Matched_Organism]'].str.split(']').apply(Series, 1).stack()
    s.index = s.index.droplevel(-1)
    s.name = 'category'
    s = s.str.split('[').str[-1]
    del cm['Matched_Epitope(s) [Matched_Antigen->Matched_Organism]']
    cm = cm.join(s)
    cm_df=pd.DataFrame(list(zip(cm.Index.tolist(),cm.category.tolist(),['CM']*len(cm),cm.Frequency.tolist(),[0]*len(cm))),columns=['sample','category','value','CM','PM'])
    cm_df=cm_df[cm_df.category!='']
    cm_df.drop_duplicates(inplace=True)
    cm_df.sort_values(by=['CM','sample'], ascending=False, inplace=True)
    s = pm['Matched_Epitope(s) [Matched_Antigen->Matched_Organism]'].str.split(']').apply(Series, 1).stack()
    s.index = s.index.droplevel(-1)
    s.name = 'category'
    s = s.str.split('[').str[-1]
    del pm['Matched_Epitope(s) [Matched_Antigen->Matched_Organism]']
    pm = pm.join(s)
    pm_df=pd.DataFrame(list(zip(pm.Index.tolist(),pm.category.tolist(),['PM']*len(pm),[0]*len(pm),pm.Frequency.tolist())),columns=['sample','category','value','CM','PM'])
    pm_df=pm_df[pm_df.category!='']
    pm_df.drop_duplicates(inplace=True)
    pm_df.sort_values(by=['PM','sample'], ascending=False, inplace=True)
    df = pd.concat([cm_df,pm_df],ignore_index=True)
    df=df[df['category'].isin(items)]
    plot_landscape(df, items, 'tcr2ag', outprefix, extension)
    
def plot_tcr2org(tcr2tcr_df, items, outprefix, extension):
    cm = tcr2tcr_df[tcr2tcr_df.Rank=='Complete_Match']
    pm = tcr2tcr_df[tcr2tcr_df['Index'].str.endswith('-1')]
    s = cm['Matched_Epitope(s) [Matched_Antigen->Matched_Organism]'].str.split(']').apply(Series, 1).stack()
    s.index = s.index.droplevel(-1)
    s.name = 'category'
    s = s.str.split('->').str[-1]
    del cm['Matched_Epitope(s) [Matched_Antigen->Matched_Organism]']
    cm = cm.join(s)
    cm_df=pd.DataFrame(list(zip(cm.Index.tolist(),cm.category.tolist(),['CM']*len(cm),cm.Frequency.tolist(),[0]*len(cm))),columns=['sample','category','value','CM','PM'])
    cm_df=cm_df[cm_df.category!='']
    cm_df.drop_duplicates(inplace=True)
    cm_df.sort_values(by=['CM','sample'], ascending=False, inplace=True)
    s = pm['Matched_Epitope(s) [Matched_Antigen->Matched_Organism]'].str.split(']').apply(Series, 1).stack()
    s.index = s.index.droplevel(-1)
    s.name = 'category'
    s = s.str.split('->').str[-1]
    del pm['Matched_Epitope(s) [Matched_Antigen->Matched_Organism]']
    pm = pm.join(s)
    pm_df=pd.DataFrame(list(zip(pm.Index.tolist(),pm.category.tolist(),['PM']*len(pm),[0]*len(pm),pm.Frequency.tolist())),columns=['sample','category','value','CM','PM'])
    pm_df=pm_df[pm_df.category!='']
    pm_df.drop_duplicates(inplace=True)
    pm_df.sort_values(by=['PM','sample'], ascending=False, inplace=True)
    df = pd.concat([cm_df,pm_df],ignore_index=True)
    df=df[df['category'].isin(items)]
    plot_landscape(df, items, 'tcr2org', outprefix, extension)
    
def parse_tcr2tcr(infile):
    tcr2tcr = pd.read_csv(infile,sep='\t',low_memory=False,skiprows=3)
    tcr2tcr_df = tcr2tcr[tcr2tcr.Rank!='Input']
    return tcr2tcr_df

def parse_tcr2x(infile):
    tcr2x = pd.read_csv(infile,sep='\t',low_memory=False,skiprows=4)
    item_lst=tcr2x[tcr2x.columns[3]].tolist()
    item_lst = [x for x in item_lst if x!='-']
    items = [x[:x.rfind(' (')] for x in item_lst if float(x[x.rfind(', ')+2:x.rfind('%')])>1]
    items = list(dict.fromkeys(items))[::-1]
    return items
    
def plot_landscape(df, items, anno_type, outprefix, extension='.pdf'):
    if anno_type == 'tcr2org':
        plot_name = 'Top Organisms'
        height_dict = {'Top Epitopes': 6, 'Match Type': 1, 'Frequency': 3}
        figsize=(18,6)
    elif anno_type == 'tcr2ag':
        plot_name = 'Top Antigens'
        height_dict = {'Top Antigens': 8, 'Match Type': 1, 'Frequency': 3}
        figsize=(18,6)
    elif anno_type == 'tcr2ept':
        plot_name = 'Top Epitopes'
        height_dict = {'Top Epitopes': 7.5, 'Match Type': 0.5, 'Frequency': 2}
        figsize=(18,10)
    else:
        print('Unknown anno_type:',plot_name);exit(0)
    example_comut = CoMut_tcranno()
    example_comut.samples = None
    bold_10 = palettable.cartocolors.qualitative.Bold_10.mpl_colors
    item_mapping = {'PM': 'forestgreen', 'CM': bold_10[6]}
    mt_mapping = {'CM': 'gold','PM': bold_10[0]}
    mt = pd.DataFrame(list(zip(df['sample'].tolist(),['Match type']*len(df),df['value'].tolist())),columns=['sample','category','value'])
    mt.drop_duplicates(inplace=True)
    freq_mapping = {'CM': 'forestgreen', 'PM': bold_10[2]}
    df2 = df[['sample','CM','PM']].drop_duplicates()
    sum_frac_pm = [sum(df[df.category==x].PM.tolist()) for x in items]
    sum_frac_cm = [sum(df[df.category==x].CM.tolist()) for x in items]
    df3 = pd.DataFrame(list(zip(items,sum_frac_cm,sum_frac_pm)),columns=['category','CM','PM'])
    side_mapping = {'PM':bold_10[3], 'CM': 'forestgreen'}
    example_comut.add_categorical_data(df, name = plot_name, category_order = items, mapping = item_mapping)
    example_comut.add_categorical_data(mt, name = 'Match Type', mapping = mt_mapping)
    example_comut.add_bar_data(df2, name = 'Frequency', mapping = freq_mapping, stacked=True, ylabel = 'Frequency')
    example_comut.add_side_bar_data(df3, paired_name = plot_name, name = 'Sum of Fraction', position = 'right', mapping = side_mapping, stacked=True, xlabel = 'Sum of Fraction (CM, PM)')
    example_comut.plot_comut(hspace = 0.03,wspace=0.02,heights = height_dict,figsize=figsize)
    example_comut.add_unified_legend(axis_name = 'Frequency',ncol=3)
    example_comut.figure.savefig(outprefix+'_'+anno_type+extension, dpi = 300, bbox_inches = 'tight')
    
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Inputs')
    parser.add_argument('--tcr2tcr', type=str, required=True, default='', help='tcr2tcr file')
    parser.add_argument('--outprefix', type=str, required=True, default='', help='prefix of the output file')
    parser.add_argument('--tcr2org', type=str, default='', help='tcr2org annotation file')
    parser.add_argument('--tcr2ag', type=str, default='', help='tcr2ag annotation file')
    parser.add_argument('--tcr2ept', type=str, default='', help='tcr2ept annotation file')
    parser.add_argument('--anno_type', type=str, default='all',help='annotation type:one of tcr2org, tcr2ag, tcr2ept, all')
    parser.add_argument('--extension', type=str, default='.png', help='extension of the plot, one of .pdf,.svg,.png')
    args = parser.parse_args()
    infile = args.tcr2tcr
    outprefix = args.outprefix
    extension = args.extension
    if args.anno_type=='all':
        if args.tcr2ag=='' or args.tcr2org=='' or args.tcr2ept=='':
            print('required tcr2org/tcr2ag/tcr2ept anno file is not given.');exit(0)
        else:
            tcr2tcr_df = parse_tcr2tcr(infile)
            epts = parse_tcr2x(args.tcr2ept)
            ags = parse_tcr2x(args.tcr2ag)
            orgs = parse_tcr2x(args.tcr2org)
            plot_tcr2ept(tcr2tcr_df, epts, outprefix, extension)
            plot_tcr2ag(tcr2tcr_df, ags, outprefix, extension)
            plot_tcr2org(tcr2tcr_df, orgs, outprefix, extension)        
    elif args.anno_type=='tcr2ept':
        if args.tcr2ept=='':
            print('required tcr2ept anno file is not given.');exit(0)
        else:
            tcr2tcr_df = parse_tcr2tcr(infile)
            epts = parse_tcr2x(args.tcr2ept)
            plot_tcr2ept(tcr2tcr_df, epts, outprefix, extension)
    elif args.anno_type=='tcr2ag':
        if args.tcr2ag=='':
            print('required tcr2ag anno file is not given.');exit(0)
        else:
            tcr2tcr_df = parse_tcr2tcr(infile)
            ags = parse_tcr2x(args.tcr2ag)
            plot_tcr2ag(tcr2tcr_df, ags, outprefix, extension)
    elif args.anno_type=='tcr2org':
        if args.tcr2org=='':
            print('required tcr2org anno file is not given.');exit(0)
        else:
            tcr2tcr_df = parse_tcr2tcr(infile)
            orgs = parse_tcr2x(args.tcr2org)
            plot_tcr2org(tcr2tcr_df, orgs, outprefix, extension)

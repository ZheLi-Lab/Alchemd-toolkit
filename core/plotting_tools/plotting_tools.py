import matplotlib.pyplot as plt
import matplotlib.font_manager
import numpy as np
import pandas as pd
import matplotlib
import seaborn as sns
matplotlib.use('agg')
from contextlib import contextmanager

"""
specify_fonts = ['DejaVu Sans', 'FreeMono', 'NotoSans', 'FreeSans']
specify_font = 'DejaVu Sans'
avail_fonts = set([f.name for f in matplotlib.font_manager.fontManager.ttflist])
for font in specify_fonts:
    if font in avail_fonts:
        specify_font = font
        break
"""
@contextmanager
def plot_settings(figsize, font_family='DejaVu Sans', font_weight='bold', font_size=20, lineweight_dict={'bottom':4, 'left': 4, 'right': 4, 'top': 4}, x_label=None, y_label=None, y_lim_tuple=None, x_lim_tuple=None, png_file_name=None, iflegend=False, ifminorticks=True):
    plt.clf()
    figure, ax = plt.subplots(figsize=figsize)
    plt.rcParams['axes.unicode_minus'] = False#Use superscript/subscript with smaller font size
    plt.rcParams['font.sans-serif']=[font_family]#Set global font

    if ifminorticks:
        # Set axis tick thickness and font size
        plt.minorticks_on()#Enable minor ticks
        plt.tick_params(which='minor',width=2, length=4)#Set minor tick size

    plt.tick_params(which='major',width=3, length=6)#Set major tick size
    plt.yticks(fontproperties=font_family, size=font_size, weight=font_weight)#Set y-axis tick font size and boldness
    plt.xticks(fontproperties=font_family, size=font_size, weight=font_weight)#Set x-axis tick font size and boldness

    # Set x-axis and y-axis labels
    if x_label is not None:
        plt.xlabel(x_label, fontsize=font_size, fontweight=font_weight)
    if y_label is not None:
        plt.ylabel(y_label, fontsize=font_size, fontweight=font_weight)

    # Set axis line thickness
    ax=plt.gca();#Get axis handle
    bottom_lineweight = lineweight_dict['bottom']
    left_lineweight = lineweight_dict['left']
    right_lineweight = lineweight_dict['right']
    top_lineweight = lineweight_dict['top']
    ax.spines['bottom'].set_linewidth(bottom_lineweight);#Set bottom axis thickness
    ax.spines['left'].set_linewidth(left_lineweight);#Set left axis thickness
    ax.spines['right'].set_linewidth(right_lineweight);#Set right axis thickness
    ax.spines['top'].set_linewidth(top_lineweight);#Set top axis thickness

    yield # Execute plotting operations here

    # Set axis ranges
    if y_lim_tuple is not None:
        plt.ylim(y_lim_tuple[0], y_lim_tuple[1])#Limit y-axis range
    if x_lim_tuple is not None:
        plt.xlim(x_lim_tuple[0], x_lim_tuple[1])#Limit x-axis range

    # Set legend
    if iflegend:
        plt.legend(loc="best",fontsize=font_size,scatterpoints=1,shadow=True,frameon=False)#Add legend,loc控制图例位置，“best”为最佳位置，“bottom”,"top"，“topringt"等，shadow adds border shadow, frameon controls whether to have border

    plt.tight_layout()#Prevent incomplete display due to overly long or large ticklabels or titles

    # Image storage related
    if png_file_name is not None:
        plt.savefig(png_file_name, format='png', bbox_inches='tight', transparent=True, dpi=600)
        # plt.show()
        plt.clf()
    else:
        plt.clf()


class PLOTTING():
    
    def __init__(self, _specify_font='DejaVu Sans'):
        specify_fonts = ['DejaVu Sans', 'FreeMono', 'NotoSans', 'FreeSans']
        _specify_font = 'DejaVu Sans'
        avail_fonts = set([f.name for f in matplotlib.font_manager.fontManager.ttflist])
        for font in specify_fonts:
            if font in avail_fonts:
                _specify_font = font
                break
        self._specify_font = _specify_font
        pass
    
    def unpack_esti_dir(self, esti_df):
        # print(f'esti_df: {esti_df}')
        frame_ratio = np.around(np.array(esti_df.index), decimals=2)
        # print('frame_ratio')
        # print(frame_ratio)
        fe = np.array(esti_df.iloc[:,0])
        # print('fe')
        # print(fe)
        std = np.array(esti_df.iloc[:,1])
        # print('fe_std')
        # print(std)
        fe_up = fe+std
        fe_down = fe-std
        # print(frame_ratio, fe, fe_up, fe_down)
        return frame_ratio, fe, fe_up, fe_down

    def plot_dU_distribution(self, d_u,png_file=None, ifplt=True, bins_=50):
        if ifplt:
            png_file_name = png_file
        else:
            png_file_name = None
        with plot_settings(figsize=(8, 6), font_family=self._specify_font, font_weight='bold', font_size=20, lineweight_dict={'bottom':4, 'left': 4, 'right': 4, 'top': 4}, x_label='dU(kBT)', y_label='P(dU)', y_lim_tuple=None, x_lim_tuple=None, png_file_name=png_file_name, iflegend=False):
            count,bin_edges = np.histogram(d_u,bins=bins_,density=True)
            xu = []
            for i in range(len(bin_edges)-1):
                xu.append((bin_edges[i]+bin_edges[i+1])/2)
            plt.plot(xu,count,"o", color = "green", markersize=4)
        return count, xu

    def plot_weighted_dU_distribution(self, d_u, weights_array, png_file=None, ifplt=False, bins_=50):
        if ifplt:
            png_file_name = png_file
        else:
            png_file_name = None
        with plot_settings(figsize=(8, 6), font_family=self._specify_font, font_weight='bold', font_size=20, lineweight_dict={'bottom':4, 'left': 4, 'right': 4, 'top': 4}, x_label='dU(kBT)', y_label='P(dU)', y_lim_tuple=None, x_lim_tuple=None, png_file_name=png_file_name, iflegend=False):
            count,bin_edges = np.histogram(d_u,bins=bins_,density=True, weights = weights_array)
            xu = []
            for i in range(len(bin_edges)-1):
                xu.append((bin_edges[i]+bin_edges[i+1])/2)
            plt.plot(xu,count,"o", color = "green", markersize=4)
        return count, xu

    def plot_resample_dU_distribution(self, ori_d_u, resample_dU, png_file=None, ifplt=True,bins=100):
        plt.clf()
        with plot_settings(figsize=(8, 6), font_family=self._specify_font, font_weight='bold', font_size=20, lineweight_dict={'bottom':4, 'left': 4, 'right': 4, 'top': 4}, x_label='dU(kBT)', y_label='P(dU)', y_lim_tuple=None, x_lim_tuple=None, png_file_name=png_file, iflegend=True):
            pre_y,pre_bin = self.plot_dU_distribution(ori_d_u,png_file=None,ifplt=False,bins_=bins)
            resample_count, xu = self.plot_dU_distribution(resample_dU,png_file=None, ifplt=False,bins_=bins)
            plt.plot(pre_bin, pre_y,"o",label='original target dU distrib',color="green",markersize=5)
            plt.plot(xu,resample_count,"x",label='resampled target dU distrib',color="blue",markersize=5)
        # ori_d_u_mean, ori_d_u_std = np.mean(ori_d_u), np.std(ori_d_u)
        # resample_d_u_mean, resample_d_u_std = np.mean(resample_dU), np.std(resample_dU)
        # plt.text(union_x[-1],0.5*max(prob_y),'mean_pre_d_u : {:.3f}\nstd_pre_d_u : {:.5f}\nmean_resample_d_u : {:.3f}\nstd_resample_d_u : {:.5f}'.format\(ori_d_u_mean,ori_d_u_std,resample_d_u_mean,resample_d_u_std),\
        #          horizontalalignment='left',verticalalignment='center')#Text annotation for mean and std values

    def plot_Boltzmann_weight_PdU_result(self, input_interpolate_bin, plot_input_interpolate_bin, B_wei_times_PdU, interpolate_PdU, plot_Boltzmann_weight, y_lim_tuple, png_file):
        with plot_settings(figsize=(8, 6), font_family=self._specify_font, font_weight='bold', font_size=20, lineweight_dict={'bottom':4, 'left': 4, 'right': 4, 'top': 4}, x_label='dU_bin(kBT)', y_label='P(dU) or Boltzman_weight or Boltzman_weight*P(dU)', y_lim_tuple=y_lim_tuple, x_lim_tuple=None, png_file_name=png_file, iflegend=True):
            plt.plot(input_interpolate_bin,B_wei_times_PdU,"o-", color = "tab:red", markersize=5,linewidth=3, label='Boltzman_weight*P(dU)')
            plt.plot(input_interpolate_bin,interpolate_PdU,"o-", color = "tab:blue", markersize=5,linewidth=3, label='P(dU)')
            plt.plot(plot_input_interpolate_bin,plot_Boltzmann_weight,"-", color = "black", linewidth=3, label='Boltzman_weight')

    def plot_cfm_checking_result(self, df, df_std, d_u_bin_filter_nan, diff_filter_nan, input_interpolate_x, png_file):
        xlabel = r'$\Delta U_{i,i+1}\/(k_{B}T)$'
        ylabel = r'$\Delta g_{i+1,i}\/(k_{B}T)$'
        with plot_settings(figsize=(8, 6), font_family=self._specify_font, font_weight='bold', font_size=20, lineweight_dict={'bottom':4, 'left': 4, 'right': 4, 'top': 4}, x_label=xlabel, y_label=ylabel, y_lim_tuple=None, x_lim_tuple=None, png_file_name=png_file, iflegend=True):
            plt.plot(d_u_bin_filter_nan,diff_filter_nan,"o-", color = "tab:red", markersize=5,linewidth=3)
            plt.fill_between(input_interpolate_x, df-df_std, df+df_std, color='tab:blue', zorder=-1)

    def plot_fe_time_serial(self, png_file_name, **fe_he_std_dir,):
        figsize=8,6
        figure, ax = plt.subplots(figsize=figsize)
        plt.rcParams['axes.unicode_minus'] = False#Use superscript/subscript with smaller font size
        plt.rcParams['font.sans-serif']=[self._specify_font]#Set global font，可选择需要的字体替换掉‘Times New Roman’
        #plt.rcParams['font.sans-serif']=['SimHei']#Use SimHei as global font to display Chinese characters
        font1={'family': self._specify_font, 'weight': 'bold', 'size': 14}#Set font template
        font2={'family': self._specify_font, 'weight': 'bold', 'size': 20}#wight为字体的粗细，可选 ‘normal\bold\light’等
        font3={'family': self._specify_font, 'weight': 'light', 'size': 12}#size is the font size
        plt.minorticks_on()#Enable minor ticks
        plt.tick_params(which='major',width=3, length=6)#Set major tick size
        plt.tick_params(which='minor',width=2, length=4)#Set minor tick size
        ax=plt.gca();#Get axis handle
        ax.spines['bottom'].set_linewidth(4);#Set bottom axis thickness
        ax.spines['left'].set_linewidth(4);#Set left axis thickness
        ax.spines['right'].set_linewidth(4);#Set right axis thickness
        ax.spines['top'].set_linewidth(4);#Set top axis thickness
        # print(fe_he_std_dir)
        ######
        if len(fe_he_std_dir) == 5:
            dir_move = fe_he_std_dir['moving']
            dir_forw = fe_he_std_dir['forward']
            dir_reve = fe_he_std_dir['reverse']
            frame_ratio_move, fe_move, fe_up_move, fe_down_move = self.unpack_esti_dir(dir_move)
            frame_ratio_forw, fe_forw, fe_up_forw, fe_down_forw = self.unpack_esti_dir(dir_forw)
            frame_ratio_reve, fe_reve, fe_up_reve, fe_down_reve = self.unpack_esti_dir(dir_reve)
            y_min = fe_up_forw[-1]-2
            y_max = fe_up_forw[-1]+2
            ###moving estimate plot
            plt.plot(frame_ratio_move, fe_move, '-', lw=2, color='#75b84f', label='moving estimate', alpha=1)#Data body
#             plt.plot(frame_ratio_move, fe_up_move, '-', lw=0.001, color='#ffffff')
#             plt.plot(frame_ratio_move, fe_down_move, '-', lw=0.001, color='#ffffff')
            plt.fill_between(frame_ratio_move, fe_down_move, fe_up_move, where=fe_down_move <= fe_up_move,
                     facecolor='#a9f971', interpolate=True,alpha=0.5)
            ###forward estimate plot
            plt.plot(frame_ratio_forw, fe_forw, '-', lw=2, color='#C11B17', label='forward estimate', alpha=1)#Data body
#             plt.plot(frame_ratio_forw, fe_up_forw, '-', lw=0.001, color='#ffffff')
#             plt.plot(frame_ratio_forw, fe_down_forw, '-', lw=0.001, color='#ffffff')
            plt.fill_between(frame_ratio_forw, fe_down_forw, fe_up_forw, where=fe_down_forw <= fe_up_forw,
                     facecolor='#ff9a8a', interpolate=True, alpha=0.5)
            ###reverse estimate plot
            plt.plot(frame_ratio_reve, fe_reve, '-', lw=2, color='#736AFF', label='reverse estimate', alpha=1)#Data body
#             plt.plot(frame_ratio_reve, fe_up_reve, '-', lw=0.001, color='#ffffff')
#             plt.plot(frame_ratio_reve, fe_down_reve, '-', lw=0.001, color='#ffffff')
            plt.fill_between(frame_ratio_reve, fe_down_reve, fe_up_reve, where=fe_down_reve <= fe_up_reve,
                     facecolor='#a2cffe', interpolate=True, alpha=0.5)

        elif len(fe_he_std_dir) == 4:
            dir_forw = fe_he_std_dir['forward']
            dir_reve = fe_he_std_dir['reverse']
            frame_ratio_forw, fe_forw, fe_up_forw, fe_down_forw = self.unpack_esti_dir(dir_forw)
            frame_ratio_reve, fe_reve, fe_up_reve, fe_down_reve = self.unpack_esti_dir(dir_reve)
            y_min = fe_up_forw[-1]-2
            y_max = fe_up_forw[-1]+2
            ###forward estimate plot
            plt.plot(frame_ratio_forw, fe_forw, '-', lw=2, color='#C11B17', label='forward estimate', alpha=1)#Data body
            plt.fill_between(frame_ratio_forw, fe_down_forw, fe_up_forw, where=fe_down_forw <= fe_up_forw,
                     facecolor='#ff9a8a', interpolate=True, alpha=0.5)
            ###reverse estimate plot
            plt.plot(frame_ratio_reve, fe_reve, '-', lw=2, color='#736AFF', label='reverse estimate', alpha=1)#Data body
            plt.fill_between(frame_ratio_reve, fe_down_reve, fe_up_reve, where=fe_down_reve <= fe_up_reve,
                     facecolor='#a2cffe', interpolate=True, alpha=0.5)

        elif len(fe_he_std_dir) == 3:
            aly_stra = list(fe_he_std_dir.keys())[0]
            dir_ = fe_he_std_dir[aly_stra]
            frame_ratio_, fe_, fe_up_, fe_down_ = self.unpack_esti_dir(dir_)
            y_min = fe_up_[-1]-2
            y_max = fe_up_[-1]+2
            ### estimate plot
            if aly_stra == 'forward':
                plt.plot(frame_ratio_, fe_, '-', lw=2, color='#C11B17', label='forward estimate', alpha=1)#Data body
                plt.fill_between(frame_ratio_, fe_down_, fe_up_, where=fe_down_ <= fe_up_,
                        facecolor='#ff9a8a', interpolate=True, alpha=0.5)
            elif aly_stra == 'reverse':
                plt.plot(frame_ratio_, fe_, '-', lw=2, color='#736AFF', label='reverse estimate', alpha=1)#Data body
                plt.fill_between(frame_ratio_, fe_down_, fe_up_, where=fe_down_ <= fe_up_,
                        facecolor='#a2cffe', interpolate=True, alpha=0.5)
            elif aly_stra == 'moving':
                plt.plot(frame_ratio_, fe_, '-', lw=2, color='#75b84f', label='moving estimate', alpha=1)#Data body
                plt.fill_between(frame_ratio_, fe_down_, fe_up_, where=fe_down_ <= fe_up_,
                        facecolor='#a9f971', interpolate=True,alpha=0.5)
        if fe_he_std_dir['fe'] is None:
            estimate_fe = None
        else:
            fe_x_plot = np.linspace(0,1,10000)
            estimate_fe = np.array([fe_he_std_dir['fe'] for i in range(0, len(fe_x_plot))])
            estimate_std = np.array([fe_he_std_dir['fe_std'] for i in range(0, len(fe_x_plot))])
            estimate_std_range_up = estimate_fe+estimate_std
            estimate_std_range_down = estimate_fe-estimate_std
            # estimate_std_range_up = [fe_he_std_dir['fe']+fe_he_std_dir['fe_std'] for i in range(0, len(fe_x_plot))]
            # estimate_std_range_down = [fe_he_std_dir['fe']-fe_he_std_dir['fe_std'] for i in range(0, len(fe_x_plot))]
            # print(fe_x_plot)
            # print(len(estimate_fe))
            # print(len(estimate_std_range_down))
            # print(estimate_std_range_up)
        ###estimate_fe_horizontal_line
        if estimate_fe is not None:
            plt.plot(fe_x_plot, estimate_fe, '-', lw=2, color='#FFE11A', label='BAR estimate result', alpha=1)
            plt.fill_between(fe_x_plot, estimate_std_range_down, estimate_std_range_up, where=estimate_std_range_down <= estimate_std_range_up,
                    facecolor='#FFE11A', interpolate=True, alpha=0.5)
        ######
        #plt.xlim(-15,-5)#Limit x-axis range
        # print('___________ymin:',y_min,"\n____________",y_max)
        plt.ylim(y_min,y_max)#Limit y-axis range
        #plt.title("Statistical analysis on the restraint strategies used in the test system ",fontdict=font2)#Title
        #plt.xlabel(r'$\Delta G_{exp} $ (kcal/mol)', fontdict=font2)#x轴标签
        #plt.ylabel(r'$\Delta G_{cal}^{MM-GBSA} $ (kcal/mol)',fontdict=font2)#y轴标签
        plt.legend(loc="best",scatterpoints=1,prop=font2,shadow=True,frameon=False)#Add legend,\
        # # loc控制图例位置，“best”为最佳位置，“bottom”,"top"，“topringt"等，\
        # # shadow adds border shadow, frameon controls whether to have border
        plt.tick_params(\
            axis='x',#Set x-axis
            direction='in',# Minor tick direction, in/out
            which='both',      # Show major and minor ticks together, major/minor/both
            bottom=True,      #Enable bottom ticks
            top=False,         #Disable top ticks
            labelbottom=True, #Enable x-axis labels
            labelsize=20) #X-axis label size
        plt.tick_params(\
            axis='y',
            direction='in',
            which='both',
            left=True,
            right=False,
            labelbottom=True,
            labelsize=20)
        plt.yticks(fontproperties=self._specify_font, size=20,weight='bold')#Set x,y axis font size and boldness
        plt.xticks(fontproperties=self._specify_font, size=20,weight='bold')#Set x,y axis font size and boldness
        plt.ticklabel_format(axis='both',style='sci')#Scientific article style
        plt.tight_layout(rect=(0,0,1,1))#rect=[left,bottom,right,top]
        plt.tight_layout()#Prevent incomplete display due to overly long or large ticklabels or titles
        plt.savefig(png_file_name,format="png",dpi=600,transparent=True)#Save png
        #plt.show()

    def plot_heatmap(self, df, vmin=-1, vmax=1, center=0, fontsize=10, png_file=None):
        
        # xLabel = df.columns
        # yLabel = df.index
        with plot_settings(figsize=(15,11), font_family=self._specify_font, font_weight='bold', font_size=20, lineweight_dict={'bottom':4, 'left': 4, 'right': 4, 'top': 4}, x_label=None, y_label=None, y_lim_tuple=None, x_lim_tuple=None, png_file_name=png_file, iflegend=False):
            cmap = 'coolwarm'
        #     im = ax.imshow(df,cmap = 'gray' , interpolation = 'nearest',origin='lower')
            df = df.astype(float)
            # print(df)
            ax = sns.heatmap(data=df, linewidth=0.8, cmap=cmap, vmin=vmin, vmax=vmax, center=center)
            # ax.set_yticks(range(len(yLabel)))
            # ax.set_yticklabels(yLabel,fontsize=fontsize, rotation=0)
            # # ax.set_xticks(range(len(xLabel)))
            # ax.set_xticklabels(xLabel,fontsize=fontsize, rotation=45) 

    def plot_heatmap_cmap(self, df, error_max=2, png_file=None, font_size=10, lineweight_dict=None):
        if lineweight_dict is None:
            print('----')
            lineweight_dict = {'bottom':8, 'left': 8, 'right': 8, 'top': 8}

        df = df.dropna(how='all')
        df = df.astype(float)
        # Define segmented boundary values 'tab20c'
        # diff_acept = np.around(error_max/df.shape[0], decimals=3)
        diff_acept = error_max/2
        thresholds = [-2*diff_acept, -diff_acept, -1/2*diff_acept, 1/2*diff_acept, diff_acept, 2*diff_acept]
        norm = plt.Normalize(vmin=thresholds[0], vmax=thresholds[-1])
        ax = sns.heatmap(data=df, linewidth=0.4, cmap='coolwarm', cbar=True, norm=norm, xticklabels=False, yticklabels=False, linecolor='black')
        # ax = sns.heatmap(data=df, linewidth=0.4, cmap='coolwarm', cbar=True, norm=norm, )
        # Set scale font size
        cbar = ax.collections[0].colorbar
        # cbar.ax.tick_params(labelsize=20)
        bottom_lineweight = lineweight_dict['bottom']
        left_lineweight = lineweight_dict['left']
        right_lineweight = lineweight_dict['right']
        top_lineweight = lineweight_dict['top']
        ax.figure.axes[-1].spines['bottom'].set_linewidth(bottom_lineweight)
        ax.figure.axes[-1].spines['left'].set_linewidth(left_lineweight)
        ax.figure.axes[-1].spines['right'].set_linewidth(right_lineweight)
        ax.figure.axes[-1].spines['top'].set_linewidth(top_lineweight)

        if png_file is not None:
            plt.savefig(png_file, format='png', bbox_inches='tight', transparent=True, dpi=600)
            # plt.show()
            plt.clf()
        else:
            plt.clf()
        # for spine_name, spine in ax.spines.items():
        #     spine.set_linewidth(lineweight_dict[spine_name])
            # Set scale tick values
            # cbar.set_ticks(thresholds)
            # cbar.set_ticklabels(thresholds)



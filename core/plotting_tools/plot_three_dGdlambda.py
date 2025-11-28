from dG_dlambda import *

lambda_col_name = "delta_A_what_to_what"
fe_col_name="free_energy(kcal/mol)"

class All_sys_data():
    def __init__(self, csv_dict):
        self.df_dict = {}
        for key, value in csv_dict.items():
            dG_dlamda_obj = ABFE_dG_dlambda(value, 'BAR')
            dG_dlamda_obj.get_dG_dlambda(lambda_col_name, fe_col_name)
            df_csv = dG_dlamda_obj.df_csv.drop(dG_dlamda_obj.df_csv[dG_dlamda_obj.df_csv['lambda_labels']=='mix'].index)
            df_crg = df_csv[df_csv['lambda_labels']=='charge']
            self.df_dict[key] = df_crg
    
    def plot_dG_dlambda(self,png_file):
        plt.clf()
        for key in self.df_dict.keys():
            df = self.df_dict[key]
            plt.plot(df.index, df['dG_dlambda'], label=key)
        # plt.xticks(df.index, [ np.around(i, 2) for i in np.arange(0, 1.0, 0.01)], rotation=90, fontsize=5)
        plt.ylabel(r'$\delta G_{cal}/\delta \lambda $ (kcal/mol)', )#y-axis label
        plt.xlabel(r'$\lambda$ info',)#x-axis label
        plt.legend() 
        # plt.minorticks_on()#Enable minor ticks
        # plt.tick_params(which='major',width=3, length=6)#Set major tick size
        # plt.tick_params(which='minor',width=2, length=4)#Set minor tick size
        ax=plt.gca();#Get axis handle
        ax.spines['bottom'].set_linewidth(4);#Set bottom axis thickness
        ax.spines['left'].set_linewidth(4);#Set left axis thickness
        ax.spines['right'].set_linewidth(4);#Set right axis thickness
        ax.spines['top'].set_linewidth(4);#Set top axis thickness
        plt.tight_layout()
            # Display figure
        if png_file is not None:
            plt.savefig(png_file, format="png",dpi=600, transparent=True)
        else:
            plt.show()

if __name__ == '__main__':
    pass
    
